/******************************************************************************
 * $Id: ogrgeojsondatasource.cpp 23367 2011-11-12 22:46:13Z rouault $
 *
 * Project:  OpenGIS Simple Features Reference Implementation
 * Purpose:  Implementation of OGRGeoJSONDataSource class (OGR GeoJSON Driver).
 * Author:   Mateusz Loskot, mateusz@loskot.net
 *
 ******************************************************************************
 * Copyright (c) 2007, Mateusz Loskot
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ****************************************************************************/
#include "ogr_geojson.h"
#include "ogrgeojsonutils.h"
#include "ogrgeojsonreader.h"
#include <cpl_http.h>
#include <jsonc/json.h> // JSON-C
#include <cstddef>
#include <cstdlib>
using namespace std;

/************************************************************************/
/*                           OGRGeoJSONDataSource()                     */
/************************************************************************/

OGRGeoJSONDataSource::OGRGeoJSONDataSource()
    : pszName_(NULL), pszGeoData_(NULL),
        papoLayers_(NULL), nLayers_(0), fpOut_(NULL),
        flTransGeom_( OGRGeoJSONDataSource::eGeometryPreserve ),
        flTransAttrs_( OGRGeoJSONDataSource::eAtributesPreserve ),
        bFpOutputIsSeekable_( FALSE ),
        nBBOXInsertLocation_(0)
{
    // I've got constructed. Lunch time!
}

/************************************************************************/
/*                           ~OGRGeoJSONDataSource()                    */
/************************************************************************/

OGRGeoJSONDataSource::~OGRGeoJSONDataSource()
{
    Clear();
    
    if( NULL != fpOut_ )
    {
        VSIFCloseL( fpOut_ );
        fpOut_ = NULL;
    }
}

/************************************************************************/
/*                           Open()                                     */
/************************************************************************/

int OGRGeoJSONDataSource::Open( const char* pszName )
{
    CPLAssert( NULL != pszName );

/* -------------------------------------------------------------------- */
/*      Release resources allocated during previous request.            */
/* -------------------------------------------------------------------- */
    if( NULL != papoLayers_ )
    {
        CPLAssert( nLayers_ > 0 );
        Clear();
    }

/* -------------------------------------------------------------------- */
/*      Determine type of data source: text file (.geojson, .json),     */
/*      Web Service or text passed directly and load data.              */
/* -------------------------------------------------------------------- */
    GeoJSONSourceType nSrcType;
    
    nSrcType = GeoJSONGetSourceType( pszName );
    if( eGeoJSONSourceService == nSrcType )
    {
        if( (strstr(pszName, "SERVICE=WFS") || strstr(pszName, "service=WFS") ||
             strstr(pszName, "service=wfs")) && !strstr(pszName, "json"))
            return FALSE;

        if( !ReadFromService( pszName ) )
            return FALSE;
    }
    else if( eGeoJSONSourceText == nSrcType )
    {
        pszGeoData_ = CPLStrdup( pszName );
    }
    else if( eGeoJSONSourceFile == nSrcType )
    {
        if( !ReadFromFile( pszName ) )
            return FALSE;
    }
    else
    {
        Clear();
        return FALSE;
    }

/* -------------------------------------------------------------------- */
/*      Construct OGR layer and feature objects from                    */
/*      GeoJSON text tree.                                              */
/* -------------------------------------------------------------------- */
    if( NULL == pszGeoData_ ||
        strncmp(pszGeoData_, "{\"couchdb\":\"Welcome\"", strlen("{\"couchdb\":\"Welcome\"")) == 0 ||
        strncmp(pszGeoData_, "{\"db_name\":\"", strlen("{\"db_name\":\"")) == 0 ||
        strncmp(pszGeoData_, "{\"total_rows\":", strlen("{\"total_rows\":")) == 0 ||
        strncmp(pszGeoData_, "{\"rows\":[", strlen("{\"rows\":[")) == 0)
    {
        Clear();
        return FALSE;
    }

    OGRGeoJSONLayer* poLayer = LoadLayer();
    if( NULL == poLayer )
    {
        Clear();
        
        CPLError( CE_Failure, CPLE_OpenFailed, 
                  "Failed to read GeoJSON data" );
        return FALSE;
    }

    poLayer->DetectGeometryType();

/* -------------------------------------------------------------------- */
/*      NOTE: Currently, the driver generates only one layer per        */
/*      single GeoJSON file, input or service request.                  */
/* -------------------------------------------------------------------- */
    const int nLayerIndex = 0;
    nLayers_ = 1;
    
    papoLayers_ =
        (OGRGeoJSONLayer**)CPLMalloc( sizeof(OGRGeoJSONLayer*) * nLayers_ );
    papoLayers_[nLayerIndex] = poLayer; 

    CPLAssert( NULL != papoLayers_ );
    CPLAssert( nLayers_ > 0 );
    return TRUE;
}

/************************************************************************/
/*                           GetName()                                  */
/************************************************************************/

const char* OGRGeoJSONDataSource::GetName()
{
    return pszName_ ? pszName_ : "";
}

/************************************************************************/
/*                           GetLayerCount()                            */
/************************************************************************/

int OGRGeoJSONDataSource::GetLayerCount()
{
    return nLayers_;
}

/************************************************************************/
/*                           GetLayer()                                 */
/************************************************************************/

OGRLayer* OGRGeoJSONDataSource::GetLayer( int nLayer )
{
    if( 0 <= nLayer || nLayer < nLayers_ )
    {
        CPLAssert( NULL != papoLayers_[nLayer] );

        OGRLayer* poLayer = papoLayers_[nLayer];

        /* Return layer in readable state. */
        poLayer->ResetReading();
        return poLayer;
    }

    return NULL;
}

/************************************************************************/
/*                            CreateLayer()                             */
/************************************************************************/

OGRLayer* OGRGeoJSONDataSource::CreateLayer( const char* pszName_,
                                             OGRSpatialReference* poSRS,
                                             OGRwkbGeometryType eGType,
                                             char** papszOptions )
{
    OGRGeoJSONLayer* poLayer = NULL;
    poLayer = new OGRGeoJSONLayer( pszName_, poSRS, eGType, papszOptions, this );

/* -------------------------------------------------------------------- */
/*      Add layer to data source layer list.                            */
/* -------------------------------------------------------------------- */
    
    // TOOD: Waiting for multi-layer support
    if ( nLayers_ != 0 )
    {
        CPLError(CE_Failure, CPLE_NotSupported,
                 "GeoJSON driver doesn't support creating more than one layer");
        return NULL;
    }

    papoLayers_ = (OGRGeoJSONLayer **)
        CPLRealloc( papoLayers_,  sizeof(OGRGeoJSONLayer*) * (nLayers_ + 1) );
    
    papoLayers_[nLayers_++] = poLayer;

    if( NULL != fpOut_ )
    {
        VSIFPrintfL( fpOut_, "{\n\"type\": \"FeatureCollection\",\n" );

        if (bFpOutputIsSeekable_)
        {
            nBBOXInsertLocation_ = (int) VSIFTellL( fpOut_ );

            char szSpaceForBBOX[SPACE_FOR_BBOX+1];
            memset(szSpaceForBBOX, ' ', SPACE_FOR_BBOX);
            szSpaceForBBOX[SPACE_FOR_BBOX] = '\0';
            VSIFPrintfL( fpOut_, "%s\n", szSpaceForBBOX);
        }

        VSIFPrintfL( fpOut_, "\"features\": [\n" );
    }

    return poLayer;
}

/************************************************************************/
/*                           TestCapability()                           */
/************************************************************************/

int OGRGeoJSONDataSource::TestCapability( const char* pszCap )
{
    if( EQUAL( pszCap, ODsCCreateLayer ) )
        return TRUE;
    else if( EQUAL( pszCap, ODsCDeleteLayer ) )
        return FALSE;
    else
        return FALSE;
}

int OGRGeoJSONDataSource::Create( const char* pszName, char** papszOptions )
{
    UNREFERENCED_PARAM(papszOptions);

    CPLAssert( NULL == fpOut_ );

    if (strcmp(pszName, "/dev/stdout") == 0)
        pszName = "/vsistdout/";

    bFpOutputIsSeekable_ =  !(strcmp(pszName,"/vsistdout/") == 0 ||
                              strncmp(pszName,"/vsigzip/", 9) == 0 ||
                              strncmp(pszName,"/vsizip/", 8) == 0);

/* -------------------------------------------------------------------- */
/*     File overwrite not supported.                                    */
/* -------------------------------------------------------------------- */
    VSIStatBufL sStatBuf;
    if( 0 == VSIStatL( pszName, &sStatBuf ) )
    {
        CPLError( CE_Failure, CPLE_NotSupported,
                  "The GeoJSON driver does not overwrite existing files." );
        return FALSE;
    }
    
/* -------------------------------------------------------------------- */
/*      Create the output file.                                         */
/* -------------------------------------------------------------------- */
    fpOut_ = VSIFOpenL( pszName, "w" );
    if( NULL == fpOut_)
    {
        CPLError( CE_Failure, CPLE_OpenFailed, 
                  "Failed to create GeoJSON datasource: %s.", 
                  pszName );
        return FALSE;
    }

    pszName_ = CPLStrdup( pszName );

    return TRUE;
}

/************************************************************************/
/*                           SetGeometryTranslation()                   */
/************************************************************************/

void
OGRGeoJSONDataSource::SetGeometryTranslation( GeometryTranslation type )
{
    flTransGeom_ = type;
}

/************************************************************************/
/*                           SetAttributesTranslation()                 */
/************************************************************************/

void OGRGeoJSONDataSource::SetAttributesTranslation( AttributesTranslation type )
{
    flTransAttrs_ = type;
}

/************************************************************************/
/*                  PRIVATE FUNCTIONS IMPLEMENTATION                    */
/************************************************************************/

void OGRGeoJSONDataSource::Clear()
{
    for( int i = 0; i < nLayers_; i++ )
    {
        CPLAssert( NULL != papoLayers_ );
        delete papoLayers_[i];
    }

    CPLFree( papoLayers_ );
    papoLayers_ = NULL;
    nLayers_ = 0;

    CPLFree( pszName_ );
    pszName_ = NULL;

    CPLFree( pszGeoData_ );
    pszGeoData_ = NULL;

    if( NULL != fpOut_ )
    {
        VSIFCloseL( fpOut_ );
    }
    fpOut_ = NULL;
}

/************************************************************************/
/*                           ReadFromFile()                             */
/************************************************************************/

int OGRGeoJSONDataSource::ReadFromFile( const char* pszSource )
{
    CPLAssert( NULL == pszGeoData_ );

    if( NULL == pszSource )
    {
        CPLDebug( "GeoJSON", "Input file path is null" );
        return FALSE;
    }

    VSILFILE* fp = NULL;
    fp = VSIFOpenL( pszSource, "rb" );
    if( NULL == fp )
    {
        CPLDebug( "GeoJSON", "Failed to open input file '%s'", pszSource );
        return FALSE;
    }

    vsi_l_offset nDataLen = 0;

    VSIFSeekL( fp, 0, SEEK_END );
    nDataLen = VSIFTellL( fp );

    // With "large" VSI I/O API we can read data chunks larger than VSIMalloc
    // could allocate. Catch it here.
    if ( nDataLen > (vsi_l_offset)(size_t)nDataLen )
    {
        CPLDebug( "GeoJSON", "Input file too large to be opened" );
        VSIFCloseL( fp );
        return FALSE;
    }

    VSIFSeekL( fp, 0, SEEK_SET );

    pszGeoData_ = (char*)VSIMalloc((size_t)(nDataLen + 1));
    if( NULL == pszGeoData_ )
    {
        VSIFCloseL(fp);
        return FALSE;
    }

    pszGeoData_[nDataLen] = '\0';
    if( ( nDataLen != VSIFReadL( pszGeoData_, 1, (size_t)nDataLen, fp ) ) )
    {
        Clear();
        VSIFCloseL( fp );
        return FALSE;
    }
    VSIFCloseL( fp );

    pszName_ = CPLStrdup( pszSource );

    CPLAssert( NULL != pszGeoData_ );
    return TRUE;
}

/************************************************************************/
/*                           ReadFromService()                          */
/************************************************************************/

int OGRGeoJSONDataSource::ReadFromService( const char* pszSource )
{
    CPLAssert( NULL == pszGeoData_ );
    CPLAssert( NULL != pszSource );

    if( eGeoJSONProtocolUnknown == GeoJSONGetProtocolType( pszSource ) )
    {
        CPLDebug( "GeoJSON", "Unknown service type (use HTTP, HTTPS, FTP)" );
        return FALSE;
    }

/* -------------------------------------------------------------------- */
/*      Fetch the GeoJSON result.                                        */
/* -------------------------------------------------------------------- */
    CPLErrorReset();

    CPLHTTPResult* pResult = NULL;
    char* papsOptions[] = { (char*) "HEADERS=Accept: text/plain Accept: application/json", NULL };

    pResult = CPLHTTPFetch( pszSource, papsOptions );

/* -------------------------------------------------------------------- */
/*      Try to handle CURL/HTTP errors.                                 */
/* -------------------------------------------------------------------- */
    if( NULL == pResult
        || 0 == pResult->nDataLen || 0 != CPLGetLastErrorNo() )
    {
        CPLHTTPDestroyResult( pResult );
        return FALSE;
    }

   if( 0 != pResult->nStatus )
    {
        CPLError( CE_Failure, CPLE_AppDefined, 
                  "Curl reports error: %d: %s",
                  pResult->nStatus, pResult->pszErrBuf );
        CPLHTTPDestroyResult( pResult );
        return FALSE;
    }

/* -------------------------------------------------------------------- */
/*      Copy returned GeoJSON data to text buffer.                      */
/* -------------------------------------------------------------------- */
    const char* pszData = reinterpret_cast<char*>(pResult->pabyData);

    if ( eGeoJSONProtocolUnknown != GeoJSONGetProtocolType( pszData ) )
    {
        CPLError( CE_Failure, CPLE_AppDefined, 
            "The data that was downloaded also starts with "
            "protocol prefix (http://, https:// or ftp://) "
            "and cannot be processed as GeoJSON data.");
        CPLHTTPDestroyResult( pResult );
        return FALSE;
    }

    // TODO: Eventually, CPLHTTPResult::pabyData could be assigned
    //       to pszGeoData_, so we will avoid copying of potentially (?) big data.
    pszGeoData_ = (char*)VSIMalloc( sizeof(char) * pResult->nDataLen + 1 );
    if( NULL == pszGeoData_ )
    {
        CPLHTTPDestroyResult( pResult );
        return FALSE;
    }

    strncpy( pszGeoData_, pszData, pResult->nDataLen );
    pszGeoData_[pResult->nDataLen] = '\0';

    pszName_ = CPLStrdup( pszSource );

/* -------------------------------------------------------------------- */
/*      Cleanup HTTP resources.                                         */
/* -------------------------------------------------------------------- */
    CPLHTTPDestroyResult( pResult );

    CPLAssert( NULL != pszGeoData_ );
    return TRUE;
}

/************************************************************************/
/*                           LoadLayer()                          */
/************************************************************************/

OGRGeoJSONLayer* OGRGeoJSONDataSource::LoadLayer()
{
    if( NULL == pszGeoData_ )
    {
        CPLError( CE_Failure, CPLE_ObjectNull,
                  "GeoJSON data buffer empty" );
        return NULL;
    }

    if ( !GeoJSONIsObject( pszGeoData_) )
    {
        CPLDebug( "GeoJSON", "No valid GeoJSON data found in source '%s'", pszName_ );
        return NULL;
    }

    OGRErr err = OGRERR_NONE;
    OGRGeoJSONLayer* poLayer = NULL;

/* -------------------------------------------------------------------- */
/*      Is it ESRI Feature Service data ?                               */
/* -------------------------------------------------------------------- */
    if ( strstr(pszGeoData_, "esriGeometry") ||
         strstr(pszGeoData_, "esriFieldTypeOID") )
    {
        OGRESRIJSONReader reader;
        err = reader.Parse( pszGeoData_ );
        if( OGRERR_NONE == err )
        {
            // TODO: Think about better name selection
            poLayer = reader.ReadLayer( OGRGeoJSONLayer::DefaultName, this );
        }

        return poLayer;
    }
    
/* -------------------------------------------------------------------- */
/*      Configure GeoJSON format translator.                            */
/* -------------------------------------------------------------------- */
    OGRGeoJSONReader reader;

    if( eGeometryAsCollection == flTransGeom_ )
    {
        reader.SetPreserveGeometryType( false );
        CPLDebug( "GeoJSON", "Geometry as OGRGeometryCollection type." );
    }
    
    if( eAtributesSkip == flTransAttrs_ )
    {
        reader.SetSkipAttributes( true );
        CPLDebug( "GeoJSON", "Skip all attributes." );
    }
    
/* -------------------------------------------------------------------- */
/*      Parse GeoJSON and build valid OGRLayer instance.                */
/* -------------------------------------------------------------------- */
    err = reader.Parse( pszGeoData_ );
    if( OGRERR_NONE == err )
    {
        // TODO: Think about better name selection
        poLayer = reader.ReadLayer( OGRGeoJSONLayer::DefaultName, this );
    }

    return poLayer;
}
