/******************************************************************************
 * $Id: ogrgeojsonutils.cpp 23293 2011-10-30 11:11:39Z rouault $
 *
 * Project:  OpenGIS Simple Features Reference Implementation
 * Purpose:  Implementation of private utilities used within OGR GeoJSON Driver.
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
#include "ogrgeojsonutils.h"
#include <cpl_port.h>
#include <cpl_conv.h>
#include <ogr_geometry.h>
#include <jsonc/json.h> // JSON-C

/************************************************************************/
/*                           GeoJSONIsObject()                          */
/************************************************************************/

int GeoJSONIsObject( const char* pszText )
{
    if( NULL == pszText )
        return FALSE;

/* -------------------------------------------------------------------- */
/*      This is a primitive test, but we need to perform it fast.       */
/* -------------------------------------------------------------------- */
    while( *pszText != '\0' && isspace( (unsigned char)*pszText ) )
        pszText++;

    if( EQUALN( pszText, "{", 1) )
        return TRUE;

    return FALSE;
}

/************************************************************************/
/*                           GeoJSONFileIsObject()                      */
/************************************************************************/

int GeoJSONFileIsObject( const char* pszSource ) 
{ 
    CPLAssert( NULL != pszSource ); 
 
    VSILFILE* fp = NULL; 
    fp = VSIFOpenL( pszSource, "rb" ); 
    if( NULL == fp ) 
    { 
        return FALSE; 
    } 
    
    // by default read first 6000 bytes 
    // 6000 was chosen as enough bytes to  
    // enable all current tests to pass 
    vsi_l_offset nToReadLen = 6000; 
    vsi_l_offset nDataLen = 0; 

    char* pszGeoData = (char*)VSIMalloc((size_t)(nToReadLen + 1)); 
    if( NULL == pszGeoData ) 
    { 
        VSIFCloseL(fp); 
        return FALSE; 
    } 

    nDataLen = VSIFReadL( pszGeoData, 1, (size_t)nToReadLen, fp );
    pszGeoData[nDataLen] = '\0'; 
    if( nDataLen == 0 ) 
    { 
        VSIFCloseL( fp ); 
        CPLFree( pszGeoData ); 
        return FALSE; 
    } 
    VSIFCloseL( fp ); 
    
    char* pszIter = pszGeoData;
    while( *pszIter != '\0' && isspace( (unsigned char)*pszIter ) )
        pszIter++;

    if( *pszIter != '{' )
    {
        CPLFree( pszGeoData ); 
        return FALSE;
    }

    int bRet = FALSE; 

    if ((strstr(pszGeoData, "\"type\"") != NULL && strstr(pszGeoData, "\"coordinates\"") != NULL) 
        || strstr(pszGeoData, "\"FeatureCollection\"") != NULL
        || (strstr(pszGeoData, "\"geometryType\"") != NULL && strstr(pszGeoData, "\"esriGeometry") != NULL)) 
        bRet = TRUE; 
     
    CPLFree( pszGeoData ); 
    return bRet; 
} 

/************************************************************************/
/*                           GeoJSONGetSourceType()                     */
/************************************************************************/

GeoJSONSourceType GeoJSONGetSourceType( const char* pszSource )
{
    GeoJSONSourceType srcType = eGeoJSONSourceUnknown;

    // NOTE: Sometimes URL ends with .geojson token, for example
    //       http://example/path/2232.geojson
    //       It's important to test beginning of source first.
    if ( eGeoJSONProtocolUnknown != GeoJSONGetProtocolType( pszSource ) )
    {
        srcType = eGeoJSONSourceService;
    }
    else if( EQUAL( CPLGetExtension( pszSource ), "geojson" )
             || EQUAL( CPLGetExtension( pszSource ), "json" )
             || ((EQUALN( pszSource, "/vsigzip/", 9) || EQUALN( pszSource, "/vsizip/", 8)) &&
                 (strstr( pszSource, ".json") || strstr( pszSource, ".JSON") ||
                  strstr( pszSource, ".geojson") || strstr( pszSource, ".GEOJSON")) ))
    {
        srcType = eGeoJSONSourceFile;
    }
    else if( GeoJSONIsObject( pszSource ) )
    {
        srcType = eGeoJSONSourceText;
    }
    else if( GeoJSONFileIsObject( pszSource ) )
    {
        srcType = eGeoJSONSourceFile;
    }

    return srcType;
}

/************************************************************************/
/*                           GeoJSONGetProtocolType()                   */
/************************************************************************/

GeoJSONProtocolType GeoJSONGetProtocolType( const char* pszSource )
{
    GeoJSONProtocolType ptclType = eGeoJSONProtocolUnknown;

    if( EQUALN( pszSource, "http:", 5 ) )
        ptclType = eGeoJSONProtocolHTTP;
    else if( EQUALN( pszSource, "https:", 6 ) )
        ptclType = eGeoJSONProtocolHTTPS;
    else if( EQUALN( pszSource, "ftp:", 4 ) )
        ptclType = eGeoJSONProtocolFTP;

    return ptclType;
}

/************************************************************************/
/*                           GeoJSONPropertyToFieldType()               */
/************************************************************************/

OGRFieldType GeoJSONPropertyToFieldType( json_object* poObject )
{
    if (poObject == NULL) { return OFTString; }

    json_type type = json_object_get_type( poObject );

    if( json_type_boolean == type )
        return OFTInteger;
    else if( json_type_double == type )
        return OFTReal;
    else if( json_type_int == type )
        return OFTInteger;
    else if( json_type_string == type )
        return OFTString;
    else if( json_type_array == type )
    {
        int nSize = json_object_array_length(poObject);
        if (nSize == 0)
            return OFTStringList; /* we don't know, so let's assume it's a string list */
        OGRFieldType eType = OFTIntegerList;
        for(int i=0;i<nSize;i++)
        {
            json_object* poRow = json_object_array_get_idx(poObject, i);
            if (poRow != NULL)
            {
                type = json_object_get_type( poRow );
                if (type == json_type_string)
                    return OFTStringList;
                else if (type == json_type_double)
                    eType = OFTRealList;
                else if (type != json_type_int &&
                         type != json_type_boolean)
                    return OFTString;
            }
        }
        return eType;
    }
    else
        return OFTString; /* null, object */
}

/************************************************************************/
/*                           OGRGeoJSONGetGeometryName()                */
/************************************************************************/

const char* OGRGeoJSONGetGeometryName( OGRGeometry const* poGeometry )
{
    CPLAssert( NULL != poGeometry );
    
    OGRwkbGeometryType eType = poGeometry->getGeometryType();

    if( wkbPoint == eType || wkbPoint25D == eType )
        return "Point";
    else if( wkbLineString == eType || wkbLineString25D == eType )
        return "LineString";
    else if( wkbPolygon == eType || wkbPolygon25D == eType )
        return "Polygon";
    else if( wkbMultiPoint == eType || wkbMultiPoint25D == eType )
        return "MultiPoint";
    else if( wkbMultiLineString == eType || wkbMultiLineString25D == eType )
        return "MultiLineString";
    else if( wkbMultiPolygon == eType || wkbMultiPolygon25D == eType )
        return "MultiPolygon";
    else if( wkbGeometryCollection == eType || wkbGeometryCollection25D == eType )
        return "GeometryCollection";
    else
        return "Unknown";
}

