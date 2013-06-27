/******************************************************************************
 *
 * Project:  KML Translator
 * Purpose:  Implements OGRLIBKMLDriver
 * Author:   Brian Case, rush at winkey dot org
 *
 ******************************************************************************
 * Copyright (c) 2010, Brian Case
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
 *****************************************************************************/

#include "ogr_libkml.h"
//#include "cpl_conv.h"
//#include "cpl_string.h"
#include "cpl_error.h"

#include <kml/dom.h>

using kmldom::KmlFactory;
using kmldom::PlacemarkPtr;
using kmldom::Placemark;
using kmldom::DocumentPtr;
using kmldom::ContainerPtr;
using kmldom::FeaturePtr;
using kmldom::GroundOverlayPtr;
using kmldom::KmlPtr;
using kmldom::Kml;
using kmlengine::KmzFile;
using kmlengine::KmlFile;
using kmlengine::Bbox;
using kmldom::ExtendedDataPtr;
using kmldom::SchemaDataPtr;
using kmldom::DataPtr;

#include "ogrlibkmlfeature.h"
#include "ogrlibkmlfield.h"
#include "ogrlibkmlstyle.h"

/******************************************************************************
 OGRLIBKMLLayer constructor

 Args:          pszLayerName    the name of the layer
                poSpatialRef    the spacial Refrance for the layer
                eGType          the layers geometry type
                poOgrDS         pointer to the datasource the layer is in
                poKmlRoot       pointer to the root kml element of the layer
                pszFileName     the filename of the layer
                bNew            true if its a new layer
                bUpdate         true if the layer is writeable
 
 Returns:       nothing
                
******************************************************************************/

OGRLIBKMLLayer::OGRLIBKMLLayer ( const char *pszLayerName,
                                 OGRSpatialReference * poSpatialRef,
                                 OGRwkbGeometryType eGType,
                                 OGRLIBKMLDataSource * poOgrDS,
                                 ElementPtr poKmlRoot,
                                 ContainerPtr poKmlContainer,
                                 const char *pszFileName,
                                 int bNew,
                                 int bUpdate )
{

    m_poStyleTable = NULL;
    iFeature = 0;
    nFeatures = 0;
    nFID = 1;

    this->bUpdate = bUpdate;
    m_pszName = CPLStrdup ( pszLayerName );
    m_pszFileName = CPLStrdup ( pszFileName );
    m_poOgrDS = poOgrDS;

    m_poOgrSRS = new OGRSpatialReference ( NULL );
    m_poOgrSRS->SetWellKnownGeogCS ( "WGS84" );

    m_poOgrFeatureDefn = new OGRFeatureDefn ( pszLayerName );
    m_poOgrFeatureDefn->Reference (  );
    m_poOgrFeatureDefn->SetGeomType ( eGType );

    /***** store the root element pointer *****/

    m_poKmlLayerRoot = poKmlRoot;
    
    /***** store the layers container *****/

    m_poKmlLayer = poKmlContainer;

    /***** was the layer created from a DS::Open *****/

    m_bReadGroundOverlay = CSLTestBoolean(CPLGetConfigOption("LIBKML_READ_GROUND_OVERLAY", "YES"));

    if ( !bNew ) {

        /***** get the number of features on the layer *****/

        nFeatures = m_poKmlLayer->get_feature_array_size (  );

        /***** get the field config *****/
        
        struct fieldconfig oFC;
        get_fieldconfig( &oFC );

        /***** name field *****/
        
        OGRFieldDefn oOgrFieldName ( oFC.namefield,OFTString );
        m_poOgrFeatureDefn->AddFieldDefn ( &oOgrFieldName );

        /***** descripton field *****/
        
        OGRFieldDefn oOgrFieldDesc ( oFC.descfield, OFTString );
        m_poOgrFeatureDefn->AddFieldDefn ( &oOgrFieldDesc );

        /***** timestamp field *****/

        OGRFieldDefn oOgrFieldTs ( oFC.tsfield, OFTDateTime );
        m_poOgrFeatureDefn->AddFieldDefn ( &oOgrFieldTs );

        /*****  timespan begin field *****/

        OGRFieldDefn oOgrFieldBegin ( oFC.beginfield, OFTDateTime );
        m_poOgrFeatureDefn->AddFieldDefn ( &oOgrFieldBegin );

        /*****  timespan end field *****/

        OGRFieldDefn oOgrFieldEnd ( oFC.endfield, OFTDateTime );
        m_poOgrFeatureDefn->AddFieldDefn ( &oOgrFieldEnd );

        /*****  altitudeMode field *****/

        OGRFieldDefn oOgrFieldAltitudeMode ( oFC.altitudeModefield, OFTString );
        m_poOgrFeatureDefn->AddFieldDefn ( &oOgrFieldAltitudeMode );

        /***** tessellate field *****/

        OGRFieldDefn oOgrFieldTessellate ( oFC.tessellatefield, OFTInteger );
        m_poOgrFeatureDefn->AddFieldDefn ( &oOgrFieldTessellate );

        /***** extrude field *****/

        OGRFieldDefn oOgrFieldExtrude ( oFC.extrudefield, OFTInteger );
        m_poOgrFeatureDefn->AddFieldDefn ( &oOgrFieldExtrude );

        /***** visibility field *****/

        OGRFieldDefn oOgrFieldVisibility ( oFC.visibilityfield, OFTInteger );
        m_poOgrFeatureDefn->AddFieldDefn ( &oOgrFieldVisibility );

        /***** draw order field *****/

        OGRFieldDefn oOgrFieldDrawOrder ( oFC.drawOrderfield, OFTInteger );
        m_poOgrFeatureDefn->AddFieldDefn ( &oOgrFieldDrawOrder );

        /***** icon field *****/

        OGRFieldDefn oOgrFieldIcon ( oFC.iconfield, OFTString );
        m_poOgrFeatureDefn->AddFieldDefn ( &oOgrFieldIcon );

        /***** get the styles *****/

        if ( m_poKmlLayer->IsA ( kmldom::Type_Document ) )
            ParseStyles ( AsDocument ( m_poKmlLayer ), &m_poStyleTable );

        /***** get the schema if the layer is a Document *****/

        m_poKmlSchema = NULL;

        if ( m_poKmlLayer->IsA ( kmldom::Type_Document ) ) {
            DocumentPtr poKmlDocument = AsDocument ( m_poKmlLayer );

            if ( poKmlDocument->get_schema_array_size (  ) ) {
                m_poKmlSchema = poKmlDocument->get_schema_array_at ( 0 );
                kml2FeatureDef ( m_poKmlSchema, m_poOgrFeatureDefn );
            }
        }

        /***** the schema is somewhere else *****/

        if (m_poKmlSchema == NULL) {

            /***** try to find the correct schema *****/

            FeaturePtr poKmlFeature;

            /***** find the first placemark *****/

            do {
                if ( iFeature >= nFeatures )
                    break;

                poKmlFeature =
                    m_poKmlLayer->get_feature_array_at ( iFeature++ );

            } while ( poKmlFeature->Type (  ) != kmldom::Type_Placemark );

            if ( iFeature <= nFeatures && poKmlFeature &&
                 poKmlFeature->Type (  ) == kmldom::Type_Placemark &&
                 poKmlFeature->has_extendeddata (  ) ) {

                ExtendedDataPtr poKmlExtendedData = poKmlFeature->
                    get_extendeddata (  );

                if ( poKmlExtendedData->get_schemadata_array_size (  ) > 0 ) {
                    SchemaDataPtr poKmlSchemaData = poKmlExtendedData->
                        get_schemadata_array_at ( 0 );

                    if ( poKmlSchemaData->has_schemaurl (  ) ) {

                        std::string oKmlSchemaUrl = poKmlSchemaData->
                            get_schemaurl (  );
                        if ( ( m_poKmlSchema =
                               m_poOgrDS->FindSchema ( oKmlSchemaUrl.
                                                       c_str (  ) ) ) ) {
                            kml2FeatureDef ( m_poKmlSchema,
                                             m_poOgrFeatureDefn );
                        }
                    }
                }
                else if ( poKmlExtendedData->get_data_array_size() > 0 )
                {
                    /* Use the <Data> of the first placemark to build the feature definition */
                    /* If others have different fields, too bad... */
                    int bLaunderFieldNames =
                        CSLTestBoolean(CPLGetConfigOption("LIBKML_LAUNDER_FIELD_NAMES", "YES"));
                    size_t nDataArraySize = poKmlExtendedData->get_data_array_size();
                    for(size_t i=0; i < nDataArraySize; i++)
                    {
                        const DataPtr& data = poKmlExtendedData->get_data_array_at(i);
                        if (data->has_name())
                        {
                            CPLString osName = data->get_name();
                            if (bLaunderFieldNames)
                                osName = LaunderFieldNames(osName);
                            OGRFieldDefn oOgrField ( osName,
                                                    OFTString );
                            m_poOgrFeatureDefn->AddFieldDefn ( &oOgrField );
                        }
                    }
                }
            }

            iFeature = 0;

        }



        /***** check if any features are another layer *****/

        m_poOgrDS->ParseLayers ( m_poKmlLayer, poSpatialRef );

    }

    /***** it was from a DS::CreateLayer *****/

    else {

        /***** mark the layer as updated *****/

        bUpdated = TRUE;

        /***** create a new schema *****/

        KmlFactory *poKmlFactory = m_poOgrDS->GetKmlFactory (  );

        m_poKmlSchema = poKmlFactory->CreateSchema (  );

        /***** set the id on the new schema *****/

        std::string oKmlSchemaID = m_pszName;
        oKmlSchemaID.append ( ".schema" );
        m_poKmlSchema->set_id ( oKmlSchemaID );
    }




}

/******************************************************************************
 OGRLIBKMLLayer Destructor

 Args:          none
 
 Returns:       nothing
                
******************************************************************************/

OGRLIBKMLLayer::~OGRLIBKMLLayer (  )
{

    CPLFree ( ( void * )m_pszName );
    CPLFree ( ( void * )m_pszFileName );
    m_poOgrSRS->Release();

    m_poOgrFeatureDefn->Release (  );


}


/******************************************************************************
 Method to get the next feature on the layer

 Args:          none
 
 Returns:       The next feature, or NULL if there is no more

 this function copyed from the sqlite driver
******************************************************************************/

OGRFeature *OGRLIBKMLLayer::GetNextFeature()

{
    for( ; TRUE; )
    {
        OGRFeature      *poFeature;

        poFeature = GetNextRawFeature();
        if( poFeature == NULL )
            return NULL;

        if( (m_poFilterGeom == NULL
            || FilterGeometry( poFeature->GetGeometryRef() ) )
            && (m_poAttrQuery == NULL
                || m_poAttrQuery->Evaluate( poFeature )) )
            return poFeature;

        delete poFeature;
    }
}

/******************************************************************************
 Method to get the next feature on the layer

 Args:          none
 
 Returns:       The next feature, or NULL if there is no more
                
******************************************************************************/

OGRFeature *OGRLIBKMLLayer::GetNextRawFeature (
     )
{
    FeaturePtr poKmlFeature;
    OGRFeature *poOgrFeature = NULL;

    /***** loop over the kml features to find the next placemark *****/

    do {
        if ( iFeature >= nFeatures )
            break;

        /***** get the next kml feature in the container *****/
        
        poKmlFeature = m_poKmlLayer->get_feature_array_at ( iFeature++ );

        /***** what type of kml feature in the container? *****/

        switch (poKmlFeature->Type (  )) {

            case kmldom::Type_Placemark:
                poOgrFeature = kml2feat ( AsPlacemark ( poKmlFeature ),
                                          m_poOgrDS, this,
                                          m_poOgrFeatureDefn, m_poOgrSRS );
                break;    

            case kmldom::Type_GroundOverlay:
                if (m_bReadGroundOverlay) {
                    poOgrFeature =
                        kmlgroundoverlay2feat ( AsGroundOverlay ( poKmlFeature ),
                                                m_poOgrDS, this,
                                                m_poOgrFeatureDefn,
                                                m_poOgrSRS );
                }
                break;
                
            default:
                break;

        }

    } while ( !poOgrFeature );

    /***** set the FID on the ogr feature *****/
    
    if (poOgrFeature)
        poOgrFeature->SetFID(nFID ++);
    
    return poOgrFeature;
}

/******************************************************************************
 method to add a feature to a layer

 Args:          poOgrFeat   pointer to the feature to add
 
 Returns:       OGRERR_NONE, or OGRERR_UNSUPPORTED_OPERATION of the layer is
                not writeable
                
******************************************************************************/

OGRErr OGRLIBKMLLayer::CreateFeature (
    OGRFeature * poOgrFeat )
{

    if ( !bUpdate )
        return OGRERR_UNSUPPORTED_OPERATION;

    PlacemarkPtr poKmlPlacemark =
        feat2kml ( m_poOgrDS, this, poOgrFeat, m_poOgrDS->GetKmlFactory (  ) );

    m_poKmlLayer->add_feature ( poKmlPlacemark );

    /***** update the layer class count of features  *****/

    nFeatures++;
    
    /***** mark the layer as updated *****/

    bUpdated = TRUE;
    m_poOgrDS->Updated (  );

    return OGRERR_NONE;
}

/******************************************************************************
 method to get the number of features on the layer

 Args:          bForce      no effect as of now
 
 Returns:       the number of features on the layer

 Note:          the result can include links, folders and other items that are
                not supported by OGR
                
******************************************************************************/

int OGRLIBKMLLayer::GetFeatureCount (
                                     int bForce )
{


    int i = 0; 
    if (m_poFilterGeom != NULL || m_poAttrQuery != NULL ) {
        i = OGRLayer::GetFeatureCount( bForce );
    }

    else {
        size_t iKmlFeature; 
        size_t nKmlFeatures = m_poKmlLayer->get_feature_array_size (  );
        FeaturePtr poKmlFeature;

        /***** loop over the kml features in the container *****/

        for ( iKmlFeature = 0; iKmlFeature < nKmlFeatures; iKmlFeature++ ) {
            poKmlFeature = m_poKmlLayer->get_feature_array_at ( iKmlFeature );

            /***** what type of kml feature? *****/

            switch (poKmlFeature->Type (  )) {

                case kmldom::Type_Placemark:
                    i++;
                    break;

                case kmldom::Type_GroundOverlay:
                    if (m_bReadGroundOverlay)
                        i++;
                    break;

                default:
                    break;

            } 
        }
    }

    return i;
}

/******************************************************************************
 GetExtent()

 Args:          psExtent    pointer to the Envelope to store the result in
                bForce      no effect as of now 
 
 Returns:       nothing
                
******************************************************************************/

OGRErr OGRLIBKMLLayer::GetExtent (
    OGREnvelope * psExtent,
    int bForce )
{
    Bbox oKmlBbox;

    if ( kmlengine::
         GetFeatureBounds ( AsFeature ( m_poKmlLayer ), &oKmlBbox ) ) {
        psExtent->MinX = oKmlBbox.get_west (  );
        psExtent->MinY = oKmlBbox.get_south (  );
        psExtent->MaxX = oKmlBbox.get_east (  );
        psExtent->MaxY = oKmlBbox.get_north (  );

        return OGRERR_NONE;
    }
    else
        return OGRLayer::GetExtent(psExtent, bForce);
}




/******************************************************************************
 Method to create a field on a layer

 Args:          poField     pointer to the Field Definition to add
                bApproxOK   no effect as of now 

 Returns:       OGRERR_NONE on success or OGRERR_UNSUPPORTED_OPERATION if the
                layer is not writeable
                
******************************************************************************/

OGRErr OGRLIBKMLLayer::CreateField (
    OGRFieldDefn * poField,
    int bApproxOK )
{

    if ( !bUpdate )
        return OGRERR_UNSUPPORTED_OPERATION;

    SimpleFieldPtr poKmlSimpleField = NULL;

    if ( (poKmlSimpleField =
         FieldDef2kml ( poField, m_poOgrDS->GetKmlFactory (  ) )) != NULL )
        m_poKmlSchema->add_simplefield ( poKmlSimpleField );

    m_poOgrFeatureDefn->AddFieldDefn ( poField );

    /***** mark the layer as updated *****/

    bUpdated = TRUE;
    m_poOgrDS->Updated (  );

    return OGRERR_NONE;
}


/******************************************************************************
 method to write the datasource to disk

 Args:      none

 Returns    nothing
                
******************************************************************************/

OGRErr OGRLIBKMLLayer::SyncToDisk (
     )
{

    return OGRERR_NONE;
}

/******************************************************************************
 method to get a layers style table
 
 Args:          none
 
 Returns:       pointer to the layers style table, or NULL if it does
                not have one
                
******************************************************************************/

OGRStyleTable *OGRLIBKMLLayer::GetStyleTable (
     )
{

    return m_poStyleTable;
}

/******************************************************************************
 method to write a style table to a layer
 
 Args:          poStyleTable    pointer to the style table to add
 
 Returns:       nothing

 note: this method assumes ownership of the style table
******************************************************************************/

void OGRLIBKMLLayer::SetStyleTableDirectly (
    OGRStyleTable * poStyleTable )
{

    if ( !bUpdate )
        return;

    KmlFactory *poKmlFactory = m_poOgrDS->GetKmlFactory (  );

    if ( m_poStyleTable )
        delete m_poStyleTable;

    m_poStyleTable = poStyleTable;

    if ( m_poKmlLayer->IsA ( kmldom::Type_Document ) ) {

        /***** delete all the styles *****/

        DocumentPtr poKmlDocument = AsDocument ( m_poKmlLayer );
        size_t nKmlStyles = poKmlDocument->get_schema_array_size (  );
        int iKmlStyle;

        for ( iKmlStyle = nKmlStyles - 1; iKmlStyle >= 0; iKmlStyle-- ) {
            poKmlDocument->DeleteStyleSelectorAt ( iKmlStyle );
        }

        /***** add the new style table to the document *****/

        styletable2kml ( poStyleTable, poKmlFactory,
                         AsContainer ( poKmlDocument ) );

    }

    /***** mark the layer as updated *****/

    bUpdated = TRUE;
    m_poOgrDS->Updated (  );

    return;
}

/******************************************************************************
 method to write a style table to a layer
 
 Args:          poStyleTable    pointer to the style table to add
 
 Returns:       nothing

 note:  this method copys the style table, and the user will still be
        responsible for its destruction
******************************************************************************/

void OGRLIBKMLLayer::SetStyleTable (
    OGRStyleTable * poStyleTable )
{

    if ( !bUpdate )
        return;

    if ( poStyleTable )
        SetStyleTableDirectly ( poStyleTable->Clone (  ) );
    else
        SetStyleTableDirectly ( NULL );
    return;
}

/******************************************************************************
 Test if capability is available.

 Args:          pszCap  layer capability name to test
 
 Returns:       True if the layer supports the capability, otherwise false

******************************************************************************/

int OGRLIBKMLLayer::TestCapability (
    const char *pszCap )
{
    int result = FALSE;

    if ( EQUAL ( pszCap, OLCRandomRead ) )
        result = FALSE;
    else if ( EQUAL ( pszCap, OLCSequentialWrite ) )
        result = bUpdate;
    else if ( EQUAL ( pszCap, OLCRandomWrite ) )
        result = FALSE;
    else if ( EQUAL ( pszCap, OLCFastFeatureCount ) )
        result = FALSE;
    else if ( EQUAL ( pszCap, OLCFastSetNextByIndex ) )
        result = FALSE;
    else if ( EQUAL ( pszCap, OLCCreateField ) )
        result = bUpdate;
    else if ( EQUAL ( pszCap, OLCDeleteFeature ) )
        result = FALSE;
    else if ( EQUAL(pszCap, OLCStringsAsUTF8) )
        result = TRUE;

    return result;
}

/************************************************************************/
/*                        LaunderFieldNames()                           */
/************************************************************************/

CPLString OGRLIBKMLLayer::LaunderFieldNames(CPLString osName)
{
    CPLString osLaunderedName;
    for(int i=0;i<(int)osName.size();i++)
    {
        char ch = osName[i];
        if ((ch >= '0' && ch <= '9') ||
            (ch >= 'a' && ch <= 'z') ||
            (ch >= 'A' && ch <= 'Z') ||
            (ch == '_'))
            osLaunderedName += ch;
        else
            osLaunderedName += "_";
    }
    return osLaunderedName;
}
