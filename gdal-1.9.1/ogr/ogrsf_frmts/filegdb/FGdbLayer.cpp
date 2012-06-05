/******************************************************************************
* $Id: FGdbLayer.cpp 23834 2012-01-31 19:02:11Z rouault $
*
* Project:  OpenGIS Simple Features Reference Implementation
* Purpose:  Implements FileGDB OGR layer.
* Author:   Ragi Yaser Burhum, ragi@burhum.com
*           Paul Ramsey, pramsey at cleverelephant.ca
*
******************************************************************************
* Copyright (c) 2010, Ragi Yaser Burhum
* Copyright (c) 2011, Paul Ramsey <pramsey at cleverelephant.ca>
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

#include "ogr_fgdb.h"
#include "ogrpgeogeometry.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include "FGdbUtils.h"
#include "cpl_minixml.h" // the only way right now to extract schema information

CPL_CVSID("$Id: FGdbLayer.cpp 23834 2012-01-31 19:02:11Z rouault $");

using std::string;
using std::wstring;

/************************************************************************/
/*                              FGdbLayer()                               */
/************************************************************************/
FGdbLayer::FGdbLayer():
    OGRLayer(), m_pDS(NULL), m_pTable(NULL), m_pFeatureDefn(NULL), 
    m_pSRS(NULL), m_wstrSubfields(L"*"), m_pOGRFilterGeometry(NULL), 
    m_pEnumRows(NULL), m_bFilterDirty(true), 
    m_supressColumnMappingError(false), m_forceMulti(false), m_bLaunderReservedKeywords(true)
{
    m_pEnumRows = new EnumRows;
}

/************************************************************************/
/*                            ~FGdbLayer()                              */
/************************************************************************/

FGdbLayer::~FGdbLayer()
{
    if (m_pFeatureDefn)
    {
        m_pFeatureDefn->Release();
        m_pFeatureDefn = NULL;
    }

    if (m_pSRS)
    {
        m_pSRS->Release();
        m_pSRS = NULL;
    }

    if (m_pEnumRows)
    {
        delete m_pEnumRows;
        m_pEnumRows = NULL;
    }

    // NOTE: never delete m_pDS - the memory doesn't belong to us
    // TODO: check if we need to close the table or if the destructor 
    // takes care of closing as it should
    if (m_pTable)
    {
        delete m_pTable;
        m_pTable = NULL;
    }

    if (m_pOGRFilterGeometry)
    {
        OGRGeometryFactory::destroyGeometry(m_pOGRFilterGeometry);
        m_pOGRFilterGeometry = NULL;
    }
    
}


/************************************************************************/
/*                            CreateFeature()                           */
/* Create an FGDB Row and populate it from an OGRFeature.               */
/*                                                                      */
/************************************************************************/

OGRErr FGdbLayer::CreateFeature( OGRFeature *poFeature )
{
    Table *fgdb_table = m_pTable;
    Row fgdb_row;
    fgdbError hr;
    ShapeBuffer shape;

    hr = fgdb_table->CreateRowObject(fgdb_row);

    /* Check the status of the Row create */
    if (FAILED(hr))
    {
        GDBErr(hr, "Failed at creating Row in CreateFeature.");
        return OGRERR_FAILURE;
    }

    OGRFeatureDefn* poFeatureDefn = m_pFeatureDefn;
    int nFieldCount = poFeatureDefn->GetFieldCount();

    /* Copy the OGR visible fields (everything except geometry and FID) */
    for( int i = 0; i < nFieldCount; i++ )
    {
        std::string field_name = poFeatureDefn->GetFieldDefn(i)->GetNameRef();
        std::wstring wfield_name = StringToWString(field_name);

        /* Set empty fields to NULL */
        if( !poFeature->IsFieldSet( i ) )
        {
            if (FAILED(hr = fgdb_row.SetNull(wfield_name)))
            {
                GDBErr(hr, "Failed setting field to NULL.");
                return OGRERR_FAILURE;
            }
            continue;
        }

        /* Set the information using the appropriate FGDB function */
        int nOGRFieldType = poFeatureDefn->GetFieldDefn(i)->GetType();

        if ( nOGRFieldType == OFTInteger )
        {
            /* Integers (we don't do FGDB Shorts) */
            int fldvalue = poFeature->GetFieldAsInteger(i);
            hr = fgdb_row.SetInteger(wfield_name, fldvalue);
        }
        else if ( nOGRFieldType == OFTReal )
        {
            /* Doubles (we don't handle FGDB Floats) */
            double fldvalue = poFeature->GetFieldAsDouble(i);
            hr = fgdb_row.SetDouble(wfield_name, fldvalue);
        }
        else if ( nOGRFieldType == OFTString )
        {
            /* Strings we convert to wstring */
            std::string fldvalue = poFeature->GetFieldAsString(i);
            std::wstring wfldvalue = StringToWString(fldvalue);
            hr = fgdb_row.SetString(wfield_name, wfldvalue);
        }
        else if ( nOGRFieldType == OFTDateTime || nOGRFieldType == OFTDate )
        {
            /* Dates we need to coerce a little */
            struct tm val;
            poFeature->GetFieldAsDateTime(i, &(val.tm_year), &(val.tm_mon), &(val.tm_mday),
                                          &(val.tm_hour), &(val.tm_min), &(val.tm_sec), NULL);
            val.tm_year -= 1900;
            val.tm_mon = val.tm_mon - 1; /* OGR months go 1-12, FGDB go 0-11 */
            hr = fgdb_row.SetDate(wfield_name, val);
        }
        else if ( nOGRFieldType == OFTBinary )
        {
            /* Binary data */
            ByteArray fgdb_bytearray;
            int bytesize;
            GByte *bytes = poFeature->GetFieldAsBinary(i, &bytesize);
            if ( bytesize )
            {
                fgdb_bytearray.Allocate(bytesize);
                memcpy(fgdb_bytearray.byteArray, bytes, bytesize);
                fgdb_bytearray.inUseLength = bytesize;
                hr = fgdb_row.SetBinary(wfield_name, fgdb_bytearray);
            }
            else
            {
                hr = fgdb_row.SetNull(wfield_name);
            }
        }
        else
        {
            /* We can't handle this type */
            CPLError( CE_Failure, CPLE_AppDefined,
                "FGDB driver does not support OGR type." );
            return OGRERR_FAILURE;
        }
    }

    /* Done with attribute fields, now do geometry */
    OGRGeometry *poGeom = poFeature->GetGeometryRef();

    /* Write geometry to a buffer */
    GByte *pabyShape = NULL;
    int nShapeSize = 0;
    OGRErr err = OGRWriteToShapeBin( poGeom, &pabyShape, &nShapeSize );
    if ( err != OGRERR_NONE )
        return err;

    /* Copy it into a ShapeBuffer */
    if ( nShapeSize > 0 )
    {
        shape.Allocate(nShapeSize);
        memcpy(shape.shapeBuffer, pabyShape, nShapeSize);
        shape.inUseLength = nShapeSize;
    }

    /* Free the shape buffer */
    CPLFree(pabyShape);

    /* Write ShapeBuffer into the Row */
    hr = fgdb_row.SetGeometry(shape);
    if (FAILED(hr))
    {
        GDBErr(hr, "Failed at writing Geometry to Row in CreateFeature.");
        return OGRERR_FAILURE;
    }

    /* Cannot write to FID field - it is managed by GDB*/
    //std::wstring wfield_name = StringToWString(m_strOIDFieldName);
    //hr = fgdb_row.SetInteger(wfield_name, poFeature->GetFID());

    /* Write the row to the table */
    hr = fgdb_table->Insert(fgdb_row);
    if (FAILED(hr))
    {
        GDBErr(hr, "Failed at writing Row to Table in CreateFeature.");
        return OGRERR_FAILURE;
    }

    return OGRERR_NONE;

}

/************************************************************************/
/*                            CreateField()                             */
/*  Build up an FGDB XML field definition and use it to create a Field  */
/*  Update the OGRFeatureDefn to reflect the new field.                 */
/*                                                                      */
/************************************************************************/

OGRErr FGdbLayer::CreateField(OGRFieldDefn* poField, int bApproxOK)
{
    OGRFieldDefn oField(poField);
    std::string fieldname = oField.GetNameRef();
    std::string fidname = std::string(GetFIDColumn());
    std::string nullable = "true";
    Table *fgdb_table = m_pTable;

    /* Try to map the OGR type to an ESRI type */
    OGRFieldType fldtype = oField.GetType();
    std::string gdbFieldType;
    if ( ! OGRToGDBFieldType(fldtype, &gdbFieldType) )
    {
        GDBErr(-1, "Failed converting field type.");
        return OGRERR_FAILURE;
    }

    /* If we don't have our FGDB Table pointer intialized, we can quit now. */
    if ( ! m_pTable )
    {
        GDBErr(-1, "FGDB Table has not been initialized.");
        return OGRERR_FAILURE;
    }

    /* Clean field names */
    std::string fieldname_clean = FGDBLaunderName(fieldname);

    if (m_bLaunderReservedKeywords)
        fieldname_clean = FGDBEscapeReservedKeywords(fieldname_clean);

    /* Truncate to 64 characters */
    if (fieldname_clean.size() > 64)
        fieldname_clean.resize(64);

    std::string temp_fieldname = fieldname_clean;

    /* Ensures uniqueness of field name */
    int numRenames = 1;
    while ((m_pFeatureDefn->GetFieldIndex(temp_fieldname.c_str()) >= 0) && (numRenames < 10))
    {
        temp_fieldname = CPLSPrintf("%s_%d", fieldname_clean.substr(0, 62).c_str(), numRenames);
        numRenames ++;
    }
    while ((m_pFeatureDefn->GetFieldIndex(temp_fieldname.c_str()) >= 0) && (numRenames < 100))
    {
        temp_fieldname = CPLSPrintf("%s_%d", fieldname_clean.substr(0, 61).c_str(), numRenames);
        numRenames ++;
    }

    if (temp_fieldname != fieldname)
    {
        if( !bApproxOK || (m_pFeatureDefn->GetFieldIndex(temp_fieldname.c_str()) >= 0) )
        {
            CPLError( CE_Failure, CPLE_NotSupported,
                "Failed to add field named '%s'",
                fieldname.c_str() );
            return OGRERR_FAILURE;
        }
        CPLError(CE_Warning, CPLE_NotSupported,
            "Normalized/laundered field name: '%s' to '%s'",
            fieldname.c_str(), temp_fieldname.c_str());

        fieldname_clean = temp_fieldname;
        oField.SetName(fieldname_clean.c_str());
    }

    /* Then the Field definition */
    CPLXMLNode *defn_xml = CPLCreateXMLNode(NULL, CXT_Element, "esri:Field");

    /* Add the XML attributes to the Field node */
    FGDB_CPLAddXMLAttribute(defn_xml, "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    FGDB_CPLAddXMLAttribute(defn_xml, "xmlns:xs", "http://www.w3.org/2001/XMLSchema");
    FGDB_CPLAddXMLAttribute(defn_xml, "xmlns:esri", "http://www.esri.com/schemas/ArcGIS/10.1");
    FGDB_CPLAddXMLAttribute(defn_xml, "xsi:type", "esri:Field");

    /* Basic field information */
    CPLCreateXMLElementAndValue(defn_xml, "Name", fieldname_clean.c_str());
    CPLCreateXMLElementAndValue(defn_xml, "Type", gdbFieldType.c_str());
    CPLCreateXMLElementAndValue(defn_xml, "IsNullable", nullable.c_str());

    /* Get the Width and Precision if we know them */
    int width = oField.GetWidth();
    int precision = oField.GetPrecision();
    if ( width <= 0 )
        GDBFieldTypeToWidthPrecision(gdbFieldType, &width, &precision);

    /* Write out the Width and Precision */
    char buf[100];
    snprintf(buf, 100, "%d", width);
    CPLCreateXMLElementAndValue(defn_xml,"Length", buf);
    snprintf(buf, 100, "%d", precision);
    CPLCreateXMLElementAndValue(defn_xml,"Precision", buf);

    /* We know nothing about Scale, so zero it out */
    CPLCreateXMLElementAndValue(defn_xml,"Scale", "0");

    /*  Attempt to preserve the original fieldname */
    if (fieldname != fieldname_clean)
    {
        CPLCreateXMLElementAndValue(defn_xml, "AliasName", fieldname.c_str());
    }

    /* Default values are discouraged in OGR API docs */
    /* <DefaultValue xsi:type="xs:string">afternoon</DefaultValue> */

    /* Convert our XML tree into a string for FGDB */
    char *defn_str = CPLSerializeXMLTree(defn_xml);
    CPLDebug("FGDB", "CreateField() generated XML for FGDB\n%s", defn_str);

    /* Add the FGDB Field to the FGDB Table. */
    fgdbError hr = fgdb_table->AddField(defn_str);

    /* Free the XML */
    CPLFree(defn_str);
    CPLDestroyXMLNode(defn_xml);

    /* Check the status of the Field add */
    if (FAILED(hr))
    {
        GDBErr(hr, "Failed at creating Field for " + fieldname);
        return OGRERR_FAILURE;
    }

    /* Now add the OGRFieldDefn to the OGRFeatureDefn */
    m_pFeatureDefn->AddFieldDefn(&oField);

    m_vOGRFieldToESRIField.push_back(StringToWString(fieldname_clean));
    m_vOGRFieldToESRIFieldType.push_back( gdbFieldType );

    /* All done and happy */
    return OGRERR_NONE;

}

/************************************************************************/
/*                      XMLSpatialReference()                           */
/*  Build up an XML representation of an OGRSpatialReference.           */
/*  Used in layer creation.                                             */
/*                                                                      */
/************************************************************************/

CPLXMLNode* XMLSpatialReference(OGRSpatialReference* poSRS, char** papszOptions)
{
    /* We always need a SpatialReference */
    CPLXMLNode *srs_xml = CPLCreateXMLNode(NULL, CXT_Element, "SpatialReference");

    /* Extract the WKID before morphing */
    char *wkid = NULL;
    if ( poSRS && poSRS->GetAuthorityCode(NULL) )
    {
        wkid = CPLStrdup(poSRS->GetAuthorityCode(NULL));
    }

    /* NULL poSRS => UnknownCoordinateSystem */
    if ( ! poSRS )
    {
        FGDB_CPLAddXMLAttribute(srs_xml, "xsi:type", "esri:UnknownCoordinateSystem");
    }
    else
    {
        /* Make a clone so we can morph it without morphing the original */
        OGRSpatialReference* poSRSClone = poSRS->Clone();

        /* Flip the WKT to ESRI form, return UnknownCoordinateSystem if we can't */
        if ( poSRSClone->morphToESRI() != OGRERR_NONE )
        {
            delete poSRSClone;
            FGDB_CPLAddXMLAttribute(srs_xml, "xsi:type", "esri:UnknownCoordinateSystem");
            return srs_xml;
        }

        /* Set the SpatialReference type attribute correctly for GEOGCS/PROJCS */
        if ( poSRSClone->IsProjected() )
            FGDB_CPLAddXMLAttribute(srs_xml, "xsi:type", "esri:ProjectedCoordinateSystem");
        else
            FGDB_CPLAddXMLAttribute(srs_xml, "xsi:type", "esri:GeographicCoordinateSystem");

        /* Add the WKT to the XML */
        char *wkt = NULL;
        poSRSClone->exportToWkt(&wkt);
        if (wkt)
        {
            CPLCreateXMLElementAndValue(srs_xml,"WKT", wkt);
            OGRFree(wkt);
        }

        /* Dispose of our close */
        delete poSRSClone;
    }
    
    /* Handle Origin/Scale/Tolerance */
    const char* grid[7] = {
      "XOrigin", "YOrigin", "XYScale",
      "ZOrigin", "ZScale",
      "XYTolerance", "ZTolerance" };
    const char* gridvalues[7] = {
      "-2147483647", "-2147483647", "1000000000",
      "-2147483647", "1000000000",
      "0.0001", "0.0001" };

    /* Convert any layer creation options available, use defaults otherwise */
    for( int i = 0; i < 7; i++ )
    {
        if ( CSLFetchNameValue( papszOptions, grid[i] ) != NULL )
            gridvalues[i] = CSLFetchNameValue( papszOptions, grid[i] );

        CPLCreateXMLElementAndValue(srs_xml, grid[i], gridvalues[i]);
    }

    /* FGDB is always High Precision */
    CPLCreateXMLElementAndValue(srs_xml, "HighPrecision", "true");     

    /* Add the WKID to the XML */
    if ( wkid ) 
    {
        CPLCreateXMLElementAndValue(srs_xml, "WKID", wkid);
        CPLFree(wkid);
    }

    return srs_xml;
}

/************************************************************************/
/*                    CreateFeatureDataset()                            */
/************************************************************************/

bool FGdbLayer::CreateFeatureDataset(FGdbDataSource* pParentDataSource, 
                                     std::string feature_dataset_name,
                                     OGRSpatialReference* poSRS,
                                     char** papszOptions )
{
    /* XML node */
    CPLXMLNode *xml_xml = CPLCreateXMLNode(NULL, CXT_Element, "?xml");
    FGDB_CPLAddXMLAttribute(xml_xml, "version", "1.0");
    FGDB_CPLAddXMLAttribute(xml_xml, "encoding", "UTF-8");

    /* First build up a bare-bones feature definition */
    CPLXMLNode *defn_xml = CPLCreateXMLNode(NULL, CXT_Element, "esri:DataElement");
    CPLAddXMLSibling(xml_xml, defn_xml);

    /* Add the attributes to the DataElement */
    FGDB_CPLAddXMLAttribute(defn_xml, "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    FGDB_CPLAddXMLAttribute(defn_xml, "xmlns:xs", "http://www.w3.org/2001/XMLSchema");
    FGDB_CPLAddXMLAttribute(defn_xml, "xmlns:esri", "http://www.esri.com/schemas/ArcGIS/10.1");

    /* Need to set this to esri:DEFeatureDataset or esri:DETable */
    FGDB_CPLAddXMLAttribute(defn_xml, "xsi:type", "esri:DEFeatureDataset");

    /* Add in more children */
    std::string catalog_page = "\\" + feature_dataset_name;
    CPLCreateXMLElementAndValue(defn_xml,"CatalogPath", catalog_page.c_str());
    CPLCreateXMLElementAndValue(defn_xml,"Name", feature_dataset_name.c_str());
    CPLCreateXMLElementAndValue(defn_xml,"ChildrenExpanded", "false");
    CPLCreateXMLElementAndValue(defn_xml,"DatasetType", "esriDTFeatureDataset");
    CPLCreateXMLElementAndValue(defn_xml,"Versioned", "false");
    CPLCreateXMLElementAndValue(defn_xml,"CanVersion", "false");

    /* Add in empty extent */
    CPLXMLNode *extent_xml = CPLCreateXMLNode(NULL, CXT_Element, "Extent");
    FGDB_CPLAddXMLAttribute(extent_xml, "xsi:nil", "true");
    CPLAddXMLChild(defn_xml, extent_xml);

    /* Add the SRS */
    if( TRUE ) // TODO: conditional on existence of SRS
    {
        CPLXMLNode *srs_xml = XMLSpatialReference(poSRS, papszOptions);
        if ( srs_xml )
            CPLAddXMLChild(defn_xml, srs_xml);
    }

    /* Convert our XML tree into a string for FGDB */
    char *defn_str = CPLSerializeXMLTree(xml_xml);
    CPLDestroyXMLNode(xml_xml);

    /* TODO, tie this to debugging levels */
    CPLDebug("FGDB", "%s", defn_str);

    /* Create the FeatureDataset. */
    Geodatabase *gdb = pParentDataSource->GetGDB();
    fgdbError hr = gdb->CreateFeatureDataset(defn_str);

    /* Free the XML */
    CPLFree(defn_str);

    /* Check table create status */
    if (FAILED(hr))
    {
        return GDBErr(hr, "Failed at creating FeatureDataset " + feature_dataset_name);
    }

    return true;
}

/************************************************************************/
/*                            Create()                                  */
/* Build up an FGDB XML layer definition and use it to create a Table   */
/* or Feature Class to work from.                                       */
/*                                                                      */
/* Layer creation options:                                              */
/*   FEATURE_DATASET, nest layer inside a FeatureDataset folder         */
/*   GEOMETRY_NAME, user-selected name for the geometry column          */
/*   OID_NAME, user-selected name for the FID column                    */
/*   XORIGIN, YORIGIN, ZORIGIN, origin of the snapping grid             */
/*   XYSCALE, ZSCALE, inverse resolution of the snapping grid           */
/*   XYTOLERANCE, ZTOLERANCE, snapping tolerance for topology/networks  */
/*                                                                      */
/************************************************************************/

bool FGdbLayer::Create(FGdbDataSource* pParentDataSource, 
                       const char* pszLayerNameIn, 
                       OGRSpatialReference* poSRS, 
                       OGRwkbGeometryType eType, 
                       char** papszOptions)
{
    std::string parent_path = "";
    std::wstring wtable_path, wparent_path;
    std::string geometry_name = FGDB_GEOMETRY_NAME;
    std::string fid_name = FGDB_OID_NAME;
    std::string esri_type;
    bool has_z = false;

    /* Launder the Layer name */
    std::string layerName = pszLayerNameIn;

    layerName = FGDBLaunderName(pszLayerNameIn);
    layerName = FGDBEscapeReservedKeywords(layerName);
    layerName = FGDBEscapeUnsupportedPrefixes(layerName);

    if (layerName.size() > 160)
        layerName.resize(160);

    /* Ensures uniqueness of layer name */
    int numRenames = 1;
    while ((pParentDataSource->GetLayerByName(layerName.c_str()) != NULL) && (numRenames < 10))
    {
        layerName = CPLSPrintf("%s_%d", layerName.substr(0, 158).c_str(), numRenames);
        numRenames ++;
    }
    while ((pParentDataSource->GetLayerByName(layerName.c_str()) != NULL) && (numRenames < 100))
    {
        layerName = CPLSPrintf("%s_%d", layerName.substr(0, 157).c_str(), numRenames);
        numRenames ++;
    }

    if (layerName != pszLayerNameIn)
    {
        CPLError(CE_Warning, CPLE_NotSupported,
                "Normalized/laundered layer name: '%s' to '%s'",
                pszLayerNameIn, layerName.c_str());
    }

    std::string table_path = "\\" + std::string(layerName);

    /* Handle the FEATURE_DATASET case */
    if (  CSLFetchNameValue( papszOptions, "FEATURE_DATASET") != NULL )
    {
        std::string feature_dataset = CSLFetchNameValue( papszOptions, "FEATURE_DATASET");

        /* Check if FEATURE_DATASET exists. Otherwise create it */
        std::vector<wstring> featuredatasets;
        Geodatabase *gdb = pParentDataSource->GetGDB();
        int bFeatureDataSetExists = FALSE;
        fgdbError hr;
        if ( !FAILED(hr = gdb->GetChildDatasets(L"\\", L"Feature Dataset", featuredatasets)) )
        {
            std::wstring feature_dataset_with_slash = L"\\" + StringToWString(feature_dataset);
            for ( unsigned int i = 0; i < featuredatasets.size(); i++ )
            {
                if (featuredatasets[i] == feature_dataset_with_slash)
                    bFeatureDataSetExists = TRUE;
            }
        }

        if (!bFeatureDataSetExists)
        {
            bool rv = CreateFeatureDataset(pParentDataSource, feature_dataset, poSRS, papszOptions);
            if ( ! rv )
                return rv;
        }

        table_path = "\\" + feature_dataset + table_path;
        parent_path = "\\" + feature_dataset;
    }

    /* Convert table_path into wstring */
    wtable_path = StringToWString(table_path);
    wparent_path = StringToWString(parent_path);

    /* Over-ride the geometry name if necessary */
    if ( CSLFetchNameValue( papszOptions, "GEOMETRY_NAME") != NULL )
        geometry_name = CSLFetchNameValue( papszOptions, "GEOMETRY_NAME");

    /* Over-ride the OID name if necessary */
    if ( CSLFetchNameValue( papszOptions, "OID_NAME") != NULL )
        fid_name = CSLFetchNameValue( papszOptions, "OID_NAME");

    /* Figure out our geometry type */
    if ( eType != wkbNone )
    {
        if ( wkbFlatten(eType) == wkbUnknown )
        {
            return GDBErr(-1, "FGDB layers cannot be created with a wkbUnknown layer geometry type.");
        }
        if ( ! OGRGeometryToGDB(eType, &esri_type, &has_z) )
            return GDBErr(-1, "Unable to map OGR type to ESRI type");
    }

    m_bLaunderReservedKeywords = CSLFetchBoolean( papszOptions, "LAUNDER_RESERVED_KEYWORDS", TRUE) == TRUE;

    /* XML node */
    CPLXMLNode *xml_xml = CPLCreateXMLNode(NULL, CXT_Element, "?xml");
    FGDB_CPLAddXMLAttribute(xml_xml, "version", "1.0");
    FGDB_CPLAddXMLAttribute(xml_xml, "encoding", "UTF-8");

    /* First build up a bare-bones feature definition */
    CPLXMLNode *defn_xml = CPLCreateXMLNode(NULL, CXT_Element, "esri:DataElement");
    CPLAddXMLSibling(xml_xml, defn_xml);

    /* Add the attributes to the DataElement */
    FGDB_CPLAddXMLAttribute(defn_xml, "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    FGDB_CPLAddXMLAttribute(defn_xml, "xmlns:xs", "http://www.w3.org/2001/XMLSchema");
    FGDB_CPLAddXMLAttribute(defn_xml, "xmlns:esri", "http://www.esri.com/schemas/ArcGIS/10.1");

    /* Need to set this to esri:DEFeatureDataset or esri:DETable */
    FGDB_CPLAddXMLAttribute(defn_xml, "xsi:type", (eType == wkbNone ? "esri:DETable" : "esri:DEFeatureClass"));

    /* Add in more children */
    CPLCreateXMLElementAndValue(defn_xml,"CatalogPath",table_path.c_str());
    CPLCreateXMLElementAndValue(defn_xml,"Name", layerName.c_str());
    CPLCreateXMLElementAndValue(defn_xml,"ChildrenExpanded", "false");

    /* WKB type of none implies this is a 'Table' otherwise it's a 'Feature Class' */
    std::string datasettype = (eType == wkbNone ? "esriDTTable" : "esriDTFeatureClass");
    CPLCreateXMLElementAndValue(defn_xml,"DatasetType", datasettype.c_str() );
    CPLCreateXMLElementAndValue(defn_xml,"Versioned", "false");
    CPLCreateXMLElementAndValue(defn_xml,"CanVersion", "false");

    /* We might need to make OID optional later, but OGR likes to have a FID */
    CPLCreateXMLElementAndValue(defn_xml,"HasOID", "true");
    CPLCreateXMLElementAndValue(defn_xml,"OIDFieldName", fid_name.c_str());

    /* Add in empty Fields */
    CPLXMLNode *fields_xml = CPLCreateXMLNode(defn_xml, CXT_Element, "Fields");
    FGDB_CPLAddXMLAttribute(fields_xml, "xsi:type", "esri:Fields");
    CPLXMLNode *fieldarray_xml = CPLCreateXMLNode(fields_xml, CXT_Element, "FieldArray");
    FGDB_CPLAddXMLAttribute(fieldarray_xml, "xsi:type", "esri:ArrayOfField");

    /* Feature Classes have an implicit geometry column, so we'll add it at creation time */
    if ( eType != wkbNone )
    {
        CPLXMLNode *shape_xml = CPLCreateXMLNode(fieldarray_xml, CXT_Element, "Field");
        FGDB_CPLAddXMLAttribute(shape_xml, "xsi:type", "esri:Field");
        CPLCreateXMLElementAndValue(shape_xml, "Name", geometry_name.c_str());
        CPLCreateXMLElementAndValue(shape_xml, "Type", "esriFieldTypeGeometry");
        CPLCreateXMLElementAndValue(shape_xml, "IsNullable", "false");
        CPLCreateXMLElementAndValue(shape_xml, "Length", "0");
        CPLCreateXMLElementAndValue(shape_xml, "Precision", "0");
        CPLCreateXMLElementAndValue(shape_xml, "Scale", "0");
        CPLCreateXMLElementAndValue(shape_xml, "Required", "true");
        CPLXMLNode *geom_xml = CPLCreateXMLNode(shape_xml, CXT_Element, "GeometryDef");
        FGDB_CPLAddXMLAttribute(geom_xml, "xsi:type", "esri:GeometryDef");
        CPLCreateXMLElementAndValue(geom_xml, "AvgNumPoints", "0");
        CPLCreateXMLElementAndValue(geom_xml, "GeometryType", esri_type.c_str());
        CPLCreateXMLElementAndValue(geom_xml,"HasM", "false");
        CPLCreateXMLElementAndValue(geom_xml,"HasZ", (has_z ? "true" : "false"));

        /* Add the SRS if we have one */
        CPLXMLNode *srs_xml = XMLSpatialReference(poSRS, papszOptions);
        if ( srs_xml )
            CPLAddXMLChild(geom_xml, srs_xml);
    }

    /* All (?) Tables and Feature Classes will have an ObjectID */
    CPLXMLNode *oid_xml = CPLCreateXMLNode(fieldarray_xml, CXT_Element, "Field");
    FGDB_CPLAddXMLAttribute(oid_xml, "xsi:type", "esri:Field");
    CPLCreateXMLElementAndValue(oid_xml, "Name", fid_name.c_str());
    CPLCreateXMLElementAndValue(oid_xml, "Type", "esriFieldTypeOID");
    CPLCreateXMLElementAndValue(oid_xml, "IsNullable", "false");
    CPLCreateXMLElementAndValue(oid_xml, "Length", "12");
    CPLCreateXMLElementAndValue(oid_xml, "Precision", "0");
    CPLCreateXMLElementAndValue(oid_xml, "Scale", "0");
    CPLCreateXMLElementAndValue(oid_xml, "Required", "true");

    /* Add in empty Indexes */
    CPLXMLNode *indexes_xml = CPLCreateXMLNode(defn_xml, CXT_Element, "Indexes");
    FGDB_CPLAddXMLAttribute(indexes_xml, "xsi:type", "esri:Indexes");
    CPLXMLNode *indexarray_xml = CPLCreateXMLNode(indexes_xml, CXT_Element, "IndexArray");
    FGDB_CPLAddXMLAttribute(indexarray_xml, "xsi:type", "esri:ArrayOfIndex");

    /* CLSID http://forums.arcgis.com/threads/34536?p=118484#post118484 */
    if ( eType == wkbNone )
    {
        CPLCreateXMLElementAndValue(defn_xml, "CLSID", "{7A566981-C114-11D2-8A28-006097AFF44E}");
        CPLCreateXMLElementAndValue(defn_xml, "EXTCLSID", "");
    }
    else
    {
        CPLCreateXMLElementAndValue(defn_xml, "CLSID", "{52353152-891A-11D0-BEC6-00805F7C4268}");
        CPLCreateXMLElementAndValue(defn_xml, "EXTCLSID", "");
    }

    /* Set the alias for the Feature Class */
    if (pszLayerNameIn != layerName)
    {
        CPLCreateXMLElementAndValue(defn_xml, "AliasName", pszLayerNameIn);
    }

    /* Map from OGR WKB type to ESRI type */
    if ( eType != wkbNone )
    {
        /* Declare our feature type */
        CPLCreateXMLElementAndValue(defn_xml,"FeatureType", "esriFTSimple");
        CPLCreateXMLElementAndValue(defn_xml,"ShapeType", esri_type.c_str());
        CPLCreateXMLElementAndValue(defn_xml,"ShapeFieldName", geometry_name.c_str());

        /* Dimensionality */
        CPLCreateXMLElementAndValue(defn_xml,"HasM", "false");
        CPLCreateXMLElementAndValue(defn_xml,"HasZ", (has_z ? "true" : "false"));

        /* TODO: Handle spatial indexes (layer creation option?) */
        CPLCreateXMLElementAndValue(defn_xml,"HasSpatialIndex", "false");

        /* These field are required for Arcmap to display aliases correctly */
        CPLCreateXMLNode(defn_xml, CXT_Element, "AreaFieldName");
        CPLCreateXMLNode(defn_xml, CXT_Element, "LengthFieldName");

        /* We can't know the extent at this point <Extent xsi:nil='true'/> */
        CPLXMLNode *extn_xml = CPLCreateXMLNode(defn_xml, CXT_Element, "Extent");
        FGDB_CPLAddXMLAttribute(extn_xml, "xsi:nil", "true");
    }

    /* Feature Class with known SRS gets an SRS entry */
    if( eType != wkbNone )
    {
        CPLXMLNode *srs_xml = XMLSpatialReference(poSRS, papszOptions);
        if ( srs_xml )
            CPLAddXMLChild(defn_xml, srs_xml);
    }

    /* Convert our XML tree into a string for FGDB */
    char *defn_str = CPLSerializeXMLTree(xml_xml);
    CPLDestroyXMLNode(xml_xml);

    /* TODO, tie this to debugging levels */
    CPLDebug("FGDB", "%s", defn_str);
    //std::cout << defn_str << std::endl;

    /* Create the table. */
    Table *table = new Table;
    Geodatabase *gdb = pParentDataSource->GetGDB();
    fgdbError hr = gdb->CreateTable(defn_str, wparent_path, *table);

    /* Free the XML */
    CPLFree(defn_str);

    /* Check table create status */
    if (FAILED(hr))
    {
        delete table;
        return GDBErr(hr, "Failed at creating table for " + table_path);
    }

    /* Store the new FGDB Table pointer and set up the OGRFeatureDefn */
    return FGdbLayer::Initialize(pParentDataSource, table, wtable_path, L"Table");
}

/*************************************************************************/
/*                            Initialize()                               */
/* Has ownership of the table as soon as it is called.                   */
/************************************************************************/

bool FGdbLayer::Initialize(FGdbDataSource* pParentDataSource, Table* pTable,
                           std::wstring wstrTablePath, std::wstring wstrType)
{
    long hr;

    m_pDS = pParentDataSource; // we never assume ownership of the parent - so our destructor should not delete

    m_pTable = pTable;

    m_wstrTablePath = wstrTablePath;
    m_wstrType = wstrType;

    wstring wstrQueryName;
    if (FAILED(hr = pParentDataSource->GetGDB()->GetQueryName(wstrTablePath, wstrQueryName)))
        return GDBErr(hr, "Failed at getting underlying table name for " +
                      WStringToString(wstrTablePath));

    m_strName = WStringToString(wstrQueryName);

    m_pFeatureDefn = new OGRFeatureDefn(m_strName.c_str()); //TODO: Should I "new" an OGR smart pointer - sample says so, but it doesn't seem right
    //as long as we use the same compiler & settings in both the ogr build and this
    //driver, we should be OK
    m_pFeatureDefn->Reference();

    string tableDef;
    if (FAILED(hr = m_pTable->GetDefinition(tableDef)))
        return GDBErr(hr, "Failed at getting table definition for " +
                      WStringToString(wstrTablePath));

    //xxx  printf("Table definition = %s", tableDef.c_str() );

    bool abort = false;

    // extract schema information from table
    CPLXMLNode *psRoot = CPLParseXMLString( tableDef.c_str() );

    if (psRoot == NULL)
    {
        CPLError( CE_Failure, CPLE_AppDefined, "%s",
                  ("Failed parsing GDB Table Schema XML for " + m_strName).c_str());
        return false;
    }

    CPLXMLNode *pDataElementNode = psRoot->psNext; // Move to next field which should be DataElement

    if( pDataElementNode != NULL
        && pDataElementNode->psChild != NULL
        && pDataElementNode->eType == CXT_Element
        && EQUAL(pDataElementNode->pszValue,"esri:DataElement") )
    {
        CPLXMLNode *psNode;

        for( psNode = pDataElementNode->psChild;
        psNode != NULL;
        psNode = psNode->psNext )
        {
            if( psNode->eType == CXT_Element && psNode->psChild != NULL )
            {
                if (EQUAL(psNode->pszValue,"OIDFieldName") )
                {
                    char* pszUnescaped = CPLUnescapeString(
                    psNode->psChild->pszValue, NULL, CPLES_XML);
                    m_strOIDFieldName = pszUnescaped;
                    CPLFree(pszUnescaped);
                }
                else if (EQUAL(psNode->pszValue,"ShapeFieldName") )
                {
                    char* pszUnescaped = CPLUnescapeString(
                    psNode->psChild->pszValue, NULL, CPLES_XML);
                    m_strShapeFieldName = pszUnescaped;
                    CPLFree(pszUnescaped);
                }
                else if (EQUAL(psNode->pszValue,"Fields") )
                {
                    if (!GDBToOGRFields(psNode))
                    {
                        abort = true;
                        break;
                    }
                }
            }
        }

        if (m_strShapeFieldName.size() == 0)
            m_pFeatureDefn->SetGeomType(wkbNone);
    }
    else
    {
        CPLError( CE_Failure, CPLE_AppDefined, "%s",
                ("Failed parsing GDB Table Schema XML (DataElement) for " + m_strName).c_str());
        return false;
    }
    CPLDestroyXMLNode( psRoot );

    if (abort)
        return false;

    return true; //AOToOGRFields(ipFields, m_pFeatureDefn, m_vOGRFieldToESRIField);
}

/************************************************************************/
/*                          ParseGeometryDef()                          */
/************************************************************************/

bool FGdbLayer::ParseGeometryDef(CPLXMLNode* psRoot)
{
    CPLXMLNode *psGeometryDefItem;

    string geometryType;
    bool hasZ = false;
    string wkt, wkid;

    for (psGeometryDefItem = psRoot->psChild;
        psGeometryDefItem != NULL;
        psGeometryDefItem = psGeometryDefItem->psNext )
    {
        //loop through all "GeometryDef" elements
        //

        if (psGeometryDefItem->eType == CXT_Element &&
            psGeometryDefItem->psChild != NULL)
        {
            if (EQUAL(psGeometryDefItem->pszValue,"GeometryType"))
            {
                char* pszUnescaped = CPLUnescapeString(
                    psGeometryDefItem->psChild->pszValue, NULL, CPLES_XML);

                geometryType = pszUnescaped;

                CPLFree(pszUnescaped);
            }
            else if (EQUAL(psGeometryDefItem->pszValue,"SpatialReference"))
            {
                ParseSpatialReference(psGeometryDefItem, &wkt, &wkid); // we don't check for success because it
                                                                // may not be there
            }
            /* No M support in OGR yet
            else if (EQUAL(psFieldNode->pszValue,"HasM")
            {
                char* pszUnescaped = CPLUnescapeString(psNode->psChild->pszValue, NULL, CPLES_XML);

                if (!strcmp(szUnescaped, "true"))
                hasM = true;

                CPLFree(pszUnescaped);
            }
            */
            else if (EQUAL(psGeometryDefItem->pszValue,"HasZ"))
            {
                char* pszUnescaped = CPLUnescapeString(
                    psGeometryDefItem->psChild->pszValue, NULL, CPLES_XML);

                if (!strcmp(pszUnescaped, "true"))
                hasZ = true;

                CPLFree(pszUnescaped);
            }
        }

    }

    OGRwkbGeometryType ogrGeoType;
    if (!GDBToOGRGeometry(geometryType, hasZ, &ogrGeoType))
        return false;

    m_pFeatureDefn->SetGeomType(ogrGeoType);

    if (wkbFlatten(ogrGeoType) == wkbMultiLineString ||
        wkbFlatten(ogrGeoType) == wkbMultiPoint)
        m_forceMulti = true;

    if (wkid.length() > 0)
    {
        m_pSRS = new OGRSpatialReference();
        if (m_pSRS->importFromEPSG(atoi(wkid.c_str())) != OGRERR_NONE)
        {
            delete m_pSRS;
            m_pSRS = NULL;
        }
        else
            return true;
    }

    if (wkt.length() > 0)
    {
        if (!GDBToOGRSpatialReference(wkt, &m_pSRS))
        {
            //report error, but be passive about it
            CPLError( CE_Warning, CPLE_AppDefined,
                      "Failed Mapping ESRI Spatial Reference");
        }
    }
    else
    {
        //report error, but be passive about it
        CPLError( CE_Warning, CPLE_AppDefined, "Empty Spatial Reference");
    }

    return true;
}

/************************************************************************/
/*                        ParseSpatialReference()                       */
/************************************************************************/

bool FGdbLayer::ParseSpatialReference(CPLXMLNode* psSpatialRefNode,
                                      string* pOutWkt, string* pOutWKID)
{
    *pOutWkt = "";
    *pOutWKID = "";

    CPLXMLNode* psSRItemNode;

    /* Loop through all the SRS elements we want to store */
    for( psSRItemNode = psSpatialRefNode->psChild;
         psSRItemNode != NULL;
         psSRItemNode = psSRItemNode->psNext )
    {
        /* The WKID maps (mostly) to an EPSG code */
        if( psSRItemNode->eType == CXT_Element &&
            psSRItemNode->psChild != NULL &&
            EQUAL(psSRItemNode->pszValue,"WKID") )
        {
            char* pszUnescaped = CPLUnescapeString(psSRItemNode->psChild->pszValue, NULL, CPLES_XML);
            *pOutWKID = pszUnescaped;
            CPLFree(pszUnescaped);
        }
        /* The WKT well-known text can be converted by OGR */
        else if( psSRItemNode->eType == CXT_Element &&
                psSRItemNode->psChild != NULL &&
                EQUAL(psSRItemNode->pszValue,"WKT") )
        {
            char* pszUnescaped = CPLUnescapeString(psSRItemNode->psChild->pszValue, NULL, CPLES_XML);
            *pOutWkt = pszUnescaped;
            CPLFree(pszUnescaped);
        }

    }
    return (*pOutWkt != "" || *pOutWKID != "");
}

/************************************************************************/
/*                          GDBToOGRFields()                           */
/************************************************************************/

bool FGdbLayer::GDBToOGRFields(CPLXMLNode* psRoot)
{
    m_vOGRFieldToESRIField.clear();

    if (psRoot->psChild == NULL || psRoot->psChild->psNext == NULL)
    {
        CPLError( CE_Failure, CPLE_AppDefined, "Unrecognized GDB XML Schema");

        return false;
    }

    psRoot = psRoot->psChild->psNext; //change root to "FieldArray"

    //CPLAssert(ogrToESRIFieldMapping.size() == pOGRFeatureDef->GetFieldCount());

    CPLXMLNode* psFieldNode;

    for( psFieldNode = psRoot->psChild;
        psFieldNode != NULL;
        psFieldNode = psFieldNode->psNext )
    {
        //loop through all "Field" elements
        //

        if( psFieldNode->eType == CXT_Element && psFieldNode->psChild != NULL &&
            EQUAL(psFieldNode->pszValue,"Field"))
        {

            CPLXMLNode* psFieldItemNode;
            std::string fieldName;
            std::string fieldType;
            int nLength = 0;
            int nPrecision = 0;

            // loop through all items in Field element
            //

            for( psFieldItemNode = psFieldNode->psChild;
                psFieldItemNode != NULL;
                psFieldItemNode = psFieldItemNode->psNext )
            {
                if (psFieldItemNode->eType == CXT_Element)
                {

                if (EQUAL(psFieldItemNode->pszValue,"Name"))
                {
                    char* pszUnescaped = CPLUnescapeString(
                        psFieldItemNode->psChild->pszValue, NULL, CPLES_XML);
                    fieldName = pszUnescaped;
                    CPLFree(pszUnescaped);
                }
                else if (EQUAL(psFieldItemNode->pszValue,"Type") )
                {
                    char* pszUnescaped = CPLUnescapeString(
                        psFieldItemNode->psChild->pszValue, NULL, CPLES_XML);
                    fieldType = pszUnescaped;
                    CPLFree(pszUnescaped);
                }
                else if (EQUAL(psFieldItemNode->pszValue,"GeometryDef") )
                {
                    if (!ParseGeometryDef(psFieldItemNode))
                        return false; // if we failed parsing the GeometryDef, we are done!
                }
                else if (EQUAL(psFieldItemNode->pszValue,"Length") )
                {
                    nLength = atoi(psFieldItemNode->psChild->pszValue);
                }
                else if (EQUAL(psFieldItemNode->pszValue,"Precision") )
                {
                    nPrecision = atoi(psFieldItemNode->psChild->pszValue);
                }
                }
            }


            ///////////////////////////////////////////////////////////////////
            // At this point we have parsed everything about the current field


            if (fieldType == "esriFieldTypeGeometry")
            {
                m_strShapeFieldName = fieldName;

                continue; // finish here for special field - don't add as OGR fielddef
            }
            else if (fieldType == "esriFieldTypeOID")
            {
                //m_strOIDFieldName = fieldName; // already set by this point

                continue; // finish here for special field - don't add as OGR fielddef
            }

            OGRFieldType ogrType;
            //CPLDebug("FGDB", "name = %s, type = %s", fieldName.c_str(), fieldType.c_str() );
            if (!GDBToOGRFieldType(fieldType, &ogrType))
            {
                // field cannot be mapped, skipping further processing
                CPLError( CE_Warning, CPLE_AppDefined, "Skipping field: [%s] type: [%s] ",
                fieldName.c_str(), fieldType.c_str() );
                continue;
            }


            //TODO: Optimization - modify m_wstrSubFields so it only fetches fields that are mapped

            OGRFieldDefn fieldTemplate( fieldName.c_str(), ogrType);
            //fieldTemplate.SetWidth(nLength);
            //fieldTemplate.SetPrecision(nPrecision);
            m_pFeatureDefn->AddFieldDefn( &fieldTemplate );

            m_vOGRFieldToESRIField.push_back(StringToWString(fieldName));
            m_vOGRFieldToESRIFieldType.push_back( fieldType );

        }
    }

    return true;
}


/************************************************************************/
/*                            ResetReading()                            */
/************************************************************************/

void FGdbLayer::ResetReading()
{
    long hr;

    if (m_pOGRFilterGeometry && !m_pOGRFilterGeometry->IsEmpty())
    {
        // Search spatial
        // As of beta1, FileGDB only supports bbox searched, if we have GEOS installed,
        // we can do the rest ourselves.

        OGREnvelope ogrEnv;

        m_pOGRFilterGeometry->getEnvelope(&ogrEnv);

        //spatial query
        FileGDBAPI::Envelope env(ogrEnv.MinX, ogrEnv.MaxX, ogrEnv.MinY, ogrEnv.MaxY);

        if FAILED(hr = m_pTable->Search(m_wstrSubfields, m_wstrWhereClause, env, true, *m_pEnumRows))
        GDBErr(hr, "Failed Searching");

        m_bFilterDirty = false;

        return;
    }

    // Search non-spatial

    if FAILED(hr = m_pTable->Search(m_wstrSubfields, m_wstrWhereClause, true, *m_pEnumRows))
        GDBErr(hr, "Failed Searching");

    m_bFilterDirty = false;
  
}

/************************************************************************/
/*                         SetSpatialFilter()                           */
/************************************************************************/

void FGdbLayer::SetSpatialFilter( OGRGeometry* pOGRGeom )
{
    if (m_pOGRFilterGeometry)
    {
        OGRGeometryFactory::destroyGeometry(m_pOGRFilterGeometry);
        m_pOGRFilterGeometry = NULL;
    }

    if (pOGRGeom == NULL || pOGRGeom->IsEmpty())
    {
        m_bFilterDirty = true;

        return;
    }

    m_pOGRFilterGeometry = pOGRGeom->clone();

    m_pOGRFilterGeometry->transformTo(m_pSRS);

    m_bFilterDirty = true;
}

/************************************************************************/
/*                         SetSpatialFilterRect()                       */
/************************************************************************/

void FGdbLayer::SetSpatialFilterRect (double dfMinX, double dfMinY, double dfMaxX, double dfMaxY)
{

    //TODO: can optimize this by changing how the filter gets generated -
    //this will work for now

    OGRGeometry* pTemp = OGRGeometryFactory::createGeometry(wkbPolygon);

    pTemp->assignSpatialReference(m_pSRS);

    OGRLinearRing ring;

    ring.addPoint( dfMinX, dfMinY );
    ring.addPoint( dfMinX, dfMaxY );
    ring.addPoint( dfMaxX, dfMaxY );
    ring.addPoint( dfMaxX, dfMinY );
    ring.addPoint( dfMinX, dfMinY );
    ((OGRPolygon *) pTemp)->addRing( &ring );

    SetSpatialFilter(pTemp);

    OGRGeometryFactory::destroyGeometry(pTemp);
}


/************************************************************************/
/*                         SetAttributeFilter()                         */
/************************************************************************/

OGRErr FGdbLayer::SetAttributeFilter( const char* pszQuery )
{
    m_wstrWhereClause = StringToWString( (pszQuery != NULL) ? pszQuery : "" );

    m_bFilterDirty = true;

    return OGRERR_NONE;
}

/************************************************************************/
/*                           OGRFeatureFromGdbRow()                      */
/************************************************************************/

bool FGdbLayer::OGRFeatureFromGdbRow(Row* pRow, OGRFeature** ppFeature)
{
    long hr;

    OGRFeature* pOutFeature = new OGRFeature(m_pFeatureDefn);

    /////////////////////////////////////////////////////////
    // Translate OID
    //

    int32 oid = -1;
    if (FAILED(hr = pRow->GetOID(oid)))
    {
        //this should never happen
        delete pOutFeature;
        return false;
    }
    pOutFeature->SetFID(oid);


    /////////////////////////////////////////////////////////
    // Translate Geometry
    //

    ShapeBuffer gdbGeometry;
    if (!FAILED(hr = pRow->GetGeometry(gdbGeometry)))
    {
        OGRGeometry* pOGRGeo = NULL;

        if ((!GDBGeometryToOGRGeometry(m_forceMulti, &gdbGeometry, m_pSRS, &pOGRGeo)) || pOGRGeo == NULL)
        {
            delete pOutFeature;
            return GDBErr(hr, "Failed to translate FileGDB Geometry to OGR Geometry for row " + string(CPLSPrintf("%d", (int)oid)));
        }

        pOutFeature->SetGeometryDirectly(pOGRGeo);
    }


    //////////////////////////////////////////////////////////
    // Map fields
    //


    size_t mappedFieldCount = m_vOGRFieldToESRIField.size();

    bool foundBadColumn = false;

    for (size_t i = 0; i < mappedFieldCount; ++i)
    {
        const wstring & wstrFieldName = m_vOGRFieldToESRIField[i];
        const std::string & strFieldType = m_vOGRFieldToESRIFieldType[i];

        bool isNull = false;

        if (FAILED(hr = pRow->IsNull(wstrFieldName, isNull)))
        {
            GDBErr(hr, "Failed to determine NULL status from column " +
                   WStringToString(wstrFieldName));
            foundBadColumn = true;
            continue;
        }

        if (isNull)
        {
            continue; //leave as unset
        }

        //
        // NOTE: This switch statement needs to be kept in sync with GDBToOGRFieldType utility function
        //       since we are only checking for types we mapped in that utility function

        switch (m_pFeatureDefn->GetFieldDefn(i)->GetType())
        {

            case OFTInteger:
            {
                int32 val;

                if (FAILED(hr = pRow->GetInteger(wstrFieldName, val)))
                {
                    int16 shortval;
                    if (FAILED(hr = pRow->GetShort(wstrFieldName, shortval)))
                    {
                        GDBErr(hr, "Failed to determine integer value for column " +
                               WStringToString(wstrFieldName));
                        foundBadColumn = true;
                        continue;
                    }
                    val = shortval;
                }

                pOutFeature->SetField(i, (int)val);
            }
            break;

            case OFTReal:
            {
                if (strFieldType == "esriFieldTypeSingle")
                {
                    float val;

                    if (FAILED(hr = pRow->GetFloat(wstrFieldName, val)))
                    {
                        GDBErr(hr, "Failed to determine float value for column " +
                               WStringToString(wstrFieldName));
                        foundBadColumn = true;
                        continue;
                    }

                    pOutFeature->SetField(i, val);
                }
                else
                {
                    double val;

                    if (FAILED(hr = pRow->GetDouble(wstrFieldName, val)))
                    {
                        GDBErr(hr, "Failed to determine real value for column " +
                               WStringToString(wstrFieldName));
                        foundBadColumn = true;
                        continue;
                    }

                    pOutFeature->SetField(i, val);
                }
            }
            break;
            case OFTString:
            {
                wstring val;

                if (FAILED(hr = pRow->GetString(wstrFieldName, val)))
                {
                    GDBErr(hr, "Failed to determine string value for column " +
                        WStringToString(wstrFieldName));
                    foundBadColumn = true;
                    continue;
                }

                pOutFeature->SetField(i, WStringToString(val).c_str());
            }
            break;

            /* TODO: Need to get test dataset to implement these leave it as NULL for now
            case OFTBinary:
            {
                ByteArray binaryBuf;

                if (FAILED(hr = pRow->GetBinary(wstrFieldName, binaryBuf)))
                {
                GDBErr(hr, "Failed to determine binary value for column " + WStringToString(wstrFieldName));
                foundBadColumn = true;
                continue;
                }

                pOutFeature->SetField(i, (int)binaryBuf.inUseLength, (GByte*)binaryBuf.byteArray);
            }
            break;
            */

            case OFTDateTime:
            {
                struct tm val;

                if (FAILED(hr = pRow->GetDate(wstrFieldName, val)))
                {
                    GDBErr(hr, "Failed to determine date value for column " +
                           WStringToString(wstrFieldName));
                    foundBadColumn = true;
                    continue;
                }

                pOutFeature->SetField(i, val.tm_year + 1900, val.tm_mon + 1,
                                      val.tm_mday, val.tm_hour, val.tm_min, val.tm_sec);
            // Examine test data to figure out how to extract that
            }
            break;

            default:
            {
                if (!m_supressColumnMappingError)
                {
                    foundBadColumn = true;
                    CPLError( CE_Warning, CPLE_AppDefined,
                            "Row id: %d col:%d has unhandled col type (%d). Setting to NULL.",
                            (int)oid, (int)i, m_pFeatureDefn->GetFieldDefn(i)->GetType());
                }
            }
        }
    }

    if (foundBadColumn)
        m_supressColumnMappingError = true;


    *ppFeature = pOutFeature;

    return true;
}


/************************************************************************/
/*                           GetNextFeature()                           */
/************************************************************************/

OGRFeature* FGdbLayer::GetNextFeature()
{
    if (m_bFilterDirty)
        ResetReading();


    while (1) //want to skip errors
    {
        if (m_pEnumRows == NULL)
            return NULL;

        long hr;

        Row row;

        if (FAILED(hr = m_pEnumRows->Next(row)))
        {
            GDBErr(hr, "Failed fetching features");
            return NULL;
        }

        if (hr != S_OK)
        {
        // It's OK, we are done fetching - failure is catched by FAILED macro
            return NULL;
        }

        OGRFeature* pOGRFeature = NULL;

        if (!OGRFeatureFromGdbRow(&row,  &pOGRFeature))
        {
            int32 oid = -1;
            row.GetOID(oid);

            GDBErr(hr, CPLSPrintf("Failed translating FGDB row [%d] to OGR Feature", oid));

            //return NULL;
            continue; //skip feature
        }

        return pOGRFeature;
    }
}

/************************************************************************/
/*                             GetFeature()                             */
/************************************************************************/

OGRFeature *FGdbLayer::GetFeature( long oid )
{
    // do query to fetch individual row

    long           hr;
    Row            row;
    EnumRows       enumRows;
    CPLString      osQuery;

    osQuery.Printf("%s = %ld", m_strOIDFieldName.c_str(), oid);

    if (FAILED(hr = m_pTable->Search(m_wstrSubfields, StringToWString(osQuery.c_str()), true, enumRows)))
    {
        GDBErr(hr, "Failed fetching row ");
        return NULL;
    }

    if (FAILED(hr = enumRows.Next(row)))
    {
        GDBErr(hr, "Failed fetching row ");
        return NULL;
    }

    if (hr != S_OK)
        return NULL; //none found - but no failure


    OGRFeature* pOGRFeature = NULL;

    if (!OGRFeatureFromGdbRow(&row,  &pOGRFeature))
    {
        GDBErr(hr, "Failed translating ArcObjects row to OGR Feature");
        return NULL;
    }

    return pOGRFeature;
}


/************************************************************************/
/*                          GetFeatureCount()                           */
/************************************************************************/

int FGdbLayer::GetFeatureCount( int bForce )
{
    long           hr;
    int32          rowCount = 0;

    if (m_pOGRFilterGeometry != NULL || m_wstrWhereClause.size() != 0)
        return OGRLayer::GetFeatureCount(bForce);

    if (FAILED(hr = m_pTable->GetRowCount(rowCount)))
    {
        GDBErr(hr, "Failed counting rows");
        return 0;
    }

#if 0
  Row            row;
  EnumRows       enumRows;

  if (FAILED(hr = m_pTable->Search(StringToWString(m_strOIDFieldName), L"", true, enumRows)))
  {
    GDBErr(hr, "Failed counting rows");
    return -1;
  }

  while (S_OK == (hr = enumRows.Next(row)))
    ++rowCount;

  if (FAILED(hr))
  {
    GDBErr(hr, "Failed counting rows (during fetch)");
    return -1;
  }
#endif

    return static_cast<int>(rowCount);
}



/************************************************************************/
/*                             GetExtent()                              */
/************************************************************************/

OGRErr FGdbLayer::GetExtent (OGREnvelope* psExtent, int bForce)
{
    if (m_pOGRFilterGeometry != NULL || m_wstrWhereClause.size() != 0 ||
        m_strShapeFieldName.size() == 0)
        return OGRLayer::GetExtent(psExtent, bForce);

    long hr;
    Envelope envelope;
    if (FAILED(hr = m_pTable->GetExtent(envelope)))
    {
        GDBErr(hr, "Failed fetching extent");
        return OGRERR_FAILURE;
    }

    psExtent->MinX = envelope.xMin;
    psExtent->MinY = envelope.yMin;
    psExtent->MaxX = envelope.xMax;
    psExtent->MaxY = envelope.yMax;

    if (CPLIsNan(psExtent->MinX) ||
        CPLIsNan(psExtent->MinY) ||
        CPLIsNan(psExtent->MaxX) ||
        CPLIsNan(psExtent->MaxY))
        return OGRERR_FAILURE;

    return OGRERR_NONE;
}

/* OGRErr FGdbLayer::StartTransaction ()
{
    if ( ! m_pTable ) 
        return OGRERR_FAILURE;
        
    m_pTable->LoadOnlyMode(true);
    m_pTable->SetWriteLock();
    return OGRERR_NONE;
    
} */


/* OGRErr FGdbLayer::CommitTransaction ()
{
    if ( ! m_pTable ) 
        return OGRERR_FAILURE;
    
    m_pTable->LoadOnlyMode(false);
    m_pTable->FreeWriteLock();
    return OGRERR_NONE;
} */

/* OGRErr FGdbLayer::RollbackTransaction ()
{
    if ( ! m_pTable ) 
        return OGRERR_FAILURE;
    
    m_pTable->LoadOnlyMode(false);
    m_pTable->FreeWriteLock();
    return OGRERR_NONE;
} */


/************************************************************************/
/*                           GetLayerXML()                              */
/* Return XML definition of the Layer as provided by FGDB. Caller must  */
/* free result.                                                         */
/* Not currently used by the driver, but can be used by external code   */
/* for specific purposes.                                               */
/************************************************************************/

OGRErr FGdbLayer::GetLayerXML (char **ppXml)
{
    long hr;
    std::string xml;

    if ( FAILED(hr = m_pTable->GetDefinition(xml)) )
    {
        GDBErr(hr, "Failed fetching XML table definition");
        return OGRERR_FAILURE;
    }

    *ppXml = CPLStrdup(xml.c_str());
    return OGRERR_NONE;
}

/************************************************************************/
/*                           GetLayerMetadataXML()                      */
/* Return XML metadata for the Layer as provided by FGDB. Caller must  */
/* free result.                                                         */
/* Not currently used by the driver, but can be used by external code   */
/* for specific purposes.                                               */
/************************************************************************/

OGRErr FGdbLayer::GetLayerMetadataXML (char **ppXml)
{
    long hr;
    std::string xml;

    if ( FAILED(hr = m_pTable->GetDocumentation(xml)) )
    {
        GDBErr(hr, "Failed fetching XML table metadata");
        return OGRERR_FAILURE;
    }

    *ppXml = CPLStrdup(xml.c_str());
    return OGRERR_NONE;
}

/************************************************************************/
/*                           TestCapability()                           */
/************************************************************************/

int FGdbLayer::TestCapability( const char* pszCap )
{

    if (EQUAL(pszCap,OLCRandomRead))
        return TRUE;

    else if (EQUAL(pszCap,OLCFastFeatureCount)) 
        return m_pOGRFilterGeometry == NULL && m_wstrWhereClause.size() == 0;

    else if (EQUAL(pszCap,OLCFastSpatialFilter))
        return TRUE;

    else if (EQUAL(pszCap,OLCFastGetExtent))
        return m_pOGRFilterGeometry == NULL && m_wstrWhereClause.size() == 0;

    else if (EQUAL(pszCap,OLCCreateField)) /* CreateField() */
        return TRUE;

    else if (EQUAL(pszCap,OLCSequentialWrite)) /* CreateFeature() */
        return TRUE;

    else if (EQUAL(pszCap,OLCStringsAsUTF8)) /* Native UTF16, converted to UTF8 */
        return TRUE;

    else if (EQUAL(pszCap,OLCReorderFields)) /* TBD ReorderFields() */
        return FALSE;

    else if (EQUAL(pszCap,OLCDeleteFeature)) /* TBD DeleteFeature() */
        return FALSE;

    else if (EQUAL(pszCap,OLCRandomWrite)) /* TBD SetFeature() */
        return FALSE;

    else if (EQUAL(pszCap,OLCDeleteField)) /* TBD DeleteField() */
        return FALSE;
        
    else if (EQUAL(pszCap,OLCFastSetNextByIndex)) /* TBD FastSetNextByIndex() */
        return FALSE;

    else if (EQUAL(pszCap,OLCTransactions)) /* TBD Start/End Transactions() */
        return FALSE;
        
    else 
        return FALSE;
}
