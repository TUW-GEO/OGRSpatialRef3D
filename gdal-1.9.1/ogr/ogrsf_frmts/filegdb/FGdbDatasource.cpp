/******************************************************************************
 * $Id: FGdbDatasource.cpp 23394 2011-11-19 19:20:12Z rouault $
 *
 * Project:  OpenGIS Simple Features Reference Implementation
 * Purpose:  Implements FileGDB OGR Datasource.
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
#include "cpl_conv.h"
#include "cpl_string.h"
#include "gdal.h"
#include "FGdbUtils.h"

CPL_CVSID("$Id: FGdbDatasource.cpp 23394 2011-11-19 19:20:12Z rouault $");

using std::vector;
using std::wstring;

/************************************************************************/
/*                          FGdbDataSource()                           */
/************************************************************************/

FGdbDataSource::FGdbDataSource():
OGRDataSource(),
m_pszName(0), m_pGeodatabase(NULL)
{
}

/************************************************************************/
/*                          ~FGdbDataSource()                          */
/************************************************************************/

FGdbDataSource::~FGdbDataSource()
{   
    CPLFree( m_pszName );

    size_t count = m_layers.size();
    for(size_t i = 0; i < count; ++i )
        delete m_layers[i];

    if (m_pGeodatabase)
    {
        ::CloseGeodatabase(*m_pGeodatabase);
        delete m_pGeodatabase;
    }
}


/************************************************************************/
/*                                Open()                                */
/************************************************************************/

int FGdbDataSource::Open(Geodatabase* pGeodatabase, const char * pszNewName, int bUpdate )
{
    m_pszName = CPLStrdup( pszNewName );
    m_pGeodatabase = pGeodatabase;

    std::vector<std::wstring> typesRequested;

	// We're only interested in Tables, Feature Datasets and Feature Classes
	typesRequested.push_back(L"Feature Class");
	typesRequested.push_back(L"Table");
	typesRequested.push_back(L"Feature Dataset");
	
    bool rv = LoadLayers(L"\\");

    return rv;
}

/************************************************************************/
/*                          OpenFGDBTables()                            */
/************************************************************************/

bool FGdbDataSource::OpenFGDBTables(const std::wstring &type,
                                    const std::vector<std::wstring> &layers)
{
    fgdbError hr;
    for ( unsigned int i = 0; i < layers.size(); i++ )
    {
        Table* pTable = new Table;
        //CPLDebug("FGDB", "Opening %s", WStringToString(layers[i]).c_str());
        if (FAILED(hr = m_pGeodatabase->OpenTable(layers[i], *pTable)))
        {
            delete pTable;
            GDBDebug(hr, "Error opening " + WStringToString(layers[i]) + ". Skipping it");
            continue;
        }
        FGdbLayer* pLayer = new FGdbLayer;
        if (!pLayer->Initialize(this, pTable, layers[i], type))
        {
            delete pLayer;
            return GDBErr(hr, "Error initializing OGRLayer for " + WStringToString(layers[i]));
        }

        m_layers.push_back(pLayer);
    }
    return true;
}

/************************************************************************/
/*                            LoadLayers()                             */
/************************************************************************/

bool FGdbDataSource::LoadLayers(const std::wstring &root) 
{
    std::vector<wstring> tables;
    std::vector<wstring> featureclasses;
    std::vector<wstring> featuredatasets;
    fgdbError hr;

    /* Find all the Tables in the root */
    if ( FAILED(hr = m_pGeodatabase->GetChildDatasets(root, L"Table", tables)) )
    {
        return GDBErr(hr, "Error reading Tables in " + WStringToString(root));
    }
    /* Open the tables we found */
    if ( tables.size() > 0 && ! OpenFGDBTables(L"Table", tables) )
        return false;

    /* Find all the Feature Classes in the root */
    if ( FAILED(hr = m_pGeodatabase->GetChildDatasets(root, L"Feature Class", featureclasses)) )
    {
        return GDBErr(hr, "Error reading Feature Classes in " + WStringToString(root));
    }
    /* Open the tables we found */
    if ( featureclasses.size() > 0 && ! OpenFGDBTables(L"Feature Class", featureclasses) )
        return false;

    /* Find all the Feature Datasets in the root */
    if ( FAILED(hr = m_pGeodatabase->GetChildDatasets(root, L"Feature Dataset", featuredatasets)) )
    {
        return GDBErr(hr, "Error reading Feature Datasets in " + WStringToString(root));
    }
    /* Look for Feature Classes inside the Feature Dataset */
    for ( unsigned int i = 0; i < featuredatasets.size(); i++ )
    {
        if ( FAILED(hr = m_pGeodatabase->GetChildDatasets(featuredatasets[i], L"Feature Class", featureclasses)) )
        {
            return GDBErr(hr, "Error reading Feature Classes in " + WStringToString(featuredatasets[i]));
        }
        if ( featureclasses.size() > 0 && ! OpenFGDBTables(L"Feature Class", featureclasses) )
            return false;
    }
    return true;
}


#if 0
/************************************************************************/
/*                            LoadLayersOld()                              */
/************************************************************************/

/* Old recursive LoadLayers. Removed in favor of simple one that only
   looks at FeatureClasses and Tables. */

// Flattens out hierarchichal GDB structure
bool FGdbDataSource::LoadLayersOld(const std::vector<wstring> & datasetTypes,
                                const wstring & parent)
{
    long hr = S_OK;

    // I didn't find an API to give me the type of the dataset based on name - I am *not*
    // parsing XML for something like this - in the meantime I can use this hack to see
    // if the dataset had any children whatsoever - if so, then I won't attempt to open it
    // otherwise, do attempt to do that

    bool childrenFound = false;
    bool errorsEncountered = false;

    for (size_t dsTypeIndex = 0; dsTypeIndex < datasetTypes.size(); dsTypeIndex++)
    {
        std::vector<wstring> childDatasets;
        m_pGeodatabase->GetChildDatasets( parent, datasetTypes[dsTypeIndex], childDatasets);

        if (childDatasets.size() > 0)
        {
            //it is a container of other datasets

            for (size_t childDatasetIndex = 0;
                 childDatasetIndex < childDatasets.size();
                 childDatasetIndex++)
            {
                childrenFound = true;

                // do something with it
                // For now, we just ignore dataset containers and only open the children
                //std::wcout << datasetTypes[dsTypeIndex] << L" " << childDatasets[childDatasetIndex] << std::endl;

                if (!LoadLayersOld(datasetTypes, childDatasets[childDatasetIndex]))
                    errorsEncountered = true;
            }
        }
    }

    //it is a full fledged dataset itself without children - open it (except the root)

    if ((!childrenFound) && parent != L"\\")
    {
        //wcout << "Opening " << parent << "...";
        Table* pTable = new Table;
        if (FAILED(hr = m_pGeodatabase->OpenTable(parent,*pTable)))
        {
            delete pTable;
            return GDBErr(hr, "Error opening " + WStringToString(parent));
        }

        FGdbLayer* pLayer = new FGdbLayer;

        //pLayer has ownership of the table pointer as soon Initialize is called
        if (!pLayer->Initialize(this, pTable, parent))
        {
            delete pLayer;

            return GDBErr(hr, "Error initializing OGRLayer for " +
                          WStringToString(parent));
        }

        m_layers.push_back(pLayer);
    }

    return !errorsEncountered;
}
#endif


/************************************************************************/
/*                            DeleteLayer()                             */
/************************************************************************/

OGRErr FGdbDataSource::DeleteLayer( int iLayer )
{
    if( iLayer < 0 || iLayer >= static_cast<int>(m_layers.size()) )
        return OGRERR_FAILURE;

    // Fetch FGDBAPI Table before deleting OGR layer object

    Table* pTable = m_layers[iLayer]->GetTable();

    std::string name = m_layers[iLayer]->GetLayerDefn()->GetName();
    std::wstring strPath = m_layers[iLayer]->GetTablePath();
    std::wstring strType = m_layers[iLayer]->GetType();

    // delete OGR layer
    delete m_layers[iLayer];

    pTable = NULL; // OGR Layer had ownership of FGDB Table

    m_layers.erase(m_layers.begin() + iLayer);

    long hr;
  
    if (FAILED(hr = m_pGeodatabase->Delete(strPath, strType)))
    {
        CPLError( CE_Warning, CPLE_AppDefined,
                 "%s was not deleted however it has been closed", name.c_str());
        GDBErr(hr, "Failed deleting dataset");
        return OGRERR_FAILURE;
    }

    return OGRERR_NONE;
}

/************************************************************************/
/*                           TestCapability()                           */
/************************************************************************/

int FGdbDataSource::TestCapability( const char * pszCap )
{
    if( EQUAL(pszCap,ODsCCreateLayer) )
        return TRUE;

    else if( EQUAL(pszCap,ODsCDeleteLayer) )
        return TRUE;

    return FALSE;
}


/************************************************************************/
/*                              GetLayer()                              */
/************************************************************************/

OGRLayer *FGdbDataSource::GetLayer( int iLayer )
{ 
    int count = static_cast<int>(m_layers.size());

    if( iLayer < 0 || iLayer >= count )
        return NULL;
    else
        return m_layers[iLayer];
}

/************************************************************************/
/*                              CreateLayer()                           */
/*                                                                      */
/* See FGdbLayer::Create for creation options                           */
/************************************************************************/

OGRLayer *
FGdbDataSource::CreateLayer( const char * pszLayerName,
                              OGRSpatialReference *poSRS,
                              OGRwkbGeometryType eType,
                              char ** papszOptions )
{
    FGdbLayer* pLayer = new FGdbLayer;
    if (!pLayer->Create(this, pszLayerName, poSRS, eType, papszOptions))
    {
        delete pLayer;
        return NULL;
    }

    m_layers.push_back(pLayer);

    return pLayer;  
}
