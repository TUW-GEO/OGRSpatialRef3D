/******************************************************************************
 * $Id: ogrvrtdatasource.cpp 24156 2012-03-23 21:48:56Z warmerdam $
 *
 * Project:  OpenGIS Simple Features Reference Implementation
 * Purpose:  Implements OGRVRTDataSource class.
 * Author:   Frank Warmerdam, warmerdam@pobox.com
 *
 ******************************************************************************
 * Copyright (c) 2003, Frank Warmerdam <warmerdam@pobox.com>
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

#include "ogr_vrt.h"
#include "cpl_conv.h"
#include "cpl_string.h"

CPL_CVSID("$Id: ogrvrtdatasource.cpp 24156 2012-03-23 21:48:56Z warmerdam $");

/************************************************************************/
/*                          OGRVRTDataSource()                          */
/************************************************************************/

OGRVRTDataSource::OGRVRTDataSource()

{
    pszName = NULL;
    papoLayers = NULL;
    nLayers = 0;
    psTree = NULL;
    nCallLevel = 0;
}

/************************************************************************/
/*                         ~OGRVRTDataSource()                         */
/************************************************************************/

OGRVRTDataSource::~OGRVRTDataSource()

{
    int         i;

    CPLFree( pszName );

    for( i = 0; i < nLayers; i++ )
        delete papoLayers[i];
    
    CPLFree( papoLayers );

    if( psTree != NULL)
        CPLDestroyXMLNode( psTree );
}

/************************************************************************/
/*                             Initialize()                             */
/************************************************************************/

int OGRVRTDataSource::Initialize( CPLXMLNode *psTree, const char *pszNewName,
                                  int bUpdate )

{
    CPLAssert( nLayers == 0 );

    this->psTree = psTree;

/* -------------------------------------------------------------------- */
/*      Set name, and capture the directory path so we can use it       */
/*      for relative datasources.                                       */
/* -------------------------------------------------------------------- */
    CPLString osVRTDirectory = CPLGetPath( pszNewName );

    pszName = CPLStrdup( pszNewName );

/* -------------------------------------------------------------------- */
/*      Look for the OGRVRTDataSource node, it might be after an        */
/*      <xml> node.                                                     */
/* -------------------------------------------------------------------- */
    CPLXMLNode *psVRTDSXML = CPLGetXMLNode( psTree, "=OGRVRTDataSource" );
    if( psVRTDSXML == NULL )
    {
        CPLError( CE_Failure, CPLE_AppDefined, 
                  "Did not find the <OGRVRTDataSource> node in the root of the document,\n"
                  "this is not really an OGR VRT." );
        return FALSE;
    }

/* -------------------------------------------------------------------- */
/*      Look for layers.                                                */
/* -------------------------------------------------------------------- */
    CPLXMLNode *psLTree;

    for( psLTree=psVRTDSXML->psChild; psLTree != NULL; psLTree=psLTree->psNext )
    {
        if( psLTree->eType != CXT_Element
            || !EQUAL(psLTree->pszValue,"OGRVRTLayer") )
            continue;

/* -------------------------------------------------------------------- */
/*      Create the layer object.                                        */
/* -------------------------------------------------------------------- */
        OGRVRTLayer  *poLayer;
        
        poLayer = new OGRVRTLayer(this);
        
        if( !poLayer->FastInitialize( psLTree, osVRTDirectory, bUpdate ) )
        {
            delete poLayer;
            return FALSE;
        }
        
/* -------------------------------------------------------------------- */
/*      Add layer to data source layer list.                            */
/* -------------------------------------------------------------------- */
        papoLayers = (OGRVRTLayer **)
            CPLRealloc( papoLayers,  sizeof(OGRVRTLayer *) * (nLayers+1) );
        papoLayers[nLayers++] = poLayer;
    }

    return TRUE;
}

/************************************************************************/
/*                           TestCapability()                           */
/************************************************************************/

int OGRVRTDataSource::TestCapability( const char * pszCap )

{
    return FALSE;
}

/************************************************************************/
/*                              GetLayer()                              */
/************************************************************************/

OGRLayer *OGRVRTDataSource::GetLayer( int iLayer )

{
    if( iLayer < 0 || iLayer >= nLayers )
        return NULL;
    else
        return papoLayers[iLayer];
}

/************************************************************************/
/*                         AddForbiddenNames()                          */
/************************************************************************/

void OGRVRTDataSource::AddForbiddenNames(const char* pszOtherDSName)
{
    aosOtherDSNameSet.insert(pszOtherDSName);
}

/************************************************************************/
/*                         IsInForbiddenNames()                         */
/************************************************************************/

int OGRVRTDataSource::IsInForbiddenNames(const char* pszOtherDSName)
{
    return aosOtherDSNameSet.find(pszOtherDSName) != aosOtherDSNameSet.end();
}
