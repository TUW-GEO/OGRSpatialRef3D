/******************************************************************************
 * $Id: test_ogrsf.cpp 23602 2011-12-19 21:35:25Z rouault $
 *
 * Project:  OpenGIS Simple Features Reference Implementation
 * Purpose:  Formal test harnass for OGRLayer implementations.
 * Author:   Frank Warmerdam, warmerdam@pobox.com
 *
 ******************************************************************************
 * Copyright (c) 1999, Frank Warmerdam
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

#include "ogrsf_frmts.h"
#include "cpl_conv.h"
#include "ogr_api.h"
#include "ogr_p.h"

CPL_CVSID("$Id: test_ogrsf.cpp 23602 2011-12-19 21:35:25Z rouault $");

int     bReadOnly = FALSE;
int     bVerbose = TRUE;

static void Usage();
static int TestOGRLayer( OGRDataSource * poDS, OGRLayer * poLayer, int bIsSQLLayer );
static int TestInterleavedReading( const char* pszDataSource, char** papszLayers );

/************************************************************************/
/*                                main()                                */
/************************************************************************/

int main( int nArgc, char ** papszArgv )

{
    const char  *pszDataSource = NULL;
    char** papszLayers = NULL;
    const char  *pszSQLStatement = NULL;
    int bRet = TRUE;

    /* Must process OGR_SKIP before OGRRegisterAll(), but we can't call */
    /* OGRGeneralCmdLineProcessor before it needs the drivers to be registered */
    /* for the --format or --formats options */
    for( int iArg = 1; iArg < nArgc; iArg++ )
    {
        if( EQUAL(papszArgv[iArg], "--config") && iArg + 2 < nArgc &&
            EQUAL(papszArgv[iArg+1], "OGR_SKIP") )
        {
            CPLSetConfigOption(papszArgv[iArg+1], papszArgv[iArg+2]);
            break;
        }
    }
    
/* -------------------------------------------------------------------- */
/*      Register format(s).                                             */
/* -------------------------------------------------------------------- */
    OGRRegisterAll();

/* -------------------------------------------------------------------- */
/*      Processing command line arguments.                              */
/* -------------------------------------------------------------------- */
    nArgc = OGRGeneralCmdLineProcessor( nArgc, &papszArgv, 0 );

    if( nArgc < 1 )
        exit( -nArgc );

/* -------------------------------------------------------------------- */
/*      Processing command line arguments.                              */
/* -------------------------------------------------------------------- */
    for( int iArg = 1; iArg < nArgc; iArg++ )
    {
        if( EQUAL(papszArgv[iArg], "--utility_version") )
        {
            printf("%s was compiled against GDAL %s and is running against GDAL %s\n",
                   papszArgv[0], GDAL_RELEASE_NAME, GDALVersionInfo("RELEASE_NAME"));
            return 0;
        }
        else if( EQUAL(papszArgv[iArg],"-ro") )
            bReadOnly = TRUE;
        else if( EQUAL(papszArgv[iArg],"-q") || EQUAL(papszArgv[iArg],"-quiet"))
            bVerbose = FALSE;
        else if( EQUAL(papszArgv[iArg],"-sql") && iArg + 1 < nArgc)
            pszSQLStatement = papszArgv[++iArg];
        else if( papszArgv[iArg][0] == '-' )
        {
            Usage();
        }
        else if (pszDataSource == NULL)
            pszDataSource = papszArgv[iArg];
        else
            papszLayers = CSLAddString(papszLayers, papszArgv[iArg]);
    }

    if( pszDataSource == NULL )
        Usage();

/* -------------------------------------------------------------------- */
/*      Open data source.                                               */
/* -------------------------------------------------------------------- */
    OGRDataSource       *poDS;
    OGRSFDriver         *poDriver;

    poDS = OGRSFDriverRegistrar::Open( pszDataSource, !bReadOnly, &poDriver );
    if( poDS == NULL && !bReadOnly )
    {
        poDS = OGRSFDriverRegistrar::Open( pszDataSource, FALSE, &poDriver );
        if( poDS != NULL && bVerbose )
        {
            printf( "Had to open data source read-only.\n" );
            bReadOnly = TRUE;
        }
    }

/* -------------------------------------------------------------------- */
/*      Report failure                                                  */
/* -------------------------------------------------------------------- */
    if( poDS == NULL )
    {
        OGRSFDriverRegistrar    *poR = OGRSFDriverRegistrar::GetRegistrar();
        
        printf( "FAILURE:\n"
                "Unable to open datasource `%s' with the following drivers.\n",
                pszDataSource );

        for( int iDriver = 0; iDriver < poR->GetDriverCount(); iDriver++ )
        {
            printf( "  -> %s\n", poR->GetDriver(iDriver)->GetName() );
        }

        exit( 1 );
    }

/* -------------------------------------------------------------------- */
/*      Some information messages.                                      */
/* -------------------------------------------------------------------- */
    if( bVerbose )
        printf( "INFO: Open of `%s' using driver `%s' successful.\n",
                pszDataSource, poDriver->GetName() );

    if( bVerbose && !EQUAL(pszDataSource,poDS->GetName()) )
    {
        printf( "INFO: Internal data source name `%s'\n"
                "      different from user name `%s'.\n",
                poDS->GetName(), pszDataSource );
    }
    
/* -------------------------------------------------------------------- */
/*      Process optionnal SQL request.                                  */
/* -------------------------------------------------------------------- */
    if (pszSQLStatement != NULL)
    {
        OGRLayer  *poResultSet = poDS->ExecuteSQL(pszSQLStatement, NULL, NULL);
        if (poResultSet == NULL)
            exit(1);
            
        printf( "INFO: Testing layer %s.\n",
                    poResultSet->GetName() );
        bRet = TestOGRLayer( poDS, poResultSet, TRUE );
        
        poDS->ReleaseResultSet(poResultSet);
    }
/* -------------------------------------------------------------------- */
/*      Process each data source layer.                                 */
/* -------------------------------------------------------------------- */
    else if (papszLayers == NULL)
    {
        for( int iLayer = 0; iLayer < poDS->GetLayerCount(); iLayer++ )
        {
            OGRLayer        *poLayer = poDS->GetLayer(iLayer);

            if( poLayer == NULL )
            {
                printf( "FAILURE: Couldn't fetch advertised layer %d!\n",
                        iLayer );
                exit( 1 );
            }

            printf( "INFO: Testing layer %s.\n",
                    poLayer->GetName() );
            bRet &= TestOGRLayer( poDS, poLayer, FALSE );
        }

        if (poDS->GetLayerCount() >= 2)
        {
            OGRDataSource::DestroyDataSource(poDS);
            poDS = NULL;
            bRet &= TestInterleavedReading( pszDataSource, NULL );
        }
    }
    else
    {
/* -------------------------------------------------------------------- */
/*      Or process layers specified by the user                         */
/* -------------------------------------------------------------------- */
        char** papszLayerIter = papszLayers;
        while (*papszLayerIter)
        {
            OGRLayer        *poLayer = poDS->GetLayerByName(*papszLayerIter);

            if( poLayer == NULL )
            {
                printf( "FAILURE: Couldn't fetch requested layer %s!\n",
                        *papszLayerIter );
                exit( 1 );
            }
            
            printf( "INFO: Testing layer %s.\n",
                    poLayer->GetName() );
            bRet &= TestOGRLayer( poDS, poLayer, FALSE );
            
            papszLayerIter ++;
        }

        if (CSLCount(papszLayers) >= 2)
        {
            OGRDataSource::DestroyDataSource(poDS);
            poDS = NULL;
            bRet &= TestInterleavedReading( pszDataSource, papszLayers );
        }
    }

/* -------------------------------------------------------------------- */
/*      Close down.                                                     */
/* -------------------------------------------------------------------- */
    OGRDataSource::DestroyDataSource(poDS);

    OGRCleanupAll();

    CSLDestroy(papszLayers);
    CSLDestroy(papszArgv);
    
#ifdef DBMALLOC
    malloc_dump(1);
#endif
    
    return (bRet) ? 0 : 1;
}

/************************************************************************/
/*                               Usage()                                */
/************************************************************************/

static void Usage()

{
    printf( "Usage: test_ogrsf [-ro] [-q] datasource_name [[layer1_name, layer2_name, ...] | [-sql statement]]\n" );
    exit( 1 );
}

/************************************************************************/
/*                      TestOGRLayerFeatureCount()                      */
/*                                                                      */
/*      Verify that the feature count matches the actual number of      */
/*      features returned during sequential reading.                    */
/************************************************************************/

static int TestOGRLayerFeatureCount( OGRDataSource* poDS, OGRLayer *poLayer, int bIsSQLLayer )

{
    int bRet = TRUE;
    int         nFC = 0, nClaimedFC = poLayer->GetFeatureCount();
    OGRFeature  *poFeature;
    OGRSpatialReference * poSRS = poLayer->GetSpatialRef();
    int         bWarnAboutSRS = FALSE;
    OGRFeatureDefn* poLayerDefn = poLayer->GetLayerDefn();

    poLayer->ResetReading();

    while( (poFeature = poLayer->GetNextFeature()) != NULL )
    {
        nFC++;

        if (poFeature->GetDefnRef() != poLayerDefn)
        {
            bRet = FALSE;
            printf( "ERROR: Feature defn differs from layer defn.\n"
                    "Feature defn = %p\n"
                    "Layer defn = %p\n",
                     poFeature->GetDefnRef(), poLayerDefn);
        }

        if( poFeature->GetGeometryRef() != NULL
            && poFeature->GetGeometryRef()->getSpatialReference() != poSRS
            && !bWarnAboutSRS )
        {
            char        *pszLayerSRSWKT, *pszFeatureSRSWKT;
            
            bWarnAboutSRS = TRUE;

            if( poSRS != NULL )
                poSRS->exportToWkt( &pszLayerSRSWKT );
            else
                pszLayerSRSWKT = CPLStrdup("(NULL)");

            if( poFeature->GetGeometryRef()->getSpatialReference() != NULL )
                poFeature->GetGeometryRef()->
                    getSpatialReference()->exportToWkt( &pszFeatureSRSWKT );
            else
                pszFeatureSRSWKT = CPLStrdup("(NULL)");

            bRet = FALSE;
            printf( "ERROR: Feature SRS differs from layer SRS.\n"
                    "Feature SRS = %s (%p)\n"
                    "Layer SRS = %s (%p)\n",
                    pszFeatureSRSWKT, poFeature->GetGeometryRef()->getSpatialReference(),
                    pszLayerSRSWKT, poSRS );
            CPLFree( pszLayerSRSWKT );
            CPLFree( pszFeatureSRSWKT );
        }
        
        OGRFeature::DestroyFeature(poFeature);
    }

    if( nFC != nClaimedFC )
    {
        bRet = FALSE;
        printf( "ERROR: Claimed feature count %d doesn't match actual, %d.\n",
                nClaimedFC, nFC );
    }
    else if( nFC != poLayer->GetFeatureCount() )
    {
        bRet = FALSE;
        printf( "ERROR: Feature count at end of layer %d differs "
                "from at start, %d.\n",
                nFC, poLayer->GetFeatureCount() );
    }
    else if( bVerbose )
        printf( "INFO: Feature count verified.\n" );
        
    if (!bIsSQLLayer)
    {
        CPLString osSQL;
        const char* pszLayerName = poLayer->GetName();
        int i;
        char ch;
        for(i=0;(ch = pszLayerName[i]) != 0;i++)
        {
            if (ch >= '0' && ch <= '9')
            {
                if (i == 0)
                    break;
            }
            else if (!((ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z')))
                break;
        }
        /* Only quote if needed. Quoting conventions depend on the driver... */
        if (ch == 0)
            osSQL.Printf("SELECT COUNT(*) FROM %s", pszLayerName);
        else
        {
            if (EQUAL(poDS->GetDriver()->GetName(), "MYSQL"))
                osSQL.Printf("SELECT COUNT(*) FROM `%s`", pszLayerName);
            else if (EQUAL(poDS->GetDriver()->GetName(), "PostgreSQL") &&
                     strchr(pszLayerName, '.'))
            {
                char** papszTokens = CSLTokenizeStringComplex(pszLayerName, ".", 0, 0);
                if (CSLCount(papszTokens) == 2)
                {
                    osSQL.Printf("SELECT COUNT(*) FROM \"%s\".\"%s\"", papszTokens[0], papszTokens[1]);
                }
                else
                    osSQL.Printf("SELECT COUNT(*) FROM \"%s\"", pszLayerName);
                CSLDestroy(papszTokens);
            }
            else if (EQUAL(poDS->GetDriver()->GetName(), "SQLAnywhere"))
                osSQL.Printf("SELECT COUNT(*) FROM %s", pszLayerName);
            else
                osSQL.Printf("SELECT COUNT(*) FROM \"%s\"", pszLayerName);
        }
        OGRLayer* poSQLLyr = poDS->ExecuteSQL(osSQL.c_str(), NULL, NULL);
        if (poSQLLyr)
        {
            OGRFeature* poFeatCount = poSQLLyr->GetNextFeature();
            if (poFeatCount == NULL)
            {
                bRet = FALSE;
                printf( "ERROR: '%s' failed.\n", osSQL.c_str() );
            }
            else if (nClaimedFC != poFeatCount->GetFieldAsInteger(0))
            {
                bRet = FALSE;
                printf( "ERROR: Claimed feature count %d doesn't match '%s' one, %d.\n",
                        nClaimedFC, osSQL.c_str(), poFeatCount->GetFieldAsInteger(0) );
            }
            OGRFeature::DestroyFeature(poFeatCount);
            poDS->ReleaseResultSet(poSQLLyr);
        }
    }

    if( bVerbose && !bWarnAboutSRS )
    {
        printf("INFO: Feature/layer spatial ref. consistency verified.\n");
    }

    return bRet;
}

/************************************************************************/
/*                       TestOGRLayerRandomRead()                       */
/*                                                                      */
/*      Read the first 5 features, and then try to use random           */
/*      reading to reread 2 and 5 and verify that this works OK.        */
/*      Don't attempt if there aren't at least 5 features.              */
/************************************************************************/

static int TestOGRLayerRandomRead( OGRLayer *poLayer )

{
    int bRet = TRUE;
    OGRFeature  *papoFeatures[5], *poFeature = NULL;
    int         iFeature;

    poLayer->SetSpatialFilter( NULL );
    
    if( poLayer->GetFeatureCount() < 5 )
    {
        if( bVerbose )
            printf( "INFO: Only %d features on layer,"
                    "skipping random read test.\n",
                    poLayer->GetFeatureCount() );
        
        return bRet;
    }

/* -------------------------------------------------------------------- */
/*      Fetch five features.                                            */
/* -------------------------------------------------------------------- */
    poLayer->ResetReading();
    
    for( iFeature = 0; iFeature < 5; iFeature++ )
    {
        papoFeatures[iFeature] = NULL;
    }
    for( iFeature = 0; iFeature < 5; iFeature++ )
    {
        papoFeatures[iFeature] = poLayer->GetNextFeature();
        if( papoFeatures[iFeature] == NULL )
        {
            if( bVerbose )
                printf( "INFO: Only %d features on layer,"
                        "skipping random read test.\n",
                        iFeature );
            goto end;
        }
    }

/* -------------------------------------------------------------------- */
/*      Test feature 2.                                                 */
/* -------------------------------------------------------------------- */
    poFeature = poLayer->GetFeature( papoFeatures[1]->GetFID() );
    if (poFeature == NULL)
    {
        printf( "ERROR: Cannot fetch feature %ld.\n",
                 papoFeatures[1]->GetFID() );
        goto end;
    }

    if( !poFeature->Equal( papoFeatures[1] ) )
    {
        bRet = FALSE;
        printf( "ERROR: Attempt to randomly read feature %ld appears to\n"
                "       have returned a different feature than sequential\n"
                "       reading indicates should have happened.\n",
                papoFeatures[1]->GetFID() );
        poFeature->DumpReadable(stdout);
        papoFeatures[1]->DumpReadable(stdout);

        goto end;
    }

    OGRFeature::DestroyFeature(poFeature);

/* -------------------------------------------------------------------- */
/*      Test feature 5.                                                 */
/* -------------------------------------------------------------------- */
    poFeature = poLayer->GetFeature( papoFeatures[4]->GetFID() );
    if( !poFeature->Equal( papoFeatures[4] ) )
    {
        bRet = FALSE;
        printf( "ERROR: Attempt to randomly read feature %ld appears to\n"
                "       have returned a different feature than sequential\n"
                "       reading indicates should have happened.\n",
                papoFeatures[4]->GetFID() );

        goto end;
    }

    if( bVerbose )
        printf( "INFO: Random read test passed.\n" );

end:
    OGRFeature::DestroyFeature(poFeature);

/* -------------------------------------------------------------------- */
/*      Cleanup.                                                        */
/* -------------------------------------------------------------------- */
    for( iFeature = 0; iFeature < 5; iFeature++ )
        OGRFeature::DestroyFeature(papoFeatures[iFeature]);

    return bRet;
}


/************************************************************************/
/*                       TestOGRLayerSetNextByIndex()                   */
/*                                                                      */
/************************************************************************/

static int TestOGRLayerSetNextByIndex( OGRLayer *poLayer )

{
    int bRet = TRUE;
    OGRFeature  *papoFeatures[5], *poFeature = NULL;
    int         iFeature;

    memset(papoFeatures, 0, sizeof(papoFeatures));

    poLayer->SetSpatialFilter( NULL );
    
    if( poLayer->GetFeatureCount() < 5 )
    {
        if( bVerbose )
            printf( "INFO: Only %d features on layer,"
                    "skipping SetNextByIndex test.\n",
                    poLayer->GetFeatureCount() );
        
        return bRet;
    }

/* -------------------------------------------------------------------- */
/*      Fetch five features.                                            */
/* -------------------------------------------------------------------- */
    poLayer->ResetReading();
    
    for( iFeature = 0; iFeature < 5; iFeature++ )
    {
        papoFeatures[iFeature] = poLayer->GetNextFeature();
        if( papoFeatures[iFeature] == NULL )
        {
            bRet = FALSE;
            printf( "ERROR: Cannot get feature %d.\n", iFeature );
            goto end;
        }
    }

/* -------------------------------------------------------------------- */
/*      Test feature at index 1.                                        */
/* -------------------------------------------------------------------- */
    if (poLayer->SetNextByIndex(1) != OGRERR_NONE)
    {
        bRet = FALSE;
        printf( "ERROR: SetNextByIndex(%d) failed.\n", 1 );
        goto end;
    }
    
    poFeature = poLayer->GetNextFeature();
    if( !poFeature->Equal( papoFeatures[1] ) )
    {
        bRet = FALSE;
        printf( "ERROR: Attempt to read feature at index %d appears to\n"
                "       have returned a different feature than sequential\n"
                "       reading indicates should have happened.\n",
                1 );

        goto end;
    }

    OGRFeature::DestroyFeature(poFeature);
    
    poFeature = poLayer->GetNextFeature();
    if( !poFeature->Equal( papoFeatures[2] ) )
    {
        bRet = FALSE;
        printf( "ERROR: Attempt to read feature after feature at index %d appears to\n"
                "       have returned a different feature than sequential\n"
                "       reading indicates should have happened.\n",
                1 );

        goto end;
    }

    OGRFeature::DestroyFeature(poFeature);
    poFeature = NULL;
    
/* -------------------------------------------------------------------- */
/*      Test feature at index 3.                                        */
/* -------------------------------------------------------------------- */
    if (poLayer->SetNextByIndex(3) != OGRERR_NONE)
    {
        bRet = FALSE;
        printf( "ERROR: SetNextByIndex(%d) failed.\n", 3 );
        goto end;
    }
    
    poFeature = poLayer->GetNextFeature();
    if( !poFeature->Equal( papoFeatures[3] ) )
    {
        bRet = FALSE;
        printf( "ERROR: Attempt to read feature at index %d appears to\n"
                "       have returned a different feature than sequential\n"
                "       reading indicates should have happened.\n",
                3 );

        goto end;
    }

    OGRFeature::DestroyFeature(poFeature);
    
    poFeature = poLayer->GetNextFeature();
    if( !poFeature->Equal( papoFeatures[4] ) )
    {
        bRet = FALSE;
        printf( "ERROR: Attempt to read feature after feature at index %d appears to\n"
                "       have returned a different feature than sequential\n"
                "       reading indicates should have happened.\n",
                3 );

        goto end;
    }


    if( bVerbose )
        printf( "INFO: SetNextByIndex() read test passed.\n" );

end:
    OGRFeature::DestroyFeature(poFeature);

/* -------------------------------------------------------------------- */
/*      Cleanup.                                                        */
/* -------------------------------------------------------------------- */
    for( iFeature = 0; iFeature < 5; iFeature++ )
        OGRFeature::DestroyFeature(papoFeatures[iFeature]);

    return bRet;
}

/************************************************************************/
/*                      TestOGRLayerRandomWrite()                       */
/*                                                                      */
/*      Test random writing by trying to switch the 2nd and 5th         */
/*      features.                                                       */
/************************************************************************/

static int TestOGRLayerRandomWrite( OGRLayer *poLayer )

{
    int bRet = TRUE;
    OGRFeature  *papoFeatures[5], *poFeature;
    int         iFeature;
    long        nFID2, nFID5;

    memset(papoFeatures, 0, sizeof(papoFeatures));

    poLayer->SetSpatialFilter( NULL );

    if( poLayer->GetFeatureCount() < 5 )
    {
        if( bVerbose )
            printf( "INFO: Only %d features on layer,"
                    "skipping random write test.\n",
                    poLayer->GetFeatureCount() );
        
        return bRet;
    }

    if( !poLayer->TestCapability( OLCRandomRead ) )
    {
        if( bVerbose )
            printf( "INFO: Skipping random write test since this layer "
                    "doesn't support random read.\n" );
        return bRet;
    }

/* -------------------------------------------------------------------- */
/*      Fetch five features.                                            */
/* -------------------------------------------------------------------- */
    poLayer->ResetReading();
    
    for( iFeature = 0; iFeature < 5; iFeature++ )
    {
        papoFeatures[iFeature] = poLayer->GetNextFeature();
        if( papoFeatures[iFeature] == NULL )
        {
            bRet = FALSE;
            printf( "ERROR: Cannot get feature %d.\n", iFeature );
            goto end;
        }
    }

/* -------------------------------------------------------------------- */
/*      Switch feature ids of feature 2 and 5.                          */
/* -------------------------------------------------------------------- */
    nFID2 = papoFeatures[1]->GetFID();
    nFID5 = papoFeatures[4]->GetFID();

    papoFeatures[1]->SetFID( nFID5 );
    papoFeatures[4]->SetFID( nFID2 );

/* -------------------------------------------------------------------- */
/*      Rewrite them.                                                   */
/* -------------------------------------------------------------------- */
    if( poLayer->SetFeature( papoFeatures[1] ) != OGRERR_NONE )
    {
        bRet = FALSE;
        printf( "ERROR: Attempt to SetFeature(1) failed.\n" );
        goto end;
    }
    if( poLayer->SetFeature( papoFeatures[4] ) != OGRERR_NONE )
    {
        bRet = FALSE;
        printf( "ERROR: Attempt to SetFeature(4) failed.\n" );
        goto end;
    }

/* -------------------------------------------------------------------- */
/*      Now re-read feature 2 to verify the effect stuck.               */
/* -------------------------------------------------------------------- */
    poFeature = poLayer->GetFeature( nFID5 );
    if(poFeature == NULL)
    {
        bRet = FALSE;
        printf( "ERROR: Attempt to GetFeature( nFID5 ) failed.\n" );
        goto end;
    }
    if( !poFeature->Equal(papoFeatures[1]) )
    {
        bRet = FALSE;
        printf( "ERROR: Written feature didn't seem to retain value.\n" );
    }
    else
    {
        printf( "INFO: Random write test passed.\n" );
    }
    OGRFeature::DestroyFeature(poFeature);

/* -------------------------------------------------------------------- */
/*      Re-invert the features to restore to original state             */
/* -------------------------------------------------------------------- */

    papoFeatures[1]->SetFID( nFID2 );
    papoFeatures[4]->SetFID( nFID5 );

    if( poLayer->SetFeature( papoFeatures[1] ) != OGRERR_NONE )
    {
        bRet = FALSE;
        printf( "ERROR: Attempt to restore SetFeature(1) failed.\n" );
    }
    if( poLayer->SetFeature( papoFeatures[4] ) != OGRERR_NONE )
    {
        bRet = FALSE;
        printf( "ERROR: Attempt to restore SetFeature(4) failed.\n" );
    }

end:
/* -------------------------------------------------------------------- */
/*      Cleanup.                                                        */
/* -------------------------------------------------------------------- */

    for( iFeature = 0; iFeature < 5; iFeature++ )
        OGRFeature::DestroyFeature(papoFeatures[iFeature]);

    return bRet;
}

/************************************************************************/
/*                         TestSpatialFilter()                          */
/*                                                                      */
/*      This is intended to be a simple test of the spatial             */
/*      filtering.  We read the first feature.  Then construct a        */
/*      spatial filter geometry which includes it, install and          */
/*      verify that we get the feature.  Next install a spatial         */
/*      filter that doesn't include this feature, and test again.       */
/************************************************************************/

static int TestSpatialFilter( OGRLayer *poLayer )

{
    int bRet = TRUE;
    OGRFeature  *poFeature, *poTargetFeature;
    OGRPolygon  oInclusiveFilter, oExclusiveFilter;
    OGRLinearRing oRing;
    OGREnvelope sEnvelope;
    int         nInclusiveCount;

/* -------------------------------------------------------------------- */
/*      Read the target feature.                                        */
/* -------------------------------------------------------------------- */
    poLayer->ResetReading();
    poTargetFeature = poLayer->GetNextFeature();

    if( poTargetFeature == NULL )
    {
        printf( "INFO: Skipping Spatial Filter test for %s.\n"
                "      No features in layer.\n",
                poLayer->GetName() );
        return bRet;
    }

    if( poTargetFeature->GetGeometryRef() == NULL )
    {
        printf( "INFO: Skipping Spatial Filter test for %s,\n"
                "      target feature has no geometry.\n",
                poTargetFeature->GetDefnRef()->GetName() );
        OGRFeature::DestroyFeature(poTargetFeature);
        return bRet;
    }

    poTargetFeature->GetGeometryRef()->getEnvelope( &sEnvelope );

/* -------------------------------------------------------------------- */
/*      Construct inclusive filter.                                     */
/* -------------------------------------------------------------------- */
    
    oRing.setPoint( 0, sEnvelope.MinX - 20.0, sEnvelope.MinY - 20.0 );
    oRing.setPoint( 1, sEnvelope.MinX - 20.0, sEnvelope.MaxY + 10.0 );
    oRing.setPoint( 2, sEnvelope.MaxX + 10.0, sEnvelope.MaxY + 10.0 );
    oRing.setPoint( 3, sEnvelope.MaxX + 10.0, sEnvelope.MinY - 20.0 );
    oRing.setPoint( 4, sEnvelope.MinX - 20.0, sEnvelope.MinY - 20.0 );
    
    oInclusiveFilter.addRing( &oRing );

    poLayer->SetSpatialFilter( &oInclusiveFilter );

/* -------------------------------------------------------------------- */
/*      Verify that we can find the target feature.                     */
/* -------------------------------------------------------------------- */
    poLayer->ResetReading();

    while( (poFeature = poLayer->GetNextFeature()) != NULL )
    {
        if( poFeature->Equal(poTargetFeature) )
        {
            OGRFeature::DestroyFeature(poFeature);
            break;
        }
        else
            OGRFeature::DestroyFeature(poFeature);
    }

    if( poFeature == NULL )
    {
        bRet = FALSE;
        printf( "ERROR: Spatial filter eliminated a feature unexpectedly!\n");
    }
    else if( bVerbose )
    {
        printf( "INFO: Spatial filter inclusion seems to work.\n" );
    }

    nInclusiveCount = poLayer->GetFeatureCount();

/* -------------------------------------------------------------------- */
/*      Construct exclusive filter.                                     */
/* -------------------------------------------------------------------- */
    oRing.setPoint( 0, sEnvelope.MinX - 20.0, sEnvelope.MinY - 20.0 );
    oRing.setPoint( 1, sEnvelope.MinX - 10.0, sEnvelope.MinY - 20.0 );
    oRing.setPoint( 2, sEnvelope.MinX - 10.0, sEnvelope.MinY - 10.0 );
    oRing.setPoint( 3, sEnvelope.MinX - 20.0, sEnvelope.MinY - 10.0 );
    oRing.setPoint( 4, sEnvelope.MinX - 20.0, sEnvelope.MinY - 20.0 );
    
    oExclusiveFilter.addRing( &oRing );

    poLayer->SetSpatialFilter( &oExclusiveFilter );

/* -------------------------------------------------------------------- */
/*      Verify that we can find the target feature.                     */
/* -------------------------------------------------------------------- */
    poLayer->ResetReading();

    while( (poFeature = poLayer->GetNextFeature()) != NULL )
    {
        if( poFeature->Equal(poTargetFeature) )
        {
            OGRFeature::DestroyFeature(poFeature);
            break;
        }
        else
            OGRFeature::DestroyFeature(poFeature);
    }

    if( poFeature != NULL )
    {
        bRet = FALSE;
        printf( "ERROR: Spatial filter failed to eliminate"
                "a feature unexpectedly!\n");
    }
    else if( poLayer->GetFeatureCount() >= nInclusiveCount )
    {
        bRet = FALSE;
        printf( "ERROR: GetFeatureCount() may not be taking spatial "
                "filter into account.\n" );
    }
    else if( bVerbose )
    {
        printf( "INFO: Spatial filter exclusion seems to work.\n" );
    }

    OGRFeature::DestroyFeature(poTargetFeature);

    poLayer->SetSpatialFilter( NULL );

    return bRet;
}


/************************************************************************/
/*                      TestAttributeFilter()                           */
/*                                                                      */
/*      This is intended to be a simple test of the attribute           */
/*      filtering.  We read the first feature.  Then construct a        */
/*      attribute filter which includes it, install and                 */
/*      verify that we get the feature.  Next install a attribute       */
/*      filter that doesn't include this feature, and test again.       */
/************************************************************************/

static int TestAttributeFilter( OGRDataSource* poDS, OGRLayer *poLayer )

{
    int bRet = TRUE;
    OGRFeature  *poFeature, *poTargetFeature;
    int         nInclusiveCount, nExclusiveCount, nTotalCount;
    CPLString osAttributeFilter;

/* -------------------------------------------------------------------- */
/*      Read the target feature.                                        */
/* -------------------------------------------------------------------- */
    poLayer->ResetReading();
    poTargetFeature = poLayer->GetNextFeature();

    if( poTargetFeature == NULL )
    {
        printf( "INFO: Skipping Attribute Filter test for %s.\n"
                "      No features in layer.\n",
                poLayer->GetName() );
        return bRet;
    }

    int i;
    OGRFieldType eType = OFTString;
    for(i=0;i<poTargetFeature->GetFieldCount();i++)
    {
        eType = poTargetFeature->GetFieldDefnRef(i)->GetType();
        if (poTargetFeature->IsFieldSet(i) &&
            (eType == OFTString || eType == OFTInteger || eType == OFTReal))
        {
            break;
        }
    }
    if( i == poTargetFeature->GetFieldCount() )
    {
        printf( "INFO: Skipping Attribute Filter test for %s.\n"
                "      Could not find non NULL field.\n",
                poLayer->GetName() );
        OGRFeature::DestroyFeature(poTargetFeature);
        return bRet;
    }

    const char* pszFieldName = poTargetFeature->GetFieldDefnRef(i)->GetNameRef();
    CPLString osValue = poTargetFeature->GetFieldAsString(i);

/* -------------------------------------------------------------------- */
/*      Construct inclusive filter.                                     */
/* -------------------------------------------------------------------- */

    if (EQUAL(poDS->GetDriver()->GetName(), "PostgreSQL") &&
        (strchr(pszFieldName, '_') || strchr(pszFieldName, ' ')))
    {
        osAttributeFilter = "\"";
        osAttributeFilter += pszFieldName;
        osAttributeFilter += "\"";
    }
    else if (strchr(pszFieldName, ' ') || pszFieldName[0] == '_')
    {
        osAttributeFilter = "'";
        osAttributeFilter += pszFieldName;
        osAttributeFilter += "'";
    }
    else
        osAttributeFilter = pszFieldName;
    osAttributeFilter += " = ";
    if (eType == OFTString)
        osAttributeFilter += "'";
    osAttributeFilter += osValue;
    if (eType == OFTString)
        osAttributeFilter += "'";
    /* Make sure that the literal will be recognized as a float value */
    /* to avoid int underflow/overflow */
    else if (eType == OFTReal && strchr(osValue, '.') == NULL)
        osAttributeFilter += ".";
    poLayer->SetAttributeFilter( osAttributeFilter );

/* -------------------------------------------------------------------- */
/*      Verify that we can find the target feature.                     */
/* -------------------------------------------------------------------- */
    poLayer->ResetReading();

    while( (poFeature = poLayer->GetNextFeature()) != NULL )
    {
        if( poFeature->Equal(poTargetFeature) )
        {
            OGRFeature::DestroyFeature(poFeature);
            break;
        }
        else
            OGRFeature::DestroyFeature(poFeature);
    }

    if( poFeature == NULL )
    {
        bRet = FALSE;
        printf( "ERROR: Attribute filter eliminated a feature unexpectedly!\n");
    }
    else if( bVerbose )
    {
        printf( "INFO: Attribute filter inclusion seems to work.\n" );
    }

    nInclusiveCount = poLayer->GetFeatureCount();

/* -------------------------------------------------------------------- */
/*      Construct exclusive filter.                                     */
/* -------------------------------------------------------------------- */
    if (EQUAL(poDS->GetDriver()->GetName(), "PostgreSQL") &&
        (strchr(pszFieldName, '_') || strchr(pszFieldName, ' ')))
    {
        osAttributeFilter = "\"";
        osAttributeFilter += pszFieldName;
        osAttributeFilter += "\"";
    }
    else if (strchr(pszFieldName, ' ') || pszFieldName[0] == '_')
    {
        osAttributeFilter = "'";
        osAttributeFilter += pszFieldName;
        osAttributeFilter += "'";
    }
    else
        osAttributeFilter = pszFieldName;
    osAttributeFilter += " <> ";
    if (eType == OFTString)
        osAttributeFilter += "'";
    osAttributeFilter += osValue;
    if (eType == OFTString)
        osAttributeFilter += "'";
    /* Make sure that the literal will be recognized as a float value */
    /* to avoid int underflow/overflow */
    else if (eType == OFTReal && strchr(osValue, '.') == NULL)
        osAttributeFilter += ".";
    poLayer->SetAttributeFilter( osAttributeFilter );

/* -------------------------------------------------------------------- */
/*      Verify that we can find the target feature.                     */
/* -------------------------------------------------------------------- */
    poLayer->ResetReading();

    int nExclusiveCountWhileIterating = 0;
    while( (poFeature = poLayer->GetNextFeature()) != NULL )
    {
        if( poFeature->Equal(poTargetFeature) )
        {
            OGRFeature::DestroyFeature(poFeature);
            break;
        }
        else
            OGRFeature::DestroyFeature(poFeature);
        nExclusiveCountWhileIterating ++;
    }

    nExclusiveCount = poLayer->GetFeatureCount();

    poLayer->SetAttributeFilter( NULL );

    nTotalCount = poLayer->GetFeatureCount();

    if( poFeature != NULL )
    {
        bRet = FALSE;
        printf( "ERROR: Attribute filter failed to eliminate "
                "a feature unexpectedly!\n");
    }
    else if( nExclusiveCountWhileIterating != nExclusiveCount ||
             nExclusiveCount >= nTotalCount ||
             nInclusiveCount > nTotalCount ||
             (nInclusiveCount == nTotalCount && nExclusiveCount != 0))
    {
        bRet = FALSE;
        printf( "ERROR: GetFeatureCount() may not be taking attribute "
                "filter into account (nInclusiveCount = %d, nExclusiveCount = %d, nExclusiveCountWhileIterating = %d, nTotalCount = %d).\n",
                 nInclusiveCount, nExclusiveCount, nExclusiveCountWhileIterating, nTotalCount);
    }
    else if( bVerbose )
    {
        printf( "INFO: Attribute filter exclusion seems to work.\n" );
    }

    OGRFeature::DestroyFeature(poTargetFeature);

    return bRet;
}

/************************************************************************/
/*                         TestOGRLayerUTF8()                           */
/************************************************************************/

static int TestOGRLayerUTF8 ( OGRLayer *poLayer )
{
    int bRet = TRUE;

    poLayer->SetSpatialFilter( NULL );
    poLayer->SetAttributeFilter( NULL );
    poLayer->ResetReading();

    int bIsAdvertizedAsUTF8 = poLayer->TestCapability( OLCStringsAsUTF8 );
    int nFields = poLayer->GetLayerDefn()->GetFieldCount();
    int bFoundString = FALSE;
    int bFoundNonASCII = FALSE;
    int bFoundUTF8 = FALSE;
    int bCanAdvertizeUTF8 = TRUE;

    OGRFeature* poFeature = NULL;
    while( bRet && (poFeature = poLayer->GetNextFeature()) != NULL )
    {
        for(int i = 0; i<nFields; i++)
        {
            if (!poFeature->IsFieldSet(i))
                continue;
            if (poFeature->GetFieldDefnRef(i)->GetType() == OFTString)
            {
                const char* pszVal = poFeature->GetFieldAsString(i);
                if (pszVal[0] != 0)
                {
                    bFoundString = TRUE;
                    const GByte* pszIter = (const GByte*) pszVal;
                    int bIsASCII = TRUE;
                    while(*pszIter)
                    {
                        if (*pszIter >= 128)
                        {
                            bFoundNonASCII = TRUE;
                            bIsASCII = FALSE;
                            break;
                        }
                        pszIter ++;
                    }
                    int bIsUTF8 = CPLIsUTF8(pszVal, -1);
                    if (bIsUTF8 && !bIsASCII)
                        bFoundUTF8 = TRUE;
                    if (bIsAdvertizedAsUTF8)
                    {
                        if (!bIsUTF8)
                        {
                            printf( "ERROR: Found non-UTF8 content at field %d of feature %ld, but layer is advertized as UTF-8.\n",
                                    i, poFeature->GetFID() );
                            bRet = FALSE;
                            break;
                        }
                    }
                    else
                    {
                        if (!bIsUTF8)
                            bCanAdvertizeUTF8 = FALSE;
                    }
                }
            }
        }
        OGRFeature::DestroyFeature(poFeature);
    }

    if (!bFoundString)
    {
    }
    else if (bCanAdvertizeUTF8)
    {
        if (bIsAdvertizedAsUTF8)
        {
            if (bFoundUTF8)
            {
                printf( "INFO: Layer has UTF-8 content and is consistently declared as having UTF-8 content.\n" );
            }
            else if (!bFoundNonASCII)
            {
                printf( "INFO: Layer has ASCII only content and is consistently declared as having UTF-8 content.\n" );
            }
        }
        else
        {
            if (bFoundUTF8)
            {
                printf( "INFO: Layer could perhaps be advertized as UTF-8 compatible (and it has non-ASCII UTF-8 content).\n" );
            }
            else if (!bFoundNonASCII)
            {
                printf( "INFO: Layer could perhaps be advertized as UTF-8 compatible (it has only ASCII content).\n" );
            }
        }
    }
    else
    {
        printf( "INFO: Layer has non UTF-8 content (and is consistently declared as not being UTF-8 compatible).\n" );
    }

    return bRet;
}

/************************************************************************/
/*                         TestGetExtent()                              */
/************************************************************************/

static int TestGetExtent ( OGRLayer *poLayer )
{
    int bRet = TRUE;

    poLayer->SetSpatialFilter( NULL );
    poLayer->SetAttributeFilter( NULL );
    poLayer->ResetReading();

    OGREnvelope sExtent;
    OGREnvelope sExtentSlow;

    OGRErr eErr = poLayer->GetExtent(&sExtent, TRUE);
    OGRErr eErr2 = poLayer->OGRLayer::GetExtent(&sExtentSlow, TRUE);

    if (eErr != eErr2)
    {
        if (eErr == OGRERR_NONE && eErr2 != OGRERR_NONE)
        {
            /* with the LIBKML driver and test_ogrsf ../autotest/ogr/data/samples.kml "Styles and Markup" */
            printf("INFO: GetExtent() succeeded but OGRLayer::GetExtent() failed.\n");
        }
        else
        {
            bRet = FALSE;
            printf("ERROR: GetExtent() failed but OGRLayer::GetExtent() succeeded.\n");
        }
    }
    else if (eErr == OGRERR_NONE)
    {
        if (fabs(sExtentSlow.MinX - sExtent.MinX) < 1e-10 &&
            fabs(sExtentSlow.MinY - sExtent.MinY) < 1e-10 &&
            fabs(sExtentSlow.MaxX - sExtent.MaxX) < 1e-10 &&
            fabs(sExtentSlow.MaxY - sExtent.MaxY) < 1e-10)
        {
            printf("INFO: GetExtent() test passed.\n");
        }
        else
        {
            if (sExtentSlow.Contains(sExtent))
            {
                printf("INFO: sExtentSlow.Contains(sExtent)\n");
            }
            else if (sExtent.Contains(sExtentSlow))
            {
                printf("INFO: sExtent.Contains(sExtentSlow)\n");
            }
            else
            {
                printf("INFO: unknown relationship between sExtent and sExentSlow.\n");
            }
            printf("INFO: sExtentSlow.MinX = %.15f\n", sExtentSlow.MinX);
            printf("INFO: sExtentSlow.MinY = %.15f\n", sExtentSlow.MinY);
            printf("INFO: sExtentSlow.MaxX = %.15f\n", sExtentSlow.MaxX);
            printf("INFO: sExtentSlow.MaxY = %.15f\n", sExtentSlow.MaxY);
            printf("INFO: sExtent.MinX = %.15f\n", sExtent.MinX);
            printf("INFO: sExtent.MinY = %.15f\n", sExtent.MinY);
            printf("INFO: sExtent.MaxX = %.15f\n", sExtent.MaxX);
            printf("INFO: sExtent.MaxY = %.15f\n", sExtent.MaxY);
        }
    }

    return bRet;
}

/*************************************************************************/
/*             TestOGRLayerDeleteAndCreateFeature()                      */
/*                                                                       */
/*      Test delete feature by trying to delete the last feature and     */
/*      recreate it.                                                     */
/*************************************************************************/

static int TestOGRLayerDeleteAndCreateFeature( OGRLayer *poLayer )

{
    int bRet = TRUE;
    OGRFeature  * poFeature = NULL;
    OGRFeature  * poFeatureTest = NULL;
    long        nFID;

    poLayer->SetSpatialFilter( NULL );
    
    if( !poLayer->TestCapability( OLCRandomRead ) )
    {
        if( bVerbose )
            printf( "INFO: Skipping delete feature test since this layer "
                    "doesn't support random read.\n" );
        return bRet;
    }

    if( poLayer->GetFeatureCount() == 0 )
    {
        if( bVerbose )
            printf( "INFO: No feature available on layer '%s',"
                    "skipping delete/create feature test.\n",
                    poLayer->GetName() );
        
        return bRet;
    }
/* -------------------------------------------------------------------- */
/*      Fetch the last feature                                          */
/* -------------------------------------------------------------------- */
    poLayer->ResetReading();

    poLayer->SetNextByIndex(poLayer->GetFeatureCount() - 1);
    poFeature = poLayer->GetNextFeature();
    if (poFeature == NULL)
    {
        bRet = FALSE;
        printf( "ERROR: Could not get last feature of layer.\n" );
        goto end;
    }

/* -------------------------------------------------------------------- */
/*      Get the feature ID of the last feature                          */
/* -------------------------------------------------------------------- */
    nFID = poFeature->GetFID();

/* -------------------------------------------------------------------- */
/*      Delete the feature.                                             */
/* -------------------------------------------------------------------- */
    if( poLayer->DeleteFeature( nFID ) != OGRERR_NONE )
    {
        bRet = FALSE;
        printf( "ERROR: Attempt to DeleteFeature() failed.\n" );
        goto end;
    }
    
/* -------------------------------------------------------------------- */
/*      Now re-read the feature to verify the delete effect worked.     */
/* -------------------------------------------------------------------- */
    CPLPushErrorHandler(CPLQuietErrorHandler); /* silent legitimate error message */
    poFeatureTest = poLayer->GetFeature( nFID );
    CPLPopErrorHandler();
    if( poFeatureTest != NULL)
    {
        bRet = FALSE;
        printf( "ERROR: The feature was not deleted.\n" );
    }
    else
    {
        printf( "INFO: Delete Feature test passed.\n" );
    }
    OGRFeature::DestroyFeature(poFeatureTest);

/* -------------------------------------------------------------------- */
/*      Re-insert the features to restore to original state             */
/* -------------------------------------------------------------------- */
    if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
    {
        bRet = FALSE;
        printf( "ERROR: Attempt to restore feature failed.\n" );
    }

    if( poFeature->GetFID() != nFID )
    {
        /* Case of shapefile driver for example that will not try to */
        /* reuse the existing FID, but will assign a new one */
        printf( "INFO: Feature was created, but with not its original FID.\n" );
        nFID = poFeature->GetFID();
    }

/* -------------------------------------------------------------------- */
/*      Now re-read the feature to verify the create effect worked.     */
/* -------------------------------------------------------------------- */
    poFeatureTest = poLayer->GetFeature( nFID );
    if( poFeatureTest == NULL)
    {
        bRet = FALSE;
        printf( "ERROR: The feature was not created.\n" );
    }
    else
    {
        printf( "INFO: Create Feature test passed.\n" );
    }
    OGRFeature::DestroyFeature(poFeatureTest);
    
end:
/* -------------------------------------------------------------------- */
/*      Cleanup.                                                        */
/* -------------------------------------------------------------------- */

    OGRFeature::DestroyFeature(poFeature);

    return bRet;
}

/************************************************************************/
/*                            TestOGRLayer()                            */
/************************************************************************/

static int TestOGRLayer( OGRDataSource* poDS, OGRLayer * poLayer, int bIsSQLLayer )

{
    int bRet = TRUE;

/* -------------------------------------------------------------------- */
/*      Verify that there is no spatial filter in place by default.     */
/* -------------------------------------------------------------------- */
    if( poLayer->GetSpatialFilter() != NULL )
    {
        printf( "WARN: Spatial filter in place by default on layer %s.\n",
                poLayer->GetName() );
        poLayer->SetSpatialFilter( NULL );
    }

/* -------------------------------------------------------------------- */
/*      Test feature count accuracy.                                    */
/* -------------------------------------------------------------------- */
    bRet &= TestOGRLayerFeatureCount( poDS, poLayer, bIsSQLLayer );

/* -------------------------------------------------------------------- */
/*      Test spatial filtering                                          */
/* -------------------------------------------------------------------- */
    bRet &= TestSpatialFilter( poLayer );

/* -------------------------------------------------------------------- */
/*      Test attribute filtering                                        */
/* -------------------------------------------------------------------- */
    bRet &= TestAttributeFilter( poDS, poLayer );

/* -------------------------------------------------------------------- */
/*      Test GetExtent()                                                */
/* -------------------------------------------------------------------- */
    bRet &= TestGetExtent( poLayer );

/* -------------------------------------------------------------------- */
/*      Test random reading.                                            */
/* -------------------------------------------------------------------- */
    if( poLayer->TestCapability( OLCRandomRead ) )
    {
        bRet &= TestOGRLayerRandomRead( poLayer );
    }
    
/* -------------------------------------------------------------------- */
/*      Test SetNextByIndex.                                            */
/* -------------------------------------------------------------------- */
    if( poLayer->TestCapability( OLCFastSetNextByIndex ) )
    {
        bRet &= TestOGRLayerSetNextByIndex( poLayer );
    }
    
/* -------------------------------------------------------------------- */
/*      Test delete feature.                                            */
/* -------------------------------------------------------------------- */
    if( poLayer->TestCapability( OLCDeleteFeature ) )
    {
        bRet &= TestOGRLayerDeleteAndCreateFeature( poLayer );
    }
    
/* -------------------------------------------------------------------- */
/*      Test random writing.                                            */
/* -------------------------------------------------------------------- */
    if( poLayer->TestCapability( OLCRandomWrite ) )
    {
        bRet &= TestOGRLayerRandomWrite( poLayer );
    }

    bRet &= TestOGRLayerUTF8( poLayer );

    return bRet;
}

/************************************************************************/
/*                        TestInterleavedReading()                      */
/************************************************************************/

static int TestInterleavedReading( const char* pszDataSource, char** papszLayers )
{
    int bRet = TRUE;
    OGRDataSource* poDS = NULL;
    OGRDataSource* poDS2 = NULL;
    OGRLayer* poLayer1 = NULL;
    OGRLayer* poLayer2 = NULL;
    OGRFeature* poFeature11_Ref = NULL;
    OGRFeature* poFeature12_Ref = NULL;
    OGRFeature* poFeature21_Ref = NULL;
    OGRFeature* poFeature22_Ref = NULL;
    OGRFeature* poFeature11 = NULL;
    OGRFeature* poFeature12 = NULL;
    OGRFeature* poFeature21 = NULL;
    OGRFeature* poFeature22 = NULL;

    /* Check that we have 2 layers with at least 2 features */
    poDS = OGRSFDriverRegistrar::Open( pszDataSource, FALSE, NULL );
    if (poDS == NULL)
    {
        printf( "INFO: Skipping TestInterleavedReading(). Cannot reopen datasource\n" );
        goto bye;
    }

    poLayer1 = papszLayers ? poDS->GetLayerByName(papszLayers[0]) : poDS->GetLayer(0);
    poLayer2 = papszLayers ? poDS->GetLayerByName(papszLayers[1]) : poDS->GetLayer(1);
    if (poLayer1 == NULL || poLayer2 == NULL ||
        poLayer1->GetFeatureCount() < 2 || poLayer2->GetFeatureCount() < 2)
    {
        printf( "INFO: Skipping TestInterleavedReading(). Test conditions are not met\n" );
        goto bye;
    }

    /* Test normal reading */
    OGRDataSource::DestroyDataSource(poDS);
    poDS = OGRSFDriverRegistrar::Open( pszDataSource, FALSE, NULL );
    poDS2 = OGRSFDriverRegistrar::Open( pszDataSource, FALSE, NULL );
    if (poDS == NULL || poDS2 == NULL)
    {
        printf( "INFO: Skipping TestInterleavedReading(). Cannot reopen datasource\n" );
        goto bye;
    }

    poLayer1 = papszLayers ? poDS->GetLayerByName(papszLayers[0]) : poDS->GetLayer(0);
    poLayer2 = papszLayers ? poDS->GetLayerByName(papszLayers[1]) : poDS->GetLayer(1);
    if (poLayer1 == NULL || poLayer2 == NULL)
    {
        printf( "ERROR: Skipping TestInterleavedReading(). Test conditions are not met\n" );
        bRet = FALSE;
        goto bye;
    }

    poFeature11_Ref = poLayer1->GetNextFeature();
    poFeature12_Ref = poLayer1->GetNextFeature();
    poFeature21_Ref = poLayer2->GetNextFeature();
    poFeature22_Ref = poLayer2->GetNextFeature();
    if (poFeature11_Ref == NULL || poFeature12_Ref == NULL || poFeature21_Ref == NULL || poFeature22_Ref == NULL)
    {
        printf( "ERROR: TestInterleavedReading() failed: poFeature11_Ref=%p, poFeature12_Ref=%p, poFeature21_Ref=%p, poFeature22_Ref=%p\n",
                poFeature11_Ref, poFeature12_Ref, poFeature21_Ref, poFeature22_Ref);
        bRet = FALSE;
        goto bye;
    }

    /* Test interleaved reading */
    poLayer1 = papszLayers ? poDS2->GetLayerByName(papszLayers[0]) : poDS2->GetLayer(0);
    poLayer2 = papszLayers ? poDS2->GetLayerByName(papszLayers[1]) : poDS2->GetLayer(1);
    if (poLayer1 == NULL || poLayer2 == NULL)
    {
        printf( "ERROR: Skipping TestInterleavedReading(). Test conditions are not met\n" );
        bRet = FALSE;
        goto bye;
    }

    poFeature11 = poLayer1->GetNextFeature();
    poFeature21 = poLayer2->GetNextFeature();
    poFeature12 = poLayer1->GetNextFeature();
    poFeature22 = poLayer2->GetNextFeature();

    if (poFeature11 == NULL || poFeature21 == NULL || poFeature12 == NULL || poFeature22 == NULL)
    {
        printf( "ERROR: TestInterleavedReading() failed: poFeature11=%p, poFeature21=%p, poFeature12=%p, poFeature22=%p\n",
                poFeature11, poFeature21, poFeature12, poFeature22);
        bRet = FALSE;
        goto bye;
    }

    if (poFeature12->Equal(poFeature11))
    {
        printf( "WARN: TestInterleavedReading() failed: poFeature12 == poFeature11. "
                "The datasource resets the layer reading when interleaved layer reading pattern is detected. Acceptable but could be improved\n" );
        goto bye;
    }

    /* We cannot directly compare the feature as they don't share */
    /* the same (pointer) layer definition, so just compare FIDs */
    if (poFeature12_Ref->GetFID() != poFeature12->GetFID())
    {
        printf( "ERROR: TestInterleavedReading() failed: poFeature12_Ref != poFeature12\n" );
        poFeature12_Ref->DumpReadable(stdout, NULL);
        poFeature12->DumpReadable(stdout, NULL);
        bRet = FALSE;
        goto bye;
    }

    if( bVerbose )
    {
        printf("INFO: TestInterleavedReading() successfull.\n");
    }

bye:
    OGRFeature::DestroyFeature(poFeature11_Ref);
    OGRFeature::DestroyFeature(poFeature12_Ref);
    OGRFeature::DestroyFeature(poFeature21_Ref);
    OGRFeature::DestroyFeature(poFeature22_Ref);
    OGRFeature::DestroyFeature(poFeature11);
    OGRFeature::DestroyFeature(poFeature21);
    OGRFeature::DestroyFeature(poFeature12);
    OGRFeature::DestroyFeature(poFeature22);
    OGRDataSource::DestroyDataSource(poDS);
    OGRDataSource::DestroyDataSource(poDS2);
    return bRet;
}
