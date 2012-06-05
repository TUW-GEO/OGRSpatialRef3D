/******************************************************************************
 * $Id: gdaltindex.c 19806 2010-06-06 11:56:19Z rouault $
 *
 * Project:  MapServer
 * Purpose:  Commandline App to build tile index for raster files.
 * Author:   Frank Warmerdam, warmerdam@pobox.com
 *
 ******************************************************************************
 * Copyright (c) 2001, Frank Warmerdam, DM Solutions Group Inc
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

#include "ogr_api.h"
#include "ogr_srs_api.h"
#include "gdal.h"
#include "cpl_port.h"
#include "cpl_conv.h"
#include "cpl_string.h"

CPL_CVSID("$Id: gdaltindex.c 19806 2010-06-06 11:56:19Z rouault $");

/************************************************************************/
/*                               Usage()                                */
/************************************************************************/

static void Usage()

{
    fprintf(stdout, "%s", 
            "\n"
            "Usage: gdaltindex [-tileindex field_name] [-write_absolute_path] \n"
            "                  [-skip_different_projection] index_file [gdal_file]*\n"
            "\n"
            "eg.\n"
            "  % gdaltindex doq_index.shp doq/*.tif\n" 
            "\n"
            "NOTES:\n"
            "  o The shapefile (index_file) will be created if it doesn't already exist.\n" 
            "  o The default tile index field is 'location'.\n"
            "  o Raster filenames will be put in the file exactly as they are specified\n"
            "    on the commandline unless the option -write_absolute_path is used.\n"
            "  o If -skip_different_projection is specified, only files with same projection ref\n"
            "    as files already inserted in the tileindex will be inserted.\n"
            "  o Simple rectangular polygons are generated in the same\n"
            "    coordinate system as the rasters.\n" );
    exit(1);
}

/************************************************************************/
/*                                main()                                */
/************************************************************************/

int main(int argc, char *argv[])
{
    const char *index_filename = NULL;
    const char *tile_index = "location";
    int		i_arg, ti_field;
    OGRDataSourceH hTileIndexDS;
    OGRLayerH hLayer = NULL;
    OGRFeatureDefnH hFDefn;
    int write_absolute_path = FALSE;
    char* current_path = NULL;
    int i;
    int nExistingFiles;
    int skip_different_projection = FALSE;
    char** existingFilesTab = NULL;
    int alreadyExistingProjectionRefValid = FALSE;
    char* alreadyExistingProjectionRef = NULL;
    char* index_filename_mod;
    int bExists;
    VSIStatBuf sStatBuf;

    /* Check that we are running against at least GDAL 1.4 */
    /* Note to developers : if we use newer API, please change the requirement */
    if (atoi(GDALVersionInfo("VERSION_NUM")) < 1400)
    {
        fprintf(stderr, "At least, GDAL >= 1.4.0 is required for this version of %s, "
                "which was compiled against GDAL %s\n", argv[0], GDAL_RELEASE_NAME);
        exit(1);
    }

    GDALAllRegister();
    OGRRegisterAll();

    argc = GDALGeneralCmdLineProcessor( argc, &argv, 0 );
    if( argc < 1 )
        exit( -argc );

/* -------------------------------------------------------------------- */
/*      Get commandline arguments other than the GDAL raster filenames. */
/* -------------------------------------------------------------------- */
    for( i_arg = 1; i_arg < argc; i_arg++ )
    {
        if( EQUAL(argv[i_arg], "--utility_version") )
        {
            printf("%s was compiled against GDAL %s and is running against GDAL %s\n",
                   argv[0], GDAL_RELEASE_NAME, GDALVersionInfo("RELEASE_NAME"));
            return 0;
        }
        else if( strcmp(argv[i_arg],"-tileindex") == 0 )
        {
            tile_index = argv[++i_arg];
        }
        else if ( strcmp(argv[i_arg],"-write_absolute_path") == 0 )
        {
            write_absolute_path = TRUE;
        }
        else if ( strcmp(argv[i_arg],"-skip_different_projection") == 0 )
        {
            skip_different_projection = TRUE;
        }
        else if( argv[i_arg][0] == '-' )
            Usage();
        else if( index_filename == NULL )
        {
            index_filename = argv[i_arg];
            i_arg++;
            break;
        }
    }

    if( index_filename == NULL || i_arg == argc )
        Usage();

/* -------------------------------------------------------------------- */
/*      Open or create the target shapefile and DBF file.               */
/* -------------------------------------------------------------------- */
    index_filename_mod = CPLStrdup(CPLResetExtension(index_filename, "shp"));

    bExists = (VSIStat(index_filename_mod, &sStatBuf) == 0);
    if (!bExists)
    {
        CPLFree(index_filename_mod);
        index_filename_mod = CPLStrdup(CPLResetExtension(index_filename, "SHP"));
        bExists = (VSIStat(index_filename_mod, &sStatBuf) == 0);
    }
    CPLFree(index_filename_mod);

    if (bExists)
    {
        hTileIndexDS = OGROpen( index_filename, TRUE, NULL );
        if (hTileIndexDS != NULL)
        {
            hLayer = OGR_DS_GetLayer(hTileIndexDS, 0);
        }
    }
    else
    {
        OGRSFDriverH hDriver;
        const char* pszDriverName = "ESRI Shapefile";

        printf( "Creating new index file...\n" );
        hDriver = OGRGetDriverByName( pszDriverName );
        if( hDriver == NULL )
        {
            printf( "%s driver not available.\n", pszDriverName );
            exit( 1 );
        }

        hTileIndexDS = OGR_Dr_CreateDataSource( hDriver, index_filename, NULL );
        if (hTileIndexDS)
        {
            char* pszLayerName = CPLStrdup(CPLGetBasename(index_filename));
            OGRSpatialReferenceH hSpatialRef = NULL;
            GDALDatasetH hDS = GDALOpen( argv[i_arg], GA_ReadOnly );
            if (hDS)
            {
                const char* pszWKT = GDALGetProjectionRef(hDS);
                if (pszWKT != NULL && pszWKT[0] != '\0')
                {
                    hSpatialRef = OSRNewSpatialReference(pszWKT);
                }
                GDALClose(hDS);
            }

            hLayer = OGR_DS_CreateLayer( hTileIndexDS, pszLayerName, hSpatialRef, wkbPolygon, NULL );
            CPLFree(pszLayerName);
            if (hSpatialRef)
                OSRRelease(hSpatialRef);

            if (hLayer)
            {
                OGRFieldDefnH hFieldDefn = OGR_Fld_Create( tile_index, OFTString );
                OGR_Fld_SetWidth( hFieldDefn, 255);
                OGR_L_CreateField( hLayer, hFieldDefn, TRUE );
                OGR_Fld_Destroy(hFieldDefn);
            }
        }
    }

    if( hTileIndexDS == NULL || hLayer == NULL )
    {
        fprintf( stderr, "Unable to open/create shapefile `%s'.\n", 
                 index_filename );
        exit(2);
    }

    hFDefn = OGR_L_GetLayerDefn(hLayer);

    for( ti_field = 0; ti_field < OGR_FD_GetFieldCount(hFDefn); ti_field++ )
    {
        OGRFieldDefnH hFieldDefn = OGR_FD_GetFieldDefn( hFDefn, ti_field );
        if( strcmp(OGR_Fld_GetNameRef(hFieldDefn), tile_index) == 0 )
            break;
    }

    if( ti_field == OGR_FD_GetFieldCount(hFDefn) )
    {
        fprintf( stderr, "Unable to find field `%s' in DBF file `%s'.\n", 
                 tile_index, index_filename );
        exit(2);
    }

    /* Load in memory existing file names in SHP */
    nExistingFiles = OGR_L_GetFeatureCount(hLayer, FALSE);
    if (nExistingFiles)
    {
        OGRFeatureH hFeature;
        existingFilesTab = (char**)CPLMalloc(nExistingFiles * sizeof(char*));
        for(i=0;i<nExistingFiles;i++)
        {
            hFeature = OGR_L_GetNextFeature(hLayer);
            existingFilesTab[i] = CPLStrdup(OGR_F_GetFieldAsString( hFeature, ti_field ));
            if (i == 0)
            {
                GDALDatasetH hDS = GDALOpen(existingFilesTab[i], GA_ReadOnly );
                if (hDS)
                {
                    alreadyExistingProjectionRefValid = TRUE;
                    alreadyExistingProjectionRef = CPLStrdup(GDALGetProjectionRef(hDS));
                    GDALClose(hDS);
                }
            }
            OGR_F_Destroy( hFeature );
        }
    }

    if (write_absolute_path)
    {
        current_path = CPLGetCurrentDir();
        if (current_path == NULL)
        {
            fprintf( stderr, "This system does not support the CPLGetCurrentDir call. "
                             "The option -write_absolute_path will have no effect\n");
            write_absolute_path = FALSE;
        }
    }

/* -------------------------------------------------------------------- */
/*      loop over GDAL files, processing.                               */
/* -------------------------------------------------------------------- */
    for( ; i_arg < argc; i_arg++ )
    {
        GDALDatasetH	hDS;
        double	        adfGeoTransform[6];
        double		adfX[5], adfY[5];
        int		nXSize, nYSize;
        char* fileNameToWrite;
        const char* projectionRef;
        VSIStatBuf sStatBuf;
        int k;
        OGRFeatureH hFeature;
        OGRGeometryH hPoly, hRing;

        /* Make sure it is a file before building absolute path name */
        if (write_absolute_path && CPLIsFilenameRelative( argv[i_arg] ) &&
            VSIStat( argv[i_arg], &sStatBuf ) == 0)
        {
            fileNameToWrite = CPLStrdup(CPLProjectRelativeFilename(current_path, argv[i_arg]));
        }
        else
        {
            fileNameToWrite = CPLStrdup(argv[i_arg]);
        }

        /* Checks that file is not already in tileindex */
        for(i=0;i<nExistingFiles;i++)
        {
            if (EQUAL(fileNameToWrite, existingFilesTab[i]))
            {
                fprintf(stderr, "File %s is already in tileindex. Skipping it.\n",
                        fileNameToWrite);
                break;
            }
        }
        if (i != nExistingFiles)
        {
            CPLFree(fileNameToWrite);
            continue;
        }

        hDS = GDALOpen( argv[i_arg], GA_ReadOnly );
        if( hDS == NULL )
        {
            fprintf( stderr, "Unable to open %s, skipping.\n", 
                     argv[i_arg] );
            CPLFree(fileNameToWrite);
            continue;
        }

        GDALGetGeoTransform( hDS, adfGeoTransform );
        if( adfGeoTransform[0] == 0.0 
            && adfGeoTransform[1] == 1.0
            && adfGeoTransform[3] == 0.0
            && ABS(adfGeoTransform[5]) == 1.0 )
        {
            fprintf( stderr, 
                     "It appears no georeferencing is available for\n"
                     "`%s', skipping.\n", 
                     argv[i_arg] );
            GDALClose( hDS );
            CPLFree(fileNameToWrite);
            continue;
        }

        projectionRef = GDALGetProjectionRef(hDS);
        if (alreadyExistingProjectionRefValid)
        {
            int projectionRefNotNull, alreadyExistingProjectionRefNotNull;
            projectionRefNotNull = projectionRef && projectionRef[0];
            alreadyExistingProjectionRefNotNull = alreadyExistingProjectionRef && alreadyExistingProjectionRef[0];
            if ((projectionRefNotNull &&
                 alreadyExistingProjectionRefNotNull &&
                 EQUAL(projectionRef, alreadyExistingProjectionRef) == 0) ||
                (projectionRefNotNull != alreadyExistingProjectionRefNotNull))
            {
                fprintf(stderr, "Warning : %s is not using the same projection system as "
                                "other files in the tileindex. This may cause problems when "
                                "using it in MapServer for example.%s\n", argv[i_arg],
                                (skip_different_projection) ? " Skipping it" : "");
                if (skip_different_projection)
                {
                    CPLFree(fileNameToWrite);
                    GDALClose( hDS );
                    continue;
                }
            }
        }
        else
        {
            alreadyExistingProjectionRefValid = TRUE;
            alreadyExistingProjectionRef = CPLStrdup(projectionRef);
        }

        nXSize = GDALGetRasterXSize( hDS );
        nYSize = GDALGetRasterYSize( hDS );
        
        adfX[0] = adfGeoTransform[0] 
            + 0 * adfGeoTransform[1] 
            + 0 * adfGeoTransform[2];
        adfY[0] = adfGeoTransform[3] 
            + 0 * adfGeoTransform[4] 
            + 0 * adfGeoTransform[5];
        
        adfX[1] = adfGeoTransform[0] 
            + nXSize * adfGeoTransform[1] 
            + 0 * adfGeoTransform[2];
        adfY[1] = adfGeoTransform[3] 
            + nXSize * adfGeoTransform[4] 
            + 0 * adfGeoTransform[5];
        
        adfX[2] = adfGeoTransform[0] 
            + nXSize * adfGeoTransform[1] 
            + nYSize * adfGeoTransform[2];
        adfY[2] = adfGeoTransform[3] 
            + nXSize * adfGeoTransform[4] 
            + nYSize * adfGeoTransform[5];
        
        adfX[3] = adfGeoTransform[0] 
            + 0 * adfGeoTransform[1] 
            + nYSize * adfGeoTransform[2];
        adfY[3] = adfGeoTransform[3] 
            + 0 * adfGeoTransform[4] 
            + nYSize * adfGeoTransform[5];
        
        adfX[4] = adfGeoTransform[0] 
            + 0 * adfGeoTransform[1] 
            + 0 * adfGeoTransform[2];
        adfY[4] = adfGeoTransform[3] 
            + 0 * adfGeoTransform[4] 
            + 0 * adfGeoTransform[5];

        hFeature = OGR_F_Create( OGR_L_GetLayerDefn( hLayer ) );
        OGR_F_SetFieldString( hFeature, ti_field, fileNameToWrite );

        hPoly = OGR_G_CreateGeometry(wkbPolygon);
        hRing = OGR_G_CreateGeometry(wkbLinearRing);
        for(k=0;k<5;k++)
            OGR_G_SetPoint_2D(hRing, k, adfX[k], adfY[k]);
        OGR_G_AddGeometryDirectly( hPoly, hRing );
        OGR_F_SetGeometryDirectly( hFeature, hPoly );

        if( OGR_L_CreateFeature( hLayer, hFeature ) != OGRERR_NONE )
        {
           printf( "Failed to create feature in shapefile.\n" );
           break;
        }

        OGR_F_Destroy( hFeature );

        
        CPLFree(fileNameToWrite);

        GDALClose( hDS );
    }
    
    CPLFree(current_path);
    
    if (nExistingFiles)
    {
        for(i=0;i<nExistingFiles;i++)
        {
            CPLFree(existingFilesTab[i]);
        }
        CPLFree(existingFilesTab);
    }
    CPLFree(alreadyExistingProjectionRef);

    OGR_DS_Destroy( hTileIndexDS );
    
    GDALDestroyDriverManager();
    OGRCleanupAll();
    CSLDestroy(argv);
    
    exit( 0 );
} 
