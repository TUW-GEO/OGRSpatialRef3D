/******************************************************************************
 *
 * Project:  OpenGIS Simple Features Reference Implementation
 * Purpose:  Classes for manipulating spatial reference systems with 
 *           vetical datum suppurt in a platform non-specific manner.
 * Authors:  Gottfried Mandlburger, Johannes Otepka, Bhargav patel, Peb Ruswono Aryan
 *
 ******************************************************************************
 * Copyright (c) 2012,  I.P.F., TU Vienna.
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
#include <iostream>

#include "ogr_spatialref3D.h"

#define RAD_TO_DEG	57.29577951308232
#define DEG_TO_RAD	.0174532925199432958
#define EPSILON 1e-5

using namespace std;
/************************************************************************/
/*                         OGRSpatialReference3D                        */
/************************************************************************/

/**
 * This class respresents a OpenGIS 3D Spatial Reference System. 
 * It is derived from OGRSpatialReference class and adds elements to 
 * define the Reference System's Vertical Datum. The class provides:
 *
 *  I  ) a (quasi) geoid raster model for transformation between 
 *       ellipsoidal and (quasi) orthometric heights
 *  II ) an additional height correction raster model to compensate the
 *       errors of historic heigth system realizations (e.g. MGI, DHHN...)
 *  III) a constant offset to deal with local height origins, like tidal 
 *       height systems, local height definitions (e.g. Wiener Null) ...
 *  IV ) a scaling factor (e.g. foot-meter-conversion)
 * 
/************************************************************************/


OGRSpatialReference3D::OGRSpatialReference3D()
{
	GDALAllRegister();
	
	bHasGeoid = false;
	bHasVCorr = false;

	poGeoid = NULL;
	poVCorr = NULL;

	dfVOffset_ = 0.0;
	dfVScale_ = 0.0;

	is_debug = false;
	dbg_geoid = NULL;
	dbg_vcorr = NULL;
}

OGRSpatialReference3D::~OGRSpatialReference3D()
{
}

OGRErr      
	OGRSpatialReference3D::importFromWkt3D( char ** pszWKT)
{
	if(OGRERR_NONE != importFromWkt(pszWKT))
	{
		//handle for import error
		std::cerr << "Error import WKT3D : " << std::endl;
	}

	const OGR_SRSNode *poNode = GetAttrNode( "GEOID" );
	if (poNode != NULL){
		const char *pszData = poNode->GetChild(1)->GetChild(0)->GetValue();
		if(OGRERR_NONE != SetGeoidModel(pszData))
		{
			//handle for loading error
			std::cerr << "Error Loading GEOID : " << pszData << std::endl;
			throw OGRERR_FAILURE;
		}
	}

	poNode = GetAttrNode( "VCORR" );
	if (poNode != NULL){
		const char *pszData = poNode->GetChild(1)->GetChild(0)->GetValue();
		if(OGRERR_NONE != SetVCorrModel(pszData))
		{
			//handle for loading error
			std::cerr << "Error Loading VCORR : " << pszData << std::endl;
			throw OGRERR_FAILURE;
		}
	}

	//currently these parameters are unused, deferred for future implementations
	poNode = GetAttrNode( "VSHIFT" );
	if (poNode != NULL){
		const char *pszData = poNode->GetChild(0)->GetValue();
		SetVOffset(atof(pszData));
	}
	
	poNode = GetAttrNode( "VSCALE" );
	if (poNode != NULL){
		const char *pszData = poNode->GetChild(0)->GetValue();
		std::cout << "Loading VSCALE : " << pszData << std::endl;
		SetVScale(atof(pszData));
	}
	return OGRERR_NONE;
}

OGRErr OGRSpatialReference3D::SetVOffset( double  dfVOffset )
{
	dfVOffset_ = dfVOffset;
	return OGRERR_NONE;
}

double OGRSpatialReference3D::GetVOffset ()
{
	return dfVOffset_;
}


OGRErr OGRSpatialReference3D::SetVScale( double  dfVScale )
{
	dfVScale_ = dfVScale;
	return OGRERR_NONE;
}

double OGRSpatialReference3D::GetVScale ()
{
	return dfVScale_;
}

OGRErr OGRSpatialReference3D::SetGeoidModel( const char * pszGeoidModel )
{
	if(poGeoid != NULL)
		delete poGeoid;

	poGeoid = new RasterResampler();
	poGeoid->Open( pszGeoidModel );
	//poGeoid = (GDALDataset *) GDALOpen( pszGeoidModel, GA_ReadOnly );
	
	if( poGeoid == NULL )
    {
        printf("gdal failed - unable to open '%s'.\n",
                 pszGeoidModel );
		return OGRERR_FAILURE;
	}
	bHasGeoid = true;
	return OGRERR_NONE;
}

OGRErr OGRSpatialReference3D::SetVCorrModel( const char * pszVCorrModel )
{
	if(poVCorr != NULL)
		delete poVCorr;

	poVCorr = new RasterResampler();
	poVCorr->Open( pszVCorrModel );

	if( poVCorr == NULL )
    {
        std::cerr << "gdal failed - unable to open '" << pszVCorrModel << "'.\n";
		return OGRERR_FAILURE;
	}
	bHasVCorr = true;
	return OGRERR_NONE;
}

OGRErr OGRSpatialReference3D::ApplyVerticalCorrection(int is_inverse, unsigned int point_count, double *x, double *y, double *z)
{
	double* dZCorr = (double*)CPLMalloc(sizeof(double)*point_count);
	double* dZTemp = (double*)CPLMalloc(sizeof(double)*point_count);

	for(unsigned int i=0; i<point_count; ++i){ 
		dZCorr[i] = dfVOffset_;
		dZTemp[i] = 0.0;
	}

	if(HasGeoidModel()){
		poGeoid->GetValueAt(point_count, x, y, dZTemp);
		for(unsigned int i=0; i<point_count; ++i){
			dZCorr[i] += dZTemp[i];

			if(is_debug && dbg_geoid != NULL)
				dbg_geoid[i] = dZTemp[i];
		}
	}

	if(HasVCorrModel()){
		poVCorr->GetValueAt(point_count, x, y, dZTemp);
		for(unsigned int i=0; i<point_count; ++i){
			dZCorr[i] += dZTemp[i];
			
			if(is_debug && dbg_vcorr != NULL)
				dbg_vcorr[i] = dZTemp[i];
		}
	}

	for(unsigned int i=0; i<point_count; ++i)
	{
		if(is_inverse)
			z[i] -= dZCorr[i];
		else
			z[i] += dZCorr[i];
	}

	CPLFree(dZCorr);
	CPLFree(dZTemp);

	return OGRERR_NONE;
}

double OGRSpatialReference3D::GetValueAt(GDALDataset* hDataset, double x, double y)
{
	/*
	 * Deprecated (no longer in use)
	 */
	cout << "get raster size" << endl;
	int nXsize = hDataset->GetRasterXSize();
    int nYsize = hDataset->GetRasterYSize();
	cout << " " << nXsize << " " << nYsize << endl;

	cout << "get geo transform" << endl;
    double geotrans[6];
    hDataset->GetGeoTransform(geotrans);

	cout << "invert geo transform" << endl;
	double inv_geotrans[6];
    if( GDALInvGeoTransform( geotrans, inv_geotrans ) == 0 )
      throw std::exception( "inversion of geo transformation failed." );

	//Calculate Raster coordinate
	cout << "calc raster coord" << endl;
	double dPixel = inv_geotrans[0]+ 
          + x * inv_geotrans[1] * RAD_TO_DEG
          + y * inv_geotrans[2] * RAD_TO_DEG;

    double dLine = inv_geotrans[3]
          + x * inv_geotrans[4] * RAD_TO_DEG
          + y * inv_geotrans[5] * RAD_TO_DEG;

	cout << " pixel,line : " << dPixel << " " << dLine << endl;
	int px = (int)floor(dPixel);
	int py = (int)floor(dLine);
	cout << " " << px << " " << py << endl;

	// Boundary checking
	if (px < 0 || py < 0 || px >= nXsize || py >= nYsize)
	{
		cout << "invalid raster coordinate";
	}
	else{
		// raster coordinate relative to nearest topleft pixel for weighting
		double dx = dPixel - px;
		double dy = dLine - py;

		// ASSUMPTION : single band raster 
		double dNoDataValue = hDataset->GetRasterBand(1)->GetNoDataValue();
		cout << " no data:  " << dNoDataValue << endl;

		//------------------------------------------------------->8 cut here

		// allocate array for height data
        double *padH = (double *) CPLMalloc(sizeof(double)*4);

		// get 2x2 pixel neighborhood
		hDataset->RasterIO( GF_Read, px, py, 2, 2, 
                          padH, 2, 2, GDT_Float64, 
                          1, NULL, 0, 0, 0 );

		cout << " "<< padH[0] << " " << padH[1] << " " << padH[2] << " "<< padH[3]<< endl;

		double dWeight = 0.0;
		double dSum = 0.0;
		double dNorm = 0.0;

		// Calculate contribution of each neighboring pixel
		// weighted by inverse rectangular area
		//
		// tl----+------tr
		// |     |      |
		// |    dy      |
		// |     |      |
		// +--dx-+------+
		// |     |      |
		// bl----+------br
		//
		// possible extension: code below can be refactored
		// by moving to separate function for computing
		// specific interpolation scheme (e.g. kriging)
		
		if (fabs(padH[0] - dNoDataValue) > EPSILON)
		{
			// top left neighbor
			dWeight = (1.0 - dx) * (1.0 - dy);
			dSum += padH[0] * dWeight;
			dNorm += dWeight;
		}
		
		if (fabs(padH[1] - dNoDataValue) > EPSILON)
		{
			// top right neighbor
			dWeight = dx * (1.0 - dy);
			dSum += padH[1] * dWeight;
			dNorm += dWeight;
		}
		
		if (fabs(padH[2] - dNoDataValue) > EPSILON)
		{
			// bottom left neighbor
			dWeight = (1.0 - dx) * (dy);
			dSum += padH[2] * dWeight;
			dNorm += dWeight;
		}
		
		if (fabs(padH[3] - dNoDataValue) > EPSILON)
		{
			// bottom right neighbor
			dWeight = (dx) * (dy);
			dSum += padH[3] * dWeight;
			dNorm += dWeight;
		}

		cout << "(Unnormalized) H_interp " << dSum << endl;

		if(dNorm < EPSILON)  // No valid data available
			dSum = 0.0; 
		else if(dNorm < 1.0) // One or more pixels is no data
			dSum /= dNorm;

		cout << "Normalization factor " << dNorm << endl;
		cout << "(Normalized) H_interp " << dSum << endl;
		
		CPLFree(padH);

		//---------------------------------------------------------------->8 cut here

		return dSum;
	}

	return 0.0;
}

bool OGRSpatialReference3D::HasGeoidModel ()
{
	return bHasGeoid;
}

bool OGRSpatialReference3D::HasVCorrModel ()
{
	return bHasVCorr;
}

bool OGRSpatialReference3D::HasVerticalModel()
{
	return HasGeoidModel() || HasVCorrModel();
}

void OGRSpatialReference3D::SetDebug(bool debug_mode)
{
	is_debug = debug_mode;
}

void OGRSpatialReference3D::SetDebugData(double* geoid_undulation, double* vert_correction)
{
	dbg_geoid = geoid_undulation;
	dbg_vcorr = vert_correction;
}


CPL_DLL OGRCoordinateTransformation3D::OGRCoordinateTransformation3D()
{

}


