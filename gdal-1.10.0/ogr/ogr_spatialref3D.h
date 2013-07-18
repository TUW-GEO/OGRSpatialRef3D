/******************************************************************************
 *
 * Project:  OpenGIS Simple Features Reference Implementation
 * Purpose:  Classes for manipulating spatial reference systems with
 *           vetical datum suppurt in a platform non-specific manner.
 * Authors:  Gottfried Mandlburger, Johannes Otepka, Bhargav patel
 *
 ******************************************************************************
 * Copyright (c) 2012,  I.P.F., TU Vieanna.
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
#ifndef _OGR_SPATIALREF3D_H_INCLUDED
#define _OGR_SPATIALREF3D_H_INCLUDED


#include "gdal_priv.h"
#include "ogr_spatialref.h"
#include "gdalwarper.h"				//needed for GDALResampleAlg

#include <vector>

//class GDALDataset;

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


class Point2D
{
	double x,y,z;
public:
	Point2D()
	{
		x=0;
		y=0;
		z=0;
	}
	void setxy(int x1, int y1)
	{
		x=x1;
		y=y1;
	}
	void setz(double z1)
	{
		z=z1;
	}
	double X()
	{
		return x;
	}
	double Y()
	{
		return y;
	}
	double Z()
	{
		return z;
	}
};

class CPL_DLL OGRSpatialReference3D:public OGRSpatialReference
{
	bool bHasGeoid;
	bool bHasVCorr;

	double dfVOffset_;
	double dfVScale_;

	GDALDataset  *poGeoid;
	GDALDataset  *poVCorr;
public:
	OGRSpatialReference3D();
	OGRSpatialReference3D(const char * pszWKT,
                          const char * pszGeoidModel,
                          const char * pszVCorrModel,
                          double dfVOffset,
                          double dfVScale);
	 virtual    ~OGRSpatialReference3D();

 	OGRErr SetGeoidModel( const char * pszGeoidModel );
    const char * GetGeoidModel ();

    OGRErr SetVCorrModel( const char * pszVCorrModel );
    const char * GetVCorrModel ();

    OGRErr SetVOffset( double  dfVOffset );
    double GetVOffset ();

    OGRErr SetVScale( double  dfVScale );
    double GetVScale ();

	bool HasGeoidModel();
	bool HasVCorrModel();

	OGRErr ApplyVerticalCorrection(int is_inverse, unsigned int point_count, double *x, double *y, double *z);

protected:
	// This interpolatez_generalize function will take multiple point as an argument for a find out interpolation
	// Also it can find out interolation at multiple band so it can use for R-G-B image
	// GDALDataset is pointer of raster, IRC_mask is for masking purpose, Irc_Band are band of raster
	//Irc_pt are point at where you want to find interpolate value it can be single or group of points
	//Xrc_z give result of interpolation in a coloumn. Result store in a coloumn if you give 4 point as an input then four column of
	//Xrc_z give interpolated value of four point. , Last one is resampling algoritham


	void interpolateZ_Generalize( const std::vector<GDALDataset*>& IrC_inputDS,
					   const std::vector<GDALDataset*>& IrC_maskDS,
					   const std::vector<std::vector<int>>& IrC_band,
		               std::vector<Point2D>& IrC_pt,
					   std::vector<std::vector<double>>& XrC_z,
					   const GDALResampleAlg resampling);


	// InterpolationZ is special case of above function which take only one point as an input and it is not contains mask. Also
	//It process only one band so input is gray scale raster for this function. It will not process multiband raster.
	void interpolateZ( const std::vector<GDALDataset*>& IrC_inputDS,
		               std::vector<Point2D>& IrC_pt,
					   std::vector<std::vector<double>>& XrC_z,
					   const GDALResampleAlg resampling);

	void interpolateZ( const std::vector<GDALDataset*>& IrC_inputDS,
													   double x,
													   double y,
													   double z,
												       const GDALResampleAlg resampling);

	void vgridshift(double x,double y, double *z);

	double GetValueAt(GDALDataset* hDataset, double x, double y);
};

class CPL_DLL OGRCoordinateTransformation3D:public OGRCoordinateTransformation
{
public:
	OGRCoordinateTransformation3D();
};

OGRCoordinateTransformation3D CPL_DLL *
OGRCreateCoordinateTransformation3D( OGRSpatialReference3D *poSource,
                                   OGRSpatialReference3D *poTarget );


#endif