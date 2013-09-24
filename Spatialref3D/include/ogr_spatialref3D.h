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

#include "cpl_port.h"
#include "gdal_priv.h"
#include "ogr_spatialref.h"
#include "res_manager.h"

CPL_C_START

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


class CPL_DLL OGRSpatialReference3D:public OGRSpatialReference
{
	bool bHasGeoid;
	bool bHasVCorr;

	double dfVOffset_;
	double dfVScale_;

	bool is_debug;
	double *dbg_geoid;
	double *dbg_vcorr;

	RasterResampler  *poGeoid;
	RasterResampler  *poVCorr;
public:
	OGRSpatialReference3D();
	virtual    ~OGRSpatialReference3D();

	OGRErr      importFromWkt3D( char ** pszWKT );

 	OGRErr SetGeoidModel( const char * pszGeoidModel );
    const char * GetGeoidModel ();

    OGRErr SetVCorrModel( const char * pszVCorrModel );
    const char * GetVCorrModel ();

    OGRErr SetVOffset( double  dfVOffset );
    double GetVOffset ();

    OGRErr SetVScale( double  dfVScale );
    double GetVScale ();

	bool HasVerticalModel();

	OGRErr ApplyVerticalCorrection(int is_inverse, unsigned int point_count, double *x, double *y, double *z);
	void SetDebug(bool debug_mode);
	void SetDebugData(double *geoid_undulation, double *vert_correction);

protected:
	bool HasGeoidModel();
	bool HasVCorrModel();

	double GetValueAt(GDALDataset* hDataset, double x, double y);
};

class CPL_DLL OGRCoordinateTransformation3D:public OGRCoordinateTransformation
{
public:
	OGRCoordinateTransformation3D();
};

CPL_DLL OGRCoordinateTransformation3D *
OGRCreateCoordinateTransformation3D( OGRSpatialReference3D *poSource,
                                   OGRSpatialReference3D *poTarget );

CPL_C_END

#endif