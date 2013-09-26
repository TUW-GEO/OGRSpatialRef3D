/******************************************************************************
 *
 * Project:  OpenGIS Simple Features Reference Implementation
 * Purpose:  Classes for manipulating spatial reference systems with
 *           vertical datum suppurt in a platform non-specific manner.
 * Authors:  Gottfried Mandlburger, Johannes Otepka, Bhargav patel
 *
 ******************************************************************************
 * Copyright (c) 2012-2013,  I.P.F., TU Vienna.
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
 *
 *  II ) an additional height correction raster model to compensate the
 *       errors of historic heigth system realizations (e.g. MGI, DHHN...)
 *
 *  III) a constant offset to deal with local height origins, like tidal
 *       height systems, local height definitions (e.g. Wiener Null) ...
 *
 *  IV ) a scaling factor (e.g. foot-meter-conversion)
 *
 ************************************************************************/


class CPL_DLL OGRSpatialReference3D:public OGRSpatialReference
{
	bool bHasGeoid;		/**< flag to indicate the spatial reference has external geoid undulation model */
	bool bHasVCorr;		/**< flag to indicate the spatial reference has external height correction model */

	double dfVOffset_;	/**< constant shift value as additional height correction */
	double dfVScale_;	/**< constant scaling value as additional height correction (not used as of now) */

	bool is_debug;		/**< flag to indicate debugging is active */
	double *dbg_geoid;	/**< pointer to store geoid values from raster for debugging & validation test */
	double *dbg_vcorr;	/**< pointer to store height correction values from raster for debugging & validation test */

	RasterResampler  *poGeoid;	/**< pointer handle to RasterResampler object responsible for lookup from geoid undulation raster */
	RasterResampler  *poVCorr;	/**< pointer handle to RasterResampler object responsible for lookup from height correction model raster */
public:
	OGRSpatialReference3D();
	virtual    ~OGRSpatialReference3D();

	//! Extended method to load additional raster parameter in WKT
    /*!
      \param pszWKT a pointer to string containing WKT text.
      \return OGRERR_NONE if successful or OGRERR_FAILURE if unable to load extra parameter
      \sa importFromWkt()
    */
	OGRErr      importFromWkt3D( char ** pszWKT );

	//! method to manually supply additional Geoid undulation model raster
    /*!
      \param pszGeoidModel a string containing filename of GDAL compatible raster.
      \return OGRERR_NONE if successful or OGRERR_FAILURE if unable to load given filename
	  \sa GetGeoidModel()
    */
 	OGRErr SetGeoidModel( const char * pszGeoidModel );
	
	//! function to retrieve geoid undulation model filename
    /*!
      \return string indicating filename of geoid undulation model
	  \sa SetGeoidModel()
    */
    const char * GetGeoidModel ();

	//! method to manually supply additional Height undulation model raster
    /*!
      \param pszVCorrModel a string containing filename of GDAL compatible raster.
      \return OGRERR_NONE if successful or OGRERR_FAILURE if unable to load given filename
	  \sa GetVCorrModel()
    */
    OGRErr SetVCorrModel( const char * pszVCorrModel );
	
	//! function to retrieve height correction model filename
    /*!
      \return string indicating filename of height correction model
	  \sa SetVCorrModel()
    */
    const char * GetVCorrModel ();

	//! method to manually supply additional Vertical shift value
    /*!
      \param dfVOffset a float value.
    */
    void SetVOffset( double  dfVOffset );
	
	//! function to get additional Vertical shift value
    double GetVOffset ();

	//! method to manually supply additional Vertical shift value (not used for now)
    /*!
      \param dfVScale a float value.
    */
    void SetVScale( double  dfVScale );
	
	//! function to get additional Vertical scaling value
    double GetVScale ();

	//! function to check for additional lookup from external vertical model (geoid/height correction)
    /*!
      \return true if the spatial reference has additional external raster model for height correction
	  \sa ApplyVerticalCorrection(), HasGeoidModel(), HasVCorrModel()
    */
	bool HasVerticalModel();

	//! used internally by implementation of OGRCoordinateTransformation3D
	/*!
      \param is_inverse a boolean value indicating the direction of coordinate transformation.
	  \param point_count an integer value indicating the number of point involved in transformation
	  \param x pointer to double or array of double for values of first coordinate axis
	  \param y pointer to double or array of double for values of second coordinate axis
	  \param z pointer to double or array of double for values of third coordinate axis
	  \return OGRERR_NONE if transformation succesful
    */
	OGRErr ApplyVerticalCorrection(int is_inverse, unsigned int point_count, double *x, double *y, double *z);

	//! method to set debug mode (retrieve raster values of the points in transformation)
    /*!
      \param debug_mode a boolean value indicating the status of debugging mode.
	  \sa SetDebugData()
    */
	void SetDebug(bool debug_mode);

	//! method to supply data buffer used in debug mode.
    /*!
      \param geoid_undulation a pointer of double or array of double having the same number of points in transformation
	  \param vert_correction a pointer of double or array of double having the same number of points in transformation
	  \sa SetDebug()
    */
	void SetDebugData(double *geoid_undulation, double *vert_correction);

protected:
	//! function to check whether or not the spatial reference has additional geoid undulation model.
	/*!
      \return true if this spatial reference has additional geoid undulation model.
	  \sa HasVerticalModel()
    */
	bool HasGeoidModel();

	//! function to check whether or not the spatial reference has additional height correction model.
	/*!
      \return true if this spatial reference has additional height correction model.
	  \sa HasVerticalModel()
    */
	bool HasVCorrModel();

	//deprecated
	//double GetValueAt(GDALDataset* hDataset, double x, double y);
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