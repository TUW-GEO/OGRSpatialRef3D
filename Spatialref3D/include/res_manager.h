/******************************************************************************
 *
 * Project: OGR SpatialRef3D
 * Purpose: helper class for managing lookup raster (Geoid, vertical correction)
 *          with ability to reduce disk seek by holding a portion of raster data
 *          in memory
 * Author: Peb Ruswono Aryan, Gottfried Mandlburger, Johannes Otepka
 *
 ******************************************************************************
 * Copyright (c) 2012-2014,  I.P.F., TU Vienna.
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
 ******************************************************************************/
#ifndef __RES_MANAGER_H__
#define __RES_MANAGER_H__

#include "gdal_priv.h"
#include "ogr_core.h"
#include "cpl_string.h"

class RasterResampler
{
	CPLString sFilename;	/**< string value containing filename of raster */

	GDALDataset *poData;	
	int nRasterWidth;		
	int nRasterHeight;		
	bool bIsSmall;			/**< is raster small enough to be kept in memory (decided when prepare) */

	//double dGeotrans[6];	// do we need to keep this here or local variable is enough ?
	double dNoDataValue;
	double dInvGeotrans[6];	/**< Inverse geotransform used for MapToRaster coordinate transform */

	double *padWindow;		/**< pointer to cached raster block of nWndWidth x nWndHeight size */

	int nWndXOffset;		/**< left position of cached raster block window */
	int nWndYOffset;		/**< top position of cached raster block window */
	int nWndWidth;			/**< width of cached raster block window */
	int nWndHeight;			/**< height of cahced raster block window */
public:
	RasterResampler();
	virtual ~RasterResampler();

	//! function to retrieve raster value at given point
	double GetValueAt(double x, double y);

	//! function to retrieve raster value at given array of points (x, y) and store in z
	void GetValueAt(int point_count, double *x, double *y, double *z);

	//! method to load GDAL compatible raster
	/*!
		\param pszFilename a string value indicating filename of raster
		\return OGRERR_NONE if loading successful or OGRERR_FAILURE otherwise
		\sa GetFilename()
	*/
	OGRErr Open(const char *pszFilename);

	//! function to retrieve raster filename
	/*!
		\return a string indicating filename of loaded raster
		\sa Open()
	*/
	const char* GetFilename();

protected:
	//! function to lookup raster value from given point in map space
	double GetValueResampled(double x, double y);

	//! caching utilities
	void Prepare();

	//! method to allocate and move raster block from file to memory
	void Request(int left, int top, int width, int height);

	//! method to release allocated buffer for storing cached raster block
	void Cleanup();

	//! function to convert coordinate from map space (lon, lat) to raster space (pixel, line)
	void MapToRaster(double *x, double *y);
};


#endif
