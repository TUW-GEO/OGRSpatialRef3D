/******************************************************************************
 *
 * Project: OGR SpatialRef3D
 * Purpose: helper class for managing lookup raster (Geoid, vertical correction)
 *          with ability to reduce disk seek by holding a portion of raster data
 *          in memory
 * Author: Peb Ruswono Aryan, Gottfried Mandlburger, Johannes Otepka
 *
 ******************************************************************************/
#ifndef __RES_MANAGER_H__
#define __RES_MANAGER_H__

#include "gdal_priv.h"
#include "ogr_core.h"
#include "cpl_string.h"

class RasterResampler
{
public:
	RasterResampler();
	virtual ~RasterResampler();

	// per-point interface
	double GetValueAt(double x, double y);

	// bulk of points interface (input x, input y, output z)
	void GetValueAt(int point_count, double *x, double *y, double *z);

	OGRErr Open(const char *pszFilename);
	const char* GetFilename();

protected:
	CPLString sFilename;	

	GDALDataset *poData;	
	int nRasterWidth;		
	int nRasterHeight;		
	bool bIsSmall;			// is raster small enough to be kept in memory (decided when prepare)

	//double dGeotrans[6];	// do we need to keep this here or local variable is enough ?
	double dNoDataValue;
	double dInvGeotrans[6];	

	double *padWindow;		// window pixels

	int nWndXOffset;		// window left
	int nWndYOffset;		// window top
	int nWndWidth;			// window width
	int nWndHeight;			// window height

	// refactor for managing interpolation
	double GetValueResampled(double x, double y);

	// caching utilities
	void Prepare();
	void Request(int left, int top, int width, int height);
	void Cleanup();
	void MapToRaster(double *x, double *y);
};


#endif
