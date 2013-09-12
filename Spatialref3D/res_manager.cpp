#include "res_manager.h"
#include "cpl_port.h"
#include "interpolation.h"
#include <iostream>

#define RAD_TO_DEG	57.29577951308232
#define MAXINT 9999999
#define MAXEXTENT 1024		// maximum window size

#define INSIDE(x,y,l,t,w,h) (x>=l && x<=l+w && y>=t && y<=t+h)

RasterResampler::RasterResampler() : bIsSmall(false), 
									nRasterWidth(0), 
									nRasterHeight(0),
									nWndXOffset(0),
									nWndYOffset(0),
									nWndWidth(0),
									nWndHeight(0)
{
	padWindow = NULL;
	poData = NULL;
}

RasterResampler::~RasterResampler()
{
	Cleanup();

	if(poData != NULL)
		GDALClose(poData);
}

double
	RasterResampler::GetValueAt(double x, double y)
{
	double dPixel = x;
	double dLine = y;
	MapToRaster(&dPixel, &dLine);

	// check if buffer not initialized or point not inside current window
	// naive caching strategy
	if (padWindow == NULL 
				|| !INSIDE((int)dPixel, (int)dLine, nWndXOffset, nWndYOffset, nWndWidth, nWndHeight)){

		int nWndLeft = MAX(0, MIN(nRasterWidth-2, (int)dPixel-MAXEXTENT/2));
		int nWndTop = MAX(0, MIN(nRasterHeight-2, (int)dLine-MAXEXTENT/2));

		//std::cout << "sample window (" << nWndLeft << " " << nWndTop << ")";
		//std::cout << " of (" << nRasterWidth << ", " << nRasterHeight << ")" << std::endl;

		Request(nWndLeft, nWndTop, MAXEXTENT, MAXEXTENT);
	}
	return GetValueResampled(dPixel, dLine);
}

void
	RasterResampler::GetValueAt(int point_count, double *x, double *y, double *z)
{
	// dummy naive approach (tested)
	if (point_count==1)
	{
		z[0] = GetValueAt(x[0], y[0]); 
		return;
	}
	//for(int i=0; i<point_count; ++i) z[i] = GetValueAt(x[i], y[i]); return;
	//
	// indexed (TESTED ON SMALL window ONLY)
	double* padX = (double*)CPLMalloc(sizeof(double)*point_count);
	double* padY = (double*)CPLMalloc(sizeof(double)*point_count);

	double dXMin = MAXINT;
	double dYMin = MAXINT;
	double dXMax = -MAXINT;
	double dYMax = -MAXINT;

	// convert all point's coordinate
	// from map coordinate to 
	// raster coordinate
	for(int i=0; i<point_count; ++i){
		double px = x[i];
		double py = y[i];

		MapToRaster(&px, &py);

		padX[i] = px;
		padY[i] = py;

		dXMin = MIN(dXMin, px);
		dYMin = MIN(dYMin, py);
		dXMax = MAX(dXMax, px);
		dYMax = MAX(dYMax, py);
	}

	double dWidth = ceil(dXMax)-floor(dXMin);
	double dHeight = ceil(dYMax)-floor(dYMin);

	if (bIsSmall || (dWidth < MAXEXTENT && dHeight < MAXEXTENT)){
		// do it in single patch
		// use naive approach embedded in single point interface
		int drwLeft = MAX(0, MIN(nRasterWidth-2, (int)dXMin));
		int drwTop = MAX(0, MIN(nRasterHeight-2,(int)dYMin));
		
		// TODO: add checking to window boundary against raster boundary
		Request(drwLeft, drwTop, (int)dWidth+1, (int)dHeight+1);

		for(int i=0; i<point_count; ++i){
			z[i] = GetValueResampled(padX[i], padY[i]);
		}
	}
	else{
		// both raster and window size not small enough to be all in memory
		// do it by window patch
		bool* panIdx = (bool*)CPLMalloc(sizeof(bool)*point_count);
		int unprocessed_count = point_count;	// # of point unprocessed

		for(int i=0; i<point_count; ++i)
				panIdx[i] = false;
		
		while(unprocessed_count < point_count){
			//initialize next cache
			for(int i=0; i<point_count; ++i){
				if(!panIdx[i]){
					if (padWindow == NULL){
				
						int nWndLeft = MAX(0, MIN(nRasterWidth-2, (int)padX[i]-MAXEXTENT/2));
						int nWndTop = MAX(0, MIN(nRasterWidth-2, (int)padY[i]-MAXEXTENT/2));

						//std::cout << nWndLeft << " " << nWndTop << std::endl;

						Request(nWndLeft, nWndTop, MAXEXTENT, MAXEXTENT);
					}

					break;
				}
			}

			// process all applicable points
			for(int i=0; i<point_count; ++i){
				// skip processed points or points outside current cache window
				if(panIdx[i]
					|| !INSIDE((int)floor(padX[i]), (int)floor(padY[i]), nWndXOffset, nWndYOffset, nWndWidth-1, nWndHeight-1)) //add floor and -1
						continue;

				z[i] = GetValueResampled(padX[i], padY[i]);

				panIdx[i] = true;	// mark processed
				unprocessed_count += 1;
			}
		}// endwhile
		
		CPLFree(panIdx);
	}

	CPLFree(padX);
	CPLFree(padY);
}

OGRErr
	RasterResampler::Open(const char *pszFilename)
{
	sFilename = pszFilename;
	poData = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
	
	if( poData == NULL )
    {
        printf("gdal failed - unable to open '%s'.\n",
                 pszFilename );
		return OGRERR_FAILURE;
	}

	Prepare();
	return OGRERR_NONE;
}

const char*
	RasterResampler::GetFilename()
{
	return sFilename;
}

double
	RasterResampler::GetValueResampled(double x, double y)
/*
 * x and y is assumed to be in raster coordinate (not buffer window coordinate)
 */
{
	int px = (int)floor(x);
	int py = (int)floor(y);

	// Boundary checking
	if (px < 0 || py < 0 || px >= nRasterWidth || py >= nRasterHeight){
		std::cerr << "point (" << px << "," << py << ") outside raster." << "(" << nRasterWidth << ", " << nRasterHeight<< ")" << std::endl;
		return 0.0;
		//throw std::exception( "point outside raster." );
	}
	else {
		double dx = x - px;
		double dy = y - py;

		// TODO: accomodate neigbor acquisition as required by other interpolation function (e.g. bicubic)
		// acquire neighbors
		int offset = (py-nWndYOffset)*nWndWidth+(px-nWndXOffset);
		double p[4];
		p[0] = padWindow[offset];
		p[1] = padWindow[offset+1];

		offset += nWndWidth;
		p[2] = padWindow[offset];
		p[3] = padWindow[offset+1];

		return bilinearInterpolation(p, dx, dy, dNoDataValue);
	}

	return 0.0;
}

void
	RasterResampler::Cleanup()
{
	if(padWindow != NULL)
		CPLFree(padWindow);
	padWindow = NULL;
}

void
	RasterResampler::Prepare()
	/*
	 * get metadata information from raster file
	 */
{
	nRasterWidth = poData->GetRasterXSize();
    nRasterHeight = poData->GetRasterYSize();
	bIsSmall = (nRasterWidth < MAXEXTENT) && (nRasterHeight < MAXEXTENT);

	dNoDataValue = poData->GetRasterBand(1)->GetNoDataValue();
	double geotrans[6];
    poData->GetGeoTransform(geotrans);

	if( GDALInvGeoTransform( geotrans, dInvGeotrans ) == 0 )
      throw std::exception( "inversion of geo transformation failed." );
}

void
	RasterResampler::Request(int left, int top, int width, int height)
	/*
	 * Access pixel data from raster file to temporary buffer (padWindow)
	 */
{
	//std::cout << "window request " << left << " " << top << " " << left+width << " " << top+height << std::endl;

	nWndXOffset = left;
	nWndYOffset = top;
	nWndWidth = MIN(nRasterWidth-left, width);
	nWndHeight = MIN(nRasterHeight-top, height);

	//std::cout << "window given " << left << " " << top << " " << left+nWndWidth << " " << top+nWndHeight << std::endl;
	
	Cleanup();

	int nWndArea = width*height;
	padWindow = (double *) CPLMalloc(sizeof(double)*nWndArea);

	poData->RasterIO( GF_Read, 
						nWndXOffset, nWndYOffset, nWndWidth, nWndHeight, 
                        padWindow, nWndWidth, nWndHeight, GDT_Float64, 
                        1, NULL, 0, 0, 0 );
}

void 
	RasterResampler::MapToRaster(double *x, double *y)
	/* 
	 * converts map coordinate (radian) to raster coordinate (pixels)
	 */
{
	*(x) = -0.5+dInvGeotrans[0]+ 
          + *(x) * dInvGeotrans[1] * RAD_TO_DEG
          + *(y) * dInvGeotrans[2] * RAD_TO_DEG;

    *(y) = -0.5+dInvGeotrans[3]
          + *(x) * dInvGeotrans[4] * RAD_TO_DEG
          + *(y) * dInvGeotrans[5] * RAD_TO_DEG;
}