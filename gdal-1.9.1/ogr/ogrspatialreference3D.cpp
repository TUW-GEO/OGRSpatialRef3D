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
#include <vector>

#include "..\gcore\gdal.h"
#include "..\alg\gdalwarper.h"
#include "..\gcore\gdal_priv.h"
#include "ogr_spatialref3D.h"


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

static int BilinearResampling( float *pf_grid, unsigned char *pi_mask, 
                               double dfNoData, int blockX, int blockY,
                               double dfSrcX, double dfSrcY,
                               double& dfValue )
{
  ////////////////////////////////////////////////////////////////////
  // Implementation based on gdalwarpkernel.cpp by f.warmadam
  ////////////////////////////////////////////////////////////////////
  double  dfValAccumulator = 0.0;
  double  dfMskAccumulator = 0.0;
  double  dfAccumulatorDivisor = 0.0;

  double dx = dfSrcX - 0.5;
  double dy = dfSrcY - 0.5;
  double dxmax = (double) (blockX-1);
  double dymax = (double) (blockY-1);

  dfValue = dfNoData;

  if( dx < 0. || dx > dxmax || dy < 0. || dy > dymax )
    return 0;

  int     iSrcX = dx == dxmax ? (int) (dxmax-1) : (int) floor(dfSrcX - 0.5);
  int     iSrcY = dy == dymax ? (int) (dymax-1) : (int) floor(dfSrcY - 0.5);
  int     iSrcOffset = iSrcX + iSrcY * blockX;
  double  dfRatioX = 1.5 - (dfSrcX - iSrcX);
  double  dfRatioY = 1.5 - (dfSrcY - iSrcY);

  float dfValUL, dfValUR, dfValLL, dfValLR;     // grid values of all four corner pixels
  float dfMskUL, dfMskUR, dfMskLL, dfMskLR;     // mask values of all four corner pixels
  dfValUL = pf_grid[iSrcOffset];
  dfValUR = pf_grid [iSrcOffset+1];
  dfValLR = pf_grid[iSrcOffset+1+blockX];
  dfValLL = pf_grid[iSrcOffset+blockX];
  if( pi_mask )
  {
    dfMskUL = (float) pi_mask[iSrcOffset];
    dfMskUR = (float) pi_mask [iSrcOffset+1];
    dfMskLR = (float) pi_mask[iSrcOffset+1+blockX];
    dfMskLL = (float) pi_mask[iSrcOffset+blockX];
  }
  else
  {
    dfMskUL = dfMskUR = dfMskLR = dfMskLL = 1.0;  // oder doch 0.0 -> wie ist die Maske definiert?
  }

  // Upper Left Pixel
  if( iSrcX >= 0 && iSrcX < blockX && iSrcY >= 0 && iSrcY < blockY && dfValUL != dfNoData )
  {
      double dfMult = dfRatioX * dfRatioY;
      dfAccumulatorDivisor += dfMult;
      dfValAccumulator += dfValUL * dfMult;
      dfMskAccumulator += dfMskUL * dfMult;
  }
      
  // Upper Right Pixel
  if( iSrcX+1 >= 0 && iSrcX+1 < blockX && iSrcY >= 0 && iSrcY < blockY && dfValUR != dfNoData )
  {
      double dfMult = (1.0-dfRatioX) * dfRatioY;
      dfAccumulatorDivisor += dfMult;
      dfValAccumulator += dfValUR * dfMult;
      dfMskAccumulator += dfMskUR * dfMult;
  }
      
  // Lower Right Pixel
  if( iSrcX+1 >= 0 && iSrcX+1 < blockX && iSrcY+1 >= 0 && iSrcY+1 < blockY && dfValLR != dfNoData )
  {
      double dfMult = (1.0-dfRatioX) * (1.0-dfRatioY);
      dfAccumulatorDivisor += dfMult;
      dfValAccumulator += dfValLR * dfMult;
      dfMskAccumulator += dfMskLR * dfMult;
  }
      
  // Lower Left Pixel
  if( iSrcX >= 0 && iSrcX < blockX && iSrcY+1 >= 0 && iSrcY+1 < blockY && dfValLL != dfNoData )
  {
      double dfMult = dfRatioX * (1.0-dfRatioY);
      dfAccumulatorDivisor += dfMult;
      dfValAccumulator += dfValLL * dfMult;
      dfMskAccumulator += dfMskLL * dfMult;
  }

	if ( dfAccumulatorDivisor == 0.0)
		return 0;

/* -------------------------------------------------------------------- */
/*      Return result.                                                  */
/* -------------------------------------------------------------------- */
  if( dfMskAccumulator / dfAccumulatorDivisor > 0.85 )
  {
    if( dfAccumulatorDivisor == 1.0 )
      dfValue = dfValAccumulator;
  //  else if( dfAccumulatorDivisor > 0.00001 )
    else if( dfAccumulatorDivisor > 0.5 )
      dfValue = dfValAccumulator / dfAccumulatorDivisor;
  }
  return dfValue != dfNoData;
}



OGRSpatialReference3D::OGRSpatialReference3D()
{
	GDALAllRegister();
}

OGRSpatialReference3D::~OGRSpatialReference3D()
{
}


OGRSpatialReference3D::OGRSpatialReference3D(const OGRSpatialReference&)
{
}

OGRSpatialReference3D::OGRSpatialReference3D(const char * pszWKT,
                          const char * pszGeoidModel,
                          const char * pszVCorrModel,
                          double dfVOffset,
                          double dfVScale)
{
}

OGRErr OGRSpatialReference3D::SetVOffset( double  dfVOffset )
{
	dfVOffset_=dfVOffset;
	return OGRERR_NONE;
}

double OGRSpatialReference3D::GetVOffset ()
{
	return dfVOffset_;
}


OGRErr OGRSpatialReference3D::SetVScale( double  dfVScale )
{
	dfVScale_=dfVScale;
	return OGRERR_NONE;
}

double OGRSpatialReference3D::GetVScale ()
{
	return dfVScale_;
}

OGRErr OGRSpatialReference3D::SetGeoidModel( const char * pszGeoidModel )
{
	poDataset = (GDALDataset *) GDALOpen( pszGeoidModel, GA_ReadOnly );
	
	if( poDataset == NULL )
    {
        printf("gdal failed - unable to open '%s'.\n",
                 pszGeoidModel );
	}
	return OGRERR_NONE;
}

OGRErr OGRSpatialReference3D::SetVCorrModel( const char * pszVCorrModel )
{
	hDatasetVCorrModel_=(GDALDataset *) GDALOpen( pszVCorrModel, GA_ReadOnly );

	if( poDataset == NULL )
    {
        printf("gdal failed - unable to open '%s'.\n",
                 pszVCorrModel );
	}
	return OGRERR_NONE;
}

void OGRSpatialReference3D::vgridshift(double x,double y,double *z)
{
	//printf("%f %f",x,y);

	vector <GDALDataset*> IrC_inputDS;

	IrC_inputDS.push_back(hDatasetVCorrModel_);
	IrC_inputDS.push_back(poDataset);

	Point2D p1;
	p1.setxy(x,y);
	
	vector <Point2D> setPoint;
	setPoint.push_back(p1);
	
	//interpolateZ_Generalize(a,Irc_mask,irc_band,setPoint,z,GDALResampleAlg::GRA_Bilinear);

	vector< vector<double>> z1;

	
	interpolateZ(IrC_inputDS,setPoint,z1,GDALResampleAlg::GRA_Bilinear);

	double Zell=150.490,Zortho;

	*z=(Zell+GetVOffset())*GetVScale()+z1[0][0]+z1[1][0];
	
	/*printf("Zortho value is %f\n",Zortho);

	p1.setz(Zortho);

	printf("X =%f Y=%f  Z=%f\n",p1.X(),p1.Y(),p1.Z());*/

	
}

void OGRSpatialReference3D::interpolateZ_Generalize( const std::vector<GDALDataset*>& IrC_inputDS,
										const std::vector<GDALDataset*>& IrC_maskDS,
										const std::vector<std::vector<int>>& IrC_band,
										std::vector<Point2D>& IrC_pt,
										std::vector<std::vector<double>>& XrC_z,
										const GDALResampleAlg resampling)
{

 // at the moment only bilinear resampling is supported
  CPLAssert( resampling == GDALResampleAlg::GRA_Bilinear || 
             resampling == GDALResampleAlg::GRA_NearestNeighbour );

  size_t sizeDS = IrC_inputDS.size();
  size_t sizeMask = IrC_maskDS.size() ? IrC_maskDS.size() : sizeDS;

  if( sizeDS != IrC_band.size() || sizeDS != sizeMask)
    throw std::exception("Number of raster datasets, grid masks and band vectors must be identical.");

   // determine minimum block size depending on specified resampling method
  int minBlockSize = resampling == GDALResampleAlg::GRA_Bilinear ? 2 : 1;

  // initialie coordinate vectors
  size_t nrPts = IrC_pt.size();
  std::vector<double> x( nrPts );
  std::vector<double> y( nrPts );

   // check existence of raster datasets and bands
  for( size_t ids=0; ids<IrC_inputDS.size(); ids++ )
  {
    GDALDataset* poDS = IrC_inputDS[ids];
    GDALDataset* poMaskDS = IrC_maskDS.size() ? IrC_maskDS[ids] : NULL;

    if( !poDS )
      throw std::exception( "Raster dataset(s) not properly initialized." );
    for( size_t ibd=0; ibd<IrC_band.size(); ibd++ )
    {
      GDALRasterBand* poBand = poDS->GetRasterBand( IrC_band[ids][ibd] );
      if( !poBand )
        throw std::exception( "Raster band(s) not properly initialized." );
      if( poMaskDS )
      {
        poBand = poMaskDS->GetRasterBand( IrC_band[ids][ibd] );
        if( !poBand )
          throw std::exception( "Mask band(s) not properly initialized." );
      }
    }
  }

   // resize and initialize result vector
  size_t nrBands = 0;
  for( size_t i=0; i < IrC_band.size(); i++ )
    nrBands += IrC_band[i].size();

  XrC_z.resize( nrBands ); 
  for( size_t ids=0; ids<IrC_inputDS.size(); ids++ )
    for( size_t ibd=0; ibd<IrC_band.size(); ibd++ )
      XrC_z[ids*sizeDS+ibd].resize( nrPts, IrC_inputDS[ids]->GetRasterBand(IrC_band[ids][ibd])->GetNoDataValue() );
   for( size_t ids=0; ids < IrC_inputDS.size(); ids++ )
  {
    GDALDataset* pC_inputDS = IrC_inputDS[ids];
    GDALDataset* pC_maskDS  = IrC_maskDS.size() ? IrC_maskDS[ids] : NULL;

    //------------------------------------------------------------
    // fetch geotransform and pixel/lines of input raster band
    int nXsize, nYsize;
    double geotrans[6];
    pC_inputDS->GetGeoTransform(geotrans);
    nXsize = pC_inputDS->GetRasterXSize();
    nYsize = pC_inputDS->GetRasterYSize();

    // check raster size
    if( nXsize < minBlockSize || nYsize < minBlockSize )
      continue;

    //------------------------------------------------------------
    // get inverse geotransformation
    double geotransInv[6];
    if( GDALInvGeoTransform( geotrans, geotransInv ) == 0 )
      throw std::exception( "inversion of geo transformation failed." );

    bool b_useMask = pC_maskDS!= 0;
    if( b_useMask )
    {
      int nx, ny;
      double geotrf[6];
      pC_maskDS->GetGeoTransform(geotrf);
      nx = pC_maskDS->GetRasterXSize();
      ny = pC_maskDS->GetRasterYSize();
      bool b_identicalGridStruct = nx == nXsize && ny == nYsize;
      if( b_identicalGridStruct )
      {
        for( int i=0; i<6; i++ )
        {
          if( geotrf[i] != geotrans[i] )
          {
            b_identicalGridStruct = false;
            break;
          }
        }
      }
      if( !b_identicalGridStruct )
        throw std::exception( "grid sturcture of surface grid and mask model not identical." );
    }

    //------------------------------------------------------------
    // get pixel/lines coordinates of all points
    double xmin =  DBL_MAX;
    double ymin =  DBL_MAX;
    double xmax = -DBL_MAX;
    double ymax = -DBL_MAX;
    for( int i=0; i<nrPts; i++ )
    {
		double dPixel = geotransInv[0]+ 
          + IrC_pt[i].X() * geotransInv[1]
          + IrC_pt[i].Y() * geotransInv[2];
      double dLine= geotransInv[3]
          + IrC_pt[i].X() * geotransInv[4]
          + IrC_pt[i].Y() * geotransInv[5];
      x[i] = dPixel;
      y[i] = dLine;
      xmin = min( xmin, dPixel );
      xmax = max( xmax, dPixel );
      ymin = min( ymin, dLine );
      ymax = max( ymax, dLine );
    }

    // abort if none of the points are within the extents
    if( xmax < 0. || xmin > (double) nXsize ||
        ymax < 0. || ymin > (double) nYsize )
      break;
    // reduce points extents to grid extents
    xmin = floor( max( xmin-0.5, 0. ) );
    ymin = floor( max( ymin-0.5, 0. ) );
    xmax = ceil( min( xmax+0.5, (double) nXsize ) );
    ymax = ceil( min( ymax+0.5, (double) nYsize ) );

    // read data from input grid/mask
    int ixmin = (int) xmin;  int ixmax = (int) xmax;
    int iymin = (int) ymin;  int iymax = (int) ymax;
    int blockX = max( ixmax-ixmin, minBlockSize );
    int blockY = max( iymax-iymin, minBlockSize );

    if( ixmin+blockX > nXsize )
      ixmin = nXsize - blockX;
    if( iymin+blockY > nYsize )
      iymin = nYsize - blockY;

    float *pf_grid = (float *) CPLMalloc(sizeof(float)*blockX*blockY);
    unsigned char *pi_mask = 0;
    if( b_useMask )
      pi_mask = (unsigned char *) CPLMalloc(sizeof(unsigned char)*blockX*blockY);

    //------------------------------------------------------------
    // process all bands
    for( unsigned ibd=0; ibd<IrC_band.size(); ibd++ )
    {
      int bandIdx = IrC_band[ids][ibd];
      GDALRasterBand* pC_inband = pC_inputDS->GetRasterBand(bandIdx);
      GDALRasterBand* pC_maskband = b_useMask ? pC_maskDS->GetRasterBand(bandIdx) : NULL;

      // fetch no data value
      double dNoData = pC_inband->GetNoDataValue();

      if( pC_inband->RasterIO( GF_Read, ixmin, iymin, blockX, blockY, 
                               pf_grid, blockX, blockY, GDT_Float32, 0, 0 ) != CE_None )
      {
        throw std::exception( "internal error: reading raster data for z-interpolation failed." );
      }
      if( b_useMask  )
      {
        if( pC_maskband->RasterIO( GF_Read, ixmin, iymin, blockX, blockY, 
                                 pi_mask, blockX, blockY, GDT_Byte, 0, 0 ) != CE_None )
        {
          throw std::exception( "internal error: reading mask data for z-interpolation failed." );
        }
        else
        {
          // replace all valid mask values (i.e. > 0 9 by "1" 
          // The mask is considered by simply interpolating the mask value
          // analogue to the grid interpolation and all grid values are 
          // considered valid, if the mask value is higher than 0.85
          for( int ii=0; ii<blockX*blockY; ii++ )
            pi_mask[ii] = pi_mask[ii] > 0 ? 1 : 0;
        }
      }

      //----------------------------------------------------------------------------
      // interpolation of z-values
      //----------------------------------------------------------------------------
	  if( resampling == GDALResampleAlg::GRA_NearestNeighbour )
      {
        for( unsigned int ii=0; ii<nrPts; ii++ )
        {
          // perform bilinear resampling
          double dfSrcX = x[ii]-xmin;
          double dfSrcY = y[ii]-ymin;
   //       NearestResampling( pf_grid, pi_mask, dNoData, blockX, blockY, dfSrcX, dfSrcY, XrC_z[ids*sizeDS+ibd][ii] );
        }
      }
      else
      {
        for( unsigned int ii=0; ii<nrPts; ii++ )
        {
          // perform bilinear resampling
          double dfSrcX = x[ii]-xmin;
          double dfSrcY = y[ii]-ymin;
          BilinearResampling( pf_grid, pi_mask, dNoData, blockX, blockY, dfSrcX, dfSrcY, XrC_z[ids*sizeDS+ibd][ii] );
        }
      }
    }
	  CPLFree(pf_grid);
	  if (pi_mask)
		  CPLFree(pi_mask);
  }
}


void OGRSpatialReference3D::interpolateZ( const std::vector<GDALDataset*>& IrC_inputDS,
													   std::vector<Point2D>& IrC_pt,
												       std::vector<std::vector<double>>& XrC_z,
													   const GDALResampleAlg resampling)
{
	// at the moment only bilinear resampling is supported
  CPLAssert( resampling == GDALResampleAlg::GRA_Bilinear || 
             resampling == GDALResampleAlg::GRA_NearestNeighbour );

  size_t sizeDS = IrC_inputDS.size();

   // determine minimum block size depending on specified resampling method
  int minBlockSize = resampling == GDALResampleAlg::GRA_Bilinear ? 2 : 1;

   // initialie coordinate vectors
  size_t nrPts = IrC_pt.size();
  std::vector<double> x( nrPts );
  std::vector<double> y( nrPts );

   // check existence of raster datasets and bands
  for( size_t ids=0; ids<IrC_inputDS.size(); ids++ )
  {
    GDALDataset* poDS = IrC_inputDS[ids];
   
    if( !poDS )
      throw std::exception( "Raster dataset(s) not properly initialized." );
    
    GDALRasterBand* poBand = poDS->GetRasterBand(1);
    
	if( !poBand )
        throw std::exception( "Raster band(s) not properly initialized." );
   }  

  XrC_z.resize( 2 );  //Bhargav: because we have two raster

  for( size_t ids=0; ids<IrC_inputDS.size(); ids++ )
      XrC_z[ids].resize( nrPts);
 

  for( size_t ids=0; ids < IrC_inputDS.size(); ids++ )
  {
    GDALDataset* pC_inputDS = IrC_inputDS[ids];

	//------------------------------------------------------------
    // fetch geotransform and pixel/lines of input raster band
    int nXsize, nYsize;
    double geotrans[6];
    pC_inputDS->GetGeoTransform(geotrans);
    nXsize = pC_inputDS->GetRasterXSize();
    nYsize = pC_inputDS->GetRasterYSize();

    // check raster size
    if( nXsize < minBlockSize || nYsize < minBlockSize )
      continue;

	//------------------------------------------------------------
    // get inverse geotransformation
    double geotransInv[6];
    if( GDALInvGeoTransform( geotrans, geotransInv ) == 0 )
      throw std::exception( "inversion of geo transformation failed." );

	//------------------------------------------------------------
    // get pixel/lines coordinates of all points
    double xmin =  DBL_MAX;
    double ymin =  DBL_MAX;
    double xmax = -DBL_MAX;
    double ymax = -DBL_MAX;
    for( int i=0; i<nrPts; i++ )
    {
      double dPixel = geotransInv[0]
          + IrC_pt[i].X() * geotransInv[1]
          + IrC_pt[i].Y() * geotransInv[2];
      double dLine= geotransInv[3]
          + IrC_pt[i].X() * geotransInv[4]
          + IrC_pt[i].Y() * geotransInv[5];
      x[i] = dPixel;
      y[i] = dLine;
      xmin = min( xmin, dPixel );
      xmax = max( xmax, dPixel );
      ymin = min( ymin, dLine );
      ymax = max( ymax, dLine );
    }
	// abort if none of the points are within the extents
    if( xmax < 0. || xmin > (double) nXsize ||
        ymax < 0. || ymin > (double) nYsize )
      break;
    // reduce points extents to grid extents
    xmin = floor( max( xmin-0.5, 0. ) );
    ymin = floor( max( ymin-0.5, 0. ) );
    xmax = ceil( min( xmax+0.5, (double) nXsize ) );
    ymax = ceil( min( ymax+0.5, (double) nYsize ) );

    // read data from input grid/mask
    int ixmin = (int) xmin;  int ixmax = (int) xmax;
    int iymin = (int) ymin;  int iymax = (int) ymax;
    int blockX = max( ixmax-ixmin, minBlockSize );
    int blockY = max( iymax-iymin, minBlockSize );

    if( ixmin+blockX > nXsize )
      ixmin = nXsize - blockX;
    if( iymin+blockY > nYsize )
      iymin = nYsize - blockY;

	float *pf_grid = (float *) CPLMalloc(sizeof(float)*blockX*blockY);
	unsigned char *pi_mask = 0;
    

	GDALRasterBand* pC_inband = pC_inputDS->GetRasterBand(1);

	// fetch no data value
    double dNoData = pC_inband->GetNoDataValue();

	if( pC_inband->RasterIO( GF_Read, ixmin, iymin, blockX, blockY, 
                               pf_grid, blockX, blockY, GDT_Float32, 0, 0 ) != CE_None )
      {
        throw std::exception( "internal error: reading raster data for z-interpolation failed." );
      }

	//----------------------------------------------------------------------------
      // interpolation of z-values
      //----------------------------------------------------------------------------
      if( resampling == GDALResampleAlg::GRA_NearestNeighbour )
      {
        for( unsigned int ii=0; ii<nrPts; ii++ )
        {
          // perform bilinear resampling
          double dfSrcX = x[ii]-xmin;
          double dfSrcY = y[ii]-ymin;
      //    NearestResampling( pf_grid, pi_mask, dNoData, blockX, blockY, dfSrcX, dfSrcY, XrC_z[ids*sizeDS+ibd][ii] );
        }
      }
      else
      {
        for( unsigned int ii=0; ii<nrPts; ii++ )
        {
          // perform bilinear resampling
          double dfSrcX = x[ii]-xmin;
          double dfSrcY = y[ii]-ymin;
          BilinearResampling( pf_grid, pi_mask, dNoData, blockX, blockY, dfSrcX, dfSrcY, XrC_z[ids][ii] );
        }
  }

 CPLFree(pf_grid);
   
  
}
}




OGRCoordinateTransformation3D::OGRCoordinateTransformation3D()
{

}

