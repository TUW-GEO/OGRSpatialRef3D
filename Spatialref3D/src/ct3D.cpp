/******************************************************************************
 *
 * Project:  OpenGIS Simple Features Reference Implementation
 * Purpose:  Classes for manipulating spatial reference systems with
 *           vertical datum support in a platform non-specific manner.
 * Authors:  Peb Ruswono Aryan, Gottfried Mandlburger, Johannes Otepka, Bhargav Patel
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
 ****************************************************************************/
#include <iostream>

#include "ogr_spatialref3D.h"
#include "cpl_port.h"
#include "cpl_error.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include "cpl_multiproc.h"


#include "..\..\proj-4.8.0\src\projects.h"
#include "..\..\proj-4.8.0\src\geocent.h"



static void *hPROJMutex = NULL;
static int pj_adjust_axis( projCtx ctx, const char *axis, int denormalize_flag,
                           long point_count, int point_offset, 
                           double *x, double *y, double *z );

static const int transient_error[50] = {
    /*             0  1  2  3  4  5  6  7  8  9   */
    /* 0 to 9 */   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    /* 10 to 19 */ 0, 0, 0, 0, 1, 1, 0, 1, 1, 1,  
    /* 20 to 29 */ 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
    /* 30 to 39 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    /* 40 to 49 */ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 };


#ifndef SRS_WGS84_SEMIMAJOR
#define SRS_WGS84_SEMIMAJOR 6378137.0
#endif

#ifndef SRS_WGS84_ESQUARED
#define SRS_WGS84_ESQUARED 0.0066943799901413165
#endif

#define Dx_BF (defn->datum_params[0])
#define Dy_BF (defn->datum_params[1])
#define Dz_BF (defn->datum_params[2])
#define Rx_BF (defn->datum_params[3])
#define Ry_BF (defn->datum_params[4])
#define Rz_BF (defn->datum_params[5])
#define M_BF  (defn->datum_params[6])

/**\file ct3D.cpp
 * This class is an implementation of Coordinate transformation. 
 * It is derived from OGRCoordinateTransformation3D class and uses PROJ.4 
 * (modified) to enable inclusion of additional vertical data in 
 * Spatial reference system input.
 * 
/************************************************************************/

class CPL_DLL OGRProj4CT3D : public OGRCoordinateTransformation3D
{
	
	OGRSpatialReference3D *poSRSSource;
	projPJ        psPJSource;
    int         bSourceLatLong;
    double      dfSourceToRadians;
    double      dfSourceFromRadians;
    int         bSourceWrap;
    double      dfSourceWrapLong;

	OGRSpatialReference3D *poSRSTarget;
    projPJ        psPJTarget;
    int         bTargetLatLong;
    double      dfTargetToRadians;
    double      dfTargetFromRadians;
    int         bTargetWrap;
    double      dfTargetWrapLong;

	 int         nErrorCount;
    
    int         bCheckWithInvertProj;
    double      dfThreshold;
    
    projCtx     pjctx;

	int         InitializeNoLock( OGRSpatialReference3D *poSource, 
                                  OGRSpatialReference3D *poTarget );
  
    int         nMaxCount;
    double     *padfOriX;
    double     *padfOriY;
    double     *padfOriZ;
    double     *padfTargetX;
    double     *padfTargetY;
    double     *padfTargetZ;
public:
	OGRProj4CT3D();
	virtual OGRSpatialReference3D *GetSourceCS();
    virtual OGRSpatialReference3D *GetTargetCS();
	int         Initialize( OGRSpatialReference3D *poSource, 
                            OGRSpatialReference3D *poTarget );

	int ct3D_pj_transform(PJ *srcdefn, PJ *dstdefn, long point_count, int point_offset,
                  double *x, double *y, double *z);

	int ct3D_pj_datum_transform( PJ *srcdefn, PJ *dstdefn, 
                        long point_count, int point_offset,
                        double *x, double *y, double *z );
	
    virtual int Transform( int nCount, 
                           double *x, double *y, double *z = NULL );
    virtual int TransformEx( int nCount, 
                             double *x, double *y, double *z = NULL,
                             int *panSuccess = NULL );
};



OGRCoordinateTransformation* creat()
{
	OGRProj4CT3D *a;
	a=new OGRProj4CT3D();
	return a;
}

OGRCoordinateTransformation3D*  
OGRCreateCoordinateTransformation3D( OGRSpatialReference3D *poSource, 
                                   OGRSpatialReference3D *poTarget )
{
	OGRProj4CT3D *a;
	a = new OGRProj4CT3D();
	if( !a->Initialize( poSource, poTarget ) )
    {
        delete a;
        return NULL;
    }
    else
    {
        return a;
    }
}

OGRProj4CT3D::OGRProj4CT3D()
{
	poSRSSource = NULL;
    poSRSTarget = NULL;
    psPJSource = NULL;
    psPJTarget = NULL;
    
    nErrorCount = 0;
	bCheckWithInvertProj = FALSE;
    dfThreshold = 0;

    nMaxCount = 0;
    padfOriX = NULL;
    padfOriY = NULL;
    padfOriZ = NULL;
    padfTargetX = NULL;
    padfTargetY = NULL;
    padfTargetZ = NULL;

	pjctx=pj_ctx_alloc();
	
}

int OGRProj4CT3D::Initialize(OGRSpatialReference3D * poSourceIn, 
                            OGRSpatialReference3D * poTargetIn )
{
	if(pjctx!=NULL)
	{
        return InitializeNoLock(poSourceIn, poTargetIn);
    }

    CPLMutexHolderD( &hPROJMutex );
    return InitializeNoLock(poSourceIn, poTargetIn);
}

int OGRProj4CT3D::InitializeNoLock( OGRSpatialReference3D * poSourceIn, 
                                  OGRSpatialReference3D * poTargetIn )

{
	 if( poSourceIn == NULL || poTargetIn == NULL )
        return FALSE;

    poSRSSource = poSourceIn;//(OGRSpatialReference3D*)poSourceIn->Clone();
    poSRSTarget = poTargetIn;//(OGRSpatialReference3D*)poTargetIn->Clone();

	/*
	poSRSSource->SetGeoidModel("geoid.tif");
	poSRSSource->SetVCorrModel("vcorr.tif");
	poSRSSource->SetVOffset(100);
	poSRSSource->SetVScale(0.15);
	*/

    bSourceLatLong = poSRSSource->IsGeographic();
    bTargetLatLong = poSRSTarget->IsGeographic();

/* -------------------------------------------------------------------- */
/*      Setup source and target translations to radians for lat/long    */
/*      systems.                                                        */
/* -------------------------------------------------------------------- */
    dfSourceToRadians = DEG_TO_RAD;
    dfSourceFromRadians = RAD_TO_DEG;
    bSourceWrap = FALSE;
    dfSourceWrapLong = 0.0;

if( bSourceLatLong )
    {
        OGR_SRSNode *poUNITS = poSRSSource->GetAttrNode( "GEOGCS|UNIT" );
        if( poUNITS && poUNITS->GetChildCount() >= 2 )
        {
            dfSourceToRadians = atof(poUNITS->GetChild(1)->GetValue());
            if( dfSourceToRadians == 0.0 )
                dfSourceToRadians = DEG_TO_RAD;
            else
                dfSourceFromRadians = 1 / dfSourceToRadians;
        }
    }

    dfTargetToRadians = DEG_TO_RAD;
    dfTargetFromRadians = RAD_TO_DEG;
    bTargetWrap = FALSE;
    dfTargetWrapLong = 0.0;

    if( bTargetLatLong )
    {
        OGR_SRSNode *poUNITS = poSRSTarget->GetAttrNode( "GEOGCS|UNIT" );
        if( poUNITS && poUNITS->GetChildCount() >= 2 )
        {
            dfTargetToRadians = atof(poUNITS->GetChild(1)->GetValue());
            if( dfTargetToRadians == 0.0 )
                dfTargetToRadians = DEG_TO_RAD;
            else
                dfTargetFromRadians = 1 / dfTargetToRadians;
        }
    }

/* -------------------------------------------------------------------- */
/*      Preliminary logic to setup wrapping.                            */
/* -------------------------------------------------------------------- */
    const char *pszCENTER_LONG;

    if( CPLGetConfigOption( "CENTER_LONG", NULL ) != NULL )
    {
        bSourceWrap = bTargetWrap = TRUE;
        dfSourceWrapLong = dfTargetWrapLong = 
            atof(CPLGetConfigOption( "CENTER_LONG", "" ));
        CPLDebug( "OGRCT", "Wrap at %g.", dfSourceWrapLong );
    }

    pszCENTER_LONG = poSRSSource->GetExtension( "GEOGCS", "CENTER_LONG" );
    if( pszCENTER_LONG != NULL )
    {
        dfSourceWrapLong = atof(pszCENTER_LONG);
        bSourceWrap = TRUE;
        CPLDebug( "OGRCT", "Wrap source at %g.", dfSourceWrapLong );
    }

    pszCENTER_LONG = poSRSTarget->GetExtension( "GEOGCS", "CENTER_LONG" );
    if( pszCENTER_LONG != NULL )
    {
        dfTargetWrapLong = atof(pszCENTER_LONG);
        bTargetWrap = TRUE;
        CPLDebug( "OGRCT", "Wrap target at %g.", dfTargetWrapLong );
    }
    
    bCheckWithInvertProj = CSLTestBoolean(CPLGetConfigOption( "CHECK_WITH_INVERT_PROJ", "NO" ));
    
    /* The threshold is rather experimental... Works well with the cases of ticket #2305 */
    if (bSourceLatLong)
        dfThreshold = atof(CPLGetConfigOption( "THRESHOLD", ".1" ));
    else
        /* 1 works well for most projections, except for +proj=aeqd that requires */
        /* a tolerance of 10000 */
        dfThreshold = atof(CPLGetConfigOption( "THRESHOLD", "10000" ));

/* -------------------------------------------------------------------- */
/*      Establish PROJ.4 handle for source if projection.               */
/* -------------------------------------------------------------------- */
    // OGRThreadSafety: The following variable is not a thread safety issue 
    // since the only issue is incrementing while accessing which at worse 
    // means debug output could be one "increment" late. 
    static int   nDebugReportCount = 0;

    char        *pszProj4Defn = NULL;

   if( poSRSSource->exportToProj4( &pszProj4Defn ) != OGRERR_NONE )
    {
        CPLFree( pszProj4Defn );
        return FALSE;
    }

    if( strlen(pszProj4Defn) == 0 )
    {
        CPLFree( pszProj4Defn );
        CPLError( CE_Failure, CPLE_AppDefined, 
                  "No PROJ.4 translation for source SRS, coordinate\n"
                  "transformation initialization has failed." );
        return FALSE;
    }


	 if (pjctx)
		 psPJSource=pj_init_plus_ctx(pjctx,pszProj4Defn);
	  else
        psPJSource = pj_init_plus( pszProj4Defn );

	  if( psPJSource == NULL )
    {
        if( pjctx != NULL)
        {
            int pj_errno = pj_ctx_get_errno(pjctx);

            /* pfn_pj_strerrno not yet thread-safe in PROJ 4.8.0 */
            CPLMutexHolderD(&hPROJMutex);
            CPLError( CE_Failure, CPLE_NotSupported, 
                      "Failed to initialize PROJ.4 with `%s'.\n%s", 
                      pszProj4Defn, pj_strerrno(pj_errno) );
        }
        else if( pj_get_errno_ref != NULL
            && pj_strerrno != NULL )
        {
            int *p_pj_errno = pj_get_errno_ref();

            CPLError( CE_Failure, CPLE_NotSupported, 
                      "Failed to initialize PROJ.4 with `%s'.\n%s", 
                      pszProj4Defn, pj_strerrno(*p_pj_errno) );
        }
        else
        {
            CPLError( CE_Failure, CPLE_NotSupported, 
                      "Failed to initialize PROJ.4 with `%s'.\n", 
                      pszProj4Defn );
        }
    }

   if( nDebugReportCount < 10 )
        CPLDebug( "OGRCT", "Source: %s", pszProj4Defn );
    
    CPLFree( pszProj4Defn );

    if( psPJSource == NULL )
        return FALSE;
/* -------------------------------------------------------------------- */
/*      Establish PROJ.4 handle for target if projection.               */
/* -------------------------------------------------------------------- */
    pszProj4Defn = NULL;

    if( poSRSTarget->exportToProj4( &pszProj4Defn ) != OGRERR_NONE )
    {
        CPLFree( pszProj4Defn );
        return FALSE;
    }

    if( strlen(pszProj4Defn) == 0 )
    {
        CPLFree( pszProj4Defn );
        CPLError( CE_Failure, CPLE_AppDefined, 
                  "No PROJ.4 translation for destination SRS, coordinate\n"
                  "transformation initialization has failed." );
        return FALSE;
    }

	if (pjctx)
        psPJTarget = pj_init_plus_ctx( pjctx, pszProj4Defn );
    else
        psPJTarget = pj_init_plus( pszProj4Defn );

    if( psPJTarget == NULL )
        CPLError( CE_Failure, CPLE_NotSupported, 
                  "Failed to initialize PROJ.4 with `%s'.", 
                  pszProj4Defn );
    
    if( nDebugReportCount < 10 )
    {
        CPLDebug( "OGRCT", "Target: %s", pszProj4Defn );
        nDebugReportCount++;
    }
	
	CPLFree( pszProj4Defn );
    
    if( psPJTarget == NULL )
        return FALSE;

    return TRUE;
}

OGRSpatialReference3D* OGRProj4CT3D::GetSourceCS()
{
    return poSRSSource;
}

OGRSpatialReference3D* OGRProj4CT3D::GetTargetCS()
{
	return poSRSTarget;
}
int OGRProj4CT3D::Transform( int nCount, double *x, double *y, double *z )
{
    int *pabSuccess = (int *) CPLMalloc(sizeof(int) * nCount );
    int bOverallSuccess, i;

    bOverallSuccess = TransformEx( nCount, x, y, z, pabSuccess );

	 for( i = 0; i < nCount; i++ )
    {
        if( !pabSuccess[i] )
        {
            bOverallSuccess = FALSE;
            break;
        }
    }

    CPLFree( pabSuccess );

    return bOverallSuccess;

}

int OGRProj4CT3D::TransformEx( int nCount, double *x, double *y, double *z,int *pabSuccess )
{
    
int   err, i;

/* -------------------------------------------------------------------- */
/*      Potentially transform to radians.                               */
/* -------------------------------------------------------------------- */
    if( bSourceLatLong )
    {
        if( bSourceWrap )
        {
            for( i = 0; i < nCount; i++ )
            {
                if( x[i] != HUGE_VAL && y[i] != HUGE_VAL )
                {
                    if( x[i] < dfSourceWrapLong - 180.0 )
                        x[i] += 360.0;
                    else if( x[i] > dfSourceWrapLong + 180 )
                        x[i] -= 360.0;
                }
            }
        }

        for( i = 0; i < nCount; i++ )
        {
            if( x[i] != HUGE_VAL )
            {
                x[i] *= dfSourceToRadians;
                y[i] *= dfSourceToRadians;
            }
        }
    }

/* -------------------------------------------------------------------- */
/*      Do the transformation using PROJ.4.                             */
/* -------------------------------------------------------------------- */
    if (pjctx == NULL)
    {
        /* The mutex has already been created */
        CPLAssert(hPROJMutex != NULL);
        CPLAcquireMutex(hPROJMutex, 1000.0);
    } if (bCheckWithInvertProj)
    {
        /* For some projections, we cannot detect if we are trying to reproject */
        /* coordinates outside the validity area of the projection. So let's do */
        /* the reverse reprojection and compare with the source coordinates */
        if (nCount > nMaxCount)
        {
            nMaxCount = nCount;
            padfOriX = (double*) CPLRealloc(padfOriX, sizeof(double)*nCount);
            padfOriY = (double*) CPLRealloc(padfOriY, sizeof(double)*nCount);
            padfOriZ = (double*) CPLRealloc(padfOriZ, sizeof(double)*nCount);
            padfTargetX = (double*) CPLRealloc(padfTargetX, sizeof(double)*nCount);
            padfTargetY = (double*) CPLRealloc(padfTargetY, sizeof(double)*nCount);
            padfTargetZ = (double*) CPLRealloc(padfTargetZ, sizeof(double)*nCount);
        }
        memcpy(padfOriX, x, sizeof(double)*nCount);
        memcpy(padfOriY, y, sizeof(double)*nCount);
        if (z)
        {
            memcpy(padfOriZ, z, sizeof(double)*nCount);
        }
        err = ct3D_pj_transform( psPJSource, psPJTarget, nCount, 1, x, y, z );
        if (err == 0)
        {
            memcpy(padfTargetX, x, sizeof(double)*nCount);
            memcpy(padfTargetY, y, sizeof(double)*nCount);
            if (z)
            {
                memcpy(padfTargetZ, z, sizeof(double)*nCount);
            }
            
            err = ct3D_pj_transform( psPJTarget, psPJSource , nCount, 1,
                                    padfTargetX, padfTargetY, (z) ? padfTargetZ : NULL);
            if (err == 0)
            {
                for( i = 0; i < nCount; i++ )
                {
                    if ( x[i] != HUGE_VAL && y[i] != HUGE_VAL &&
                        (fabs(padfTargetX[i] - padfOriX[i]) > dfThreshold ||
                         fabs(padfTargetY[i] - padfOriY[i]) > dfThreshold) )
                    {
                        x[i] = HUGE_VAL;
                        y[i] = HUGE_VAL;
                    }
                }
            }
        }
    }

	 else
     {
        err = ct3D_pj_transform( psPJSource, psPJTarget, nCount, 1, x, y, z );
     }

	/* -------------------------------------------------------------------- */
/*      Try to report an error through CPL.  Get proj.4 error string    */
/*      if possible.  Try to avoid reporting thousands of error         */
/*      ... supress further error reporting on this OGRProj4CT if we    */
/*      have already reported 20 errors.                                */
/* -------------------------------------------------------------------- */
    if( err != 0 )
    {
        if( pabSuccess )
            memset( pabSuccess, 0, sizeof(int) * nCount );

        if( ++nErrorCount < 20 )
        {
            if (pjctx != NULL)
                /* pfn_pj_strerrno not yet thread-safe in PROJ 4.8.0 */
                CPLAcquireMutex(hPROJMutex, 1000.0);

            const char *pszError = NULL;
            if( pj_strerrno != NULL )
                pszError = pj_strerrno( err );
            
            if( pszError == NULL )
                CPLError( CE_Failure, CPLE_AppDefined, 
                          "Reprojection failed, err = %d", 
                          err );
            else
                CPLError( CE_Failure, CPLE_AppDefined, "%s", pszError );

            if (pjctx != NULL)
                /* pfn_pj_strerrno not yet thread-safe in PROJ 4.8.0 */
                CPLReleaseMutex(hPROJMutex);
        }
        else if( nErrorCount == 20 )
        {
            CPLError( CE_Failure, CPLE_AppDefined, 
                      "Reprojection failed, err = %d, further errors will be supressed on the transform object.", 
                      err );
        }

        if (pjctx == NULL)
            CPLReleaseMutex(hPROJMutex);
        return FALSE;
    }

    if (pjctx == NULL)
        CPLReleaseMutex(hPROJMutex);

/* -------------------------------------------------------------------- */
/*      Potentially transform back to degrees.                          */
/* -------------------------------------------------------------------- */
    if( bTargetLatLong )
    {
        for( i = 0; i < nCount; i++ )
        {
            if( x[i] != HUGE_VAL && y[i] != HUGE_VAL )
            {
                x[i] *= dfTargetFromRadians;
                y[i] *= dfTargetFromRadians;
            }
        }

        if( bTargetWrap )
        {
            for( i = 0; i < nCount; i++ )
            {
                if( x[i] != HUGE_VAL && y[i] != HUGE_VAL )
                {
                    if( x[i] < dfTargetWrapLong - 180.0 )
                        x[i] += 360.0;
                    else if( x[i] > dfTargetWrapLong + 180 )
                        x[i] -= 360.0;
                }
            }
        }
    }

/* -------------------------------------------------------------------- */
/*      Establish error information if pabSuccess provided.             */
/* -------------------------------------------------------------------- */
    if( pabSuccess )
    {
        for( i = 0; i < nCount; i++ )
        {
            if( x[i] == HUGE_VAL || y[i] == HUGE_VAL )
                pabSuccess[i] = FALSE;
            else
                pabSuccess[i] = TRUE;
        }
    }

    return TRUE;
}


/************************************************************************/
/*                       pj_geocentic_to_wgs84()                        */
/************************************************************************/

int pj_geocentric_to_wgs84( PJ *defn, 
                            long point_count, int point_offset,
                            double *x, double *y, double *z )

{
    int       i;

    if( defn->datum_type == PJD_3PARAM )
    {
        for( i = 0; i < point_count; i++ )
        {
            long io = i * point_offset;
            
            if( x[io] == HUGE_VAL )
                continue;

            x[io] = x[io] + Dx_BF;
            y[io] = y[io] + Dy_BF;
            z[io] = z[io] + Dz_BF;
        }
    }
    else if( defn->datum_type == PJD_7PARAM )
    {
        for( i = 0; i < point_count; i++ )
        {
            long io = i * point_offset;
            double x_out, y_out, z_out;

            if( x[io] == HUGE_VAL )
                continue;

            x_out = M_BF*(       x[io] - Rz_BF*y[io] + Ry_BF*z[io]) + Dx_BF;
            y_out = M_BF*( Rz_BF*x[io] +       y[io] - Rx_BF*z[io]) + Dy_BF;
            z_out = M_BF*(-Ry_BF*x[io] + Rx_BF*y[io] +       z[io]) + Dz_BF;

            x[io] = x_out;
            y[io] = y_out;
            z[io] = z_out;
        }
    }

    return 0;
}

/************************************************************************/
/*                      pj_geocentic_from_wgs84()                       */
/************************************************************************/

int pj_geocentric_from_wgs84( PJ *defn, 
                              long point_count, int point_offset,
                              double *x, double *y, double *z )

{
    int       i;

    if( defn->datum_type == PJD_3PARAM )
    {
        for( i = 0; i < point_count; i++ )
        {
            long io = i * point_offset;

            if( x[io] == HUGE_VAL )
                continue;
            
            x[io] = x[io] - Dx_BF;
            y[io] = y[io] - Dy_BF;
            z[io] = z[io] - Dz_BF;
        }
    }
    else if( defn->datum_type == PJD_7PARAM )
    {
        for( i = 0; i < point_count; i++ )
        {
            long io = i * point_offset;
            double x_tmp, y_tmp, z_tmp;

            if( x[io] == HUGE_VAL )
                continue;

            x_tmp = (x[io] - Dx_BF) / M_BF;
            y_tmp = (y[io] - Dy_BF) / M_BF;
            z_tmp = (z[io] - Dz_BF) / M_BF;

            x[io] =        x_tmp + Rz_BF*y_tmp - Ry_BF*z_tmp;
            y[io] = -Rz_BF*x_tmp +       y_tmp + Rx_BF*z_tmp;
            z[io] =  Ry_BF*x_tmp - Rx_BF*y_tmp +       z_tmp;
        }
    }

    return 0;
}



/************************************************************************/
/*                           pj_adjust_axis()                           */
/*                                                                      */
/*      Normalize or de-normalized the x/y/z axes.  The normal form     */
/*      is "enu" (easting, northing, up).                               */
/************************************************************************/
static int pj_adjust_axis( projCtx ctx, 
                           const char *axis, int denormalize_flag,
                           long point_count, int point_offset, 
                           double *x, double *y, double *z )

{
    double x_in, y_in, z_in = 0.0;
    int i, i_axis;

    if( !denormalize_flag )
    {
        for( i = 0; i < point_count; i++ )
        {
            x_in = x[point_offset*i];
            y_in = y[point_offset*i];
            if( z )
                z_in = z[point_offset*i];
     
            for( i_axis = 0; i_axis < 3; i_axis++ )
            {
                double value;

                if( i_axis == 0 )
                    value = x_in;
                else if( i_axis == 1 )
                    value = y_in;
                else
                    value = z_in;
                
                switch( axis[i_axis] )
                {
                  case 'e':
                    x[point_offset*i] = value; break;
                  case 'w':
                    x[point_offset*i] = -value; break;
                  case 'n':
                    y[point_offset*i] = value; break;
                  case 's':
                    y[point_offset*i] = -value; break;
                  case 'u':
                    if( z ) z[point_offset*i] = value; break;
                  case 'd':
                    if( z ) z[point_offset*i] = -value; break;
                  default:
                    pj_ctx_set_errno( ctx, PJD_ERR_AXIS );
                    return PJD_ERR_AXIS;
                }
            } /* i_axis */
        } /* i (point) */
    }

    else /* denormalize */
    {
        for( i = 0; i < point_count; i++ )
        {
            x_in = x[point_offset*i];
            y_in = y[point_offset*i];
            if( z )
                z_in = z[point_offset*i];
     
            for( i_axis = 0; i_axis < 3; i_axis++ )
            {
                double *target;

                if( i_axis == 2 && z == NULL )
                    continue;

                if( i_axis == 0 )
                    target = x;
                else if( i_axis == 1 )
                    target = y;
                else
                    target = z;
                
                switch( axis[i_axis] )
                {
                  case 'e':
                    target[point_offset*i] = x_in; break;
                  case 'w':
                    target[point_offset*i] = -x_in; break;
                  case 'n':
                    target[point_offset*i] = y_in; break;
                  case 's':
                    target[point_offset*i] = -y_in; break;
                  case 'u':
                    target[point_offset*i] = z_in; break;
                  case 'd':
                    target[point_offset*i] = -z_in; break;
                  default:
                    pj_ctx_set_errno( ctx, PJD_ERR_AXIS );
                    return PJD_ERR_AXIS;
                }
            } /* i_axis */
        } /* i (point) */
    }
    
    return 0;
}

//pj_apply_gridshift_2 and pj_apply_gridshift_3 are internal function not exposed in DLL
//so it is copied here.

# define assert(exp)	((void)0)
void ct3D_pj_acquire_lock()
{
}
void ct3D_pj_release_lock()
{
}

LP
ct3D_nad_intr(LP t, struct CTABLE *ct) {
	LP val, frct;
	ILP indx;
	double m00, m10, m01, m11;
	FLP *f00, *f10, *f01, *f11;
	long index;
	int in;

	indx.lam = floor(t.u /= ct->del.u);
	indx.phi = floor(t.v /= ct->del.v);
	frct.u = t.u - indx.lam;
	frct.v = t.v - indx.phi;
	val.u = val.v = HUGE_VAL;
	if (indx.lam < 0) {
		if (indx.lam == -1 && frct.u > 0.99999999999) {
			++indx.lam;
			frct.u = 0.;
		} else
			return val;
	} else if ((in = indx.lam + 1) >= ct->lim.lam) {
		if (in == ct->lim.lam && frct.u < 1e-11) {
			--indx.lam;
			frct.u = 1.;
		} else
			return val;
	}
	if (indx.phi < 0) {
		if (indx.phi == -1 && frct.v > 0.99999999999) {
			++indx.phi;
			frct.v = 0.;
		} else
			return val;
	} else if ((in = indx.phi + 1) >= ct->lim.phi) {
		if (in == ct->lim.phi && frct.v < 1e-11) {
			--indx.phi;
			frct.v = 1.;
		} else
			return val;
	}
	index = indx.phi * ct->lim.lam + indx.lam;
	f00 = ct->cvs + index++;
	f10 = ct->cvs + index;
	index += ct->lim.lam;
	f11 = ct->cvs + index--;
	f01 = ct->cvs + index;
	m11 = m10 = frct.u;
	m00 = m01 = 1. - frct.u;
	m11 *= frct.v;
	m01 *= frct.v;
	frct.v = 1. - frct.v;
	m00 *= frct.v;
	m10 *= frct.v;
	val.u = m00 * f00->lam + m10 * f10->lam +
			  m01 * f01->lam + m11 * f11->lam;
	val.v = m00 * f00->phi + m10 * f10->phi +
			  m01 * f01->phi + m11 * f11->phi;
	return val;
}

#define MAX_TRY 9
#define TOL 1e-12
LP
ct3D_nad_cvt(LP in, int inverse, struct CTABLE *ct) {
	LP t, tb;

	if (in.u == HUGE_VAL)
		return in;
	/* normalize input to ll origin */
	tb = in;
	tb.u -= ct->ll.u;
	tb.v -= ct->ll.v;
	tb.u = adjlon(tb.u - PI) + PI;
	t = ct3D_nad_intr(tb, ct);
	if (inverse) {
		LP del, dif;
		int i = MAX_TRY;

		if (t.u == HUGE_VAL) return t;
		t.u = tb.u + t.u;
		t.v = tb.v - t.v;

		do {
			del = ct3D_nad_intr(t, ct);

                        /* This case used to return failure, but I have
                           changed it to return the first order approximation
                           of the inverse shift.  This avoids cases where the
                           grid shift *into* this grid came from another grid.
                           While we aren't returning optimally correct results
                           I feel a close result in this case is better than
                           no result.  NFW
                           To demonstrate use -112.5839956 49.4914451 against
                           the NTv2 grid shift file from Canada. */
			if (del.u == HUGE_VAL) 
                        {
                            if( getenv( "PROJ_DEBUG" ) != NULL )
                                fprintf( stderr, 
                                         "Inverse grid shift iteration failed, presumably at grid edge.\n"
                                         "Using first approximation.\n" );
                            /* return del */;
                            break;
                        }

			t.u -= dif.u = t.u - del.u - tb.u;
			t.v -= dif.v = t.v + del.v - tb.v;
		} while (i-- && fabs(dif.u) > TOL && fabs(dif.v) > TOL);
		if (i < 0) {
                    if( getenv( "PROJ_DEBUG" ) != NULL )
                        fprintf( stderr, 
                                 "Inverse grid shift iterator failed to converge.\n" );
                    t.u = t.v = HUGE_VAL;
                    return t;
		}
		in.u = adjlon(t.u + ct->ll.u);
		in.v = t.v + ct->ll.v;
	} else {
		if (t.u == HUGE_VAL)
			in = t;
		else {
			in.u -= t.u;
			in.v += t.v;
		}
	}
	return in;
}

/************************************************************************/
/*                          pj_gridinfo_load()                          */
/*                                                                      */
/*      This function is intended to implement delayed loading of       */
/*      the data contents of a grid file.  The header and related       */
/*      stuff are loaded by pj_gridinfo_init().                         */
/************************************************************************/
static int  byte_order_test = 1;
#define IS_LSB	(((unsigned char *) (&byte_order_test))[0] == 1)
static void ct3D_swap_words( unsigned char *data, int word_size, int word_count )

{
    int	word;

    for( word = 0; word < word_count; word++ )
    {
        int	i;
        
        for( i = 0; i < word_size/2; i++ )
        {
            int	t;
            
            t = data[i];
            data[i] = data[word_size-i-1];
            data[word_size-i-1] = t;
        }
        
        data += word_size;
    }
}

static const char *(*pj_finder)(const char *) = NULL;
static int path_count = 0;
static char **search_path = NULL;
static char * proj_lib_name =
#ifdef PROJ_LIB
PROJ_LIB;
#else
0;
#endif

/************************************************************************/
/*                            pj_open_lib()                             */
/************************************************************************/

FILE *
ct3D_pj_open_lib(projCtx ctx, char *name, char *mode) {
    char fname[MAX_PATH_FILENAME+1];
    const char *sysname;
    FILE *fid;
    int n = 0;
    int i;
#ifdef WIN32
    static const char dir_chars[] = "/\\";
#else
    static const char dir_chars[] = "/";
#endif

#ifndef _WIN32_WCE

    /* check if ~/name */
    if (*name == '~' && strchr(dir_chars,name[1]) )
        if ((sysname = getenv("HOME")) != NULL) {
            (void)strcpy(fname, sysname);
            fname[n = strlen(fname)] = DIR_CHAR;
            fname[++n] = '\0';
            (void)strcpy(fname+n, name + 1);
            sysname = fname;
        } else
            return NULL;

    /* or fixed path: /name, ./name or ../name */
    else if (strchr(dir_chars,*name)
             || (*name == '.' && strchr(dir_chars,name[1])) 
             || (!strncmp(name, "..", 2) && strchr(dir_chars,name[2]))
             || (name[1] == ':' && strchr(dir_chars,name[2])) )
        sysname = name;

    /* or try to use application provided file finder */
    else if( pj_finder != NULL && pj_finder( name ) != NULL )
        sysname = pj_finder( name );

    /* or is environment PROJ_LIB defined */
    else if ((sysname = getenv("PROJ_LIB")) || (sysname = proj_lib_name)) {
        (void)strcpy(fname, sysname);
        fname[n = strlen(fname)] = DIR_CHAR;
        fname[++n] = '\0';
        (void)strcpy(fname+n, name);
        sysname = fname;
    } else /* just try it bare bones */
        sysname = name;

    if ((fid = fopen(sysname, mode)) != NULL)
        errno = 0;

    /* If none of those work and we have a search path, try it */
    if (!fid && path_count > 0)
    {
        for (i = 0; fid == NULL && i < path_count; i++)
        {
            sprintf(fname, "%s%c%s", search_path[i], DIR_CHAR, name);
            sysname = fname;
            fid = fopen (sysname, mode);
        }
        if (fid)
            errno = 0;
    }

    if( ctx->last_errno == 0 && errno != 0 )
        pj_ctx_set_errno( ctx, errno );

    pj_log( ctx, PJ_LOG_DEBUG_MAJOR, 
            "pj_open_lib(%s): call fopen(%s) - %s\n",
            name, sysname,
            fid == NULL ? "failed" : "succeeded" );

    return(fid);
#else
    return NULL;
#endif /* _WIN32_WCE */
}

/************************************************************************/
/*                          nad_ctable_load()                           */
/*                                                                      */
/*      Load the data portion of a ctable formatted grid.               */
/************************************************************************/

int ct3D_nad_ctable_load( projCtx ctx, struct CTABLE *ct, FILE *fid )

{
    int  a_size;

    fseek( fid, sizeof(struct CTABLE), SEEK_SET );

    /* read all the actual shift values */
    a_size = ct->lim.lam * ct->lim.phi;
    ct->cvs = (FLP *) pj_malloc(sizeof(FLP) * a_size);
    if( ct->cvs == NULL 
        || fread(ct->cvs, sizeof(FLP), a_size, fid) != a_size )
    {
        pj_dalloc( ct->cvs );
        ct->cvs = NULL;

        pj_log( ctx, PJ_LOG_ERROR, 
                "ctable loading failed on fread() - binary incompatible?\n" );
        pj_ctx_set_errno( ctx, -38 );
        return 0;
    }

    return 1;
} 

/************************************************************************/
/*                          nad_ctable_init()                           */
/*                                                                      */
/*      Read the header portion of a "ctable" format grid.              */
/************************************************************************/

struct CTABLE *ct3D_nad_ctable_init( projCtx ctx, FILE * fid )
{
    struct CTABLE *ct;
    int		id_end;

    /* read the table header */
    ct = (struct CTABLE *) pj_malloc(sizeof(struct CTABLE));
    if( ct == NULL 
        || fread( ct, sizeof(struct CTABLE), 1, fid ) != 1 )
    {
        pj_ctx_set_errno( ctx, -38 );
        return NULL;
    }

    /* do some minimal validation to ensure the structure isn't corrupt */
    if( ct->lim.lam < 1 || ct->lim.lam > 100000 
        || ct->lim.phi < 1 || ct->lim.phi > 100000 )
    {
        pj_ctx_set_errno( ctx, -38 );
        return NULL;
    }
    
    /* trim white space and newlines off id */
    for( id_end = strlen(ct->id)-1; id_end > 0; id_end-- )
    {
        if( ct->id[id_end] == '\n' || ct->id[id_end] == ' ' )
            ct->id[id_end] = '\0';
        else
            break;
    }

    ct->cvs = NULL;

    return ct;
}

/************************************************************************/
/*                          nad_ctable2_load()                          */
/*                                                                      */
/*      Load the data portion of a ctable2 formatted grid.              */
/************************************************************************/

int ct3D_nad_ctable2_load( projCtx ctx, struct CTABLE *ct, FILE *fid )

{
    int  a_size;

    fseek( fid, 160, SEEK_SET );

    /* read all the actual shift values */
    a_size = ct->lim.lam * ct->lim.phi;
    ct->cvs = (FLP *) pj_malloc(sizeof(FLP) * a_size);
    if( ct->cvs == NULL 
        || fread(ct->cvs, sizeof(FLP), a_size, fid) != a_size )
    {
        pj_dalloc( ct->cvs );
        ct->cvs = NULL;

        if( getenv("PROJ_DEBUG") != NULL )
        {
            fprintf( stderr,
            "ctable2 loading failed on fread() - binary incompatible?\n" );
        }

        pj_ctx_set_errno( ctx, -38 );
        return 0;
    }

    if( !IS_LSB )
    {
        ct3D_swap_words( (unsigned char*)ct->cvs, 4, a_size * 2 );
    }

    return 1;
} 

/************************************************************************/
/*                          nad_ctable2_init()                          */
/*                                                                      */
/*      Read the header portion of a "ctable2" format grid.             */
/************************************************************************/

struct CTABLE *ct3D_nad_ctable2_init( projCtx ctx, FILE * fid )
{
    struct CTABLE *ct;
    int		id_end;
    char        header[160];

    if( fread( header, sizeof(header), 1, fid ) != 1 )
    {
        pj_ctx_set_errno( ctx, -38 );
        return NULL;
    }

    if( !IS_LSB )
    {
        ct3D_swap_words( (unsigned char*)&(header[0]) +  96, 8, 4 );
        ct3D_swap_words( (unsigned char*)&(header[0]) + 128, 4, 2 );
    }

    if( strncmp(header,"CTABLE V2",9) != 0 )
    {
        pj_log( ctx, PJ_LOG_ERROR, "ctable2 - wrong header!" );
        pj_ctx_set_errno( ctx, -38 );
        return NULL;
    }

    /* read the table header */
    ct = (struct CTABLE *) pj_malloc(sizeof(struct CTABLE));
    if( ct == NULL )
    {
        pj_ctx_set_errno( ctx, -38 );
        return NULL;
    }

    memcpy( ct->id,       header +  16, 80 );
    memcpy( &ct->ll.u,  header +  96, 8 );
    memcpy( &ct->ll.v,  header + 104, 8 );
    memcpy( &ct->del.u, header + 112, 8 );
    memcpy( &ct->del.v, header + 120, 8 );
    memcpy( &ct->lim.lam, header + 128, 4 );
    memcpy( &ct->lim.phi, header + 132, 4 );

    /* do some minimal validation to ensure the structure isn't corrupt */
    if( ct->lim.lam < 1 || ct->lim.lam > 100000 
        || ct->lim.phi < 1 || ct->lim.phi > 100000 )
    {
        pj_ctx_set_errno( ctx, -38 );
        return NULL;
    }
    
    /* trim white space and newlines off id */
    for( id_end = strlen(ct->id)-1; id_end > 0; id_end-- )
    {
        if( ct->id[id_end] == '\n' || ct->id[id_end] == ' ' )
            ct->id[id_end] = '\0';
        else
            break;
    }

    ct->cvs = NULL;

    return ct;
}

int ct3D_pj_gridinfo_load( projCtx ctx, PJ_GRIDINFO *gi )

{
    if( gi == NULL || gi->ct == NULL )
        return 0;

/* -------------------------------------------------------------------- */
/*      Original platform specific CTable format.                       */
/* -------------------------------------------------------------------- */
    if( strcmp(gi->format,"ctable") == 0 )
    {
        FILE *fid;
        int result;

        fid = ct3D_pj_open_lib( ctx, gi->filename, "rb" );
        
        if( fid == NULL )
        {
            pj_ctx_set_errno( ctx, -38 );
            return 0;
        }

        result = ct3D_nad_ctable_load( ctx, gi->ct, fid );

        fclose( fid );

        return result;
    }

/* -------------------------------------------------------------------- */
/*      CTable2 format.                                                 */
/* -------------------------------------------------------------------- */
    else if( strcmp(gi->format,"ctable2") == 0 )
    {
        FILE *fid;
        int result;

        fid = ct3D_pj_open_lib( ctx, gi->filename, "rb" );
        
        if( fid == NULL )
        {
            pj_ctx_set_errno( ctx, -38 );
            return 0;
        }

        result = ct3D_nad_ctable2_load( ctx, gi->ct, fid );

        fclose( fid );

        return result;
    }

/* -------------------------------------------------------------------- */
/*      NTv1 format.                                                    */
/*      We process one line at a time.  Note that the array storage     */
/*      direction (e-w) is different in the NTv1 file and what          */
/*      the CTABLE is supposed to have.  The phi/lam are also           */
/*      reversed, and we have to be aware of byte swapping.             */
/* -------------------------------------------------------------------- */
    else if( strcmp(gi->format,"ntv1") == 0 )
    {
        double	*row_buf;
        int	row;
        FILE *fid;

        fid = ct3D_pj_open_lib( ctx, gi->filename, "rb" );
        
        if( fid == NULL )
        {
            pj_ctx_set_errno( ctx, -38 );
            return 0;
        }

        fseek( fid, gi->grid_offset, SEEK_SET );

        row_buf = (double *) pj_malloc(gi->ct->lim.lam * sizeof(double) * 2);
        gi->ct->cvs = (FLP *) pj_malloc(gi->ct->lim.lam*gi->ct->lim.phi*sizeof(FLP));
        if( row_buf == NULL || gi->ct->cvs == NULL )
        {
            pj_ctx_set_errno( ctx, -38 );
            return 0;
        }
        
        for( row = 0; row < gi->ct->lim.phi; row++ )
        {
            int	    i;
            FLP     *cvs;
            double  *diff_seconds;

            if( fread( row_buf, sizeof(double), gi->ct->lim.lam * 2, fid ) 
                != 2 * gi->ct->lim.lam )
            {
                pj_dalloc( row_buf );
                pj_dalloc( gi->ct->cvs );
                pj_ctx_set_errno( ctx, -38 );
                return 0;
            }

            if( IS_LSB )
                ct3D_swap_words( (unsigned char *) row_buf, 8, gi->ct->lim.lam*2 );

            /* convert seconds to radians */
            diff_seconds = row_buf;

            for( i = 0; i < gi->ct->lim.lam; i++ )
            {
                cvs = gi->ct->cvs + (row) * gi->ct->lim.lam
                    + (gi->ct->lim.lam - i - 1);

                cvs->phi = *(diff_seconds++) * ((PI/180.0) / 3600.0);
                cvs->lam = *(diff_seconds++) * ((PI/180.0) / 3600.0);
            }
        }

        pj_dalloc( row_buf );

        fclose( fid );

        return 1;
    }

/* -------------------------------------------------------------------- */
/*      NTv2 format.                                                    */
/*      We process one line at a time.  Note that the array storage     */
/*      direction (e-w) is different in the NTv2 file and what          */
/*      the CTABLE is supposed to have.  The phi/lam are also           */
/*      reversed, and we have to be aware of byte swapping.             */
/* -------------------------------------------------------------------- */
    else if( strcmp(gi->format,"ntv2") == 0 )
    {
        float	*row_buf;
        int	row;
        FILE *fid;

        pj_log( ctx, PJ_LOG_DEBUG_MINOR, 
                "NTv2 - loading grid %s", gi->ct->id );

        fid = ct3D_pj_open_lib( ctx, gi->filename, "rb" );
        
        if( fid == NULL )
        {
            pj_ctx_set_errno( ctx, -38 );
            return 0;
        }

        fseek( fid, gi->grid_offset, SEEK_SET );

        row_buf = (float *) pj_malloc(gi->ct->lim.lam * sizeof(float) * 4);
        gi->ct->cvs = (FLP *) pj_malloc(gi->ct->lim.lam*gi->ct->lim.phi*sizeof(FLP));
        if( row_buf == NULL || gi->ct->cvs == NULL )
        {
            pj_ctx_set_errno( ctx, -38 );
            return 0;
        }
        
        for( row = 0; row < gi->ct->lim.phi; row++ )
        {
            int	    i;
            FLP     *cvs;
            float   *diff_seconds;

            if( fread( row_buf, sizeof(float), gi->ct->lim.lam*4, fid ) 
                != 4 * gi->ct->lim.lam )
            {
                pj_dalloc( row_buf );
                pj_dalloc( gi->ct->cvs );
                gi->ct->cvs = NULL;
                pj_ctx_set_errno( ctx, -38 );
                return 0;
            }

            if( !IS_LSB )
                ct3D_swap_words( (unsigned char *) row_buf, 4, 
                            gi->ct->lim.lam*4 );

            /* convert seconds to radians */
            diff_seconds = row_buf;

            for( i = 0; i < gi->ct->lim.lam; i++ )
            {
                cvs = gi->ct->cvs + (row) * gi->ct->lim.lam
                    + (gi->ct->lim.lam - i - 1);

                cvs->phi = *(diff_seconds++) * ((PI/180.0) / 3600.0);
                cvs->lam = *(diff_seconds++) * ((PI/180.0) / 3600.0);
                diff_seconds += 2; /* skip accuracy values */
            }
        }

        pj_dalloc( row_buf );

        fclose( fid );

        return 1;
    }

/* -------------------------------------------------------------------- */
/*      GTX format.                                                     */
/* -------------------------------------------------------------------- */
    else if( strcmp(gi->format,"gtx") == 0 )
    {
        int   words = gi->ct->lim.lam * gi->ct->lim.phi;
        FILE *fid;

        fid = ct3D_pj_open_lib( ctx, gi->filename, "rb" );
        
        if( fid == NULL )
        {
            pj_ctx_set_errno( ctx, -38 );
            return 0;
        }

        fseek( fid, gi->grid_offset, SEEK_SET );

        gi->ct->cvs = (FLP *) pj_malloc(words*sizeof(float));
        if( gi->ct->cvs == NULL )
        {
            pj_ctx_set_errno( ctx, -38 );
            return 0;
        }
        
        if( fread( gi->ct->cvs, sizeof(float), words, fid ) != words )
        {
            pj_dalloc( gi->ct->cvs );
            gi->ct->cvs = NULL;
            return 0;
        }

        if( IS_LSB )
            ct3D_swap_words( (unsigned char *) gi->ct->cvs, 4, words );

        fclose( fid );
        return 1;
    }

    else
    {
        return 0;
    }
}

/************************************************************************/
/*                        pj_apply_gridshift_3()                        */
/*                                                                      */
/*      This is the real workhorse, given a gridlist.                   */
/************************************************************************/

int ct3D_pj_apply_gridshift_3( projCtx ctx, PJ_GRIDINFO **tables, int grid_count,
                          int inverse, long point_count, int point_offset,
                          double *x, double *y, double *z )

{
    int  i;
    static int debug_count = 0;

    if( tables == NULL || grid_count == 0 )
    {
        pj_ctx_set_errno( ctx, -38);
        return -38;
    }

    ctx->last_errno = 0;

    for( i = 0; i < point_count; i++ )
    {
        long io = i * point_offset;
        LP   input, output;
        int  itable;

        input.v = y[io];
        input.u = x[io];
        output.v = HUGE_VAL;
        output.u = HUGE_VAL;

        /* keep trying till we find a table that works */
        for( itable = 0; itable < grid_count; itable++ )
        {
            PJ_GRIDINFO *gi = tables[itable];
            struct CTABLE *ct = gi->ct;
            double epsilon = (fabs(ct->del.v)+fabs(ct->del.u))/10000.0;

            /* skip tables that don't match our point at all.  */
            if( ct->ll.v - epsilon > input.v 
                || ct->ll.u - epsilon > input.u
                || (ct->ll.v + (ct->lim.phi-1) * ct->del.v + epsilon 
                    < input.v)
                || (ct->ll.u + (ct->lim.lam-1) * ct->del.u + epsilon 
                    < input.u) )
                continue;

            /* If we have child nodes, check to see if any of them apply. */
            if( gi->child != NULL )
            {
                PJ_GRIDINFO *child;

                for( child = gi->child; child != NULL; child = child->next )
                {
                    struct CTABLE *ct1 = child->ct;
                    double epsilon = 
                        (fabs(ct1->del.v)+fabs(ct1->del.u))/10000.0;

                    if( ct1->ll.v - epsilon > input.v 
                        || ct1->ll.u - epsilon > input.u
                        || (ct1->ll.v+(ct1->lim.phi-1)*ct1->del.v + epsilon 
                            < input.v)
                        || (ct1->ll.u+(ct1->lim.lam-1)*ct1->del.u + epsilon 
                            < input.u) )
                        continue;

                    break;
                }

                /* we found a more refined child node to use */
                if( child != NULL )
                {
                    gi = child;
                    ct = child->ct;
                }
            }

            /* load the grid shift info if we don't have it. */
            if( ct->cvs == NULL && !ct3D_pj_gridinfo_load( ctx, gi ) )
            {
                pj_ctx_set_errno( ctx, -38 );
                return -38;
            }
            
            output = ct3D_nad_cvt( input, inverse, ct );
            if( output.u != HUGE_VAL )
            {
                if( debug_count++ < 20 )
                    pj_log( ctx, PJ_LOG_DEBUG_MINOR,
                            "ct3D_pj_apply_gridshift(): used %s", ct->id );
                break;
            }
        }

        if( output.u == HUGE_VAL )
        {
            if( ctx->debug_level >= PJ_LOG_DEBUG_MAJOR )
            {
                pj_log( ctx, PJ_LOG_DEBUG_MAJOR,
                    "pj_apply_gridshift(): failed to find a grid shift table for\n"
                    "                      location (%.7fdW,%.7fdN)",
                    x[io] * RAD_TO_DEG, 
                    y[io] * RAD_TO_DEG );
                for( itable = 0; itable < grid_count; itable++ )
                {
                    PJ_GRIDINFO *gi = tables[itable];
                    if( itable == 0 )
                        pj_log( ctx, PJ_LOG_DEBUG_MAJOR,
                                "   tried: %s", gi->gridname );
                    else
                        pj_log( ctx, PJ_LOG_DEBUG_MAJOR,
                                ",%s", gi->gridname );
                }
            }

            /* 
             * We don't actually have any machinery currently to set the 
             * following macro, so this is mostly kept here to make it clear 
             * how we ought to operate if we wanted to make it super clear 
             * that an error has occured when points are outside our available
             * datum shift areas.  But if this is on, we will find that "low 
             * value" points on the fringes of some datasets will completely 
             * fail causing lots of problems when it is more or less ok to 
             * just not apply a datum shift.  So rather than deal with
             * that we just fallback to no shift. (see also bug #45).
             */
#ifdef ERR_GRID_AREA_TRANSIENT_SEVERE
            y[io] = HUGE_VAL;
            x[io] = HUGE_VAL;
#else
            /* leave x/y unshifted. */
#endif
        }
        else
        {
            y[io] = output.v;
            x[io] = output.u;
        }
    }

    return 0;
}

/************************************************************************/
/*                       pj_gridinfo_init_ntv1()                        */
/*                                                                      */
/*      Load an NTv1 style Canadian grid shift file.                    */
/************************************************************************/

static int ct3D_pj_gridinfo_init_ntv1( projCtx ctx, FILE * fid, PJ_GRIDINFO *gi )

{
    unsigned char header[176];
    struct CTABLE *ct;
    LP		ur;
    
    assert( sizeof(int) == 4 );
    assert( sizeof(double) == 8 );
    if( sizeof(int) != 4 || sizeof(double) != 8 )
    {
        pj_log( ctx, PJ_LOG_ERROR,
                 "basic types of inappropraiate size in nad_load_ntv1()" );
        pj_ctx_set_errno( ctx, -38 );
        return 0;
    }

/* -------------------------------------------------------------------- */
/*      Read the header.                                                */
/* -------------------------------------------------------------------- */
    if( fread( header, sizeof(header), 1, fid ) != 1 )
    {
        pj_ctx_set_errno( ctx, -38 );
        return 0;
    }

/* -------------------------------------------------------------------- */
/*      Regularize fields of interest.                                  */
/* -------------------------------------------------------------------- */
    if( IS_LSB )
    {
        ct3D_swap_words( header+8, 4, 1 );
        ct3D_swap_words( header+24, 8, 1 );
        ct3D_swap_words( header+40, 8, 1 );
        ct3D_swap_words( header+56, 8, 1 );
        ct3D_swap_words( header+72, 8, 1 );
        ct3D_swap_words( header+88, 8, 1 );
        ct3D_swap_words( header+104, 8, 1 );
    }

    if( *((int *) (header+8)) != 12 )
    {
        pj_log( ctx, PJ_LOG_ERROR, 
                "NTv1 grid shift file has wrong record count, corrupt?" );
        pj_ctx_set_errno( ctx, -38 );
        return 0;
    }

/* -------------------------------------------------------------------- */
/*      Fill in CTABLE structure.                                       */
/* -------------------------------------------------------------------- */
    ct = (struct CTABLE *) pj_malloc(sizeof(struct CTABLE));
    strcpy( ct->id, "NTv1 Grid Shift File" );

    ct->ll.u = - *((double *) (header+72));
    ct->ll.v = *((double *) (header+24));
    ur.u = - *((double *) (header+56));
    ur.v = *((double *) (header+40));
    ct->del.u = *((double *) (header+104));
    ct->del.v = *((double *) (header+88));
    ct->lim.lam = (int) (fabs(ur.u-ct->ll.u)/ct->del.u + 0.5) + 1;
    ct->lim.phi = (int) (fabs(ur.v-ct->ll.v)/ct->del.v + 0.5) + 1;

    pj_log( ctx, PJ_LOG_DEBUG_MINOR,
            "NTv1 %dx%d: LL=(%.9g,%.9g) UR=(%.9g,%.9g)",
            ct->lim.lam, ct->lim.phi,
            ct->ll.u, ct->ll.v, ur.u, ur.v );

    ct->ll.u *= DEG_TO_RAD;
    ct->ll.v *= DEG_TO_RAD;
    ct->del.u *= DEG_TO_RAD;
    ct->del.v *= DEG_TO_RAD;
    ct->cvs = NULL;

    gi->ct = ct;
    gi->grid_offset = ftell( fid );
    gi->format = "ntv1";

    return 1;
}

/************************************************************************/
/*                       pj_gridinfo_init_ntv2()                        */
/*                                                                      */
/*      Load a ntv2 (.gsb) file.                                        */
/************************************************************************/

static int ct3D_pj_gridinfo_init_ntv2( projCtx ctx, FILE *fid, PJ_GRIDINFO *gilist )

{
    unsigned char header[11*16];
    int num_subfiles, subfile;

    assert( sizeof(int) == 4 );
    assert( sizeof(double) == 8 );
    if( sizeof(int) != 4 || sizeof(double) != 8 )
    {
        pj_log( ctx, PJ_LOG_ERROR,
             "basic types of inappropraiate size in pj_gridinfo_init_ntv2()" );
        pj_ctx_set_errno( ctx, -38 );
        return 0;
    }

/* -------------------------------------------------------------------- */
/*      Read the overview header.                                       */
/* -------------------------------------------------------------------- */
    if( fread( header, sizeof(header), 1, fid ) != 1 )
    {
        pj_ctx_set_errno( ctx, -38 );
        return 0;
    }

/* -------------------------------------------------------------------- */
/*      Byte swap interesting fields if needed.                         */
/* -------------------------------------------------------------------- */
    if( !IS_LSB )
    {
        ct3D_swap_words( header+8, 4, 1 );
        ct3D_swap_words( header+8+16, 4, 1 );
        ct3D_swap_words( header+8+32, 4, 1 );
        ct3D_swap_words( header+8+7*16, 8, 1 );
        ct3D_swap_words( header+8+8*16, 8, 1 );
        ct3D_swap_words( header+8+9*16, 8, 1 );
        ct3D_swap_words( header+8+10*16, 8, 1 );
    }

/* -------------------------------------------------------------------- */
/*      Get the subfile count out ... all we really use for now.        */
/* -------------------------------------------------------------------- */
    memcpy( &num_subfiles, header+8+32, 4 );

/* ==================================================================== */
/*      Step through the subfiles, creating a PJ_GRIDINFO for each.     */
/* ==================================================================== */
    for( subfile = 0; subfile < num_subfiles; subfile++ )
    {
        struct CTABLE *ct;
        LP ur;
        int gs_count;
        PJ_GRIDINFO *gi;

/* -------------------------------------------------------------------- */
/*      Read header.                                                    */
/* -------------------------------------------------------------------- */
        if( fread( header, sizeof(header), 1, fid ) != 1 )
        {
            pj_ctx_set_errno( ctx, -38 );
            return 0;
        }

        if( strncmp((const char *) header,"SUB_NAME",8) != 0 )
        {
            pj_ctx_set_errno( ctx, -38 );
            return 0;
        }
        
/* -------------------------------------------------------------------- */
/*      Byte swap interesting fields if needed.                         */
/* -------------------------------------------------------------------- */
        if( !IS_LSB )
        {
            ct3D_swap_words( header+8+16*4, 8, 1 );
            ct3D_swap_words( header+8+16*5, 8, 1 );
            ct3D_swap_words( header+8+16*6, 8, 1 );
            ct3D_swap_words( header+8+16*7, 8, 1 );
            ct3D_swap_words( header+8+16*8, 8, 1 );
            ct3D_swap_words( header+8+16*9, 8, 1 );
            ct3D_swap_words( header+8+16*10, 4, 1 );
        }
        
/* -------------------------------------------------------------------- */
/*      Initialize a corresponding "ct" structure.                      */
/* -------------------------------------------------------------------- */
        ct = (struct CTABLE *) pj_malloc(sizeof(struct CTABLE));
        strncpy( ct->id, (const char *) header + 8, 8 );
        ct->id[8] = '\0';

        ct->ll.u = - *((double *) (header+7*16+8)); /* W_LONG */
        ct->ll.v = *((double *) (header+4*16+8));   /* S_LAT */

        ur.u = - *((double *) (header+6*16+8));     /* E_LONG */
        ur.v = *((double *) (header+5*16+8));       /* N_LAT */

        ct->del.u = *((double *) (header+9*16+8));
        ct->del.v = *((double *) (header+8*16+8));

        ct->lim.lam = (int) (fabs(ur.u-ct->ll.u)/ct->del.u + 0.5) + 1;
        ct->lim.phi = (int) (fabs(ur.v-ct->ll.v)/ct->del.v + 0.5) + 1;

        pj_log( ctx, PJ_LOG_DEBUG_MINOR,
                "NTv2 %s %dx%d: LL=(%.9g,%.9g) UR=(%.9g,%.9g)\n",
                ct->id, 
                ct->lim.lam, ct->lim.phi,
                ct->ll.u/3600.0, ct->ll.v/3600.0,
                ur.u/3600.0, ur.v/3600.0 );
        
        ct->ll.u *= DEG_TO_RAD/3600.0;
        ct->ll.v *= DEG_TO_RAD/3600.0;
        ct->del.u *= DEG_TO_RAD/3600.0;
        ct->del.v *= DEG_TO_RAD/3600.0;

        memcpy( &gs_count, header + 8 + 16*10, 4 );
        if( gs_count != ct->lim.lam * ct->lim.phi )
        {
            pj_log( ctx, PJ_LOG_ERROR,
                    "GS_COUNT(%d) does not match expected cells (%dx%d=%d)\n",
                    gs_count, ct->lim.lam, ct->lim.phi, 
                    ct->lim.lam * ct->lim.phi );
            pj_ctx_set_errno( ctx, -38 );
            return 0;
        }

        ct->cvs = NULL;

/* -------------------------------------------------------------------- */
/*      Create a new gridinfo for this if we aren't processing the      */
/*      1st subfile, and initialize our grid info.                      */
/* -------------------------------------------------------------------- */
        if( subfile == 0 )
            gi = gilist;
        else
        {
            gi = (PJ_GRIDINFO *) pj_malloc(sizeof(PJ_GRIDINFO));
            memset( gi, 0, sizeof(PJ_GRIDINFO) );
    
            gi->gridname = strdup( gilist->gridname );
            gi->filename = strdup( gilist->filename );
            gi->next = NULL;
        }

        gi->ct = ct;
        gi->format = "ntv2";
        gi->grid_offset = ftell( fid );

/* -------------------------------------------------------------------- */
/*      Attach to the correct list or sublist.                          */
/* -------------------------------------------------------------------- */
        if( strncmp((const char *)header+24,"NONE",4) == 0 )
        {
            if( gi != gilist )
            {
                PJ_GRIDINFO *lnk;

                for( lnk = gilist; lnk->next != NULL; lnk = lnk->next ) {}
                lnk->next = gi;
            }
        }

        else
        {
            PJ_GRIDINFO *lnk;
            PJ_GRIDINFO *gp = gilist;
            
            while( gp != NULL 
                   && strncmp(gp->ct->id,(const char*)header+24,8) != 0 )
                gp = gp->next;

            if( gp == NULL )
            {
                pj_log( ctx, PJ_LOG_ERROR,
                        "pj_gridinfo_init_ntv2(): "
                        "failed to find parent %8.8s for %s.\n", 
                        (const char *) header+24, gi->ct->id );

                for( lnk = gp; lnk->next != NULL; lnk = lnk->next ) {}
                lnk->next = gi;
            }
            else if( gp->child == NULL )
            {
                gp->child = gi;
            }
            else
            {
                for( lnk = gp->child; lnk->next != NULL; lnk = lnk->next ) {}
                lnk->next = gi;
            }
        }

/* -------------------------------------------------------------------- */
/*      Seek past the data.                                             */
/* -------------------------------------------------------------------- */
        fseek( fid, gs_count * 16, SEEK_CUR );
    }

    return 1;
}

/************************************************************************/
/*                       pj_gridinfo_init_gtx()                         */
/*                                                                      */
/*      Load a NOAA .gtx vertical datum shift file.                     */
/************************************************************************/

static int ct3D_pj_gridinfo_init_gtx( projCtx ctx, FILE * fid, PJ_GRIDINFO *gi )

{
    unsigned char header[40];
    struct CTABLE *ct;
    double      xorigin,yorigin,xstep,ystep;
    int         rows, columns;

    assert( sizeof(int) == 4 );
    assert( sizeof(double) == 8 );
    if( sizeof(int) != 4 || sizeof(double) != 8 )
    {
        pj_log( ctx, PJ_LOG_ERROR,
                "basic types of inappropraiate size in nad_load_gtx()" );
        pj_ctx_set_errno( ctx, -38 );
        return 0;
    }

/* -------------------------------------------------------------------- */
/*      Read the header.                                                */
/* -------------------------------------------------------------------- */
    if( fread( header, sizeof(header), 1, fid ) != 1 )
    {
        pj_ctx_set_errno( ctx, -38 );
        return 0;
    }

/* -------------------------------------------------------------------- */
/*      Regularize fields of interest and extract.                      */
/* -------------------------------------------------------------------- */
    if( IS_LSB )
    {
        ct3D_swap_words( header+0, 8, 4 );
        ct3D_swap_words( header+32, 4, 2 );
    }

    memcpy( &yorigin, header+0, 8 );
    memcpy( &xorigin, header+8, 8 );
    memcpy( &ystep, header+16, 8 );
    memcpy( &xstep, header+24, 8 );

    memcpy( &rows, header+32, 4 );
    memcpy( &columns, header+36, 4 );

    if( xorigin < -360 || xorigin > 360 
        || yorigin < -90 || yorigin > 90 )
    {
        pj_log( ctx, PJ_LOG_ERROR, 
                "gtx file header has invalid extents, corrupt?");
        pj_ctx_set_errno( ctx, -38 );
        return 0;
    }

/* -------------------------------------------------------------------- */
/*      Fill in CTABLE structure.                                       */
/* -------------------------------------------------------------------- */
    ct = (struct CTABLE *) pj_malloc(sizeof(struct CTABLE));
    strcpy( ct->id, "GTX Vertical Grid Shift File" );

    ct->ll.u = xorigin;
    ct->ll.v = yorigin;
    ct->del.u = xstep;
    ct->del.v = ystep;
    ct->lim.lam = columns;
    ct->lim.phi = rows;

    /* some GTX files come in 0-360 and we shift them back into the
       expected -180 to 180 range if possible.  This does not solve 
       problems with grids spanning the dateline. */
    if( ct->ll.u >= 180.0 )
        ct->ll.u -= 360.0;

    if( ct->ll.u >= 0.0 && ct->ll.u + ct->del.u * ct->lim.lam > 180.0 )
    {
        pj_log( ctx, PJ_LOG_DEBUG_MAJOR,
                "This GTX spans the dateline!  This will cause problems." );
    }

    pj_log( ctx, PJ_LOG_DEBUG_MINOR,
            "GTX %dx%d: LL=(%.9g,%.9g) UR=(%.9g,%.9g)",
            ct->lim.lam, ct->lim.phi,
            ct->ll.u, ct->ll.v, 
            ct->ll.u + (columns-1)*xstep, ct->ll.v + (rows-1)*ystep);
    
    ct->ll.u *= DEG_TO_RAD;
    ct->ll.v *= DEG_TO_RAD;
    ct->del.u *= DEG_TO_RAD;
    ct->del.v *= DEG_TO_RAD;
    ct->cvs = NULL;

    gi->ct = ct;
    gi->grid_offset = 40;
    gi->format = "gtx";

    return 1;
}

/************************************************************************/
/*                          pj_gridinfo_init()                          */
/*                                                                      */
/*      Open and parse header details from a datum gridshift file       */
/*      returning a list of PJ_GRIDINFOs for the grids in that          */
/*      file.  This superceeds use of nad_init() for modern             */
/*      applications.                                                   */
/************************************************************************/

PJ_GRIDINFO *ct3D_pj_gridinfo_init( projCtx ctx, const char *gridname )

{
    char 	fname[MAX_PATH_FILENAME+1];
    PJ_GRIDINFO *gilist;
    FILE 	*fp;
    char	header[160];

    errno = pj_errno = 0;
    ctx->last_errno = 0;

/* -------------------------------------------------------------------- */
/*      Initialize a GRIDINFO with stub info we would use if it         */
/*      cannot be loaded.                                               */
/* -------------------------------------------------------------------- */
    gilist = (PJ_GRIDINFO *) pj_malloc(sizeof(PJ_GRIDINFO));
    memset( gilist, 0, sizeof(PJ_GRIDINFO) );
    
    gilist->gridname = strdup( gridname );
    gilist->filename = NULL;
    gilist->format = "missing";
    gilist->grid_offset = 0;
    gilist->ct = NULL;
    gilist->next = NULL;

/* -------------------------------------------------------------------- */
/*      Open the file using the usual search rules.                     */
/* -------------------------------------------------------------------- */
    strcpy(fname, gridname);
    if (!(fp = ct3D_pj_open_lib(ctx, fname, "rb"))) {
        ctx->last_errno = 0; /* don't treat as a persistent error */
        return gilist;
    }

    gilist->filename = strdup(fname);
    
/* -------------------------------------------------------------------- */
/*      Load a header, to determine the file type.                      */
/* -------------------------------------------------------------------- */
    if( fread( header, sizeof(header), 1, fp ) != 1 )
    {
        fclose( fp );
        pj_ctx_set_errno( ctx, -38 );
        return gilist;
    }

    fseek( fp, SEEK_SET, 0 );

/* -------------------------------------------------------------------- */
/*      Determine file type.                                            */
/* -------------------------------------------------------------------- */
    if( strncmp(header + 0, "HEADER", 6) == 0 
        && strncmp(header + 96, "W GRID", 6) == 0 
        && strncmp(header + 144, "TO      NAD83   ", 16) == 0 )
    {
        ct3D_pj_gridinfo_init_ntv1( ctx, fp, gilist );
    }
    
    else if( strncmp(header + 0, "NUM_OREC", 8) == 0 
             && strncmp(header + 48, "GS_TYPE", 7) == 0 )
    {
        ct3D_pj_gridinfo_init_ntv2( ctx, fp, gilist );
    }

    else if( strlen(gridname) > 4 
             && (strcmp(gridname+strlen(gridname)-3,"gtx") == 0 
                 || strcmp(gridname+strlen(gridname)-3,"GTX") == 0) )
    {
        ct3D_pj_gridinfo_init_gtx( ctx, fp, gilist );
    }

    else if( strncmp(header + 0,"CTABLE V2",9) == 0 )
    {
        struct CTABLE *ct = ct3D_nad_ctable2_init( ctx, fp );

        gilist->format = "ctable2";
        gilist->ct = ct;

        pj_log( ctx, PJ_LOG_DEBUG_MAJOR, 
                "Ctable2 %s %dx%d: LL=(%.9g,%.9g) UR=(%.9g,%.9g)\n",
                ct->id, 
                ct->lim.lam, ct->lim.phi,
                ct->ll.u * RAD_TO_DEG, ct->ll.v * RAD_TO_DEG,
                (ct->ll.u + (ct->lim.lam-1)*ct->del.u) * RAD_TO_DEG, 
                (ct->ll.v + (ct->lim.phi-1)*ct->del.v) * RAD_TO_DEG );
    }

    else
    {
        struct CTABLE *ct = ct3D_nad_ctable_init( ctx, fp );

        gilist->format = "ctable";
        gilist->ct = ct;

        pj_log( ctx, PJ_LOG_DEBUG_MAJOR, 
                "Ctable %s %dx%d: LL=(%.9g,%.9g) UR=(%.9g,%.9g)\n",
                ct->id, 
                ct->lim.lam, ct->lim.phi,
                ct->ll.u * RAD_TO_DEG, ct->ll.v * RAD_TO_DEG,
                (ct->ll.u + (ct->lim.lam-1)*ct->del.u) * RAD_TO_DEG, 
                (ct->ll.v + (ct->lim.phi-1)*ct->del.v) * RAD_TO_DEG );
    }

    fclose(fp);

    return gilist;
}

static PJ_GRIDINFO *grid_list = NULL;

/************************************************************************/
/*                       pj_gridlist_merge_grid()                       */
/*                                                                      */
/*      Find/load the named gridfile and merge it into the              */
/*      last_nadgrids_list.                                             */
/************************************************************************/

static int ct3D_pj_gridlist_merge_gridfile( projCtx ctx, 
                                       const char *gridname,
                                       PJ_GRIDINFO ***p_gridlist,
                                       int *p_gridcount, 
                                       int *p_gridmax )

{
    int got_match=0;
    PJ_GRIDINFO *this_grid, *tail = NULL;

/* -------------------------------------------------------------------- */
/*      Try to find in the existing list of loaded grids.  Add all      */
/*      matching grids as with NTv2 we can get many grids from one      */
/*      file (one shared gridname).                                     */
/* -------------------------------------------------------------------- */
    for( this_grid = grid_list; this_grid != NULL; this_grid = this_grid->next)
    {
        if( strcmp(this_grid->gridname,gridname) == 0 )
        {
            got_match = 1;

            /* dont add to the list if it is invalid. */
            if( this_grid->ct == NULL )
                return 0;

            /* do we need to grow the list? */
            if( *p_gridcount >= *p_gridmax - 2 )
            {
                PJ_GRIDINFO **new_list;
                int new_max = *p_gridmax + 20;

                new_list = (PJ_GRIDINFO **) pj_malloc(sizeof(void*) * new_max);
                if( *p_gridlist != NULL )
                {
                    memcpy( new_list, *p_gridlist,
                            sizeof(void*) * (*p_gridmax) );
                    pj_dalloc( *p_gridlist );
                }

                *p_gridlist = new_list;
                *p_gridmax = new_max;
            }

            /* add to the list */
            (*p_gridlist)[(*p_gridcount)++] = this_grid;
            (*p_gridlist)[*p_gridcount] = NULL;
        }

        tail = this_grid;
    }

    if( got_match )
        return 1;

/* -------------------------------------------------------------------- */
/*      Try to load the named grid.                                     */
/* -------------------------------------------------------------------- */
    this_grid = ct3D_pj_gridinfo_init( ctx, gridname );

    if( this_grid == NULL )
    {
        /* we should get at least a stub grid with a missing "ct" member */
        assert( FALSE );
        return 0;
    }
    
    if( tail != NULL )
        tail->next = this_grid;
    else
        grid_list = this_grid;

/* -------------------------------------------------------------------- */
/*      Recurse to add the grid now that it is loaded.                  */
/* -------------------------------------------------------------------- */
    return ct3D_pj_gridlist_merge_gridfile( ctx, gridname, p_gridlist, 
                                       p_gridcount, p_gridmax );
}

/************************************************************************/
/*                     pj_gridlist_from_nadgrids()                      */
/*                                                                      */
/*      This functions loads the list of grids corresponding to a       */
/*      particular nadgrids string into a list, and returns it.  The    */
/*      list is kept around till a request is made with a different     */
/*      string in order to cut down on the string parsing cost, and     */
/*      the cost of building the list of tables each time.              */
/************************************************************************/

PJ_GRIDINFO **ct3D_pj_gridlist_from_nadgrids( projCtx ctx, const char *nadgrids, 
                                         int *grid_count)

{
    const char *s;
    PJ_GRIDINFO **gridlist = NULL;
    int grid_max = 0;

    pj_errno = 0;
    *grid_count = 0;

    ct3D_pj_acquire_lock();

/* -------------------------------------------------------------------- */
/*      Loop processing names out of nadgrids one at a time.            */
/* -------------------------------------------------------------------- */
    for( s = nadgrids; *s != '\0'; )
    {
        int   end_char;
        int   required = 1;
        char  name[128];

        if( *s == '@' )
        {
            required = 0;
            s++;
        }

        for( end_char = 0; 
             s[end_char] != '\0' && s[end_char] != ','; 
             end_char++ ) {}

        if( end_char >= sizeof(name) )
        {
            pj_ctx_set_errno( ctx, -38 );
            ct3D_pj_release_lock();
            return NULL;
        }
        
        strncpy( name, s, end_char );
        name[end_char] = '\0';

        s += end_char;
        if( *s == ',' )
            s++;

        if( !ct3D_pj_gridlist_merge_gridfile( ctx, name, &gridlist, grid_count, 
                                         &grid_max) 
            && required )
        {
            pj_ctx_set_errno( ctx, -38 );
            ct3D_pj_release_lock();
            return NULL;
        }
        else
            pj_errno = 0;
    }

    ct3D_pj_release_lock();

    return gridlist;
}


/************************************************************************/
/*                        pj_apply_gridshift_2()                        */
/*                                                                      */
/*      This implmentation takes uses the gridlist from a coordinate    */
/*      system definition.  If the gridlist has not yet been            */
/*      populated in the coordinate system definition we set it up      */
/*      now.                                                            */
/************************************************************************/

int ct3D_pj_apply_gridshift_2( PJ *defn, int inverse, 
                          long point_count, int point_offset,
                          double *x, double *y, double *z )

{
    if( defn->gridlist == NULL )
    {
        defn->gridlist = 
            ct3D_pj_gridlist_from_nadgrids( pj_get_ctx( defn ),
                                       pj_param(defn->ctx, defn->params,"snadgrids").s,
                                       &(defn->gridlist_count) );

        if( defn->gridlist == NULL || defn->gridlist_count == 0 )
            return defn->ctx->last_errno;
    }
     
    return ct3D_pj_apply_gridshift_3( pj_get_ctx( defn ),
                                 defn->gridlist, defn->gridlist_count, inverse, 
                                 point_count, point_offset, x, y, z );
}

int OGRProj4CT3D::ct3D_pj_datum_transform( PJ *srcdefn, PJ *dstdefn, 
                        long point_count, int point_offset,
                        double *x, double *y, double *z )
{
    double      src_a, src_es, dst_a, dst_es;
    int         z_is_temp = FALSE;

/* -------------------------------------------------------------------- */
/*      We cannot do any meaningful datum transformation if either      */
/*      the source or destination are of an unknown datum type          */
/*      (ie. only a +ellps declaration, no +datum).  This is new        */
/*      behavior for PROJ 4.6.0.                                        */
/* -------------------------------------------------------------------- */
    if( srcdefn->datum_type == PJD_UNKNOWN
        || dstdefn->datum_type == PJD_UNKNOWN )
        return 0;

/* -------------------------------------------------------------------- */
/*      Short cut if the datums are identical.                          */
/* -------------------------------------------------------------------- */
    if( pj_compare_datums( srcdefn, dstdefn ) )
        return 0;

    src_a = srcdefn->a_orig;
    src_es = srcdefn->es_orig;

    dst_a = dstdefn->a_orig;
    dst_es = dstdefn->es_orig;

/* -------------------------------------------------------------------- */
/*      Create a temporary Z array if one is not provided.              */
/* -------------------------------------------------------------------- */
    if( z == NULL )
    {
        int	bytes = sizeof(double) * point_count * point_offset;
        z = (double *) pj_malloc(bytes);
        memset( z, 0, bytes );
        z_is_temp = TRUE;
    }

#define CHECK_RETURN(defn) {if( defn->ctx->last_errno != 0 && (defn->ctx->last_errno > 0 || transient_error[-defn->ctx->last_errno] == 0) ) { if( z_is_temp ) pj_dalloc(z); return defn->ctx->last_errno; }}

/* -------------------------------------------------------------------- */
/*	If this datum requires grid shifts, then apply it to geodetic   */
/*      coordinates.                                                    */
/* -------------------------------------------------------------------- */
    if( srcdefn->datum_type == PJD_GRIDSHIFT )
    {
				//PEB:gsoc2014 remove gridshift in datum transform
        //ct3D_pj_apply_gridshift_2( srcdefn, 0, point_count, point_offset, x, y, z );
        CHECK_RETURN(srcdefn);

        src_a = SRS_WGS84_SEMIMAJOR;
        src_es = SRS_WGS84_ESQUARED;
    }

    if( dstdefn->datum_type == PJD_GRIDSHIFT )
    {
        dst_a = SRS_WGS84_SEMIMAJOR;
        dst_es = SRS_WGS84_ESQUARED;
    }

/* ==================================================================== */
/*      Do we need to go through geocentric coordinates?                */
/* ==================================================================== */
    if( src_es != dst_es || src_a != dst_a
        || srcdefn->datum_type == PJD_3PARAM 
        || srcdefn->datum_type == PJD_7PARAM
        || dstdefn->datum_type == PJD_3PARAM 
        || dstdefn->datum_type == PJD_7PARAM)
    {
/* -------------------------------------------------------------------- */
/*      Convert to geocentric coordinates.                              */
/* -------------------------------------------------------------------- */
        srcdefn->ctx->last_errno = 
            pj_geodetic_to_geocentric( src_a, src_es,
                                       point_count, point_offset, x, y, z );
        CHECK_RETURN(srcdefn);

/* -------------------------------------------------------------------- */
/*      Convert between datums.                                         */
/* -------------------------------------------------------------------- */
        if( srcdefn->datum_type == PJD_3PARAM 
            || srcdefn->datum_type == PJD_7PARAM )
        {
            pj_geocentric_to_wgs84( srcdefn, point_count, point_offset,x,y,z);
            CHECK_RETURN(srcdefn);
        }

        if( dstdefn->datum_type == PJD_3PARAM 
            || dstdefn->datum_type == PJD_7PARAM )
        {
            pj_geocentric_from_wgs84( dstdefn, point_count,point_offset,x,y,z);
            CHECK_RETURN(dstdefn);
        }

/* -------------------------------------------------------------------- */
/*      Convert back to geodetic coordinates.                           */
/* -------------------------------------------------------------------- */
        dstdefn->ctx->last_errno = 
            pj_geocentric_to_geodetic( dst_a, dst_es,
                                       point_count, point_offset, x, y, z );
        CHECK_RETURN(dstdefn);
    }

/* -------------------------------------------------------------------- */
/*      Apply grid shift to destination if required.                    */
/* -------------------------------------------------------------------- */
    if( dstdefn->datum_type == PJD_GRIDSHIFT )
    {
				//PEB:gsoc2014 remove datum transform
        //ct3D_pj_apply_gridshift_2( dstdefn, 1, point_count, point_offset, x, y, z );
        //CHECK_RETURN(dstdefn);
    }

    if( z_is_temp )
        pj_dalloc( z );

    return 0;
}

int OGRProj4CT3D::ct3D_pj_transform(PJ *srcdefn, PJ *dstdefn, long point_count, int point_offset,
                  double *x, double *y, double *z)
{
    long      i;
    int       err;

	//double x1,y1;
	
	//x1=*x;
	//y1=*y;
	int         z_is_temp = FALSE;

    srcdefn->ctx->last_errno = 0;
    dstdefn->ctx->last_errno = 0;


    if( point_offset == 0 )
        point_offset = 1;

/* -------------------------------------------------------------------- */
/*      Transform unusual input coordinate axis orientation to          */
/*      standard form if needed.                                        */
/* -------------------------------------------------------------------- */
    if( strcmp(srcdefn->axis,"enu") != 0 )
    {
        int err;

        err = pj_adjust_axis( srcdefn->ctx, srcdefn->axis, 
                              0, point_count, point_offset, x, y, z );
        if( err != 0 )
            return err;
    }

/* -------------------------------------------------------------------- */
/*      Transform Z to meters if it isn't already.                      */
/* -------------------------------------------------------------------- */
    if( srcdefn->vto_meter != 1.0 && z != NULL )
    {
        for( i = 0; i < point_count; i++ )
            z[point_offset*i] *= srcdefn->vto_meter;
    }

/* -------------------------------------------------------------------- */
/*      Transform geocentric source coordinates to lat/long.            */
/* -------------------------------------------------------------------- */
    if( srcdefn->is_geocent )
    {
        if( z == NULL )
        {
            pj_ctx_set_errno( pj_get_ctx(srcdefn), PJD_ERR_GEOCENTRIC);
            return PJD_ERR_GEOCENTRIC;
        }

        if( srcdefn->to_meter != 1.0 )
        {
            for( i = 0; i < point_count; i++ )
            {
                if( x[point_offset*i] != HUGE_VAL )
                {
                    x[point_offset*i] *= srcdefn->to_meter;
                    y[point_offset*i] *= srcdefn->to_meter;
                }
            }
        }

        err = pj_geocentric_to_geodetic( srcdefn->a_orig, srcdefn->es_orig,
                                         point_count, point_offset, 
                                         x, y, z );
        if( err != 0 )
            return err;
    }

/* -------------------------------------------------------------------- */
/*      Transform source points to lat/long, if they aren't             */
/*      already.                                                        */
/* -------------------------------------------------------------------- */
    else if( !srcdefn->is_latlong )
    {
        if( srcdefn->inv == NULL )
        {
            pj_ctx_set_errno( pj_get_ctx(srcdefn), -17 );
            pj_log( pj_get_ctx(srcdefn), PJ_LOG_ERROR, 
                    "pj_transform(): source projection not invertable" );
            return -17;
        }

        for( i = 0; i < point_count; i++ )
        {
            XY         projected_loc;
            LP	       geodetic_loc;

            projected_loc.u = x[point_offset*i];
            projected_loc.v = y[point_offset*i];

            if( projected_loc.u == HUGE_VAL )
                continue;

            geodetic_loc = pj_inv( projected_loc, srcdefn );
            if( srcdefn->ctx->last_errno != 0 )
            {
                if( (srcdefn->ctx->last_errno != 33 /*EDOM*/ 
                     && srcdefn->ctx->last_errno != 34 /*ERANGE*/ )
                    && (srcdefn->ctx->last_errno > 0 
                        || srcdefn->ctx->last_errno < -44 || point_count == 1
                        || transient_error[-srcdefn->ctx->last_errno] == 0 ) )
                    return srcdefn->ctx->last_errno;
                else
                {
                    geodetic_loc.u = HUGE_VAL;
                    geodetic_loc.v = HUGE_VAL;
                }
            }

            x[point_offset*i] = geodetic_loc.u;
            y[point_offset*i] = geodetic_loc.v;
        }
    }

/* -------------------------------------------------------------------- */
/*      But if they are already lat long, adjust for the prime          */
/*      meridian if there is one in effect.                             */
/* -------------------------------------------------------------------- */
    if( srcdefn->from_greenwich != 0.0 )
    {
        for( i = 0; i < point_count; i++ )
        {
            if( x[point_offset*i] != HUGE_VAL )
                x[point_offset*i] += srcdefn->from_greenwich;
        }
    }

/* -------------------------------------------------------------------- */
/*      Do we need to translate from geoid to ellipsoidal vertical      */
/*      datum?                                                          */
/* -------------------------------------------------------------------- */
	/*
    if( srcdefn->has_geoid_vgrids )	
    {
        if( pj_apply_vgridshift( srcdefn, "sgeoidgrids", 
                                 &(srcdefn->vgridlist_geoid), 
                                 &(srcdefn->vgridlist_geoid_count),
                                 0, point_count, point_offset, x, y, z ) != 0 )
            return pj_ctx_get_errno(srcdefn->ctx);
    }
	*/
		
	//PEB:gsoc2014 - do gridshift on horizontal coordinates
	//TODO: the checking is not correct, deduce info from WKT
	if( srcdefn->datum_type == PJD_GRIDSHIFT )
  {
      ct3D_pj_apply_gridshift_2( srcdefn, 0, point_count, point_offset, x, y, z );
      CHECK_RETURN(srcdefn);
  }

	//PEB:gsoc2013 - use GDAL raster for vertical shift
	if(z!=NULL && poSRSSource->HasVerticalModel())
	{
		poSRSSource->ApplyVerticalCorrection(0, point_count, x, y, z);
	}

/* -------------------------------------------------------------------- */
/*      Convert datums if needed, and possible.                         */
/* -------------------------------------------------------------------- */
    if( ct3D_pj_datum_transform( srcdefn, dstdefn, point_count, point_offset, 
                            x, y, z ) != 0 )
    {
        if( srcdefn->ctx->last_errno != 0 )
            return srcdefn->ctx->last_errno;
        else
            return dstdefn->ctx->last_errno;
    }
/* -------------------------------------------------------------------- */
/*      Do we need to translate from geoid to ellipsoidal vertical      */
/*      datum?                                                          */
/* -------------------------------------------------------------------- */
	//PEB:gsoc2013
	if(z != NULL && poSRSTarget->HasVerticalModel())
	{
		//x y z coordinates are in radian
		poSRSTarget->ApplyVerticalCorrection(1, point_count, x, y, z);
	}
	//PEB:gsoc2014
	//TODO: the checking may not be correct
	if( dstdefn->datum_type == PJD_GRIDSHIFT )
  {
      ct3D_pj_apply_gridshift_2( dstdefn, 1, point_count, point_offset, x, y, z );
      CHECK_RETURN(dstdefn);
  }
	/*
    if( dstdefn->has_geoid_vgrids )
    {
        if( pj_apply_vgridshift( dstdefn, "sgeoidgrids", 
                                 &(dstdefn->vgridlist_geoid), 
                                 &(dstdefn->vgridlist_geoid_count),
                                 1, point_count, point_offset, x, y, z ) != 0 )
            return dstdefn->ctx->last_errno;
    }
      */  
/* -------------------------------------------------------------------- */
/*      But if they are staying lat long, adjust for the prime          */
/*      meridian if there is one in effect.                             */
/* -------------------------------------------------------------------- */
    if( dstdefn->from_greenwich != 0.0 )
    {
        for( i = 0; i < point_count; i++ )
        {
            if( x[point_offset*i] != HUGE_VAL )
                x[point_offset*i] -= dstdefn->from_greenwich;
        }
    }

/* -------------------------------------------------------------------- */
/*      Transform destination latlong to geocentric if required.        */
/* -------------------------------------------------------------------- */
    if( dstdefn->is_geocent )
    {
        if( z == NULL )
        {
            pj_ctx_set_errno( dstdefn->ctx, PJD_ERR_GEOCENTRIC );
            return PJD_ERR_GEOCENTRIC;
        }

        pj_geodetic_to_geocentric( dstdefn->a_orig, dstdefn->es_orig,
                                   point_count, point_offset, x, y, z );

        if( dstdefn->fr_meter != 1.0 )
        {
            for( i = 0; i < point_count; i++ )
            {
                if( x[point_offset*i] != HUGE_VAL )
                {
                    x[point_offset*i] *= dstdefn->fr_meter;
                    y[point_offset*i] *= dstdefn->fr_meter;
                }
            }
        }
    }

/* -------------------------------------------------------------------- */
/*      Transform destination points to projection coordinates, if      */
/*      desired.                                                        */
/* -------------------------------------------------------------------- */
    else if( !dstdefn->is_latlong )
    {
        for( i = 0; i < point_count; i++ )
        {
            XY         projected_loc;
            LP	       geodetic_loc;

            geodetic_loc.u = x[point_offset*i];
            geodetic_loc.v = y[point_offset*i];

            if( geodetic_loc.u == HUGE_VAL )
                continue;

            projected_loc = pj_fwd( geodetic_loc, dstdefn );
            if( dstdefn->ctx->last_errno != 0 )
            {
                if( (dstdefn->ctx->last_errno != 33 /*EDOM*/ 
                     && dstdefn->ctx->last_errno != 34 /*ERANGE*/ )
                    && (dstdefn->ctx->last_errno > 0 
                        || dstdefn->ctx->last_errno < -44 || point_count == 1
                        || transient_error[-dstdefn->ctx->last_errno] == 0 ) )
                    return dstdefn->ctx->last_errno;
                else
                {
                    projected_loc.u = HUGE_VAL;
                    projected_loc.v = HUGE_VAL;
                }
            }

            x[point_offset*i] = projected_loc.u;
            y[point_offset*i] = projected_loc.v;
        }
    }

/* -------------------------------------------------------------------- */
/*      If a wrapping center other than 0 is provided, rewrap around    */
/*      the suggested center (for latlong coordinate systems only).     */
/* -------------------------------------------------------------------- */
    else if( dstdefn->is_latlong && dstdefn->is_long_wrap_set )
    {
        for( i = 0; i < point_count; i++ )
        {
            if( x[point_offset*i] == HUGE_VAL )
                continue;

            while( x[point_offset*i] < dstdefn->long_wrap_center - PI )
                x[point_offset*i] += TWOPI;
            while( x[point_offset*i] > dstdefn->long_wrap_center + PI )
                x[point_offset*i] -= TWOPI;
        }
    }

/* -------------------------------------------------------------------- */
/*      Transform Z from meters if needed.                              */
/* -------------------------------------------------------------------- */
    if( dstdefn->vto_meter != 1.0 && z != NULL )
    {
        for( i = 0; i < point_count; i++ )
            z[point_offset*i] *= dstdefn->vfr_meter;
    }

/* -------------------------------------------------------------------- */
/*      Transform normalized axes into unusual output coordinate axis   */
/*      orientation if needed.                                          */
/* -------------------------------------------------------------------- */
    if( strcmp(dstdefn->axis,"enu") != 0 )
    {
        int err;

        err = pj_adjust_axis( dstdefn->ctx, dstdefn->axis, 
                              1, point_count, point_offset, x, y, z );
        if( err != 0 )
            return err;
	}

 return 0;
}