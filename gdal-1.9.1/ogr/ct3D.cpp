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

class OGRProj4CT3D : public OGRCoordinateTransformation3D
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
	a=new OGRProj4CT3D();
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

    poSRSSource = (OGRSpatialReference3D*)poSourceIn->Clone();
    poSRSTarget = (OGRSpatialReference3D*)poTargetIn->Clone();


	poSRSSource->SetGeoidModel("geoid.tif");
	poSRSSource->SetVCorrModel("vcorr.tif");
	poSRSSource->SetVOffset(100);
	poSRSSource->SetVScale(0.15);

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
        err =ct3D_pj_transform( psPJSource, psPJTarget, nCount, 1, x, y, z );
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


int OGRProj4CT3D::ct3D_pj_transform(PJ *srcdefn, PJ *dstdefn, long point_count, int point_offset,
                  double *x, double *y, double *z)
{
    long      i;
    int       err;

	double x1,y1,*z1;
	
	x1=*x;
	y1=*y;
	

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
    if( srcdefn->has_geoid_vgrids )	
    {
        if( pj_apply_vgridshift( srcdefn, "sgeoidgrids", 
                                 &(srcdefn->vgridlist_geoid), 
                                 &(srcdefn->vgridlist_geoid_count),
                                 0, point_count, point_offset, x, y, z ) != 0 )
            return pj_ctx_get_errno(srcdefn->ctx);
    }

/* -------------------------------------------------------------------- */
/*      Convert datums if needed, and possible.                         */
/* -------------------------------------------------------------------- */
    if( pj_datum_transform( srcdefn, dstdefn, point_count, point_offset, 
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
    if( dstdefn->has_geoid_vgrids )
    {
        if( pj_apply_vgridshift( dstdefn, "sgeoidgrids", 
                                 &(dstdefn->vgridlist_geoid), 
                                 &(dstdefn->vgridlist_geoid_count),
                                 1, point_count, point_offset, x, y, z ) != 0 )
            return dstdefn->ctx->last_errno;
    }
        

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
/*      Do we need to translate from geoid to ellipsoidal vertical      */
/*      datum?                                                          */
/* -------------------------------------------------------------------- */


// Bhargav::Replace this function with my interpolation function which will take x , y and Z as an argument and give us height
   /* if( srcdefn->has_geoid_vgrids )
    {
        if( pj_apply_vgridshift( srcdefn, "sgeoidgrids", 
                                 &(srcdefn->vgridlist_geoid), 
                                 &(srcdefn->vgridlist_geoid_count),
                                 0, point_count, point_offset, x, y, z ) != 0 )
            return pj_ctx_get_errno(srcdefn->ctx);
    }*/

	poSRSSource->vgridshift(x1,y1,z1);
/* -------------------------------------------------------------------- */
/*      Convert datums if needed, and possible.                         */
/* -------------------------------------------------------------------- */
    if( pj_datum_transform( srcdefn, dstdefn, point_count, point_offset, 
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
    if( dstdefn->has_geoid_vgrids )
    {
        if( pj_apply_vgridshift( dstdefn, "sgeoidgrids", 
                                 &(dstdefn->vgridlist_geoid), 
                                 &(dstdefn->vgridlist_geoid_count),
                                 1, point_count, point_offset, x, y, z ) != 0 )
            return dstdefn->ctx->last_errno;
    }

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

