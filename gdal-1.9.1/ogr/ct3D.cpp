#include <iostream>

#include "ogr_spatialref.h"
#include "cpl_port.h"
#include "cpl_error.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include "cpl_multiproc.h"


#include "..\..\proj-4.8.0\src\projects.h"

class transformationabc : public OGRCoordinateTransformation
{
	
	 OGRSpatialReference *poSRSSource;
public:

	virtual OGRSpatialReference *GetSourceCS();
    virtual OGRSpatialReference *GetTargetCS();
    virtual int Transform( int nCount, 
                           double *x, double *y, double *z = NULL );
    virtual int TransformEx( int nCount, 
                             double *x, double *y, double *z = NULL,
                             int *panSuccess = NULL );
};

OGRCoordinateTransformation*
creat()
{
	transformationabc *a;
	a=new transformationabc();
	return a;
}

OGRSpatialReference* transformationabc::GetSourceCS()

{
    return poSRSSource;
}

OGRSpatialReference* transformationabc::GetTargetCS()
{
	return poSRSSource;
}
int transformationabc::Transform( int nCount, double *x, double *y, double *z )
{
	return 1;
}

int transformationabc::TransformEx( int nCount, double *x, double *y, double *z,int *pabSuccess )
{
	return 1;
}