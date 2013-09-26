#include "math.h"
#include "assert.h"
#include "interpolation.h"

#define EPSILON 1e-5

double bilinearInterpolation(double p[4], double dx, double dy, double dNoDataValue)
{
	double tl = p[0];
	double tr = p[1];

	double bl = p[2];
	double br = p[3];

	double dWeight = 0.0;
	double dSum = 0.0;
	double dNorm = 0.0;

	// Calculate contribution of each neighboring pixel
	// weighted by inverse rectangular area
	//
	// tl----+------tr
	// |     |      |
	// |    dy      |
	// |     |      |
	// +--dx-+------+
	// |     |      |
	// bl----+------br
	//
	// possible extension: code below can be refactored
	// by moving to separate function for computing
	// specific interpolation scheme (e.g. kriging)
		
	if (fabs(tl - dNoDataValue) > EPSILON)
	{
		// top left neighbor
		dWeight = (1.0 - dx) * (1.0 - dy);
		dSum += tl * dWeight;
		dNorm += dWeight;
	}
		
	if (fabs(tr - dNoDataValue) > EPSILON)
	{
		// top right neighbor
		dWeight = dx * (1.0 - dy);
		dSum += tr * dWeight;
		dNorm += dWeight;
	}
		
	if (fabs(bl - dNoDataValue) > EPSILON)
	{
		// bottom left neighbor
		dWeight = (1.0 - dx) * (dy);
		dSum += bl * dWeight;
		dNorm += dWeight;
	}
		
	if (fabs(br - dNoDataValue) > EPSILON)
	{
		// bottom right neighbor
		dWeight = (dx) * (dy);
		dSum += br * dWeight;
		dNorm += dWeight;
	}

	//cout << "(Unnormalized) H_interp " << dSum << endl;

	if(dNorm < EPSILON)  // No valid data available
		dSum = 0.0; 
	else if(dNorm < 1.0) // One or more pixels is no data
		dSum /= dNorm;

	return dSum;
}

/*
 * implementation from http://www.paulinternet.nl/?page=bicubic
 * TODO: handle NoDataValue or look in gdalwarp
 */

double cubicInterpolate (double p[4], double x) {
	return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}

double bicubicInterpolate (double p[4][4], double x, double y) {
	double arr[4];
	arr[0] = cubicInterpolate(p[0], y);
	arr[1] = cubicInterpolate(p[1], y);
	arr[2] = cubicInterpolate(p[2], y);
	arr[3] = cubicInterpolate(p[3], y);
	return cubicInterpolate(arr, x);
}

double tricubicInterpolate (double p[4][4][4], double x, double y, double z) {
	double arr[4];
	arr[0] = bicubicInterpolate(p[0], y, z);
	arr[1] = bicubicInterpolate(p[1], y, z);
	arr[2] = bicubicInterpolate(p[2], y, z);
	arr[3] = bicubicInterpolate(p[3], y, z);
	return cubicInterpolate(arr, x);
}

double nCubicInterpolate (int n, double* p, double coordinates[]) {
	assert(n > 0);
	if (n == 1) {
		return cubicInterpolate(p, *coordinates);
	}
	else {
		double arr[4];
		int skip = 1 << (n - 1) * 2;
		arr[0] = nCubicInterpolate(n - 1, p, coordinates + 1);
		arr[1] = nCubicInterpolate(n - 1, p + skip, coordinates + 1);
		arr[2] = nCubicInterpolate(n - 1, p + 2*skip, coordinates + 1);
		arr[3] = nCubicInterpolate(n - 1, p + 3*skip, coordinates + 1);
		return cubicInterpolate(arr, *coordinates);
	}
}