#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>

#include "cpl_conv.h"
#include "ogr_spatialref3D.h"
#include "OptionParser.h"

#define MAX_DATA 10000

//given values
double x_etrs[MAX_DATA];
double y_etrs[MAX_DATA];
double z_etrs[MAX_DATA];

double lat_grs[MAX_DATA];
double lon_grs[MAX_DATA];

double hell_grs[MAX_DATA];
double h_orth[MAX_DATA];

double x_gebr[MAX_DATA];
double y_gebr[MAX_DATA];
double h_gebr[MAX_DATA];

double und_bess[MAX_DATA];
double und_grs[MAX_DATA];

double lat_mgi[MAX_DATA];
double lon_mgi[MAX_DATA];

double ras_val[MAX_DATA];
double h_grid[MAX_DATA];

//computed values
double hell_mgi[MAX_DATA];
double x_mgi[MAX_DATA];
double y_mgi[MAX_DATA];
double z_mgi[MAX_DATA];

int ms[MAX_DATA]; // meridian strip

// data from text
int num_data = 0;

using namespace std;

char buffer[1024];

#define GEOG_ETRS "etrs89.prj" // geo[G]raphic coord, [G]RS80
#define GEOC_ETRS "etrs89_geoc.prj" // [G]eo[C]entric coord, [G]RS80
#define GEOG_ETRS_ORTH "etrs89_ortho.prj" //orthometric height

#define GEOG_MGI "mgi.prj" // geo[G]raphic coord, [M]GI
#define GEOC_MGI "mgi_geocen.prj" // [G]eo[C]entric coord, [M]GI
#define GEOG_MGI_ORTH "etrs89_ortho.prj" //orthometric height

#define PROJ_MGI_ORTH "gk-m34_orthoH.prj" //orthometric height
#define PROJ_MGI_USEH "mgi_gk_ortho.prj" //orthometric height

char *loadWktFile(const char* sWktFilename){
	ifstream inFile;

	inFile.open(sWktFilename, ios::in);
	if (!inFile) {
	  cerr << "Can't open input file " << sWktFilename << endl;
	  exit(1);
	}
	memset(buffer, 0, 1024);

	inFile.seekg(0,ios::end);
    streampos	length = inFile.tellg();
    inFile.seekg(0,ios::beg);

	inFile.read(buffer,length);
	return buffer;
}

void compute_values()
{
	//compute ellipsoidal height of mgi from geocentric coord
	for (int row_number=0; row_number<num_data; ++row_number)
	{
		//compute per point
	}
}

/**********************************************************************
 * one-way transformation
 **********************************************************************/
void geog_etrs_ellh_to_geog_mgi_ellh()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOG) + Ellipsoidal Height" << endl;
	cout << "Target coord.: MGI (GEOG)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_ETRS);
	oSourceSRS.importFromWkt(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_MGI);
	oTargetSRS.importFromWkt(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lat_grs[row_number];
		r1[row_number] = lon_grs[row_number];
		r2[row_number] = hell_grs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		double err_min0 = HUGE_VAL;
		double err_max0 = -HUGE_VAL;
		double err_sum0 = 0.0;
		double err_ssq0 = 0.0;

		double err_min1 = HUGE_VAL;
		double err_max1 = -HUGE_VAL;
		double err_sum1 = 0.0;
		double err_ssq1 = 0.0;

		double err_min2 = HUGE_VAL;
		double err_max2 = -HUGE_VAL;
		double err_sum2 = 0.0;
		double err_ssq2 = 0.0;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			double err0 = fabs(r0[row_number]-lat_mgi[row_number]);
			double err1 = fabs(r1[row_number]-lon_mgi[row_number]);

			err_min0 = MIN(err_min0, err0);
			err_max0 = MAX(err_max0, err0);
			err_sum0 += err0;
			err_ssq0 += err0*err0;

			err_min1 = MIN(err_min1, err1);
			err_max1 = MAX(err_max1, err1);
			err_sum1 += err1;
			err_ssq1 += err1*err1;
		}

		cout << fixed;
		cout << setprecision(8);
		cout << "Error (axis 0) : " << endl; 
		cout << "\tmin : " << err_min0 << endl; 
		cout << "\tmax : " << err_max0 << endl; 
		cout << "\tavg : " << err_sum0 / (double)num_data << endl; 
		cout << setprecision(18);
		cout << "\tvar : " << ((double)num_data*err_ssq0 - err_sum0*err_sum0) / ((double)num_data*((double)num_data-1.0)) << endl; 

		cout << setprecision(8);
		cout << "Error (axis 1) : " << endl; 
		cout << "\tmin : " << err_min1 << endl; 
		cout << "\tmax : " << err_max1 << endl; 
		cout << "\tavg : " << err_sum1 / (double)num_data << endl; 
		cout << setprecision(18);
		cout << "\tvar : " << ((double)num_data*err_ssq1 - err_sum1*err_sum1) / ((double)num_data*((double)num_data-1.0)) << endl; 
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
}

void geog_etrs_ellh_to_geoc_etrs()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOG) + Ellipsoidal Height" << endl;
	cout << "Target coord.: ETRS89 (GEOC)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_ETRS);
	oSourceSRS.importFromWkt(&(wkt1));

	char *wkt2 = loadWktFile(GEOC_ETRS);
	oTargetSRS.importFromWkt(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lat_grs[row_number];
		r1[row_number] = lon_grs[row_number];
		r2[row_number] = hell_grs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		double err_min0 = HUGE_VAL;
		double err_max0 = -HUGE_VAL;
		double err_sum0 = 0.0;
		double err_ssq0 = 0.0;

		double err_min1 = HUGE_VAL;
		double err_max1 = -HUGE_VAL;
		double err_sum1 = 0.0;
		double err_ssq1 = 0.0;

		double err_min2 = HUGE_VAL;
		double err_max2 = -HUGE_VAL;
		double err_sum2 = 0.0;
		double err_ssq2 = 0.0;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			double err0 = fabs(r0[row_number]-x_etrs[row_number]);
			double err1 = fabs(r1[row_number]-y_etrs[row_number]);
			double err2 = fabs(r2[row_number]-z_etrs[row_number]);

			err_min0 = MIN(err_min0, err0);
			err_max0 = MAX(err_max0, err0);
			err_sum0 += err0;
			err_ssq0 += err0*err0;

			err_min1 = MIN(err_min1, err1);
			err_max1 = MAX(err_max1, err1);
			err_sum1 += err1;
			err_ssq1 += err1*err1;

			err_min2 = MIN(err_min2, err2);
			err_max2 = MAX(err_max2, err2);
			err_sum2 += err2;
			err_ssq2 += err2*err2;
		}

		cout << fixed;
		cout << setprecision(8);
		cout << "Error (axis 0) : " << endl; 
		cout << "\tmin : " << err_min0 << endl; 
		cout << "\tmax : " << err_max0 << endl; 
		cout << "\tavg : " << err_sum0 / (double)num_data << endl; 
		cout << setprecision(18);
		cout << "\tvar : " << ((double)num_data*err_ssq0 - err_sum0*err_sum0) / ((double)num_data*((double)num_data-1.0)) << endl; 

		cout << setprecision(8);
		cout << "Error (axis 1) : " << endl; 
		cout << "\tmin : " << err_min1 << endl; 
		cout << "\tmax : " << err_max1 << endl; 
		cout << "\tavg : " << err_sum1 / (double)num_data << endl; 
		cout << setprecision(18);
		cout << "\tvar : " << ((double)num_data*err_ssq1 - err_sum1*err_sum1) / ((double)num_data*((double)num_data-1.0)) << endl; 

		cout << setprecision(8);
		cout << "Error (axis 2) : " << endl; 
		cout << "\tmin : " << err_min2 << endl; 
		cout << "\tmax : " << err_max2 << endl; 
		cout << "\tavg : " << err_sum2 / (double)num_data << endl; 
		cout << setprecision(18);
		cout << "\tvar : " << ((double)num_data*err_ssq2 - err_sum2*err_sum2) / ((double)num_data*((double)num_data-1.0)) << endl; 
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
}

void geoc_etrs_to_geog_etrs_ellh()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOC)" << endl;
	cout << "Target coord.: ETRS89 (GEOG) + Ellipsoidal Height" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOC_ETRS);
	oSourceSRS.importFromWkt(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_ETRS);
	oTargetSRS.importFromWkt(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = x_etrs[row_number];
		r1[row_number] = y_etrs[row_number];
		r2[row_number] = z_etrs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		double err_min0 = HUGE_VAL;
		double err_max0 = -HUGE_VAL;
		double err_sum0 = 0.0;
		double err_ssq0 = 0.0;

		double err_min1 = HUGE_VAL;
		double err_max1 = -HUGE_VAL;
		double err_sum1 = 0.0;
		double err_ssq1 = 0.0;

		double err_min2 = HUGE_VAL;
		double err_max2 = -HUGE_VAL;
		double err_sum2 = 0.0;
		double err_ssq2 = 0.0;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			double err0 = fabs(r0[row_number]-lat_grs[row_number]);
			double err1 = fabs(r1[row_number]-lon_grs[row_number]);
			double err2 = fabs(r2[row_number]-hell_grs[row_number]);
			
			err_min0 = MIN(err_min0, err0);
			err_max0 = MAX(err_max0, err0);
			err_sum0 += err0;
			err_ssq0 += err0*err0;

			err_min1 = MIN(err_min1, err1);
			err_max1 = MAX(err_max1, err1);
			err_sum1 += err1;
			err_ssq1 += err1*err1;

			err_min2 = MIN(err_min2, err2);
			err_max2 = MAX(err_max2, err2);
			err_sum2 += err2;
			err_ssq2 += err2*err2;
		}

		cout << fixed;
		cout << setprecision(8);
		cout << "Error (axis 0) : " << endl; 
		cout << "\tmin : " << err_min0 << endl; 
		cout << "\tmax : " << err_max0 << endl; 
		cout << "\tavg : " << err_sum0 / (double)num_data << endl; 
		cout << setprecision(18);
		cout << "\tvar : " << ((double)num_data*err_ssq0 - err_sum0*err_sum0) / ((double)num_data*((double)num_data-1.0)) << endl; 

		cout << setprecision(8);
		cout << "Error (axis 1) : " << endl; 
		cout << "\tmin : " << err_min1 << endl; 
		cout << "\tmax : " << err_max1 << endl; 
		cout << "\tavg : " << err_sum1 / (double)num_data << endl; 
		cout << setprecision(18);
		cout << "\tvar : " << ((double)num_data*err_ssq1 - err_sum1*err_sum1) / ((double)num_data*((double)num_data-1.0)) << endl; 

		cout << setprecision(8);
		cout << "Error (axis 2) : " << endl; 
		cout << "\tmin : " << err_min2 << endl; 
		cout << "\tmax : " << err_max2 << endl; 
		cout << "\tavg : " << err_sum2 / (double)num_data << endl; 
		cout << setprecision(18);
		cout << "\tvar : " << ((double)num_data*err_ssq2 - err_sum2*err_sum2) / ((double)num_data*((double)num_data-1.0)) << endl; 
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
}

void geoc_etrs_to_geoc_mgi()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOC)" << endl;
	cout << "Target coord.: MGI (GEOC)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOC_ETRS);
	oSourceSRS.importFromWkt(&(wkt1));

	char *wkt2 = loadWktFile(GEOC_MGI);
	oTargetSRS.importFromWkt(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = x_etrs[row_number];
		r1[row_number] = y_etrs[row_number];
		r2[row_number] = z_etrs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		double err_min0 = HUGE_VAL;
		double err_max0 = -HUGE_VAL;
		double err_sum0 = 0.0;
		double err_ssq0 = 0.0;

		double err_min1 = HUGE_VAL;
		double err_max1 = -HUGE_VAL;
		double err_sum1 = 0.0;
		double err_ssq1 = 0.0;

		double err_min2 = HUGE_VAL;
		double err_max2 = -HUGE_VAL;
		double err_sum2 = 0.0;
		double err_ssq2 = 0.0;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			double err0 = fabs(r0[row_number]-x_mgi[row_number]);
			double err1 = fabs(r1[row_number]-y_mgi[row_number]);
			double err2 = fabs(r2[row_number]-z_mgi[row_number]);

			err_min0 = MIN(err_min0, err0);
			err_max0 = MAX(err_max0, err0);
			err_sum0 += err0;
			err_ssq0 += err0*err0;

			err_min1 = MIN(err_min1, err1);
			err_max1 = MAX(err_max1, err1);
			err_sum1 += err1;
			err_ssq1 += err1*err1;

			err_min2 = MIN(err_min2, err2);
			err_max2 = MAX(err_max2, err2);
			err_sum2 += err2;
			err_ssq2 += err2*err2;
		}

		cout << fixed;
		cout << setprecision(8);
		cout << "Error (axis 0) : " << endl; 
		cout << "\tmin : " << err_min0 << endl; 
		cout << "\tmax : " << err_max0 << endl; 
		cout << "\tavg : " << err_sum0 / (double)num_data << endl; 
		cout << setprecision(18);
		cout << "\tvar : " << ((double)num_data*err_ssq0 - err_sum0*err_sum0) / ((double)num_data*((double)num_data-1.0)) << endl; 

		cout << setprecision(8);
		cout << "Error (axis 1) : " << endl; 
		cout << "\tmin : " << err_min1 << endl; 
		cout << "\tmax : " << err_max1 << endl; 
		cout << "\tavg : " << err_sum1 / (double)num_data << endl; 
		cout << setprecision(18);
		cout << "\tvar : " << ((double)num_data*err_ssq1 - err_sum1*err_sum1) / ((double)num_data*((double)num_data-1.0)) << endl; 

		cout << setprecision(8);
		cout << "Error (axis 2) : " << endl; 
		cout << "\tmin : " << err_min2 << endl; 
		cout << "\tmax : " << err_max2 << endl; 
		cout << "\tavg : " << err_sum2 / (double)num_data << endl; 
		cout << setprecision(18);
		cout << "\tvar : " << ((double)num_data*err_ssq2 - err_sum2*err_sum2) / ((double)num_data*((double)num_data-1.0)) << endl; 
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
}

void geog_mgi_ellh_to_geog_etrs_ellh()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	//double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: MGI (GEOG)" << endl;
	cout << "Target coord.: ETRS89 (GEOG)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_MGI);
	oSourceSRS.importFromWkt(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_ETRS);
	oTargetSRS.importFromWkt(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lat_mgi[row_number];
		r1[row_number] = lon_mgi[row_number];
		//r2[row_number] = hell_grs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,0) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		double err_min0 = HUGE_VAL;
		double err_max0 = -HUGE_VAL;
		double err_sum0 = 0.0;
		double err_ssq0 = 0.0;

		double err_min1 = HUGE_VAL;
		double err_max1 = -HUGE_VAL;
		double err_sum1 = 0.0;
		double err_ssq1 = 0.0;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			double err0 = fabs(r0[row_number]-lat_grs[row_number]);
			double err1 = fabs(r1[row_number]-lon_grs[row_number]);

			err_min0 = MIN(err_min0, err0);
			err_max0 = MAX(err_max0, err0);
			err_sum0 += err0;
			err_ssq0 += err0*err0;

			err_min1 = MIN(err_min1, err1);
			err_max1 = MAX(err_max1, err1);
			err_sum1 += err1;
			err_ssq1 += err1*err1;
		}

		cout << fixed;
		cout << setprecision(8);
		cout << "Error (axis 0) : " << endl; 
		cout << "\tmin : " << err_min0 << endl; 
		cout << "\tmax : " << err_max0 << endl; 
		cout << "\tavg : " << err_sum0 / (double)num_data << endl; 
		cout << setprecision(18);
		cout << "\tvar : " << ((double)num_data*err_ssq0 - err_sum0*err_sum0) / ((double)num_data*((double)num_data-1.0)) << endl; 

		cout << setprecision(8);
		cout << "Error (axis 1) : " << endl; 
		cout << "\tmin : " << err_min1 << endl; 
		cout << "\tmax : " << err_max1 << endl; 
		cout << "\tavg : " << err_sum1 / (double)num_data << endl; 
		cout << setprecision(18);
		cout << "\tvar : " << ((double)num_data*err_ssq1 - err_sum1*err_sum1) / ((double)num_data*((double)num_data-1.0)) << endl; 
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	//CPLFree(r2);
}

void geoc_mgi_ellh_to_geog_etrs_ellh()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	//double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: MGI (GEOG)" << endl;
	cout << "Target coord.: ETRS89 (GEOG)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_MGI);
	oSourceSRS.importFromWkt(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_ETRS);
	oTargetSRS.importFromWkt(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lat_mgi[row_number];
		r1[row_number] = lon_mgi[row_number];
		//r2[row_number] = hell_grs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,0) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		double err_min0 = HUGE_VAL;
		double err_max0 = -HUGE_VAL;
		double err_sum0 = 0.0;
		double err_ssq0 = 0.0;

		double err_min1 = HUGE_VAL;
		double err_max1 = -HUGE_VAL;
		double err_sum1 = 0.0;
		double err_ssq1 = 0.0;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			double err0 = fabs(r0[row_number]-lat_grs[row_number]);
			double err1 = fabs(r1[row_number]-lon_grs[row_number]);

			err_min0 = MIN(err_min0, err0);
			err_max0 = MAX(err_max0, err0);
			err_sum0 += err0;
			err_ssq0 += err0*err0;

			err_min1 = MIN(err_min1, err1);
			err_max1 = MAX(err_max1, err1);
			err_sum1 += err1;
			err_ssq1 += err1*err1;
		}

		cout << fixed;
		cout << setprecision(8);
		cout << "Error (axis 0) : " << endl; 
		cout << "\tmin : " << err_min0 << endl; 
		cout << "\tmax : " << err_max0 << endl; 
		cout << "\tavg : " << err_sum0 / (double)num_data << endl; 
		cout << setprecision(18);
		cout << "\tvar : " << ((double)num_data*err_ssq0 - err_sum0*err_sum0) / ((double)num_data*((double)num_data-1.0)) << endl; 

		cout << setprecision(8);
		cout << "Error (axis 1) : " << endl; 
		cout << "\tmin : " << err_min1 << endl; 
		cout << "\tmax : " << err_max1 << endl; 
		cout << "\tavg : " << err_sum1 / (double)num_data << endl; 
		cout << setprecision(18);
		cout << "\tvar : " << ((double)num_data*err_ssq1 - err_sum1*err_sum1) / ((double)num_data*((double)num_data-1.0)) << endl; 
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	//CPLFree(r2);
}
/**********************************************************************
 * two-way transformation
 **********************************************************************/
void geog_etrs_ellh_to_geog_mgi_ellh_2()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T -> S ]------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOG) + Ellipsoidal Height" << endl;
	cout << "Target coord.: MGI (GEOG)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_ETRS);
	oSourceSRS.importFromWkt(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_MGI);
	oTargetSRS.importFromWkt(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	OGRCoordinateTransformation3D *poCT_inv = OGRCreateCoordinateTransformation3D( 
													&oTargetSRS, &oSourceSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lat_grs[row_number];
		r1[row_number] = lon_grs[row_number];
		r2[row_number] = hell_grs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "forward Transformation failed.\n" );
	else
	{
		printf( "forward Transformation successful.\n" );

		if( poCT_inv == NULL || !poCT_inv->Transform( num_data, r0, r1 ,r2) )
			printf( "inverse Transformation failed.\n" );
		else
		{
			printf( "inverse Transformation successful.\n" );

			double err_min0 = HUGE_VAL;
			double err_max0 = -HUGE_VAL;
			double err_sum0 = 0.0;
			double err_ssq0 = 0.0;

			double err_min1 = HUGE_VAL;
			double err_max1 = -HUGE_VAL;
			double err_sum1 = 0.0;
			double err_ssq1 = 0.0;

			double err_min2 = HUGE_VAL;
			double err_max2 = -HUGE_VAL;
			double err_sum2 = 0.0;
			double err_ssq2 = 0.0;

			for(int row_number=0; row_number < num_data; row_number++)
			{
				double err0 = fabs(r0[row_number]-lat_grs[row_number]);
				double err1 = fabs(r1[row_number]-lon_grs[row_number]);
				double err2 = fabs(r2[row_number]-hell_grs[row_number]);

				err_min0 = MIN(err_min0, err0);
				err_max0 = MAX(err_max0, err0);
				err_sum0 += err0;
				err_ssq0 += err0*err0;

				err_min1 = MIN(err_min1, err1);
				err_max1 = MAX(err_max1, err1);
				err_sum1 += err1;
				err_ssq1 += err1*err1;

				err_min2 = MIN(err_min2, err2);
				err_max2 = MAX(err_max2, err2);
				err_sum2 += err2;
				err_ssq2 += err2*err2;
			}

			cout << fixed;
			cout << setprecision(18);
			cout << "Error (axis 0) : " << endl; 
			cout << "\tmin : " << err_min0 << endl; 
			cout << "\tmax : " << err_max0 << endl; 
			cout << "\tavg : " << err_sum0 / (double)num_data << endl; 
			cout << "\tvar : " << ((double)num_data*err_ssq0 - err_sum0*err_sum0) / ((double)num_data*((double)num_data-1.0)) << endl; 

			cout << "Error (axis 1) : " << endl; 
			cout << "\tmin : " << err_min1 << endl; 
			cout << "\tmax : " << err_max1 << endl; 
			cout << "\tavg : " << err_sum1 / (double)num_data << endl; 
			cout << "\tvar : " << ((double)num_data*err_ssq1 - err_sum1*err_sum1) / ((double)num_data*((double)num_data-1.0)) << endl; 

			cout << "Error (axis 2) : " << endl; 
			cout << "\tmin : " << err_min2 << endl; 
			cout << "\tmax : " << err_max2 << endl; 
			cout << "\tavg : " << err_sum2 / (double)num_data << endl; 
			cout << "\tvar : " << ((double)num_data*err_ssq2 - err_sum2*err_sum2) / ((double)num_data*((double)num_data-1.0)) << endl; 
		}
	}

	delete poCT;
	delete poCT_inv;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
}

void geoc_etrs_to_geoc_mgi_2()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T -> S ]------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOC)" << endl;
	cout << "Target coord.: MGI (GEOC)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOC_ETRS);
	oSourceSRS.importFromWkt(&(wkt1));

	char *wkt2 = loadWktFile(GEOC_MGI);
	oTargetSRS.importFromWkt(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	OGRCoordinateTransformation3D *poCT_inv = OGRCreateCoordinateTransformation3D( 
													&oTargetSRS, &oSourceSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = x_etrs[row_number];
		r1[row_number] = y_etrs[row_number];
		r2[row_number] = z_etrs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "forward Transformation failed.\n" );
	else
	{
		printf( "forward Transformation successful.\n" );

		if( poCT_inv == NULL || !poCT_inv->Transform( num_data, r0, r1 ,r2) )
			printf( "inverse Transformation failed.\n" );
		else
		{
			printf( "inverse Transformation successful.\n" );

			double err_min0 = HUGE_VAL;
			double err_max0 = -HUGE_VAL;
			double err_sum0 = 0.0;
			double err_ssq0 = 0.0;

			double err_min1 = HUGE_VAL;
			double err_max1 = -HUGE_VAL;
			double err_sum1 = 0.0;
			double err_ssq1 = 0.0;

			double err_min2 = HUGE_VAL;
			double err_max2 = -HUGE_VAL;
			double err_sum2 = 0.0;
			double err_ssq2 = 0.0;

			for(int row_number=0; row_number < num_data; row_number++)
			{
				double err0 = fabs(r0[row_number]-x_etrs[row_number]);
				double err1 = fabs(r1[row_number]-y_etrs[row_number]);
				double err2 = fabs(r2[row_number]-z_etrs[row_number]);

				err_min0 = MIN(err_min0, err0);
				err_max0 = MAX(err_max0, err0);
				err_sum0 += err0;
				err_ssq0 += err0*err0;

				err_min1 = MIN(err_min1, err1);
				err_max1 = MAX(err_max1, err1);
				err_sum1 += err1;
				err_ssq1 += err1*err1;

				err_min2 = MIN(err_min2, err2);
				err_max2 = MAX(err_max2, err2);
				err_sum2 += err2;
				err_ssq2 += err2*err2;
			}

			cout << fixed;
			cout << setprecision(18);
			cout << "Error (axis 0) : " << endl; 
			cout << "\tmin : " << err_min0 << endl; 
			cout << "\tmax : " << err_max0 << endl; 
			cout << "\tavg : " << err_sum0 / (double)num_data << endl; 
			cout << "\tvar : " << ((double)num_data*err_ssq0 - err_sum0*err_sum0) / ((double)num_data*((double)num_data-1.0)) << endl; 

			cout << "Error (axis 1) : " << endl; 
			cout << "\tmin : " << err_min1 << endl; 
			cout << "\tmax : " << err_max1 << endl; 
			cout << "\tavg : " << err_sum1 / (double)num_data << endl; 
			cout << "\tvar : " << ((double)num_data*err_ssq1 - err_sum1*err_sum1) / ((double)num_data*((double)num_data-1.0)) << endl; 

			cout << "Error (axis 2) : " << endl; 
			cout << "\tmin : " << err_min2 << endl; 
			cout << "\tmax : " << err_max2 << endl; 
			cout << "\tavg : " << err_sum2 / (double)num_data << endl; 
			cout << "\tvar : " << ((double)num_data*err_ssq2 - err_sum2*err_sum2) / ((double)num_data*((double)num_data-1.0)) << endl; 
		}
	}

	delete poCT;
	delete poCT_inv;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
}

/**********************************************************************
 * main
 **********************************************************************/
int main(int argc, char *argv[])
{
	cout << "validation" << endl;
	
	optparse::OptionParser parser = optparse::OptionParser().description("OGRSpatialRef3D validation program");
	
	parser.add_option("-i", "--input-file").dest("input_file").help("set input reference coordinate data FILE").metavar("FILE");
	parser.add_option("-n", "--num-input").dest("num_input").help("number of input data N taken from sample file (-1 means all data in file)").metavar("N").set_default(-1);
	
	optparse::Values options = parser.parse_args(argc, argv);
	vector<string> args = parser.args();
	
	if ((options["input_file"].length() == 0)){
			cerr << "Input Reference Coordinate is not set." << endl;
			exit(1);
	}
	
	bool first_row = false;
	
	int max_input = atoi(options["num_input"].c_str());

	ifstream inFile;
	inFile.open(options["input_file"], ios::in);
	if (!inFile) {
		cerr << "Can't open input file " << options["input_file"] << endl;
		exit(1);
	}
	else{
		cout << "reading file " << options["input_file"] << endl;
	}
	
	string line="";
	while(!inFile.eof()){
		getline(inFile, line);
		if (!first_row){
			first_row = true;
			continue;// ignore first row
		}

		stringstream ss(line);

		double col = 0.0;  
		int icol = 0;

		//cout << line << endl;
		
		ss >> col; x_etrs[num_data] = col;
		ss >> col; y_etrs[num_data] = col;
		ss >> col; z_etrs[num_data] = col;

		ss >> col; lat_grs[num_data] = col;
		ss >> col; lon_grs[num_data] = col;

		ss >> col; hell_grs[num_data] = col;

		ss >> col; y_gebr[num_data] = col;
		ss >> col; x_gebr[num_data] = col;
		ss >> col; h_gebr[num_data] = col;

		ss >> col; und_bess[num_data] = col;
		ss >> col; und_grs[num_data] = col;

		ss >> icol; ms[num_data] = icol;

		ss >> col; h_orth[num_data] = col;

		ss >> col; lat_mgi[num_data] = col;
		ss >> col; lon_mgi[num_data] = col;

		ss >> col; ras_val[num_data] = col;
		ss >> col; h_grid[num_data] = col;

		/*
		cout << setprecision(18) << x_etrs[num_data] << " " << y_etrs[num_data] << " " << z_etrs[num_data] << endl;
		cout << lat_grs[num_data] << " " << lon_grs[num_data] << " " << hell_grs[num_data] << endl;
		cout << endl;
		*/

		num_data += 1;

		if(max_input != -1 && num_data == max_input)
			break;
	}

	cout << "# data point(s) : " << num_data << endl;

	compute_values();

	//- from etrs
	geog_etrs_ellh_to_geoc_etrs();
	geoc_etrs_to_geog_etrs_ellh();

	geog_etrs_ellh_to_geog_mgi_ellh();

	geog_etrs_ellh_to_geog_mgi_orth();

	geog_etrs_ellh_to_geog_mgi_ellh_2();

	geoc_etrs_to_geoc_mgi_2();

	//- from mgi
	//geoc_mgi_to_geog_etrs_ellh();

	geog_mgi_ellh_to_geog_etrs_ellh();
	
	//geog_mgi_ellh_to_geog_etrs_ellh_2();
	

	cout << "Press ENTER to exit.";
	cin.get();
	
	return 0;
}
