#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>

#include "cpl_conv.h"
#include "ogr_spatialref3D.h"
#include "OptionParser.h"

#define MAX_DATA 10000

class SummStat
{
public:
	double min;
	double max;
	double sum;
	double ssq;
	int n;
	SummStat() : min(HUGE_VAL), max(-HUGE_VAL), sum(0.0), ssq(0.0), n(0) {}
	void add(double value);
	void printout(int shortprec=8, int longprec=18);
};

void
	SummStat::add(double value)
{
	n += 1;
	min = MIN(min, value);
	max = MAX(max, value);
	sum += value;
	ssq += value*value;
}

void
	SummStat::printout(int shortprec, int longprec)
{
	std::cout << std::setprecision(shortprec);
	std::cout << "\tmin : " << min << std::endl; 
	std::cout << "\tmax : " << max << std::endl; 
	std::cout << "\tavg : " << sum / (double)n << std::endl; 
	std::cout << std::setprecision(longprec);
	std::cout << "\tvar : " << ((double)n*ssq - sum*sum) / (double)(n*(n-1)) << std::endl; 
}


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

#define GEOG_ETRS "etrs89.prj" 
#define GEOC_ETRS "etrs89_geoc.prj" 
#define GEOG_ETRS_ORTH "etrs89_ortho.prj" 

#define GEOG_MGI "mgi.prj" 
#define GEOC_MGI "mgi_geocen.prj" 
#define GEOG_MGI_ORTH "mgi_ortho.prj" 

#define PROJ_MGI_ORTH "gk-m34_orthoH.prj" //orthometric height
#define PROJ_MGI_USEH "gk-m34_orthoH_WienerNull.prj" //orthometric height with offset

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
//- geoc etrs

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
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOC_MGI);
	oTargetSRS.importFromWkt3D(&(wkt2));

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
		SummStat err0, err1, err2;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-x_mgi[row_number]));
			err1.add(fabs(r1[row_number]-y_mgi[row_number]));
			err2.add(fabs(r2[row_number]-z_mgi[row_number]));
		}

		cout << "Error (axis 0) : " << endl; 
		err0.printout();

		cout << "Error (axis 1) : " << endl; 
		err1.printout();

		cout << "Error (axis 2) : " << endl; 
		err2.printout();
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
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_ETRS);
	oTargetSRS.importFromWkt3D(&(wkt2));

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
		SummStat err0, err1, err2;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-lat_grs[row_number]));
			err1.add(fabs(r1[row_number]-lon_grs[row_number]));
			err2.add(fabs(r2[row_number]-hell_grs[row_number]));
		}

		cout << "Error (axis 0) : " << endl; 
		err0.printout();

		cout << "Error (axis 1) : " << endl; 
		err1.printout();

		cout << "Error (axis 2) : " << endl; 
		err2.printout();
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
}

//- geoc mgi

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
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_ETRS);
	oTargetSRS.importFromWkt3D(&(wkt2));

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
		SummStat err0, err1;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-lat_grs[row_number]));
			err1.add(fabs(r1[row_number]-lon_grs[row_number]));
		}

		cout << "Error (axis 0) : " << endl; 
		err0.printout();
		cout << "Error (axis 1) : " << endl; 
		err1.printout();
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	//CPLFree(r2);
}

//-geog etrs

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
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOC_ETRS);
	oTargetSRS.importFromWkt3D(&(wkt2));

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
		SummStat err0, err1, err2;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-x_etrs[row_number]));
			err1.add(fabs(r1[row_number]-y_etrs[row_number]));
			err2.add(fabs(r2[row_number]-z_etrs[row_number]));
		}

		cout << "Error (axis 0) : " << endl; 
		err0.printout();

		cout << "Error (axis 1) : " << endl; 
		err1.printout();

		cout << "Error (axis 2) : " << endl; 
		err2.printout();
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
}

void geog_etrs_ellh_to_geog_etrs_ortho()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r3 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOG) + Ellipsoidal Height" << endl;
	cout << "Target coord.: ETRS89 (GEOC)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_ETRS);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_ETRS_ORTH);
	oTargetSRS.importFromWkt3D(&(wkt2));

	oTargetSRS.SetDebug(true);
	oTargetSRS.SetDebugData(r3, 0);

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lon_grs[row_number];
		r1[row_number] = lat_grs[row_number];
		r2[row_number] = hell_grs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2, err3;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-lon_grs[row_number]));
			err1.add(fabs(r1[row_number]-lat_grs[row_number]));
			err2.add(fabs(r2[row_number]-h_orth[row_number]));
			err3.add(fabs(r3[row_number]-und_grs[row_number]));
		}

		cout << "Error (axis 0) : " << endl; err0.printout();
		cout << "Error (axis 1) : " << endl; err1.printout();
		cout << "Error (axis 2) : " << endl; err2.printout();
		cout << "Error (geoid undulation) : " << endl; err3.printout();
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
	CPLFree(r3);
}

void geog_etrs_orth_to_geog_mgi_orth()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOG) + Orthometric height" << endl;
	cout << "Target coord.: MGI (GEOG) + Orthometric height" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_ETRS_ORTH);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_MGI_ORTH);
	oTargetSRS.importFromWkt3D(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lon_grs[row_number];
		r1[row_number] = lat_grs[row_number];
		r2[row_number] = h_orth[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-lon_mgi[row_number]));
			err1.add(fabs(r1[row_number]-lat_mgi[row_number]));
			err2.add(fabs(r2[row_number]-h_orth[row_number]));
		}

		cout << "Error (axis 0) : " << endl; err0.printout();
		cout << "Error (axis 1) : " << endl; err1.printout();
		cout << "Error (axis 2) : " << endl; err2.printout();
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
}

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
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_MGI);
	oTargetSRS.importFromWkt3D(&(wkt2));

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
		SummStat err0, err1;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-lat_mgi[row_number]));
			err1.add(fabs(r1[row_number]-lon_mgi[row_number]));
		}

		cout << setprecision(8);
		cout << "Error (axis 0) : " << endl; 
		err0.printout();

		cout << setprecision(8);
		cout << "Error (axis 1) : " << endl; 
		err1.printout();
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
}

//- geog mgi

//- proj mgi
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
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_MGI);
	oTargetSRS.importFromWkt3D(&(wkt2));

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
			SummStat err0, err1, err2;

			for(int row_number=0; row_number < num_data; row_number++)
			{
				err0.add(fabs(r0[row_number]-lat_grs[row_number]));
				err1.add(fabs(r1[row_number]-lon_grs[row_number]));
				err2.add(fabs(r2[row_number]-hell_grs[row_number]));
			}

			cout << "Error (axis 0) : " << endl; 
			err0.printout();
			cout << "Error (axis 1) : " << endl; 
			err1.printout();
			cout << "Error (axis 2) : " << endl; 
			err2.printout();
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
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOC_MGI);
	oTargetSRS.importFromWkt3D(&(wkt2));

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
			SummStat err0, err1, err2;

			for(int row_number=0; row_number < num_data; row_number++)
			{
				err0.add(fabs(r0[row_number]-x_etrs[row_number]));
				err1.add(fabs(r1[row_number]-y_etrs[row_number]));
				err2.add(fabs(r2[row_number]-z_etrs[row_number]));
			}

			cout << "Error (axis 0) : " << endl; 
			err0.printout();
			cout << "Error (axis 1) : " << endl; 
			err1.printout();
			cout << "Error (axis 2) : " << endl; 
			err2.printout();
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
void validation_from_etrs()
{
	//geoc_etrs_to_geoc_mgi_2();

	//geoc_etrs_to_geog_etrs_ellh();

	//geog_etrs_ellh_to_geoc_etrs();

	geog_etrs_ellh_to_geog_etrs_ortho();

	//geog_etrs_ellh_to_geog_mgi_ellh();

	//geog_etrs_orth_to_geog_mgi_orth();

	//geog_etrs_ellh_to_geog_mgi_ellh_2();
	
}

void validation_from_mgi()
{
	//geoc_mgi_to_geog_etrs_ellh();

	//geog_mgi_ellh_to_geog_etrs_ellh();
	
	//geog_mgi_ellh_to_geog_etrs_ellh_2();
}

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

	cout << fixed;
	cout << "# data point(s) : " << num_data << endl;

	compute_values();

	validation_from_etrs();
	validation_from_mgi();
	
	cout << "Press ENTER to exit.";
	cin.get();
	
	return 0;
}
