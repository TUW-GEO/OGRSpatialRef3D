#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

#include "validate.h"
#include "cpl_conv.h"
#include "ogr_spatialref3D.h"

double *x_etrs, *y_etrs, *z_etrs;
double *lon_grs, *lat_grs;
double *hell_grs, *h_orth;
double *x_gebr, *y_gebr, *h_gebr;
double *und_bess, *und_grs;
double *lat_mgi, *lon_mgi;
double *ras_val, *h_grid;
int *ms; // meridian strip
int num_data;

double *hell_mgi;
double *x_mgi, *y_mgi, *z_mgi;

char buffer[1024];

using namespace std;

void compute_ellh_mgi();
void compute_geoc_mgi();

void val_init()
{
	x_etrs = (double*)CPLMalloc(sizeof(double)*MAX_DATA);
	y_etrs = (double*)CPLMalloc(sizeof(double)*MAX_DATA);
	z_etrs = (double*)CPLMalloc(sizeof(double)*MAX_DATA);

	lon_grs = (double*)CPLMalloc(sizeof(double)*MAX_DATA);
	lat_grs = (double*)CPLMalloc(sizeof(double)*MAX_DATA);

	hell_grs = (double*)CPLMalloc(sizeof(double)*MAX_DATA);
	h_orth = (double*)CPLMalloc(sizeof(double)*MAX_DATA);

	x_gebr = (double*)CPLMalloc(sizeof(double)*MAX_DATA);
	y_gebr = (double*)CPLMalloc(sizeof(double)*MAX_DATA);
	h_gebr = (double*)CPLMalloc(sizeof(double)*MAX_DATA);

	und_bess = (double*)CPLMalloc(sizeof(double)*MAX_DATA);
	und_grs = (double*)CPLMalloc(sizeof(double)*MAX_DATA);

	lon_mgi = (double*)CPLMalloc(sizeof(double)*MAX_DATA);
	lat_mgi = (double*)CPLMalloc(sizeof(double)*MAX_DATA);

	ras_val = (double*)CPLMalloc(sizeof(double)*MAX_DATA);
	h_grid = (double*)CPLMalloc(sizeof(double)*MAX_DATA);

	ms = (int*)CPLMalloc(sizeof(int)*MAX_DATA);

	hell_mgi = (double*)CPLMalloc(sizeof(double)*MAX_DATA);

	x_mgi = (double*)CPLMalloc(sizeof(double)*MAX_DATA);
	y_mgi = (double*)CPLMalloc(sizeof(double)*MAX_DATA);
	z_mgi = (double*)CPLMalloc(sizeof(double)*MAX_DATA);
}

void val_cleanup()
{
	CPLFree(x_etrs);
	CPLFree(y_etrs);
	CPLFree(z_etrs);

	CPLFree(lon_grs);
	CPLFree(lat_grs);

	CPLFree(hell_grs);
	CPLFree(h_orth);

	CPLFree(x_gebr);
	CPLFree(y_gebr);
	CPLFree(h_gebr);

	CPLFree(und_bess);
	CPLFree(und_grs);

	CPLFree(lon_mgi);
	CPLFree(lat_mgi);

	CPLFree(ras_val);
	CPLFree(h_grid);

	CPLFree(ms);
	
	CPLFree(hell_mgi);

	CPLFree(x_mgi);
	CPLFree(y_mgi);
	CPLFree(z_mgi);
}

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

void loadRefFile(string filename, int max_input)
{
	bool first_row = false;
	ifstream inFile;
	inFile.open(filename, ios::in);
	if (!inFile) {
		cerr << "Can't open input file " << filename << endl;
		exit(1);
	}
	else{
		cout << "reading file " << filename << endl;
	}
	num_data = 0;
	
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
		//cout << col << endl;

		ss >> col; x_gebr[num_data] = col;
		ss >> col; y_gebr[num_data] = col;
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
	inFile.close();

	filename = "mgi_LonLatEllH.xyz";
	inFile.open(filename, ios::in);
	if (!inFile) {
		cerr << "Can't open input file " << filename << endl;
		exit(1);
	}
	else{
		cout << "reading file " << filename << endl;
	}

	for(int i=0; i<num_data; ++i){
		getline(inFile, line);

		stringstream ss(line);
		ss >> lon_mgi[i] >> lat_mgi[i] >> hell_mgi[i];

		double col = 0.0;
	}
	inFile.close();

	compute_ellh_mgi();
	compute_geoc_mgi();
}

void compute_ellh_mgi()
{
	for(int row_number=0; row_number < num_data; row_number++)
	//{
		hell_mgi[row_number] = h_orth[row_number] + und_bess[row_number];
	//}
	//return;

	//compute X,Y in MGI GK
	OGRSpatialReference3D oSourceSRS, oTargetSRS_28, oTargetSRS_31, oTargetSRS_34;

	char *wkt1 = loadWktFile(GEOG_MGI);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(PROJ_MGI_28);
	oTargetSRS_28.importFromWkt3D(&(wkt2));

	wkt2 = loadWktFile(PROJ_MGI_31);
	oTargetSRS_31.importFromWkt3D(&(wkt2));

	wkt2 = loadWktFile(PROJ_MGI_34);
	oTargetSRS_34.importFromWkt3D(&(wkt2));

	OGRCoordinateTransformation3D *poCT_28 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS_28 );
	OGRCoordinateTransformation3D *poCT_31 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS_31 );
	OGRCoordinateTransformation3D *poCT_34 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS_34 );

	for(int row_number=0; row_number < num_data; row_number++)
		{
			x_gebr[row_number] = lon_mgi[row_number];
			y_gebr[row_number] = lat_mgi[row_number];
			//r2[row_number] = z_etrs[row_number];

			switch(ms[row_number])
			{
			case 28:
				poCT_28->Transform( 1, &(x_gebr[row_number]), &(y_gebr[row_number]), 0);//&(r2[row_number]));
				break;
			case 31:
				poCT_31->Transform( 1, &(x_gebr[row_number]), &(y_gebr[row_number]), 0);//&(r2[row_number]));
				break;
			case 34:
				poCT_34->Transform( 1, &(x_gebr[row_number]), &(y_gebr[row_number]), 0);//&(r2[row_number]));
				break;
			default:
				cerr << "invalid meridianstrip value" << ms[row_number] << endl;
			}

		}

	delete poCT_28;
	delete poCT_31;
	delete poCT_34;
}

void compute_geoc_mgi()
{
	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	char *wkt1 = loadWktFile(GEOG_MGI);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOC_MGI);
	oTargetSRS.importFromWkt3D(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		x_mgi[row_number] = lon_mgi[row_number];
		y_mgi[row_number] = lat_mgi[row_number];
		z_mgi[row_number] = hell_mgi[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, x_mgi, y_mgi ,z_mgi) )
	{
		//printf( "Transformation failed.\n" );
	}
	else
	{
		//printf( "Transformation successful.\n" );
		
	}

	delete poCT;
}