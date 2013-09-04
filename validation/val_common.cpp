#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

#include "validate.h"
#include "cpl_conv.h"

double *x_etrs, *y_etrs, *z_etrs;
double *lon_grs, *lat_grs;
double *hell_grs, *h_orth;
double *x_gebr, *y_gebr, *h_gebr;
double *und_bess, *und_grs;
double *lat_mgi, *lon_mgi;
double *ras_val, *h_grid;
int *ms; // meridian strip
int num_data;

char buffer[1024];

using namespace std;

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
}