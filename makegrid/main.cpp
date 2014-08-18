/******************************************************************************
 *
 * Project:  NTv2 Grid generator test
 * Purpose:  creating gridshift from inhomogeneous to homogeneous coordinate
 * Authors:  Peb Ruswono Aryan, Gottfried Mandlburger, Johannes Otepka
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
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "OptionParser.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include "ogr_spatialref.h"

#define STORAD 4.84813681e-6
#define RADTOS 206264.806

/**\file main.cpp
 * This file is an intra-datum gridshift generator program 
 * it take a NTv2 gridshift file as input and creates another NTv2 gridshift where
 * it stores coordinate shift in the same datum but using another reference
 * coordinate to calculate shift because of inhomogenity in the source CRs
 * 
 * The command-line options for running this program are:
 *	
 *		-r | --ref-coord=WKT_FILE		: set WKT_FILE as reference coordinate system
 *										  description
 *	
 *		-i | --input-grid=GSB_FILE			: set GSB_FILE as input NTv2 gridshift file
 * 
 *		-o | --output-grid=GSB_FILE			: set GSB_FILE as output NTv2 gridshift file
 *	
 *		-s | --source-coord=FILE		: set WKT_FILE as source coordinate system
 *										  description
 *	
 *  
 *
 *
 ************************************************************************/

extern "C" {

//! subfile header structure (key name omitted)
struct t_subfile
{
	//! Name of gridshift subfile (space-padded char string of length 8)
	char name[9];
	//! Name of parent gridshift subfile (space-padded char string of length 8)
	char parent[9];
	//! Date created (space-padded char string of length 8 e.g. (DDMMYYYY))
	char created[9];
	//! Date last modified (space-padded char string of length 8 e.g. (DDMMYYYY))
	char updated[9];

	//! Lower Latitude Limit (64-bit double)
	double s_lat;
	//! Upper Latitude Limit (64-bit double)
	double n_lat;
	//! Lower Longitude Limit (64-bit double)
	double e_long;
	//! Upper Longitude Limit (64-bit double)
	double w_long;

	//! Latitude Increment (64-bit double)
	double lat_inc;
	//! Longitude Increment (64-bit double)
	double long_inc;

	//! Number of Gridshift point (64-bit double)
	int gs_count;
};

//! NTv2 header (overview record)
struct t_header
{
	//! Number of Overview Record (int32)
	int num_orec;
	//! Number of Subfile Record (int32)
	int num_srec;
	//! Number of Subfiles (int32)
	int num_file;

	//! string indicating Unit for storing the gridshift and lat-lon (space-padded char string of length 8 e.g. "SECONDS ")
	char gs_type[9];

	//! Version of gridshift file (space-padded char string of length 8 e.g. "NTv2.0  ")
	char version[9];
	//! Source Coordinate System  (space-padded char string of length 8 e.g. " MGI    ")
	char system_f[9];
	//! Target Coordinate System  (space-padded char string of length 8 e.g. " ETRS89 ")
	char system_t[9];

	//! Length of Semi major axis in the source CRS (64-bit double)
	double major_f;
	//! Length of Semi minor axis in the source CRS (64-bit double)
	double minor_f;
	//! Length of Semi major axis in the target CRS (64-bit double)
	double major_t;
	//! Length of Semi minor axis in the target CRS (64-bit double)
	double minor_t;

	t_subfile *subfiles;
};

char hbuf[9] = "";

void read_key(FILE *f, char *key)
{
	fread(key, sizeof(char), 8, f);
	key[8] = '\0';
}

void swap_word(void *data, size_t bytecount)
{
	for(unsigned i=0; i<bytecount/2; ++i)
	{
		unsigned char tmp = *((unsigned char*)data+i);
		*((unsigned char*)data+i) = *((unsigned char*)data+bytecount-1-i);
		*((unsigned char*)data+bytecount-1-i) = tmp;
	}
}

void read_int(FILE *f, int *value)
{
	fread(value, sizeof(int), 1, f);
}

void read_double(FILE *f, double *value)
{
	fread(value, sizeof(double), 1, f);
}

//! procedure to read NTv2 header from a file handle
/*!
    \param f a pointer to FILE structure.
		\param orec a pointer to t_header structure (NTv2 header) to store the header info.
*/
void read_header(FILE *f, t_header *orec)
{
	int itmp;
	
	read_key(f, hbuf); read_int(f, &(orec->num_orec)); read_int(f, &itmp);
	read_key(f, hbuf); read_int(f, &(orec->num_srec)); read_int(f, &itmp);
	read_key(f, hbuf); read_int(f, &(orec->num_file)); read_int(f, &itmp);

	read_key(f, hbuf); read_key(f, orec->gs_type);
	read_key(f, hbuf); read_key(f, orec->version);
	read_key(f, hbuf); read_key(f, orec->system_f);
	read_key(f, hbuf); read_key(f, orec->system_t);

	read_key(f, hbuf); read_double(f, &(orec->major_f));
	read_key(f, hbuf); read_double(f, &(orec->minor_f));
	read_key(f, hbuf); read_double(f, &(orec->major_t));
	read_key(f, hbuf); read_double(f, &(orec->minor_t));

	orec->subfiles = (t_subfile*)CPLMalloc(orec->num_file*sizeof(t_subfile));
	for(int fi=0; fi<orec->num_file; ++fi)
	{
		read_key(f, hbuf); read_key(f, orec->subfiles[fi].name);
		read_key(f, hbuf); read_key(f, orec->subfiles[fi].parent);
		read_key(f, hbuf); read_key(f, orec->subfiles[fi].created);
		read_key(f, hbuf); read_key(f, orec->subfiles[fi].updated);

		read_key(f, hbuf); read_double(f, &(orec->subfiles[fi].s_lat));
		read_key(f, hbuf); read_double(f, &(orec->subfiles[fi].n_lat));
		read_key(f, hbuf); read_double(f, &(orec->subfiles[fi].e_long));
		read_key(f, hbuf); read_double(f, &(orec->subfiles[fi].w_long));

		read_key(f, hbuf); read_double(f, &(orec->subfiles[fi].lat_inc));
		read_key(f, hbuf); read_double(f, &(orec->subfiles[fi].long_inc));

		read_key(f, hbuf); read_int(f, &(orec->subfiles[fi].gs_count));
	}
}

void write_str(FILE *f, char* buf, size_t count)
{
	fwrite(buf, sizeof(char), count, f);
}

void write_int(FILE *f, int value)
{
	int tmp = value;
	fwrite(&tmp, sizeof(int), 1, f);
	tmp = 0;
	fwrite(&tmp, sizeof(int), 1, f);
}

void write_double(FILE *f, double value)
{
	double tmp = value;
	fwrite(&tmp, sizeof(double), 1, f);
}

void write_float(FILE *f, double value)
{
	float tmp = (float)value;
	fwrite(&tmp, sizeof(float), 1, f);
}

//! procedure to write NTv2 header to a file handle
/*!
    \param f a pointer to FILE structure.
		\param orec a pointer to t_header structure (NTv2 header).
*/
void write_header(FILE *f, t_header *orec)
{
	CPLStrlcpy(hbuf, "NUM_OREC", 8); write_str(f, hbuf, 8); write_int(f, orec->num_orec);
	CPLStrlcpy(hbuf, "NUM_SREC", 8); write_str(f, hbuf, 8); write_int(f, orec->num_srec);
	CPLStrlcpy(hbuf, "NUM_FILE", 8); write_str(f, hbuf, 8); write_int(f, orec->num_file);

	CPLStrlcpy(hbuf, "GS_TYPE ", 8); write_str(f, hbuf, 8); write_str(f, orec->gs_type, 8);
	CPLStrlcpy(hbuf, "VERSION ", 8); write_str(f, hbuf, 8); write_str(f, orec->version, 8);
	CPLStrlcpy(hbuf, "SYSTEM_F", 8); write_str(f, hbuf, 8); write_str(f, orec->system_f, 8);
	CPLStrlcpy(hbuf, "SYSTEM_T", 8); write_str(f, hbuf, 8); write_str(f, orec->system_t, 8);

	CPLStrlcpy(hbuf, "MAJOR_F ", 8); write_str(f, hbuf, 8); write_double(f, orec->major_f);
	CPLStrlcpy(hbuf, "MINOR_F ", 8); write_str(f, hbuf, 8); write_double(f, orec->minor_f);
	CPLStrlcpy(hbuf, "MAJOR_T ", 8); write_str(f, hbuf, 8); write_double(f, orec->major_t);
	CPLStrlcpy(hbuf, "MINOR_T ", 8); write_str(f, hbuf, 8); write_double(f, orec->minor_t);

	for(int si=0; si<orec->num_file; ++si)
	{
		CPLStrlcpy(hbuf, "SUB_NAME", 8); write_str(f, hbuf, 8); write_str(f, orec->subfiles[si].name, 8);
		CPLStrlcpy(hbuf, "PARENT  ", 8); write_str(f, hbuf, 8); write_str(f, orec->subfiles[si].parent, 8);
		CPLStrlcpy(hbuf, "CREATED ", 8); write_str(f, hbuf, 8); write_str(f, orec->subfiles[si].created, 8);
		CPLStrlcpy(hbuf, "UPDATED ", 8); write_str(f, hbuf, 8); write_str(f, orec->subfiles[si].updated, 8);

		CPLStrlcpy(hbuf, "S_LAT   ", 8); write_str(f, hbuf, 8); write_double(f, orec->subfiles[si].s_lat);
		CPLStrlcpy(hbuf, "N_LAT   ", 8); write_str(f, hbuf, 8); write_double(f, orec->subfiles[si].n_lat);
		CPLStrlcpy(hbuf, "E_LONG  ", 8); write_str(f, hbuf, 8); write_double(f, orec->subfiles[si].e_long);
		CPLStrlcpy(hbuf, "W_LONG  ", 8); write_str(f, hbuf, 8); write_double(f, orec->subfiles[si].w_long);

		CPLStrlcpy(hbuf, "LAT_INC ", 8); write_str(f, hbuf, 8); write_double(f, orec->subfiles[si].lat_inc);
		CPLStrlcpy(hbuf, "LONG_INC", 8); write_str(f, hbuf, 8); write_double(f, orec->subfiles[si].long_inc);

		CPLStrlcpy(hbuf, "GS_COUNT", 8); write_str(f, hbuf, 8); write_int(f, orec->subfiles[si].gs_count);
	}
}

};
using namespace std;

//! global variable to hold temporary string read from file
char buffer[1024];

//! function to load WKT file into a string
/*!
    \param sWktFilename a pointer to string containing WKT filename.
    \return pointer to char array if successful or halt the program
*/
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

int main(int argc, char*argv[])
{
	optparse::OptionParser parser = optparse::OptionParser().description("NTv2 GridShift Generator for intra-datum Grid-Shift (Inhom.-hom.)");

	parser.add_option("-s", "--source-coord").dest("src_coord").help("set source coordinate system WKT_FILE").metavar("WKT_FILE");
	parser.add_option("-r", "--ref-coord").dest("ref_coord").help("set destination coordinate system WKT_FILE").metavar("WKT_FILE");

	parser.add_option("-f", "--source-wkt").dest("src_wkt").help("set source coordinate system WKT").metavar("WKT");
	parser.add_option("-t", "--ref-wkt").dest("ref_wkt").help("set destination coordinate system WKT").metavar("WKT");

	parser.add_option("-i", "--input-grid").dest("input_file").help("set input grid shift filename GSB_FILE").metavar("GSB_FILE");
	parser.add_option("-o", "--input-grid").dest("output_file").help("set output grid shift filename GSB_FILE").metavar("GSB_FILE");

	parser.add_option("-v", "--verbose").dest("verbose").action("store_true").set_default("0").help("display debug information");

	optparse::Values options = parser.parse_args(argc, argv);
	vector<string> args = parser.args();

	if ((options["src_coord"].length() == 0 && options["src_wkt"].length() == 0) ||
		(options["ref_coord"].length() == 0 && options["ref_wkt"].length() == 0)){
			cerr << "Source and Reference Coordinate system (-s/-f and -r/-t) is not set." << endl;
			exit(1);
	}

	if ((options["input_file"].length() == 0 || options["output_file"].length() == 0)){
			cerr << "Input and Output gridshift file (-i and -o) is not set." << endl;
			exit(1);
	}

	bool is_verbose = options.get("verbose");

	//source CR in NTv2 grid
	OGRSpatialReference oSourceSRS, oTargetSRS;
	if(options["src_coord"].length() != 0){
		char *wkt1 = loadWktFile(options["src_coord"].c_str());
		oSourceSRS.importFromWkt(&(wkt1));
		oSourceSRS.exportToProj4(&(wkt1));
		if(is_verbose)
		cout << "SOURCE SRS: " << wkt1 << endl;
	}
	else{
		int slen = options["src_wkt"].length() + 1;
		char *cstr = new char[slen];
		CPLStrlcpy(cstr, options["src_wkt"].c_str(), slen);
		oSourceSRS.importFromWkt(&(cstr));
		delete [] cstr;
	}

	//reference CR (ETRS89)
	if(options["ref_coord"].length() != 0){
		char *wkt2 = loadWktFile(options["ref_coord"].c_str());
		oTargetSRS.importFromWkt(&(wkt2));
		oTargetSRS.exportToProj4(&(wkt2));
		if(is_verbose)
		cout << "REFERENCE SRS: "<< wkt2 << endl;
	}
	else{
		int slen = options["ref_wkt"].length() + 1;
		char *cstr = new char[slen];
		CPLStrlcpy(cstr, options["ref_wkt"].c_str(), slen);
		oTargetSRS.importFromWkt(&(cstr));
		delete [] cstr;
	}

	//forward
	OGRCoordinateTransformation *poCT = OGRCreateCoordinateTransformation( &oSourceSRS,
                                               &oTargetSRS );
	//reverse transformation
	OGRCoordinateTransformation *poCT_r = OGRCreateCoordinateTransformation( &oTargetSRS,
                                               &oSourceSRS );

	if(poCT == NULL){
		printf("ERR: cannot initialize coordinate transformation");
		return 1;
	}

	if(poCT_r == NULL){
		printf("ERR: cannot initialize rev. coordinate transformation");
		return 1;
	}

	char * input_file = (char*)&(*(options["input_file"].c_str()));
	if(!CPLCheckForFile(input_file, NULL)){
		printf("ERR: input file not exists");
		return 1;
	}
	
	//input grid shift coordinate
	FILE * fg;
	fopen_s(&fg, options["input_file"].c_str(), "rb");
	//fg = fopen(options["input_file"].c_str(), "rb");
	//output filename
	FILE * fo;
	fopen_s(&fo, options["output_file"].c_str(), "wb");
	//fo = fopen(options["output_file"].c_str(), "wb");

	if(fg != NULL && fo != NULL)
	{
		t_header orec;
		read_header(fg, &orec);
		
		if(is_verbose){
			printf("NUM_OREC : %d\n", orec.num_orec);
			printf("NUM_SREC : %d\n", orec.num_srec);
			printf("NUM_FILE : %d\n", orec.num_file);
			printf("GS_TYPE  : \"%s\"\n", orec.gs_type);
			printf("VERSION  : \"%s\"\n", orec.version);
			printf("SYSTEM_F : \"%s\"\n", orec.system_f);
			printf("SYSTEM_T : \"%s\"\n", orec.system_t);
			printf("MAJOR-MINOR_F : %f %f\n", orec.major_f, orec.minor_f);
			printf("MAJOR-MINOR_T : %f %f\n", orec.major_t, orec.minor_t);

			for(int si=0; si<orec.num_file; ++si)
			{
				printf("NAME (PARENT)  : \"%s\" (\"%s\")\n", orec.subfiles[si].name, orec.subfiles[si].parent);
				printf("CREATED (UPDATED)  : \"%s\" (\"%s\")\n", orec.subfiles[si].created, orec.subfiles[si].updated);
				printf("BBOX (S,W) (N,E)  : \n\t(%f, %f)\n\t(%f, %f)\n", orec.subfiles[si].s_lat, orec.subfiles[si].w_long, orec.subfiles[si].n_lat, orec.subfiles[si].e_long);
				printf("INCREMENT (lat, long)  : (%f, %f)\n", orec.subfiles[si].lat_inc, orec.subfiles[si].long_inc);
			}
		}
		
		CPLStrlcpy(orec.system_t, orec.system_f, 8);
		orec.major_t = orec.major_f;
		orec.minor_t = orec.minor_f;

		write_header(fo, &orec);

		double min_dx, min_dy, max_dx, max_dy;

		for(int si=0; si<orec.num_file; ++si)
		{
			int n_cols = (int)(1+(orec.subfiles[si].w_long-orec.subfiles[si].e_long)/orec.subfiles[si].long_inc);
			int n_rows = (int)(1+(orec.subfiles[si].n_lat-orec.subfiles[si].s_lat)/orec.subfiles[si].lat_inc);

			min_dx = 99999999.;
			min_dy = 99999999.;
			max_dx = -min_dx;
			max_dy = -min_dy;

			for(int ilat=0; ilat<n_rows; ++ilat)
			{
				for(int ilong=0; ilong<n_cols; ++ilong)
				{
					double lat = ilat * orec.subfiles[si].lat_inc + orec.subfiles[si].s_lat;
					double lon = ilong *orec.subfiles[si].long_inc + orec.subfiles[si].e_long;

					lat *= STORAD;
					lon *= STORAD;

					double x = lon;
					double y = lat;
					
					if(is_verbose)
					printf("lon: %f\tlat: %f ", lon, lat);
					

					if( poCT == NULL || !poCT->Transform( 1, &x, &y))
					{
						if(is_verbose)
							printf("Transformation Failed %f %f \n", x, y);
						
						continue;
					}
					if( poCT_r == NULL || !poCT_r->Transform(1, &x, &y))
					{
						if(is_verbose)
						printf("Reverse Transformation Failed %f %f\n", x, y);
						
						continue;
					}

					double dx = (x-lon)*RADTOS;
					double dy = (y-lat)*RADTOS;

					if(is_verbose)
					printf("\tdelta-lon: %f\tdelta-lat: %f\n", dx, dy);

					min_dx = MIN(min_dx, dx);
					min_dy = MIN(min_dy, dy);
					max_dx = MAX(max_dx, dx);
					max_dy = MAX(max_dy, dy);

					write_float(fo, dy);
					write_float(fo, dx);
					write_float(fo, -1.0);
					write_float(fo, -1.0);
				}
			}
		}

		//print stats
		if(is_verbose){
			printf("MIN %f %f\n", min_dx, min_dy);
			printf("MAX %f %f\n", max_dx, max_dy);
		}
		
		//write closing header
		CPLStrlcpy(hbuf, "END     ", 8); write_str(fo, hbuf, 8); write_int(fo, 0);

		CPLFree(orec.subfiles);
		fclose(fo);
		fclose(fg);
	}
	else
	{
		fprintf(stderr, "ERROR opening input/output file(s)\n");
		return 1;
	}
	
	return 0;
}
