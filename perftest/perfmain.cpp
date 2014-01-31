/******************************************************************************
 *
 * Project:  SpatialRef3D performance test
 * Purpose:  program to test the performance transformation speed in 
 *           different size of data chunk.
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
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdio>
#include <sstream>
#include <ctime>
#include <random>

#include "cpl_conv.h"
#include "ogr_spatialref3D.h"
#include "OptionParser.h"


/************************************************************************/
/*                         OGRSpatialReference3D                        */
/************************************************************************/

/**\file perfmain.cpp
 * This file is an program to test the speed of implemented code 
 * using reference data supplied in CSV. 
 * 
 * The command-line options for running this program are:
 *	
 *	
 *		-c | --chunk-size=CHUNK			: number of point taken for each measurement
 *										DEFAULT = 10
 *	
 *		-d | --dest-coord=WKT_FILE		: set WKT_FILE as target coordinate system
 *										  description
 *	
 *		-i | --input-file=FILE			: set FILE as input reference coordinate data
 *	
 *		-m | --max-input=N				: limit N maximum number of input data taken from sample file
 *										DEFAULT = MAX_DATA (16M points)
 *	
 *		-n | --num-input=N				: number of input data N taken from sample file
 *										  (value of -1 means all data in file will be used)
 *										DEFAULT = -1
 *	
 *		-r | --repeat=N					: number of repetition for each transformation to be done
 *										DEFAULT = 3 (MINIMUM = 2)
 *	
 *		-s | --source-coord=FILE		: set FILE as source coordinate system
 *										  description
 *	
 *  
 *
 *
 ************************************************************************/

//! temporary variable for storing coordinate data read from file
double *x_in;
//! temporary variable for storing coordinate data read from file
double *y_in;
//! temporary variable for storing coordinate data read from file
double *z_in;

//! temporary variable for storing coordinate data used in transformation
double *x_out;
//! temporary variable for storing coordinate data used in transformation
double *y_out;
//! temporary variable for storing coordinate data used in transformation
double *z_out;

//#define TEST_FILE "Line13.xyz"
//! maximum number of rows from data file
#define MAX_DATA 16000000	//12877662

//! macro to retrieve timer value in seconds (double)
#define GET_TIMER(x) x = (double)(clock())/CLOCKS_PER_SEC; //in [s]

//! macro to get time difference
#define DIFF_TIME(a,b) (a-b)

//! temporary storage for data read from file
char buffer[1024];

/*
#define SOURCE_SRS "utm33-etrs89.prj"
#define TARGET_SRS "utm33-etrs89-orthoH.prj"		//geoid
#define TARGET_SRS_HEIGHT "gkm34-mgi-gkm34.prj"		//datum change + geoid
#define TARGET_SRS_ALLH "gkm34-mgi-gkm34_h.prj"		//datum change + geoid + height correction
#define TARGET_SRS_ELLH "gkm34-mgi-gkm34_ell.prj"	//datum change
*/

using namespace std;

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

//! program's entry point
int main(int argc, char *argv[])
{
	optparse::OptionParser parser = optparse::OptionParser().description("OGRSpatialRef3D performance testing program");

	parser.add_option("-i", "--input-file").dest("input_file").help("set input data FILE").metavar("FILE");

	parser.add_option("-s", "--source-coord").dest("src_coord").help("set source coordinate system WKT_FILE").metavar("WKT_FILE");
	parser.add_option("-d", "--dest-coord").dest("dst_coord").help("set destination coordinate system WKT_FILE").metavar("WKT_FILE");

	parser.add_option("-m", "--max-input").dest("max_input").help("maximum number of points read from input file (default 16M)").set_default(MAX_DATA);
	parser.add_option("-c", "--chunk-size").dest("chunk_size").help("number of points per chunk in one transformation call (default 10 pts per chunk)").metavar("CHUNK").set_default(10);
	parser.add_option("-r", "--repeat").dest("num_repeat").help("number of repetition for each transformation should be done (minimum 2, default 3)").metavar("N").set_default(3);

	optparse::Values options = parser.parse_args(argc, argv);
	vector<string> args = parser.args();

	if ((options["input_file"].length() == 0)){
			cerr << "Input Reference Coordinate is not set." << endl;
			exit(1);
	}

	if ((options["src_coord"].length() == 0)){
			cerr << "Source Spatial Reference is not set." << endl;
			exit(1);
	}

	if ((options["dst_coord"].length() == 0)){
			cerr << "Target Spatial Reference is not set." << endl;
			exit(1);
	}

	int max_input = atoi(options["max_input"].c_str());
	const char *TEST_FILE = options["dst_coord"].c_str();

	double start_time, end_time;
	GET_TIMER(start_time);

	FILE *fi = fopen(options["input_file"].c_str(), "r");
	if (!fi) {
		cerr << "Can't open input file " << options["input_file"] << endl;
		exit(1);
	}
	else{
		cout << "reading file " << options["input_file"] << endl;
	}

	x_in = (double*)CPLMalloc(sizeof(double)*max_input);
	y_in = (double*)CPLMalloc(sizeof(double)*max_input);
	z_in = (double*)CPLMalloc(sizeof(double)*max_input);
	
	int last_num_data = 1;
	int num_data = 0;

	while(!feof(fi)){
		fscanf(fi, "%lf %lf %lf", &(x_in[num_data]), &(y_in[num_data]), &(z_in[num_data]));

		num_data += 1;
		if ((num_data/last_num_data)==10){
			cout << num_data << "...";
			last_num_data = num_data;
		}

		if(max_input != -1 && num_data == max_input)
			break;
	}
	fclose(fi);
	
	GET_TIMER(end_time);
	cout << fixed;
	cout << num_data << endl;
	cout << DIFF_TIME(end_time, start_time)<< " s" << endl;

	//--

	//int num_samples[] = {1,10,100,1000,10000,100000,1000000,10000000};
	int num_samples = atoi(options["chunk_size"].c_str());
	int sample_count = 0;
	double delta_time = 0.0;
	int data_offset = 0;
	int num_run = atoi(options["num_repeat"].c_str());

	double sum = 0.0;
	double sumsq = 0.0;

	x_out = (double*)CPLMalloc(sizeof(double)*num_samples);
	y_out = (double*)CPLMalloc(sizeof(double)*num_samples);
	z_out = (double*)CPLMalloc(sizeof(double)*num_samples);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;
	OGRCoordinateTransformation3D *poCT;

	char *wkt;

	//--

	wkt = loadWktFile(options["src_coord"].c_str());
	oSourceSRS.importFromWkt3D(&(wkt));

	wkt = loadWktFile(options["dst_coord"].c_str());
	oTargetSRS.importFromWkt3D(&(wkt));

	poCT = OGRCreateCoordinateTransformation3D(&oSourceSRS, &oTargetSRS );
	
	sum = 0.0;
	sumsq = 0.0;
	for(int run=0; run<num_run; ++run){
		GET_TIMER(start_time);

		data_offset = 0;
		while(data_offset<num_data){
			
			sample_count = MIN(num_samples, num_data-data_offset+1);
			for (int sample=0; sample<sample_count; ++sample)
			{
				x_out[sample] = x_in[data_offset+sample];
				y_out[sample] = y_in[data_offset+sample];
				z_out[sample] = z_in[data_offset+sample];
			}
			data_offset += sample_count; 

			if( poCT == NULL || !poCT->Transform( sample_count, x_out, y_out ,z_out) )
			{
				cout << "Transformation failed.\n";
			}
		}//process next chunk until all data used

		GET_TIMER(end_time);
		delta_time = DIFF_TIME(end_time, start_time);

		sum += delta_time;
		sumsq += delta_time*delta_time;

		//cout << SOURCE_SRS << " to " << TARGET_SRS << " : ";
		//cout << delta_time << endl;
	}
	cout << fixed;
	cout << setprecision(3);
	cout << " num samples : " << num_samples << endl;
	cout << "running " << num_run << "times. AVG: " << sum/num_run << " STDDEV: " << ((double)num_run*sumsq-sum*sum)/(double)(num_run*(num_run-1)) << endl;
	
	CPLFree(x_in);
	CPLFree(y_in);
	CPLFree(z_in);
	CPLFree(x_out);
	CPLFree(y_out);
	CPLFree(z_out);
	return 0;
}
