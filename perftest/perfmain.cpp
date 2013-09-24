/******************************************************************************
 *
 * Project:  SpatialRef3D performance test
 * Purpose:  program to test the performance transformation speed in 
 *           different size of data chunk.
 * Authors:  Peb Ruswono Aryan, Gottfried Mandlburger, Johannes Otepka
 *
 ******************************************************************************
 * Copyright (c) 2012,  I.P.F., TU Vienna.
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

double *x_in, *y_in, *z_in;
double *x_out, *y_out, *z_out;

#define TEST_FILE "Line13.xyz"
#define MAX_DATA 13000000	//12877662

#define GET_TIMER(x) x = (double)(clock())/CLOCKS_PER_SEC; //in [s]
#define DIFF_TIME(a,b) (a-b)

char buffer[1024];

#define SOURCE_SRS "utm33-etrs89.prj"
#define TARGET_SRS "utm33-etrs89-orthoH.prj"		//geoid
#define TARGET_SRS_HEIGHT "gkm34-mgi-gkm34.prj"		//datum change + geoid
#define TARGET_SRS_ALLH "gkm34-mgi-gkm34_h.prj"		//datum change + geoid + height correction
#define TARGET_SRS_ELLH "gkm34-mgi-gkm34_ell.prj"	//datum change

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

int main(int argc, char *argv[])
{
	double start_time, end_time;
	GET_TIMER(start_time);

	FILE *fi = fopen(TEST_FILE, "r");
	if (!fi) {
		cerr << "Can't open input file " << TEST_FILE << endl;
		exit(1);
	}
	else{
		cout << "reading file " << TEST_FILE << endl;
	}

	x_in = (double*)CPLMalloc(sizeof(double)*MAX_DATA);
	y_in = (double*)CPLMalloc(sizeof(double)*MAX_DATA);
	z_in = (double*)CPLMalloc(sizeof(double)*MAX_DATA);
	x_out = (double*)CPLMalloc(sizeof(double)*MAX_DATA);
	y_out = (double*)CPLMalloc(sizeof(double)*MAX_DATA);
	z_out = (double*)CPLMalloc(sizeof(double)*MAX_DATA);

	cout << fixed;
	int last_num_data = 1;
	int num_data = 0;
	int max_input = -1;

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

	int num_samples[] = {1,10,100,1000,10000,100000,1000000,10000000};
	int sample_count = 0;
	double delta_time = 0.0;
	int data_offset = 0;
	int num_run = 10;

	double sum = 0.0;
	double sumsq = 0.0;

	OGRSpatialReference3D oSourceSRS, oTargetSRS, oTargetSRS_h, oTargetSRS_hh, oTargetSRS_eh;
	OGRCoordinateTransformation3D *poCT[4];

	char *wkt;

	//--

	wkt = loadWktFile(SOURCE_SRS);
	oSourceSRS.importFromWkt3D(&(wkt));

	wkt = loadWktFile(TARGET_SRS);
	oTargetSRS.importFromWkt3D(&(wkt));

	wkt = loadWktFile(TARGET_SRS_HEIGHT);
	oTargetSRS_h.importFromWkt3D(&(wkt));

	wkt = loadWktFile(TARGET_SRS_ALLH);
	oTargetSRS_hh.importFromWkt3D(&(wkt));

	wkt = loadWktFile(TARGET_SRS_ELLH);
	oTargetSRS_eh.importFromWkt3D(&(wkt));

	poCT[0] = OGRCreateCoordinateTransformation3D(&oSourceSRS, &oTargetSRS );
	poCT[1] = OGRCreateCoordinateTransformation3D(&oSourceSRS, &oTargetSRS_h );
	poCT[2] = OGRCreateCoordinateTransformation3D(&oSourceSRS, &oTargetSRS_hh );
	poCT[3] = OGRCreateCoordinateTransformation3D(&oSourceSRS, &oTargetSRS_eh );
	
	std::cout << std::setprecision(3);

	for(int trafo=0; trafo<4; ++trafo)
	{
		for(int cur_step=0; cur_step<8; ++cur_step)
		{
			if (num_data<num_samples[cur_step]) break;
			cout << " num samples : " << num_samples[cur_step] << endl;
		
			sum = 0.0;
			sumsq = 0.0;
			for(int run=0; run<num_run; ++run){
				GET_TIMER(start_time);

				//just geoid
				data_offset = 0;
				while(data_offset<num_data){
			
					sample_count = MIN(num_samples[cur_step], num_data-data_offset+1);
					for (int sample=0; sample<sample_count; ++sample)
					{
						x_out[sample] = x_in[data_offset+sample];
						y_out[sample] = y_in[data_offset+sample];
						z_out[sample] = z_in[data_offset+sample];
					}
					data_offset += sample_count; 

					//insert transform here
					if( poCT == NULL || !poCT[trafo]->Transform( sample_count, x_out, y_out ,z_out) )
					{
						cout << "Transformation failed.\n";
					}
				}//retry

				GET_TIMER(end_time);
				delta_time = DIFF_TIME(end_time, start_time);

				sum += delta_time;
				sumsq += delta_time*delta_time;

				switch(trafo){
					case 0:cout << SOURCE_SRS << " to " << TARGET_SRS << " : "; break;
					case 1:cout << SOURCE_SRS << " to " << TARGET_SRS_HEIGHT << " : "; break;
					case 2:cout << SOURCE_SRS << " to " << TARGET_SRS_ALLH << " : "; break;
					case 3:cout << SOURCE_SRS << " to " << TARGET_SRS_ELLH << " : "; break;
				}
				cout << delta_time << endl;
			}
			cout << "running " << num_run << "times. AVG: " << sum/num_run << " STDDEV: " << ((double)num_run*sumsq-sum*sum)/(double)(num_run*(num_run-1)) << endl;
		}//step

	}

	CPLFree(x_in);
	CPLFree(y_in);
	CPLFree(z_in);
	CPLFree(x_out);
	CPLFree(y_out);
	CPLFree(z_out);
	return 0;
}
