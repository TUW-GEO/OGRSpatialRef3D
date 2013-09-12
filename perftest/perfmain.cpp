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

char buffer[1024];

#define SOURCE_SRS "utm33-etrs89.prj"
#define TARGET_SRS "utm33-etrs89-orthoH.prj"	//geoid
#define TARGET_SRS_HEIGHT "gkm34-mgi-gkm34.prj" //datum change + geoid
#define TARGET_SRS_ALLH "gkm34-mgi-gkm34_h.prj"	//datum change + geoid + height correction

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
	time_t start_time, end_time;
	time(&start_time);

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
	
	time(&end_time);
	cout << fixed;
	cout << num_data << endl;
	cout << (int)difftime(end_time, start_time)<< " s" << endl;

	//return 0;

	int num_samples[] = {1,10,100,1000,10000,100000,1000000,10000000};
	double sample_sum[8];
	double sample_sumsq[8];
	double sample_sum_h[8];
	double sample_sumsq_h[8];
	double sample_sum_hh[8];
	double sample_sumsq_hh[8];
	int num_sampling = 10;
	double delta_time = 0.0;

	OGRSpatialReference3D oSourceSRS, oTargetSRS, oTargetSRS_h, oTargetSRS_hh;
	char *wkt;
	wkt = loadWktFile(SOURCE_SRS);
	oSourceSRS.importFromWkt3D(&(wkt));

	wkt = loadWktFile(TARGET_SRS);
	oTargetSRS.importFromWkt3D(&(wkt));

	wkt = loadWktFile(TARGET_SRS_HEIGHT);
	oTargetSRS_h.importFromWkt3D(&(wkt));

	wkt = loadWktFile(TARGET_SRS_ALLH);
	oTargetSRS_hh.importFromWkt3D(&(wkt));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	OGRCoordinateTransformation3D *poCT_h = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS_h );

	OGRCoordinateTransformation3D *poCT_hh = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS_hh );
	default_random_engine generator;
	int cur_step;

	for(cur_step=0; cur_step<=7; ++cur_step)
	{
		if (num_data<num_samples[cur_step]) break;
		cout << cur_step << " num samples : " << num_samples[cur_step] << endl;
		uniform_int_distribution<int> distribution(0,num_data-num_samples[cur_step]);
		int data_offset = 0;

		sample_sum[cur_step] = 0.0;
		sample_sumsq[cur_step] = 0.0;
		sample_sum_h[cur_step] = 0.0;
		sample_sumsq_h[cur_step] = 0.0;
		sample_sum_hh[cur_step] = 0.0;
		sample_sumsq_hh[cur_step] = 0.0;
		
		for(int retry=0; retry<num_sampling; ++retry){
			//generate experiment data
			data_offset = distribution(generator);

			//just geoid

			for (int sample=0; sample<num_samples[cur_step]; ++sample)
			{
				x_out[sample] = x_in[data_offset+sample];
				y_out[sample] = y_in[data_offset+sample];
				z_out[sample] = z_in[data_offset+sample];
			}

			time(&start_time);
			//insert transform here
			if( poCT == NULL || !poCT->Transform( num_samples[cur_step], x_out, y_out ,z_out) )
			{
				cout << "Transformation failed.\n";
			}
			time(&end_time);
			delta_time = difftime(end_time, start_time)*1000;

			sample_sum[cur_step] += delta_time;
			sample_sumsq[cur_step] += (delta_time*delta_time);


			//datum shift + geoid

			for (int sample=0; sample<num_samples[cur_step]; ++sample)
			{
				x_out[sample] = x_in[data_offset+sample];
				y_out[sample] = y_in[data_offset+sample];
				z_out[sample] = z_in[data_offset+sample];
			}

			time(&start_time);
			//insert transform here
			if( poCT == NULL || !poCT_h->Transform( num_samples[cur_step], x_out, y_out ,z_out) )
			{
				cout << "Transformation failed.\n";
			}
			time(&end_time);
			delta_time = difftime(end_time, start_time)*1000;

			sample_sum_h[cur_step] += delta_time;
			sample_sumsq_h[cur_step] += (delta_time*delta_time);

			//datum shift + geoid + height correction

			for (int sample=0; sample<num_samples[cur_step]; ++sample)
			{
				x_out[sample] = x_in[data_offset+sample];
				y_out[sample] = y_in[data_offset+sample];
				z_out[sample] = z_in[data_offset+sample];
			}

			time(&start_time);
			//insert transform here
			if( poCT == NULL || !poCT_hh->Transform( num_samples[cur_step], x_out, y_out ,z_out) )
			{
				cout << "Transformation failed.\n";
			}
			time(&end_time);
			delta_time = difftime(end_time, start_time)*1000;

			sample_sum_hh[cur_step] += delta_time;
			sample_sumsq_hh[cur_step] += (delta_time*delta_time);
		}//retry

		std::cout << std::setprecision(3);
		double n = (double)num_sampling;
		cout << "num data : " << n << " retries : " << num_sampling << endl;
		cout << SOURCE_SRS << " to " << TARGET_SRS << endl;
		cout << "avg : " << sample_sum[cur_step] / n << endl;
		cout << "var : " << ((double)n*sample_sumsq[cur_step] - sample_sum[cur_step]*sample_sum[cur_step]) / (n*(n-1)) << std::endl; 
		cout << SOURCE_SRS << " to " << TARGET_SRS_HEIGHT << endl;
		cout << "avg : " << sample_sum_h[cur_step] / n << endl;
		cout << "var : " << ((double)n*sample_sumsq_h[cur_step] - sample_sum_h[cur_step]*sample_sum_h[cur_step]) / (n*(n-1)) << std::endl; 
		cout << SOURCE_SRS << " to " << TARGET_SRS_ALLH << endl;
		cout << "avg : " << sample_sum_hh[cur_step] / n << endl;
		cout << "var : " << ((double)n*sample_sumsq_hh[cur_step] - sample_sum_hh[cur_step]*sample_sum_hh[cur_step]) / (n*(n-1)) << std::endl; 
	}//step

	CPLFree(x_in);
	CPLFree(y_in);
	CPLFree(z_in);
	CPLFree(x_out);
	CPLFree(y_out);
	CPLFree(z_out);
	return 0;
}
