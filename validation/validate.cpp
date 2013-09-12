#include <iostream>
#include <iomanip>
#include <string>

#include "cpl_conv.h"
#include "ogr_spatialref3D.h"
#include "OptionParser.h"
#include "validate.h"

using namespace std;

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
	
	val_init();

	loadRefFile(options["input_file"], atoi(options["num_input"].c_str()));

	cout << fixed; cout << "# data point(s) : " << num_data << endl;

	val_geoc_etrs();
	//val_geog_etrs();
	//val_geog_etrs_ortho();

	//val_geog_mgi(); //hell_mgi not available
	//val_geog_mgi_ortho();
	//val_proj_mgi(); //error rasterio out of bound

	cout << "cleaning up.." << endl;
	val_cleanup();
	
	cout << "Press ENTER to exit.";
	cin.get();
	
	return 0;
}
