#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "cpl_conv.h"
#include "ogr_spatialref3D.h"
#include "OptionParser.h"

using namespace std;

char buffer[1024];

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

int main(int argc, char* argv[])
{
  //init gdal/proj.4 data directory path 
	optparse::OptionParser parser = optparse::OptionParser().description("OGRSpatialRef3D test program");

	parser.add_option("-g", "--gdal-data").dest("gdal_data").help("set path to gdal data").set_default("..\\gdal-1.10.0\\data");
	parser.add_option("-s", "--source-coord").dest("src_coord").help("set source coordinate system WKT_FILE").metavar("WKT_FILE");
	parser.add_option("-d", "--dest-coord").dest("dst_coord").help("set destination coordinate system WKT_FILE").metavar("WKT_FILE");

	parser.add_option("-v", "--source-wkt").dest("src_wkt").help("set source coordinate system WKT").metavar("WKT");
	parser.add_option("-w", "--destination-wkt").dest("dst_wkt").help("set destination coordinate system WKT").metavar("WKT");

	parser.add_option("-i", "--input-coord").dest("input_file").help("set input coordinate data FILE").metavar("FILE");

	optparse::Values options = parser.parse_args(argc, argv);
	vector<string> args = parser.args();

	if ((options["src_coord"].length() == 0 && options["src_wkt"].length() == 0) ||
		(options["dst_coord"].length() == 0 && options["dst_wkt"].length() == 0)){
			cerr << "Source or Destination Coordinate is not set." << endl;
			exit(1);
	}

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	if(options["src_coord"].length() != 0){
		char *wkt1 = loadWktFile(options["src_coord"].c_str());
		oSourceSRS.importFromWkt(&(wkt1));
		//oSourceSRS.exportToProj4(&(wkt1));
		//cout << wkt1 << endl;
	}
	else{
		int slen = options["src_wkt"].length() + 1;
		char *cstr = new char[slen];
		CPLStrlcpy(cstr, options["src_wkt"].c_str(), slen);
		oSourceSRS.importFromWkt(&(cstr));
		delete [] cstr;
	}

	if(options["dst_coord"].length() != 0){
		char *wkt2 = loadWktFile(options["dst_coord"].c_str());
		oTargetSRS.importFromWkt(&(wkt2));
		//oTargetSRS.exportToProj4(&(wkt2));
		//cout << wkt2 << endl;
	}
	else{
		int slen = options["dst_wkt"].length() + 1;
		char *cstr = new char[slen];
		CPLStrlcpy(cstr, options["dst_wkt"].c_str(), slen);
		oTargetSRS.importFromWkt(&(cstr));
		delete [] cstr;
	}

	const OGR_SRSNode *poGEOID = oSourceSRS.GetAttrNode( "GEOID" );
	if (poGEOID != NULL){
		const char *pszGEOID = poGEOID->GetChild(1)->GetChild(0)->GetValue();
		std::cout << "Loading source GEOID : " << pszGEOID << std::endl;
		oSourceSRS.SetGeoidModel(pszGEOID);
	}
	else{
		std::cout << "NO GEOID in SOURCE" << std::endl;
	}

	poGEOID = oTargetSRS.GetAttrNode( "GEOID" );
	if (poGEOID != NULL){
		const char *pszGEOID = poGEOID->GetChild(1)->GetChild(0)->GetValue();
		std::cout << "Loading target GEOID : " << pszGEOID << std::endl;
		oTargetSRS.SetGeoidModel(pszGEOID);
	}
	else{
		std::cout << "NO GEOID in TARGET" << std::endl;
	}

	CPLSetConfigOption("GDAL_DATA", options["gdal_data"].c_str());
	if(options["input_file"].length() == 0){
		cerr << "no data FILE given" << endl;
		exit(1);
	}

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( &oSourceSRS,
                                               &oTargetSRS );
	
	cout << "coordinate transform created" << endl;

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
		stringstream ss(line);

		double sourcex = 0.0;  
		double sourcey = 0.0;
		double sourcez = 0.0;

		ss >> sourcex >> sourcey >> sourcez;

		double targetx = sourcex;
		double targety = sourcey;
		double targetz = sourcez;

		cout << sourcex << sourcey << sourcez << endl;
		//do actual transformation
		if( poCT == NULL || !poCT->Transform( 1, &targetx, &targety ,&targetz) )
			printf( "Transformation failed.\n" );
		else
			printf( "(%f,%f,%f) -> (%f,%f,%f)\n", sourcex, sourcey, sourcez, targetx, targety, targetz );
	}

	cout << "Press ENTER to exit.";
	cin.get();

	return 0;
}