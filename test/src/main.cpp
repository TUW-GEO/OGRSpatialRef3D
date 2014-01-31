/******************************************************************************
 *
 * Project:  Spatialref3D coordinate transformation test
 * Purpose:  program for testing coordinate transformation function
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

#include "cpl_conv.h"
#include "ogr_spatialref3D.h"
#include "OptionParser.h"

/************************************************************************/
/*                         OGRSpatialReference3D                        */
/************************************************************************/

/**\file main.cpp
 * This file is an example commandline program that utilize 
 * OGRSpatialReference3D class to transform point(s) from one Spatial 
 * Reference system to another.
 * 
 * The command-line options for running this program are:
 *  
 *	
 *		-d | --dest-coord=WKT_FILE		: set WKT_FILE as target coordinate system
 *										  description
 *	
 *		-g | --gdal-data=TEXT			: set path to global data used by GDAL
 *										DEFAULT = ..\\gdal-1.10.0\\data
 *	
 *		-i | --input-coord=FILE			: set FILE as input coordinate data
 *	
 *		-s | --source-coord=FILE		: set FILE as source coordinate system
 *										  description
 *	
 *  
 *
 *
 ************************************************************************/

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

//! program's entry point
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
		oSourceSRS.importFromWkt3D(&(wkt1));
		oSourceSRS.exportToProj4(&(wkt1));
		cout << "SOURCE SRS: " << wkt1 << endl;
	}
	else{
		int slen = options["src_wkt"].length() + 1;
		char *cstr = new char[slen];
		CPLStrlcpy(cstr, options["src_wkt"].c_str(), slen);
		oSourceSRS.importFromWkt3D(&(cstr));
		delete [] cstr;
	}

	if(options["dst_coord"].length() != 0){
		char *wkt2 = loadWktFile(options["dst_coord"].c_str());
		oTargetSRS.importFromWkt3D(&(wkt2));
		oTargetSRS.exportToProj4(&(wkt2));
		cout << "TARGET SRS: "<< wkt2 << endl;
	}
	else{
		int slen = options["dst_wkt"].length() + 1;
		char *cstr = new char[slen];
		CPLStrlcpy(cstr, options["dst_wkt"].c_str(), slen);
		oTargetSRS.importFromWkt3D(&(cstr));
		delete [] cstr;
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

		cout << sourcex << " " << sourcey << " " << sourcez << endl;
		//do actual transformation
		if( poCT == NULL || !poCT->Transform( 1, &targetx, &targety ,&targetz) )
			printf( "Transformation failed.\n" );
		else
			printf( "(%f, %f, %f) -> (%f, %f, %f)\n", sourcex, sourcey, sourcez, targetx, targety, targetz );
	}

	cout << "Press ENTER to exit.";
	cin.get();

	return 0;
}