/******************************************************************************
 *
 * Project:  Spatialref3D correctness validation test
 * Purpose:  program to test the correctness of transformation output
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
#include <sstream>
#include <iterator>

#include "cpl_conv.h"
#include "ogr_spatialref3D.h"
#include "OptionParser.h"
#include "validate.h"

/************************************************************************/
/*                         OGRSpatialReference3D                        */
/************************************************************************/

/**\file validate.cpp
 * This file is an program to test the correctness of implemented code 
 * using reference data supplied in CSV. 
 * 
 * The command-line options for running this program are:
 *	
 *		-i | --input-file=FILE	: set FILE as input reference coordinate data
 *	
 *		-n | --num-input=N		: number of input data N taken from sample file
 *								  (value of -1 means all data in file will be used)
 *								DEFAULT = -1
 *	
 *  
 *
 *
 ************************************************************************/

using namespace std;

//! program's entry point
int main(int argc, char *argv[])
{
	cout << "validation" << endl;
	
	optparse::OptionParser parser = optparse::OptionParser().description("OGRSpatialRef3D validation program");
	
	parser.add_option("-i", "--input-file").dest("input_file").help("set input reference coordinate data FILE").metavar("FILE");
	parser.add_option("-n", "--num-input").dest("num_input").help("number of input data N taken from sample file (-1 means all data in file)").metavar("N").set_default(-1);
	parser.add_option("-s", "--source-col").dest("source_col").help("set column names COLS which contains source coordinate data (comma separated names without spaces)").metavar("COLS");
	parser.add_option("-t", "--target-col").dest("target_col").help("set column names COLS which contains target coordinate data (comma separated names without spaces)").metavar("COLS");
	
	optparse::Values options = parser.parse_args(argc, argv);
	vector<string> args = parser.args();
	
	if ((options["input_file"].length() == 0)){
			cerr << "Input Reference Coordinate is not set." << endl;
			exit(1);
	}
	else if ((options["source_col"].length() == 0)){
		cerr << "source column names is not set" << endl;
		exit(1);
	}
	else if ((options["target_col"].length() == 0)){
		cerr << "target column names is not set" << endl;
		exit(1);
	}

	string delimiter = ","; 
	split(options["source_col"], delimiter, src_cols);
	split(options["target_col"], delimiter, tgt_cols);

	for(vector<string>::iterator it = src_cols.begin(); it != src_cols.end(); ++it){
		cout << *it << endl;
	}
	for(vector<string>::iterator it = tgt_cols.begin(); it != tgt_cols.end(); ++it){
		cout << *it << endl;
	}
	
	
	val_init();

	loadRefFile(options["input_file"], atoi(options["num_input"].c_str()));

	cout << fixed; cout << "# data point(s) : " << num_data << endl;
	return 0;

	val_geoc_etrs();
	val_geog_etrs();
	val_geog_etrs_ortho();

	val_geog_mgi();
	val_geog_mgi_ortho();
	val_proj_mgi(); 

	cout << "cleaning up.." << endl;
	val_cleanup();
	
	cout << "Press ENTER to exit.";
	cin.get();
	
	return 0;
}
