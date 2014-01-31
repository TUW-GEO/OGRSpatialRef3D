/******************************************************************************
 *
 * Project:  OpenGIS Simple Features Reference Implementation
 * Purpose:  Classes for manipulating spatial reference systems with
 *           vertical datum suppurt in a platform non-specific manner.
 * Authors:  Gottfried Mandlburger, Johannes Otepka, Bhargav patel
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
