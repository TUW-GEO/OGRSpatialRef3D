/******************************************************************************
 *
 * Project:  OpenGIS Simple Features Reference Implementation
 * Purpose:  library for doing interpolation
 * Authors:  Gottfried Mandlburger, Johannes Otepka, Peb Ruswono Aryan
 *
 ******************************************************************************
 * Copyright (c) 2012-2013,  I.P.F., TU Vienna.
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
#ifndef __VALIDATE_H__
#define __VALIDATE_H__

#include <string>
#include <vector>
#include <map>
#include "cpl_conv.h"

#define MAX_DATA 2500

#define GEOG_ETRS "etrs89.prj" 
#define GEOC_ETRS "etrs89_geoc.prj" 
#define GEOG_ETRS_ORTH "etrs89_ortho.prj" 

#define GEOG_MGI "mgi.prj" 
#define GEOC_MGI "mgi_geocen.prj" 
#define GEOG_MGI_ORTH "mgi_ortho.prj" 

#define PROJ_MGI "mgi_proj.prj" //orthometric height with offset
#define PROJ_MGI_28 "mgi_proj_28.prj" //orthometric height with offset
#define PROJ_MGI_31 "mgi_proj_31.prj" //orthometric height with offset
#define PROJ_MGI_34 "mgi_proj_34.prj" //orthometric height with offset

using namespace std;
/*
 *
 */
// - given values
extern double *x_etrs, *y_etrs, *z_etrs;
extern double *lon_grs, *lat_grs;
extern double *hell_grs, *h_orth;
extern double *x_gebr, *y_gebr, *h_gebr;
extern double *und_bess, *und_grs;
extern double *lat_mgi, *lon_mgi;
extern double *ras_val, *h_grid;
extern int *ms; // meridian strip
extern int num_data;

extern double *x_src, *y_src, *z_src;
extern double *x_tgt, *y_tgt, *z_tgt;
extern double *und_src, *vcorr_src;
extern double *und_tgt, *vcorr_tgt;

//computed values
extern double *hell_mgi, *x_mgi, *y_mgi, *z_mgi;
extern vector<string> src_cols, tgt_cols;
extern map<string, int> columns;

class SummStat
{
public:
	double min;
	double max;
	double sum;
	double ssq;
	int n;
	SummStat() : min(HUGE_VAL), max(-HUGE_VAL), sum(0.0), ssq(0.0), n(0) {}
	void add(double value);
	void printout(int shortprec=8, int longprec=18);
};

void val_init();
void val_cleanup();

char *loadWktFile(const char* sWktFilename);
void loadRefFile(std::string, int);

void split(string source, string delimiter, vector<string> &tokens);

/*
 * combinations : 
 *		- geog_etrs_ortho
 *		- geog_etrs	(ellipsoidal)
 *		- geoc_etrs
 *
 *		- geoc_mgi
 *		- geog_mgi (ellipsoidal)
 *		- geog_mgi_ortho
 *		- proj_mgi (vertically corrected + shifted)
 */

//-
// 1-way transformation
//-

// geoc_etrs

void geoc_etrs_to_geog_etrs();
void geoc_etrs_to_geog_etrs_ortho();

void geoc_etrs_to_geoc_mgi();
void geoc_etrs_to_geog_mgi();
void geoc_etrs_to_geog_mgi_ortho();
void geoc_etrs_to_proj_mgi();

void val_geoc_etrs();
// geog_etrs

void geog_etrs_to_geoc_etrs();
void geog_etrs_to_geog_etrs_ortho();

void geog_etrs_to_geoc_mgi();
void geog_etrs_to_geog_mgi();
void geog_etrs_to_geog_mgi_ortho();
void geog_etrs_to_proj_mgi();

void val_geog_etrs();
// geog_etrs_ortho

void geog_etrs_ortho_to_geoc_etrs();
void geog_etrs_ortho_to_geog_etrs();

void geog_etrs_ortho_to_geoc_mgi();
void geog_etrs_ortho_to_geog_mgi();
void geog_etrs_ortho_to_geog_mgi_ortho();
void geog_etrs_ortho_to_proj_mgi();

void val_geog_etrs_ortho();
// geoc_mgi : no data available

//void geoc_mgi_to_geoc_etrs();
//void geoc_mgi_to_geog_etrs();
//void geoc_mgi_to_geog_etrs_ortho();

//void geoc_mgi_to_geog_mgi();
//void geoc_mgi_to_geoc_mgi_ortho();
//void geoc_mgi_to_proj_mgi();

//void val_geoc_mgi();
// geog_mgi

void geog_mgi_to_geoc_etrs();
void geog_mgi_to_geog_etrs();
void geog_mgi_to_geog_etrs_ortho();

void geog_mgi_to_geoc_mgi();
void geog_mgi_to_geog_mgi_ortho();
void geog_mgi_to_proj_mgi_gebr();

void val_geog_mgi();
// geog_mgi_ortho

void geog_mgi_ortho_to_geoc_etrs();
void geog_mgi_ortho_to_geog_etrs_ortho();
void geog_mgi_ortho_to_geog_etrs();

void geog_mgi_ortho_to_geoc_mgi();
void geog_mgi_ortho_to_geog_mgi();
void geog_mgi_ortho_to_proj_mgi();

void val_geog_mgi_ortho();
// proj_mgi

void proj_mgi_to_geoc_etrs();
void proj_mgi_to_geog_etrs();
void proj_mgi_to_geog_etrs_ortho();

void proj_mgi_to_geoc_mgi();
void proj_mgi_to_geog_mgi();
void proj_mgi_to_geog_mgi_ortho();

void val_proj_mgi();
//-
// 2-way transformation
//-

#endif
