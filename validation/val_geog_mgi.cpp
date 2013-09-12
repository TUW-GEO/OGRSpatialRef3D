#include <iostream>

#include "validate.h"

#include "ogr_spatialref3D.h"

using namespace std;

void geog_mgi_to_geoc_etrs()
{

}

void geog_mgi_to_geog_etrs(){}
void geog_mgi_to_geog_etrs_ortho(){}

void geog_mgi_to_geoc_mgi(){}

void geog_mgi_to_geog_mgi_ortho(){}
void geog_mgi_to_proj_mgi_gebr(){}

void val_geog_mgi()
{
	geog_mgi_to_geoc_etrs();

	geog_mgi_to_geog_etrs();
	geog_mgi_to_geog_etrs_ortho();

	geog_mgi_to_geoc_mgi();

	geog_mgi_to_geog_mgi_ortho();
	geog_mgi_to_proj_mgi_gebr();
}