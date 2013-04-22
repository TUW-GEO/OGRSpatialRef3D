#include <stdio.h>
#include <iostream>

#include "cpl_conv.h"
#include "ogr_spatialref3D.h"

int main()
{
  //init gdal/proj.4 data directory path 
	CPLSetConfigOption("GDAL_DATA", "..\\gdal-1.9.1\\distro\\data");

  OGRSpatialReference3D oSourceSRS, oTargetSRS;
  OGRCoordinateTransformation3D *poCT1;
           
  //init coordinate system from epsg code
	oSourceSRS.importFromEPSG( 31491 );	
  oTargetSRS.importFromEPSG( 31492 );

  oSourceSRS.SetGeoidModel("geoid.tif");	//set geoid
  oSourceSRS.SetVCorrModel("vcorr.tif");	//set vertical correction model

  oSourceSRS.SetVScale(0.15); //setting vertical scale
  oSourceSRS.SetVOffset(100); 

	//create coordinate transformation object
  poCT1 = OGRCreateCoordinateTransformation3D( &oSourceSRS,
                                               &oTargetSRS );

	double sourcex = 35.630;
  double sourcey = 47.950;
	double sourcez = 0;
            
  double targetx = sourcex;
  double targety = sourcey;
  double targetz = sourcez;
 
	//do actual transformation
  if( poCT1 == NULL || !poCT1->Transform( 1, &targetx, &targety ,&targetz) )
      printf( "Transformation failed.\n" );
  else
      printf( "(%f,%f,%f) -> (%f,%f,%f)\n", sourcex, sourcey, sourcez, targetx, targety, targetz );
	std::cout << std::endl << "Press <Enter> to end program" << std::endl;
  std::cin.get();
}