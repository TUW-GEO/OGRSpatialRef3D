#include <stdio.h>
#include <iostream>

#include "ogr_spatialref.h"
#include "ogr_spatialref3D.h"
#include "cpl_conv.h"

int main()
{
  CPLSetConfigOption("GDAL_DATA", "..\\gdal-1.9.1\\distro\\data");

  /*OGRSpatialReference oSourceSRS, oTargetSRS;
  OGRCoordinateTransformation *poCT;
            
  oSourceSRS.importFromEPSG( 31491 );
  oTargetSRS.importFromEPSG( 31492 );

  double sourcex = 0;
  double sourcey = 0;
            
  poCT = OGRCreateCoordinateTransformation( &oSourceSRS,
                                            &oTargetSRS );

  double targetx = sourcex;
  double targety = sourcey;
            
  if( poCT == NULL || !poCT->Transform( 1, &targetx, &targety ) )
      printf( "Transformation failed.\n" );
  else
      printf( "(%f,%f) -> (%f,%f)\n", sourcex, sourcey, targetx, targety );

	std::cout << std::endl << "Press <Enter> to end program" << std::endl;
  std::cin.get();*/

  OGRSpatialReference3D oSourceSRS, oTargetSRS;
  OGRCoordinateTransformation3D *poCT1;
            


  oSourceSRS.importFromEPSG( 31491 );
  oTargetSRS.importFromEPSG( 31492 );

  double sourcex = 0;
  double sourcey = 0;

  poCT1=OGRCreateCoordinateTransformation3D( &oSourceSRS,
                                            &oTargetSRS );
            

  double targetx = sourcex;
  double targety = sourcey;
            

   if( poCT1 == NULL || !poCT1->Transform( 1, &targetx, &targety ) )
      printf( "Transformation failed.\n" );
  else
      printf( "(%f,%f) -> (%f,%f)\n", sourcex, sourcey, targetx, targety );
	std::cout << std::endl << "Press <Enter> to end program" << std::endl;
  std::cin.get();

  oSourceSRS.SetGeoidModel("geoid.tif");

  oSourceSRS.SetVCorrModel("vcorr.tif");

  oSourceSRS.SetVScale(0.15); //setting vertical scale
  oSourceSRS.SetVOffset(100); //setting vertical offset
  oSourceSRS.transform();
}