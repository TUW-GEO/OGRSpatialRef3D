#include <stdio.h>
#include <iostream>

#include "ogr_spatialref.h"
#include "cpl_conv.h"

int main()
{
  CPLSetConfigOption("GDAL_DATA", "..\\gdal-1.9.1\\distro\\data");

  OGRSpatialReference oSourceSRS, oTargetSRS;
  OGRCoordinateTransformation *poCT,*poCT1;
            


  oSourceSRS.importFromEPSG( 31491 );
  oTargetSRS.importFromEPSG( 31492 );

  double sourcex = 0;
  double sourcey = 0;

  poCT1=OGRCreateCoordinateTransformation3D( &oSourceSRS,
                                            &oTargetSRS );
            
  poCT = OGRCreateCoordinateTransformation( &oSourceSRS,
                                            &oTargetSRS );



  double targetx = sourcex;
  double targety = sourcey;
            

 
 /* if( poCT == NULL || !poCT->Transform( 1, &targetx, &targety ) )
      printf( "Transformation failed.\n" );
  else
      printf( "(%f,%f) -> (%f,%f)\n", sourcex, sourcey, targetx, targety );*/

   if( poCT1 == NULL || !poCT1->Transform( 1, &targetx, &targety ) )
      printf( "Transformation failed.\n" );
  else
      printf( "(%f,%f) -> (%f,%f)\n", sourcex, sourcey, targetx, targety );
	std::cout << std::endl << "Press <Enter> to end program" << std::endl;
  std::cin.get();
}