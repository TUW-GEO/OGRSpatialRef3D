#include <iostream>
#include <iomanip>

#include "validate.h"

#include "ogr_spatialref3D.h"

using namespace std;

void geog_etrs_to_geoc_etrs()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOG) + Ellipsoidal Height" << endl;
	cout << "Target coord.: ETRS89 (GEOC)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_ETRS);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOC_ETRS);
	oTargetSRS.importFromWkt3D(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lon_grs[row_number];
		r1[row_number] = lat_grs[row_number];
		r2[row_number] = hell_grs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-x_etrs[row_number]));
			err1.add(fabs(r1[row_number]-y_etrs[row_number]));
			err2.add(fabs(r2[row_number]-z_etrs[row_number]));
		}

		cout << "Error (axis 0) : " << endl; err0.printout();
		cout << "Error (axis 1) : " << endl; err1.printout();
		cout << "Error (axis 2) : " << endl; err2.printout();
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
}

void geog_etrs_to_geog_etrs_ortho()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r3 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOG) + Ellipsoidal Height" << endl;
	cout << "Target coord.: ETRS89 (GEOG) + Orthometric Height" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_ETRS);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_ETRS_ORTH);
	oTargetSRS.importFromWkt3D(&(wkt2));

	oTargetSRS.SetDebug(true);
	oTargetSRS.SetDebugData(r3, 0);

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lon_grs[row_number];
		r1[row_number] = lat_grs[row_number];
		r2[row_number] = hell_grs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2, err3;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-lon_grs[row_number]));
			err1.add(fabs(r1[row_number]-lat_grs[row_number]));
			err2.add(fabs(r2[row_number]-h_orth[row_number]));
			err3.add(fabs(r3[row_number]-und_grs[row_number]));
		}

		cout << "Error (axis 0) : " << endl; err0.printout();
		cout << "Error (axis 1) : " << endl; err1.printout();
		cout << "Error (axis 2) : " << endl; err2.printout();
		cout << "Error (geoid undulation) : " << endl; err3.printout();
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
	CPLFree(r3);
}

void geog_etrs_to_geoc_mgi()
{
	// two way transform since reference geoc_mgi data is not available
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T -> S ]------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOG) + Ellipsoidal Height" << endl;
	cout << "Target coord.: MGI (GEOC)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_ETRS);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOC_MGI);
	oTargetSRS.importFromWkt3D(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	OGRCoordinateTransformation3D *poCT_inv = OGRCreateCoordinateTransformation3D( 
													&oTargetSRS, &oSourceSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lon_grs[row_number];
		r1[row_number] = lat_grs[row_number];
		r2[row_number] = hell_grs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "forward Transformation failed.\n" );
	else
	{
		printf( "forward Transformation successful.\n" );

		if( poCT_inv == NULL || !poCT_inv->Transform( num_data, r0, r1 ,r2) )
			printf( "inverse Transformation failed.\n" );
		else
		{
			printf( "inverse Transformation successful.\n" );
			SummStat err0, err1, err2;

			for(int row_number=0; row_number < num_data; row_number++)
			{
				err0.add(fabs(r0[row_number]-lon_grs[row_number]));
				err1.add(fabs(r1[row_number]-lat_grs[row_number]));
				err2.add(fabs(r2[row_number]-hell_grs[row_number]));
			}

			cout << "Error (axis 0) : " << endl; 
			err0.printout();
			cout << "Error (axis 1) : " << endl; 
			err1.printout();
			cout << "Error (axis 2) : " << endl; 
			err2.printout();
		}
	}

	delete poCT;
	delete poCT_inv;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
}

void geog_etrs_to_geog_mgi_2();

void geog_etrs_to_geog_mgi()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOG) + Ellipsoidal Height" << endl;
	cout << "Target coord.: MGI (GEOG)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_ETRS);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_MGI);
	oTargetSRS.importFromWkt3D(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lon_grs[row_number];
		r1[row_number] = lat_grs[row_number];
		r2[row_number] = hell_grs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-lon_mgi[row_number]));
			err1.add(fabs(r1[row_number]-lat_mgi[row_number]));
		}

		cout << setprecision(8);
		cout << "Error (axis 0) : " << endl; 
		err0.printout();

		cout << setprecision(8);
		cout << "Error (axis 1) : " << endl; 
		err1.printout();
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);

	geog_etrs_to_geog_mgi_2();
}

void geog_etrs_to_geog_mgi_ortho()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r3 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOG) + Ellipsoidal Height" << endl;
	cout << "Target coord.: MGI (GEOG) + Orthometric Height" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_ETRS);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_MGI_ORTH);
	oTargetSRS.importFromWkt3D(&(wkt2));

	oTargetSRS.SetDebug(true);
	oTargetSRS.SetDebugData(r3, 0);

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lon_grs[row_number];
		r1[row_number] = lat_grs[row_number];
		r2[row_number] = hell_grs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2, err3;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-lon_mgi[row_number]));
			err1.add(fabs(r1[row_number]-lat_mgi[row_number]));
			err2.add(fabs(r2[row_number]-h_orth[row_number]));
			err2.add(fabs(r3[row_number]-und_bess[row_number]));
		}

		cout << "Error (axis 0) : " << endl; err0.printout();
		cout << "Error (axis 1) : " << endl; err1.printout();
		cout << "Error (axis 2) : " << endl; err2.printout();
		cout << "Error (geoid undulation) : " << endl; err3.printout();
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
	CPLFree(r3);
}

void geog_etrs_to_proj_mgi()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r3 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r4 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOG) + Ellipsoidal Height" << endl;
	cout << "Target coord.: MGI (PROJ) + In-use Height" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_ETRS);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(PROJ_MGI);
	oTargetSRS.importFromWkt3D(&(wkt2));

	oTargetSRS.SetDebug(true);
	oTargetSRS.SetDebugData(r3, r4);

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lon_grs[row_number];
		r1[row_number] = lat_grs[row_number];
		r2[row_number] = hell_grs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2, err3, err4;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-x_gebr[row_number]));
			err1.add(fabs(r1[row_number]-y_gebr[row_number]));
			err2.add(fabs(r2[row_number]-h_gebr[row_number]));
			err2.add(fabs(r3[row_number]-und_bess[row_number]));
			err2.add(fabs(r4[row_number]-ras_val[row_number]));
		}

		cout << "Error (axis 0) : " << endl; err0.printout();
		cout << "Error (axis 1) : " << endl; err1.printout();
		cout << "Error (axis 2) : " << endl; err2.printout();
		cout << "Error (geoid undulation) : " << endl; err3.printout();
		cout << "Error (height correction) : " << endl; err4.printout();
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
	CPLFree(r3);
	CPLFree(r4);
}

void geog_etrs_to_geog_mgi_2()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T -> S ]------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOG) + Ellipsoidal Height" << endl;
	cout << "Target coord.: MGI (GEOG)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_ETRS);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_MGI);
	oTargetSRS.importFromWkt3D(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	OGRCoordinateTransformation3D *poCT_inv = OGRCreateCoordinateTransformation3D( 
													&oTargetSRS, &oSourceSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lon_grs[row_number];
		r1[row_number] = lat_grs[row_number];
		r2[row_number] = hell_grs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "forward Transformation failed.\n" );
	else
	{
		printf( "forward Transformation successful.\n" );

		if( poCT_inv == NULL || !poCT_inv->Transform( num_data, r0, r1 ,r2) )
			printf( "inverse Transformation failed.\n" );
		else
		{
			printf( "inverse Transformation successful.\n" );
			SummStat err0, err1, err2;

			for(int row_number=0; row_number < num_data; row_number++)
			{
				err0.add(fabs(r0[row_number]-lon_grs[row_number]));
				err1.add(fabs(r1[row_number]-lat_grs[row_number]));
				err2.add(fabs(r2[row_number]-hell_grs[row_number]));
			}

			cout << "Error (axis 0) : " << endl; 
			err0.printout();
			cout << "Error (axis 1) : " << endl; 
			err1.printout();
			cout << "Error (axis 2) : " << endl; 
			err2.printout();
		}
	}

	delete poCT;
	delete poCT_inv;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
}

void val_geog_etrs()
{
	geog_etrs_to_geoc_etrs();
	geog_etrs_to_geog_etrs_ortho();

	geog_etrs_to_geoc_mgi();//not available
	geog_etrs_to_geog_mgi();
	geog_etrs_to_geog_mgi_ortho();//error
	geog_etrs_to_proj_mgi();//error
}