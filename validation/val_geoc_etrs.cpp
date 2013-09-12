#include <iostream>

#include "validate.h"

#include "ogr_spatialref3D.h"

using namespace std;

void geoc_etrs_to_geog_etrs()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOC)" << endl;
	cout << "Target coord.: ETRS89 (GEOG) + Ellipsoidal Height" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOC_ETRS);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_ETRS);
	oTargetSRS.importFromWkt3D(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = x_etrs[row_number];
		r1[row_number] = y_etrs[row_number];
		r2[row_number] = z_etrs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-lon_grs[row_number]));
			err1.add(fabs(r1[row_number]-lat_grs[row_number]));
			err2.add(fabs(r2[row_number]-hell_grs[row_number]));
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

void geoc_etrs_to_geog_etrs_ortho()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r3 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOC)" << endl;
	cout << "Target coord.: ETRS89 (GEOG) + Orthometric Height" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOC_ETRS);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_ETRS_ORTH);
	oTargetSRS.importFromWkt3D(&(wkt2));

	oTargetSRS.SetDebug(true);
	oTargetSRS.SetDebugData(r3, 0);

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = x_etrs[row_number];
		r1[row_number] = y_etrs[row_number];
		r2[row_number] = z_etrs[row_number];
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

void geoc_etrs_to_geoc_mgi_2()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T -> S ]------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOC)" << endl;
	cout << "Target coord.: MGI (GEOC)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOC_ETRS);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOC_MGI);
	oTargetSRS.importFromWkt3D(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	OGRCoordinateTransformation3D *poCT_inv = OGRCreateCoordinateTransformation3D( 
													&oTargetSRS, &oSourceSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = x_etrs[row_number];
		r1[row_number] = y_etrs[row_number];
		r2[row_number] = z_etrs[row_number];
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
				err0.add(fabs(r0[row_number]-x_etrs[row_number]));
				err1.add(fabs(r1[row_number]-y_etrs[row_number]));
				err2.add(fabs(r2[row_number]-z_etrs[row_number]));
			}

			cout << "Error (axis 0) : " << endl; err0.printout();
			cout << "Error (axis 1) : " << endl; err1.printout();
			cout << "Error (axis 2) : " << endl; err2.printout();
		}
	}

	delete poCT;
	delete poCT_inv;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
}

void geoc_etrs_to_geoc_mgi()
{
	//geoc_etrs_to_geoc_mgi_2();
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOC)" << endl;
	cout << "Target coord.: MGI (GEOC)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOC_ETRS);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOC_MGI);
	oTargetSRS.importFromWkt3D(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = x_etrs[row_number];
		r1[row_number] = y_etrs[row_number];
		r2[row_number] = z_etrs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-x_mgi[row_number]));
			err1.add(fabs(r1[row_number]-y_mgi[row_number]));
			err2.add(fabs(r2[row_number]-z_mgi[row_number]));
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

void geoc_etrs_to_geog_mgi()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOC)" << endl;
	cout << "Target coord.: MGI (GEOG) + Ellipsoidal Height" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOC_ETRS);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_MGI);
	oTargetSRS.importFromWkt3D(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = x_etrs[row_number];
		r1[row_number] = y_etrs[row_number];
		r2[row_number] = z_etrs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-lon_mgi[row_number]));
			err1.add(fabs(r1[row_number]-lat_mgi[row_number]));
			err2.add(fabs(r2[row_number]-hell_mgi[row_number]));
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

void geoc_etrs_to_geog_mgi_ortho()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOC)" << endl;
	cout << "Target coord.: MGI (GEOG) + Orthometric Height" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOC_ETRS);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_MGI_ORTH);
	oTargetSRS.importFromWkt3D(&(wkt2));

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = x_etrs[row_number];
		r1[row_number] = y_etrs[row_number];
		r2[row_number] = z_etrs[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-lon_mgi[row_number]));
			err1.add(fabs(r1[row_number]-lat_mgi[row_number]));
			err2.add(fabs(r2[row_number]-h_orth[row_number]));
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

void geoc_etrs_to_proj_mgi()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r3 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r4 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS_28, oTargetSRS_31, oTargetSRS_34;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: ETRS89 (GEOC)" << endl;
	cout << "Target coord.: MGI (PROJ) + In-use Height" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOC_ETRS);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(PROJ_MGI_28);
	oTargetSRS_28.importFromWkt3D(&(wkt2));

	wkt2 = loadWktFile(PROJ_MGI_31);
	oTargetSRS_31.importFromWkt3D(&(wkt2));

	wkt2 = loadWktFile(PROJ_MGI_34);
	oTargetSRS_34.importFromWkt3D(&(wkt2));

	OGRCoordinateTransformation3D *poCT_28 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS_28 );
	OGRCoordinateTransformation3D *poCT_31 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS_31 );
	OGRCoordinateTransformation3D *poCT_34 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS_34 );

	oTargetSRS_28.SetDebug(true);
	oTargetSRS_31.SetDebug(true);
	oTargetSRS_34.SetDebug(true);

	if( poCT_28 == NULL || poCT_31 == NULL || poCT_28 == NULL )
	{
		cerr << "Transformation instance creation failed." << endl;

	}
	else
	{
		SummStat err0, err1, err2, err3, err4;
		for(int row_number=0; row_number < num_data; row_number++)
		{
			r0[row_number] = x_etrs[row_number];
			r1[row_number] = y_etrs[row_number];
			r2[row_number] = z_etrs[row_number];

			switch(ms[row_number])
			{
			case 28:
				oTargetSRS_28.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_28->Transform( 1, &(r0[row_number]), &(r1[row_number]), &(r2[row_number]));
				break;
			case 31:
				oTargetSRS_31.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_31->Transform( 1, &(r0[row_number]), &(r1[row_number]), &(r2[row_number]));
				break;
			case 34:
				oTargetSRS_34.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_34->Transform( 1, &(r0[row_number]), &(r1[row_number]), &(r2[row_number]));
				break;
			default:
				cerr << "invalid meridianstrip value" << ms[row_number] << endl;
			}

			err0.add(fabs(r0[row_number]-x_gebr[row_number]));
			err1.add(fabs(r1[row_number]-y_gebr[row_number]));
			err2.add(fabs(r2[row_number]-h_gebr[row_number]));
			err3.add(fabs(r3[row_number]-und_bess[row_number]));
			err4.add(fabs(r4[row_number]-ras_val[row_number]));
		}

		cout << "Error (axis 0) : " << endl; err0.printout();
		cout << "Error (axis 1) : " << endl; err1.printout();
		cout << "Error (axis 2) : " << endl; err2.printout();
		cout << "Error (geoid undulation) : " << endl; err3.printout();
		cout << "Error (height correction) : " << endl; err4.printout();
	}

	delete poCT_28;
	delete poCT_31;
	delete poCT_34;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
	CPLFree(r3);
	CPLFree(r4);
}

void val_geoc_etrs()
{
	geoc_etrs_to_geog_etrs();
	geoc_etrs_to_geog_etrs_ortho();

	geoc_etrs_to_geoc_mgi();
	geoc_etrs_to_geog_mgi();
	geoc_etrs_to_geog_mgi_ortho();
	geoc_etrs_to_proj_mgi();
}