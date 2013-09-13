#include <iostream>

#include "validate.h"

#include "ogr_spatialref3D.h"

using namespace std;

void geog_mgi_ortho_to_geoc_etrs()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r3 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: MGI (GEOG) + Orthometric Height" << endl;
	cout << "Target coord.: ETRS89 (GEOC)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_MGI_ORTH);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOC_ETRS);
	oTargetSRS.importFromWkt3D(&(wkt2));

	oSourceSRS.SetDebug(true);
	oSourceSRS.SetDebugData(r3, 0);

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lon_mgi[row_number];
		r1[row_number] = lat_mgi[row_number];
		r2[row_number] = h_orth[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2, err3;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-x_etrs[row_number]));
			err1.add(fabs(r1[row_number]-y_etrs[row_number]));
			err2.add(fabs(r2[row_number]-z_etrs[row_number]));
			err3.add(fabs(r3[row_number]-und_bess[row_number]));
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

void geog_mgi_ortho_to_geog_etrs_ortho()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r3 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r4 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: MGI (GEOG) + Orthometric Height" << endl;
	cout << "Target coord.: ETRS89 (GEOC) + Orthometric Height" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_MGI_ORTH);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_ETRS_ORTH);
	oTargetSRS.importFromWkt3D(&(wkt2));

	oSourceSRS.SetDebug(true);
	oSourceSRS.SetDebugData(r3, 0);

	oTargetSRS.SetDebug(true);
	oTargetSRS.SetDebugData(r4, 0);

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lon_mgi[row_number];
		r1[row_number] = lat_mgi[row_number];
		r2[row_number] = h_orth[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2, err3, err4;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-lon_grs[row_number]));
			err1.add(fabs(r1[row_number]-lat_grs[row_number]));
			err2.add(fabs(r2[row_number]-h_orth[row_number]));
			err3.add(fabs(r3[row_number]-und_bess[row_number]));
			err4.add(fabs(r4[row_number]-und_grs[row_number]));
		}

		cout << "Error (axis 0) : " << endl; err0.printout();
		cout << "Error (axis 1) : " << endl; err1.printout();
		cout << "Error (axis 2) : " << endl; err2.printout();
		cout << "Error (source geoid undulation) : " << endl; err3.printout();
		cout << "Error (target geoid undulation) : " << endl; err4.printout();
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
	CPLFree(r3);
	CPLFree(r4);
}

void geog_mgi_ortho_to_geog_etrs()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r3 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: MGI (GEOG) + Orthometric Height" << endl;
	cout << "Target coord.: ETRS89 (GEOC) + Ellipsoidal Height" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_MGI_ORTH);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_ETRS);
	oTargetSRS.importFromWkt3D(&(wkt2));

	oSourceSRS.SetDebug(true);
	oSourceSRS.SetDebugData(r3, 0);

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lon_mgi[row_number];
		r1[row_number] = lat_mgi[row_number];
		r2[row_number] = h_orth[row_number];
	}

	if( poCT == NULL || !poCT->Transform( num_data, r0, r1 ,r2) )
		printf( "Transformation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2, err3, err4;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			err0.add(fabs(r0[row_number]-lon_grs[row_number]));
			err1.add(fabs(r1[row_number]-lat_grs[row_number]));
			err2.add(fabs(r2[row_number]-hell_grs[row_number]));
			err3.add(fabs(r3[row_number]-und_bess[row_number]));
		}

		cout << "Error (axis 0) : " << endl; err0.printout();
		cout << "Error (axis 1) : " << endl; err1.printout();
		cout << "Error (axis 2) : " << endl; err2.printout();
		cout << "Error (source geoid undulation) : " << endl; err3.printout();
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
	CPLFree(r3);
}

void geog_mgi_ortho_to_geoc_mgi()
{
	//not-available : geoc_mgi reference data is calculated from geog_mgi
}

void geog_mgi_ortho_to_geog_mgi()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r3 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: MGI (GEOG) + Orthometric Height" << endl;
	cout << "Target coord.: MGI (GEOG)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_MGI_ORTH);
	oSourceSRS.importFromWkt3D(&(wkt1));

	char *wkt2 = loadWktFile(GEOG_MGI);
	oTargetSRS.importFromWkt3D(&(wkt2));

	oSourceSRS.SetDebug(true);
	oSourceSRS.SetDebugData(r3, 0);

	OGRCoordinateTransformation3D *poCT = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS, &oTargetSRS );

	for(int row_number=0; row_number < num_data; row_number++)
	{
		r0[row_number] = lon_mgi[row_number];
		r1[row_number] = lat_mgi[row_number];
		r2[row_number] = h_orth[row_number];
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
			//err2.add(fabs(r2[row_number]-h_gebr[row_number]));
			err3.add(fabs(r3[row_number]-und_bess[row_number]));
		}

		cout << "Error (axis 0) : " << endl; err0.printout();
		cout << "Error (axis 1) : " << endl; err1.printout();
		//cout << "Error (axis 2) : " << endl; err2.printout();
		cout << "Error (source geoid undulation) : " << endl; err3.printout();
	}

	delete poCT;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
	CPLFree(r3);
}

void geog_mgi_ortho_to_proj_mgi()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r3 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r4 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r5 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS, oTargetSRS_28, oTargetSRS_31, oTargetSRS_34;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: MGI (GEOG) + Orthometric Height" << endl;
	cout << "Target coord.: MGI (PROJ) + In-use Height" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt1 = loadWktFile(GEOG_MGI_ORTH);
	oSourceSRS.importFromWkt3D(&(wkt1));

	oSourceSRS.SetDebug(true);

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
	

	if( poCT_28 == NULL || poCT_31 == NULL || poCT_34 == NULL )
	{
		cerr << "Transformation instance creation failed." << endl;

	}
	else
	{
		SummStat err0, err1, err2, err3, err4, err5;
		for(int row_number=0; row_number < num_data; row_number++)
		{
			r0[row_number] = lon_mgi[row_number];
			r1[row_number] = lat_mgi[row_number];
			r2[row_number] = h_orth[row_number];

			oSourceSRS.SetDebugData(&(r5[row_number]), 0);
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
			err3.add(fabs(r5[row_number]-und_bess[row_number]));
			err4.add(fabs(r3[row_number]-und_bess[row_number]));
			err5.add(fabs(r4[row_number]-ras_val[row_number]));
		}

		cout << "Error (axis 0) : " << endl; err0.printout();
		cout << "Error (axis 1) : " << endl; err1.printout();
		cout << "Error (axis 2) : " << endl; err2.printout();
		cout << "Error (source geoid undulation) : " << endl; err5.printout();
		cout << "Error (target geoid undulation) : " << endl; err3.printout();
		cout << "Error (target height correction) : " << endl; err4.printout();
	}

	delete poCT_28;
	delete poCT_31;
	delete poCT_34;
	CPLFree(r0);
	CPLFree(r1);
	CPLFree(r2);
	CPLFree(r3);
	CPLFree(r4);
	CPLFree(r5);
}

void val_geog_mgi_ortho()
{
	geog_mgi_ortho_to_geoc_etrs();
	geog_mgi_ortho_to_geog_etrs();
	geog_mgi_ortho_to_geog_etrs_ortho();

	geog_mgi_ortho_to_geoc_mgi();
	geog_mgi_ortho_to_geog_mgi();
	geog_mgi_ortho_to_proj_mgi();
}