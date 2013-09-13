#include <iostream>

#include "validate.h"

#include "ogr_spatialref3D.h"

using namespace std;

void proj_mgi_to_geoc_etrs()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r3 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r4 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS_28, oSourceSRS_31, oSourceSRS_34, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: MGI (PROJ) + In-use Height" << endl;
	cout << "Target coord.: ETRS89 (GEOC)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt2 = loadWktFile(GEOC_ETRS);
	oTargetSRS.importFromWkt3D(&(wkt2));

	char *wkt1 = loadWktFile(PROJ_MGI_28);
	oSourceSRS_28.importFromWkt3D(&(wkt1));

	wkt1 = loadWktFile(PROJ_MGI_31);
	oSourceSRS_31.importFromWkt3D(&(wkt1));

	wkt1 = loadWktFile(PROJ_MGI_34);
	oSourceSRS_34.importFromWkt3D(&(wkt1));

	oSourceSRS_28.SetDebug(true);
	oSourceSRS_31.SetDebug(true);
	oSourceSRS_34.SetDebug(true);

	OGRCoordinateTransformation3D *poCT_28 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS_28, &oTargetSRS );
	OGRCoordinateTransformation3D *poCT_31 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS_31, &oTargetSRS );
	OGRCoordinateTransformation3D *poCT_34 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS_34, &oTargetSRS );


	if( poCT_28 == NULL || poCT_31 == NULL || poCT_34 == NULL )
		printf( "Transformaer creation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2, err3, err4;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			r0[row_number] = x_gebr[row_number];
			r1[row_number] = y_gebr[row_number];
			r2[row_number] = h_gebr[row_number];

			switch(ms[row_number])
			{
			case 28:
				oSourceSRS_28.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_28->Transform( 1, &(r0[row_number]), &(r1[row_number]) ,&(r2[row_number]));
				break;
			case 31:
				oSourceSRS_31.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_31->Transform( 1, &(r0[row_number]), &(r1[row_number]) ,&(r2[row_number]));
				break;
			case 34:
				oSourceSRS_34.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_34->Transform( 1, &(r0[row_number]), &(r1[row_number]) ,&(r2[row_number]));
				break;
			default:
				cerr << "invalid meridian strip value" << endl;
			}

			err0.add(fabs(r0[row_number]-x_etrs[row_number]));
			err1.add(fabs(r1[row_number]-y_etrs[row_number]));
			err2.add(fabs(r2[row_number]-z_etrs[row_number]));
			err3.add(fabs(r3[row_number]-und_bess[row_number]));
			err4.add(fabs(r4[row_number]-ras_val[row_number]));
		}

		cout << "Error (axis 0) : " << endl; err0.printout();
		cout << "Error (axis 1) : " << endl; err1.printout();
		cout << "Error (axis 2) : " << endl; err2.printout();
		cout << "Error (source geoid undulation) : " << endl; err3.printout();
		cout << "Error (source height correction) : " << endl; err4.printout();
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

void proj_mgi_to_geog_etrs()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r3 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r4 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS_28, oSourceSRS_31, oSourceSRS_34, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: MGI (PROJ) + In-use Height" << endl;
	cout << "Target coord.: ETRS89 (GEOG) + Ellipsoidal Height" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt2 = loadWktFile(GEOG_ETRS);
	oTargetSRS.importFromWkt3D(&(wkt2));

	char *wkt1 = loadWktFile(PROJ_MGI_28);
	oSourceSRS_28.importFromWkt3D(&(wkt1));

	wkt1 = loadWktFile(PROJ_MGI_31);
	oSourceSRS_31.importFromWkt3D(&(wkt1));

	wkt1 = loadWktFile(PROJ_MGI_34);
	oSourceSRS_34.importFromWkt3D(&(wkt1));

	oSourceSRS_28.SetDebug(true);
	oSourceSRS_31.SetDebug(true);
	oSourceSRS_34.SetDebug(true);

	OGRCoordinateTransformation3D *poCT_28 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS_28, &oTargetSRS );
	OGRCoordinateTransformation3D *poCT_31 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS_31, &oTargetSRS );
	OGRCoordinateTransformation3D *poCT_34 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS_34, &oTargetSRS );


	if( poCT_28 == NULL || poCT_31 == NULL || poCT_34 == NULL )
		printf( "Transformaer creation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2, err3, err4;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			r0[row_number] = x_gebr[row_number];
			r1[row_number] = y_gebr[row_number];
			r2[row_number] = h_gebr[row_number];

			switch(ms[row_number])
			{
			case 28:
				oSourceSRS_28.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_28->Transform( 1, &(r0[row_number]), &(r1[row_number]) ,&(r2[row_number]));
				break;
			case 31:
				oSourceSRS_31.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_31->Transform( 1, &(r0[row_number]), &(r1[row_number]) ,&(r2[row_number]));
				break;
			case 34:
				oSourceSRS_34.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_34->Transform( 1, &(r0[row_number]), &(r1[row_number]) ,&(r2[row_number]));
				break;
			default:
				cerr << "invalid meridian strip value" << endl;
			}

			err0.add(fabs(r0[row_number]-lon_grs[row_number]));
			err1.add(fabs(r1[row_number]-lat_grs[row_number]));
			err2.add(fabs(r2[row_number]-hell_grs[row_number]));
			err3.add(fabs(r3[row_number]-und_bess[row_number]));
			err4.add(fabs(r4[row_number]-ras_val[row_number]));
		}

		cout << "Error (axis 0) : " << endl; err0.printout();
		cout << "Error (axis 1) : " << endl; err1.printout();
		cout << "Error (axis 2) : " << endl; err2.printout();
		cout << "Error (source geoid undulation) : " << endl; err3.printout();
		cout << "Error (source height correction) : " << endl; err4.printout();
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

void proj_mgi_to_geog_etrs_ortho()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r3 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r4 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r5 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS_28, oSourceSRS_31, oSourceSRS_34, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: MGI (PROJ) + In-use Height" << endl;
	cout << "Target coord.: ETRS89 (GEOG) + Orthometric Height" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt2 = loadWktFile(GEOG_ETRS_ORTH);
	oTargetSRS.importFromWkt3D(&(wkt2));

	char *wkt1 = loadWktFile(PROJ_MGI_28);
	oSourceSRS_28.importFromWkt3D(&(wkt1));

	wkt1 = loadWktFile(PROJ_MGI_31);
	oSourceSRS_31.importFromWkt3D(&(wkt1));

	wkt1 = loadWktFile(PROJ_MGI_34);
	oSourceSRS_34.importFromWkt3D(&(wkt1));

	oSourceSRS_28.SetDebug(true);
	oSourceSRS_31.SetDebug(true);
	oSourceSRS_34.SetDebug(true);
	oTargetSRS.SetDebug(true);

	OGRCoordinateTransformation3D *poCT_28 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS_28, &oTargetSRS );
	OGRCoordinateTransformation3D *poCT_31 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS_31, &oTargetSRS );
	OGRCoordinateTransformation3D *poCT_34 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS_34, &oTargetSRS );


	if( poCT_28 == NULL || poCT_31 == NULL || poCT_34 == NULL )
		printf( "Transformaer creation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2, err3, err4, err5;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			r0[row_number] = x_gebr[row_number];
			r1[row_number] = y_gebr[row_number];
			r2[row_number] = h_gebr[row_number];

			oTargetSRS.SetDebugData(&(r5[row_number]), 0);
			switch(ms[row_number])
			{
			case 28:
				oSourceSRS_28.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_28->Transform( 1, &(r0[row_number]), &(r1[row_number]) ,&(r2[row_number]));
				break;
			case 31:
				oSourceSRS_31.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_31->Transform( 1, &(r0[row_number]), &(r1[row_number]) ,&(r2[row_number]));
				break;
			case 34:
				oSourceSRS_34.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_34->Transform( 1, &(r0[row_number]), &(r1[row_number]) ,&(r2[row_number]));
				break;
			default:
				cerr << "invalid meridian strip value" << endl;
			}

			err0.add(fabs(r0[row_number]-lon_grs[row_number]));
			err1.add(fabs(r1[row_number]-lat_grs[row_number]));
			err2.add(fabs(r2[row_number]-h_orth[row_number]));
			err3.add(fabs(r3[row_number]-und_bess[row_number]));
			err4.add(fabs(r4[row_number]-ras_val[row_number]));
			err5.add(fabs(r5[row_number]-und_grs[row_number]));
		}

		cout << "Error (axis 0) : " << endl; err0.printout();
		cout << "Error (axis 1) : " << endl; err1.printout();
		cout << "Error (axis 2) : " << endl; err2.printout();
		cout << "Error (source geoid undulation) : " << endl; err3.printout();
		cout << "Error (source height correction) : " << endl; err4.printout();
		cout << "Error (target geoid undulation) : " << endl; err5.printout();
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

void proj_mgi_to_geoc_mgi()
{

}

void proj_mgi_to_geog_mgi()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r3 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r4 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS_28, oSourceSRS_31, oSourceSRS_34, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: MGI (PROJ) + In-use Height" << endl;
	cout << "Target coord.: MGI (GEOG)" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt2 = loadWktFile(GEOG_MGI);
	oTargetSRS.importFromWkt3D(&(wkt2));

	char *wkt1 = loadWktFile(PROJ_MGI_28);
	oSourceSRS_28.importFromWkt3D(&(wkt1));

	wkt1 = loadWktFile(PROJ_MGI_31);
	oSourceSRS_31.importFromWkt3D(&(wkt1));

	wkt1 = loadWktFile(PROJ_MGI_34);
	oSourceSRS_34.importFromWkt3D(&(wkt1));

	oSourceSRS_28.SetDebug(true);
	oSourceSRS_31.SetDebug(true);
	oSourceSRS_34.SetDebug(true);

	OGRCoordinateTransformation3D *poCT_28 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS_28, &oTargetSRS );
	OGRCoordinateTransformation3D *poCT_31 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS_31, &oTargetSRS );
	OGRCoordinateTransformation3D *poCT_34 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS_34, &oTargetSRS );


	if( poCT_28 == NULL || poCT_31 == NULL || poCT_34 == NULL )
		printf( "Transformaer creation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2, err3, err4;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			r0[row_number] = x_gebr[row_number];
			r1[row_number] = y_gebr[row_number];
			r2[row_number] = h_gebr[row_number];

			switch(ms[row_number])
			{
			case 28:
				oSourceSRS_28.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_28->Transform( 1, &(r0[row_number]), &(r1[row_number]) ,&(r2[row_number]));
				break;
			case 31:
				oSourceSRS_31.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_31->Transform( 1, &(r0[row_number]), &(r1[row_number]) ,&(r2[row_number]));
				break;
			case 34:
				oSourceSRS_34.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_34->Transform( 1, &(r0[row_number]), &(r1[row_number]) ,&(r2[row_number]));
				break;
			default:
				cerr << "invalid meridian strip value" << endl;
			}

			err0.add(fabs(r0[row_number]-lon_mgi[row_number]));
			err1.add(fabs(r1[row_number]-lat_mgi[row_number]));
			err2.add(fabs(r2[row_number]-hell_mgi[row_number]));
			err3.add(fabs(r3[row_number]-und_bess[row_number]));
			err4.add(fabs(r4[row_number]-ras_val[row_number]));
		}

		cout << "Error (axis 0) : " << endl; err0.printout();
		cout << "Error (axis 1) : " << endl; err1.printout();
		cout << "Error (axis 2) : " << endl; err2.printout();
		cout << "Error (source geoid undulation) : " << endl; err3.printout();
		cout << "Error (source height correction) : " << endl; err4.printout();
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

void proj_mgi_to_geog_mgi_ortho()
{
	double *r0 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r1 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r2 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r3 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r4 = (double*)CPLMalloc(sizeof(double)*num_data);
	double *r5 = (double*)CPLMalloc(sizeof(double)*num_data);

	OGRSpatialReference3D oSourceSRS_28, oSourceSRS_31, oSourceSRS_34, oTargetSRS;

	cout << "----------------[ S -> T ]-----------------------" << endl;
	cout << "Source coord.: MGI (PROJ) + In-use Height" << endl;
	cout << "Target coord.: MGI (GEOG) + Orthometric Height" << endl;
	cout << "-------------------------------------------------" << endl;

	char *wkt2 = loadWktFile(GEOG_MGI_ORTH);
	oTargetSRS.importFromWkt3D(&(wkt2));

	oTargetSRS.SetDebug(true);
	
	char *wkt1 = loadWktFile(PROJ_MGI_28);
	oSourceSRS_28.importFromWkt3D(&(wkt1));

	wkt1 = loadWktFile(PROJ_MGI_31);
	oSourceSRS_31.importFromWkt3D(&(wkt1));

	wkt1 = loadWktFile(PROJ_MGI_34);
	oSourceSRS_34.importFromWkt3D(&(wkt1));

	oSourceSRS_28.SetDebug(true);
	oSourceSRS_31.SetDebug(true);
	oSourceSRS_34.SetDebug(true);

	OGRCoordinateTransformation3D *poCT_28 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS_28, &oTargetSRS );
	OGRCoordinateTransformation3D *poCT_31 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS_31, &oTargetSRS );
	OGRCoordinateTransformation3D *poCT_34 = OGRCreateCoordinateTransformation3D( 
													&oSourceSRS_34, &oTargetSRS );


	if( poCT_28 == NULL || poCT_31 == NULL || poCT_34 == NULL )
		printf( "Transformaer creation failed.\n" );
	else
	{
		printf( "Transformation successful.\n" );
		SummStat err0, err1, err2, err3, err4, err5;

		for(int row_number=0; row_number < num_data; row_number++)
		{
			r0[row_number] = x_gebr[row_number];
			r1[row_number] = y_gebr[row_number];
			r2[row_number] = h_gebr[row_number];

			oTargetSRS.SetDebugData(&(r5[row_number]), 0);
			switch(ms[row_number])
			{
			case 28:
				oSourceSRS_28.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_28->Transform( 1, &(r0[row_number]), &(r1[row_number]) ,&(r2[row_number]));
				break;
			case 31:
				oSourceSRS_31.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_31->Transform( 1, &(r0[row_number]), &(r1[row_number]) ,&(r2[row_number]));
				break;
			case 34:
				oSourceSRS_34.SetDebugData(&(r3[row_number]), &(r4[row_number]));
				poCT_34->Transform( 1, &(r0[row_number]), &(r1[row_number]) ,&(r2[row_number]));
				break;
			default:
				cerr << "invalid meridian strip value" << endl;
			}

			err0.add(fabs(r0[row_number]-lon_mgi[row_number]));
			err1.add(fabs(r1[row_number]-lat_mgi[row_number]));
			err2.add(fabs(r2[row_number]-h_orth[row_number]));
			err3.add(fabs(r3[row_number]-und_bess[row_number]));
			err4.add(fabs(r4[row_number]-ras_val[row_number]));
			err5.add(fabs(r5[row_number]-und_bess[row_number]));
		}

		cout << "Error (axis 0) : " << endl; err0.printout();
		cout << "Error (axis 1) : " << endl; err1.printout();
		cout << "Error (axis 2) : " << endl; err2.printout();
		cout << "Error (source geoid undulation) : " << endl; err3.printout();
		cout << "Error (source height correction) : " << endl; err4.printout();
		cout << "Error (target geoid undulation) : " << endl; err5.printout();
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

void val_proj_mgi()
{
	proj_mgi_to_geoc_etrs();
	proj_mgi_to_geog_etrs();
	proj_mgi_to_geog_etrs_ortho();

	proj_mgi_to_geoc_mgi();
	proj_mgi_to_geog_mgi();
	proj_mgi_to_geog_mgi_ortho();	
}