/******************************************************************************
 * $Id: ogr_mssqlspatial.h 24331 2012-04-28 11:49:05Z tamas $
 *
 * Project:  MSSQL Spatial driver
 * Purpose:  Definition of classes for OGR MSSQL Spatial driver.
 * Author:   Tamas Szekeres, szekerest at gmail.com
 *
 ******************************************************************************
 * Copyright (c) 2010, Tamas Szekeres
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

#ifndef _OGR_MSSQLSPATIAL_H_INCLUDED
#define _OGR_MSSQLSPATIAL_H_INCLUDED

#include "ogrsf_frmts.h"
#include "cpl_odbc.h"
#include "cpl_error.h"

class OGRMSSQLSpatialDataSource;

/* geometry format to transfer geometry column */
#define MSSQLGEOMETRY_NATIVE 0
#define MSSQLGEOMETRY_WKB 1
#define MSSQLGEOMETRY_WKT 2

/* geometry column types */
#define MSSQLCOLTYPE_GEOMETRY  0
#define MSSQLCOLTYPE_GEOGRAPHY 1
#define MSSQLCOLTYPE_BINARY 2
#define MSSQLCOLTYPE_TEXT 3

/************************************************************************/
/*                         OGRMSSQLAppendEscaped( )                     */
/************************************************************************/

void OGRMSSQLAppendEscaped( CPLODBCStatement* poStatement, const char* pszStrValue);

/************************************************************************/
/*                           OGRMSSQLGeometryParser                     */
/************************************************************************/

class OGRMSSQLGeometryValidator
{
protected:
    int bIsValid;
    OGRGeometry*    poValidGeometry;
    OGRGeometry*    poOriginalGeometry;

public:
                    OGRMSSQLGeometryValidator(OGRGeometry* poGeom);
                    ~OGRMSSQLGeometryValidator();

    int             ValidatePoint(OGRPoint * poGeom);
    int             ValidateMultiPoint(OGRMultiPoint * poGeom);
    int             ValidateLineString(OGRLineString * poGeom);
    int             ValidateLinearRing(OGRLinearRing * poGeom);
    int             ValidateMultiLineString(OGRMultiLineString * poGeom);
    int             ValidatePolygon(OGRPolygon* poGeom);
    int             ValidateMultiPolygon(OGRMultiPolygon* poGeom);
    int             ValidateGeometryCollection(OGRGeometryCollection* poGeom);
    int             ValidateGeometry(OGRGeometry* poGeom);

    OGRGeometry*    GetValidGeometryRef();
    int             IsValid() { return bIsValid; };
};

/************************************************************************/
/*                           OGRMSSQLGeometryParser                     */
/************************************************************************/

class OGRMSSQLGeometryParser
{
protected:    
    unsigned char* pszData;
    /* serialization propeties */
    char chProps;
    /* point array */
    int nPointSize;
    int nPointPos;
    int nNumPoints;
    /* figure array */
    int nFigurePos;
    int nNumFigures;
    /* shape array */
    int nShapePos;
    int nNumShapes;
    int nSRSId;
    /* geometry or geography */
    int nColType;

protected:
    OGRPoint*           ReadPoint(int iShape);
    OGRMultiPoint*      ReadMultiPoint(int iShape);
    OGRLineString*      ReadLineString(int iShape);
    OGRMultiLineString* ReadMultiLineString(int iShape);
    OGRPolygon*         ReadPolygon(int iShape);
    OGRMultiPolygon*    ReadMultiPolygon(int iShape);
    OGRGeometryCollection* ReadGeometryCollection(int iShape);

public:
                        OGRMSSQLGeometryParser( int nGeomColumnType );
    OGRErr              ParseSqlGeometry(unsigned char* pszInput, int nLen,
                                                        OGRGeometry **poGeom);
    int                 GetSRSId() { return nSRSId; };
};


/************************************************************************/
/*                             OGRMSSQLSpatialLayer                     */
/************************************************************************/

class OGRMSSQLSpatialLayer : public OGRLayer
{
    protected:
    OGRFeatureDefn     *poFeatureDefn;

    CPLODBCStatement   *poStmt;

    // Layer spatial reference system, and srid.
    OGRSpatialReference *poSRS;
    int                 nSRSId;

    int                 iNextShapeId;

    OGRMSSQLSpatialDataSource    *poDS;

    int                nGeomColumnType;
    char               *pszGeomColumn;
    char               *pszFIDColumn;

    int                *panFieldOrdinals;

    CPLErr              BuildFeatureDefn( const char *pszLayerName,
                                          CPLODBCStatement *poStmt );

    virtual CPLODBCStatement *  GetStatement() { return poStmt; }

  public:
                        OGRMSSQLSpatialLayer();
    virtual             ~OGRMSSQLSpatialLayer();

    virtual void        ResetReading();
    virtual OGRFeature *GetNextRawFeature();
    virtual OGRFeature *GetNextFeature();

    virtual OGRFeature *GetFeature( long nFeatureId );
    
    OGRFeatureDefn *    GetLayerDefn() { return poFeatureDefn; }

    virtual OGRSpatialReference *GetSpatialRef();

    virtual OGRErr     StartTransaction();
    virtual OGRErr     CommitTransaction();
    virtual OGRErr     RollbackTransaction();

    virtual const char *GetFIDColumn();
    virtual const char *GetGeometryColumn();

    virtual int         TestCapability( const char * );
    char*               GByteArrayToHexString( const GByte* pabyData, int nLen);
};

/************************************************************************/
/*                       OGRMSSQLSpatialTableLayer                      */
/************************************************************************/

class OGRMSSQLSpatialTableLayer : public OGRMSSQLSpatialLayer
{
    int                 bUpdateAccess;
    int                 bLaunderColumnNames;
    int                 bPreservePrecision;

    char                *pszQuery;

    void		ClearStatement();
    CPLODBCStatement* BuildStatement(const char* pszColumns);

    CPLString BuildFields();

    virtual CPLODBCStatement *  GetStatement();

    char               *pszTableName;
    char               *pszSchemaName;

  public:
                        OGRMSSQLSpatialTableLayer( OGRMSSQLSpatialDataSource * );
                        ~OGRMSSQLSpatialTableLayer();

    CPLErr              Initialize( const char *pszSchema,
                                    const char *pszTableName, 
                                    const char *pszGeomCol, 
                                    int nCoordDimension, 
                                    int nSRId,
                                    OGRwkbGeometryType eType);

    OGRErr              CreateSpatialIndex();
    void                DropSpatialIndex();

    virtual void        ResetReading();
    virtual int         GetFeatureCount( int );

    virtual OGRErr      SetAttributeFilter( const char * );

    virtual OGRErr      SetFeature( OGRFeature *poFeature );
    virtual OGRErr      DeleteFeature( long nFID );
    virtual OGRErr      CreateFeature( OGRFeature *poFeature );

    const char*         GetTableName() { return pszTableName; }
    const char*         GetSchemaName() { return pszSchemaName; }

    virtual OGRErr      CreateField( OGRFieldDefn *poField,
                                     int bApproxOK = TRUE );
   
    virtual OGRFeature *GetFeature( long nFeatureId );

    virtual int         TestCapability( const char * );

    void                SetLaunderFlag( int bFlag ) 
                                { bLaunderColumnNames = bFlag; }
    void                SetPrecisionFlag( int bFlag )
                                { bPreservePrecision = bFlag; }
    void                AppendFieldValue(CPLODBCStatement *poStatement,
                                       OGRFeature* poFeature, int i);

    int                 FetchSRSId();
};

/************************************************************************/
/*                      OGRMSSQLSpatialSelectLayer                      */
/************************************************************************/

class OGRMSSQLSpatialSelectLayer : public OGRMSSQLSpatialLayer
{
    char                *pszBaseStatement;

    void		ClearStatement();
    OGRErr              ResetStatement();

    virtual CPLODBCStatement *  GetStatement();

  public:
                        OGRMSSQLSpatialSelectLayer( OGRMSSQLSpatialDataSource *, 
                                           CPLODBCStatement * );
                        ~OGRMSSQLSpatialSelectLayer();

    virtual void        ResetReading();
    virtual int         GetFeatureCount( int );

    virtual OGRFeature *GetFeature( long nFeatureId );
    
    virtual OGRErr      GetExtent(OGREnvelope *psExtent, int bForce = TRUE);

    virtual int         TestCapability( const char * );
};

/************************************************************************/
/*                           OGRODBCDataSource                          */
/************************************************************************/

class OGRMSSQLSpatialDataSource : public OGRDataSource
{
    OGRMSSQLSpatialTableLayer    **papoLayers;
    int                 nLayers;
    
    char               *pszName;

    char               *pszCatalog;

    int                 bDSUpdate;
    CPLODBCSession      oSession;

    int                 nGeometryFormat;

    // We maintain a list of known SRID to reduce the number of trips to
    // the database to get SRSes. 
    int                 nKnownSRID;
    int                *panSRID;
    OGRSpatialReference **papoSRS;
    
  public:
                        OGRMSSQLSpatialDataSource();
                        ~OGRMSSQLSpatialDataSource();

    const char          *GetCatalog() { return pszCatalog; }

    int                 ParseValue(char** pszValue, char* pszSource, const char* pszKey,
                                  int nStart, int nNext, int nTerm, int bRemove);

    int                 Open( const char *, int bUpdate, int bTestOpen );
    int                 OpenTable( const char *pszSchemaName, const char *pszTableName, 
                                   const char *pszGeomCol,int nCoordDimension,
                                   int nSRID, OGRwkbGeometryType eType, int bUpdate );

    const char          *GetName() { return pszName; }
    int                 GetLayerCount();
    OGRLayer            *GetLayer( int );

    int                 GetGeometryFormat() { return nGeometryFormat; }

    virtual int         DeleteLayer( int iLayer );
    virtual OGRLayer    *CreateLayer( const char *,
                                      OGRSpatialReference * = NULL,
                                      OGRwkbGeometryType = wkbUnknown,
                                      char ** = NULL );

    int                 TestCapability( const char * );

    virtual OGRLayer *  ExecuteSQL( const char *pszSQLCommand,
                                    OGRGeometry *poSpatialFilter,
                                    const char *pszDialect );
    virtual void        ReleaseResultSet( OGRLayer * poLayer );

    char                *LaunderName( const char *pszSrcName );
    OGRErr              InitializeMetadataTables();

    OGRSpatialReference* FetchSRS( int nId );
    int                 FetchSRSId( OGRSpatialReference * poSRS );

    // Internal use
    CPLODBCSession     *GetSession() { return &oSession; }

};

/************************************************************************/
/*                             OGRODBCDriver                            */
/************************************************************************/

class OGRMSSQLSpatialDriver : public OGRSFDriver
{
  public:
    ~OGRMSSQLSpatialDriver();
                
    const char *GetName();
    OGRDataSource *Open( const char *, int );

    virtual OGRDataSource *CreateDataSource( const char *pszName,
                                             char ** = NULL );
    
    int                 TestCapability( const char * );
};

#endif /* ndef _OGR_MSSQLSPATIAL_H_INCLUDED */
