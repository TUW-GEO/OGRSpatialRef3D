/******************************************************************************
 * $Id: ogr_ili2.h 13906 2008-03-01 13:08:28Z rouault $
 *
 * Project:  Interlis 2 Translator
 * Purpose:   Definition of classes for OGR Interlis 2 driver.
 * Author:   Markus Schnider, Sourcepole AG
 *
 ******************************************************************************
 * Copyright (c) 2004, Pirmin Kalberer, Sourcepole AG
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

#ifndef _OGR_ILI2_H_INCLUDED
#define _OGR_ILI2_H_INCLUDED

#include "ogrsf_frmts.h"
#include "ili2reader.h"
#include "iom/iom.h"

#include <string>
#include <list>

class OGRILI2DataSource;

/************************************************************************/
/*                           OGRILI2Layer                               */
/************************************************************************/

class OGRILI2Layer : public OGRLayer
{
private:
    OGRSpatialReference *poSRS;
    OGRFeatureDefn     *poFeatureDefn;
    std::list<OGRFeature *>    listFeature;
    std::list<OGRFeature *>::const_iterator listFeatureIt;

    int                 bWriter;

    OGRILI2DataSource   *poDS;

  public:
                        OGRILI2Layer( const char * pszName, 
                                     OGRSpatialReference *poSRS, 
                                     int bWriter,
                                     OGRwkbGeometryType eType,
                                     OGRILI2DataSource *poDS );

                       ~OGRILI2Layer();

    OGRErr              SetFeature(OGRFeature *poFeature);
    
    void                ResetReading();
    OGRFeature *        GetNextFeature();

    int                 GetFeatureCount( int bForce = TRUE );

    OGRErr              CreateFeature( OGRFeature *poFeature );
    
    OGRFeatureDefn *    GetLayerDefn() { return poFeatureDefn; }

    OGRErr              CreateField( OGRFieldDefn *poField, int bApproxOK = TRUE );

    OGRSpatialReference *GetSpatialRef();
    
    int                 TestCapability( const char * );
};

/************************************************************************/
/*                          OGRILI2DataSource                           */
/************************************************************************/

class OGRILI2DataSource : public OGRDataSource
{
  private:
    std::list<OGRLayer *> listLayer;
    
    char        *pszName;
    IILI2Reader *poReader;
    IOM_FILE    fpTransfer;  //for writing
    IOM_BASKET  basket;

    int         nLayers;
    OGRILI2Layer** papoLayers;

  public:
                OGRILI2DataSource();
               ~OGRILI2DataSource();

    int         Open( const char *, int bTestOpen );
    int         Create( const char *pszFile, char **papszOptions );

    const char *GetName() { return pszName; }
    IOM_BASKET  GetBasket() { return basket; }
    int         GetLayerCount() { return listLayer.size(); }
    OGRLayer   *GetLayer( int );

    virtual OGRLayer *CreateLayer( const char *, 
                                      OGRSpatialReference * = NULL,
                                      OGRwkbGeometryType = wkbUnknown,
                                      char ** = NULL );

    int         TestCapability( const char * );
};

/************************************************************************/
/*                            OGRILI2Driver                             */
/************************************************************************/

class OGRILI2Driver : public OGRSFDriver
{
  public:
                ~OGRILI2Driver();

    const char *GetName();
    OGRDataSource *Open( const char *, int );

    virtual OGRDataSource *CreateDataSource( const char *pszName,
                                             char ** = NULL );

    int                 TestCapability( const char * );
};

#endif /* _OGR_ILI2_H_INCLUDED */
