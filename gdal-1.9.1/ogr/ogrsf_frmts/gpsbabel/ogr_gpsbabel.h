/******************************************************************************
 * $Id: ogr_gpsbabel.h 20996 2010-10-28 18:38:15Z rouault $
 *
 * Project:  OpenGIS Simple Features Reference Implementation
 * Purpose:  Private definitions for OGR/GPSBabel driver.
 * Author:   Even Rouault, <even dot rouault at mines dash paris dot org>
 *
 ******************************************************************************
 * Copyright (c) 2010, Even Rouault
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

#ifndef _OGR_GPSBABEL_H_INCLUDED
#define _OGR_GPSBABEL_H_INCLUDED

#include "ogrsf_frmts.h"
#include "cpl_string.h"

int ForkAndPipe(const char * const argv[], VSILFILE* fin, VSILFILE* fout);

/************************************************************************/
/*                        OGRGPSBabelDataSource                         */
/************************************************************************/

class OGRGPSBabelDataSource : public OGRDataSource
{
    int                 nLayers;
    OGRLayer*           apoLayers[5];
    char               *pszName;
    char               *pszGPSBabelDriverName;
    char               *pszFilename;
    CPLString           osTmpFileName;
    OGRDataSource      *poGPXDS;

  public:
                        OGRGPSBabelDataSource();
                        ~OGRGPSBabelDataSource();

    virtual const char  *GetName() { return pszName; }
    virtual int         GetLayerCount() { return nLayers; }
    virtual OGRLayer   *GetLayer( int );

    virtual int         TestCapability( const char * );

    int                 Open ( const char* pszFilename, int bUpdateIn );

    static int          IsSpecialFile(const char* pszFilename);
    static int          IsValidDriverName(const char* pszGPSBabelDriverName);
};


/************************************************************************/
/*                   OGRGPSBabelWriteDataSource                         */
/************************************************************************/

class OGRGPSBabelWriteDataSource : public OGRDataSource
{
    char               *pszName;
    char               *pszGPSBabelDriverName;
    char               *pszFilename;
    CPLString           osTmpFileName;
    OGRDataSource      *poGPXDS;

    int                 Convert();

  public:
                        OGRGPSBabelWriteDataSource();
                        ~OGRGPSBabelWriteDataSource();

    virtual const char  *GetName() { return pszName; }
    virtual int         GetLayerCount();
    virtual OGRLayer   *GetLayer( int );

    virtual int         TestCapability( const char * );

    virtual OGRLayer   *CreateLayer( const char * pszLayerName,
                                     OGRSpatialReference *poSRS,
                                     OGRwkbGeometryType eType,
                                     char ** papszOptions );

    int                 Create ( const char* pszFilename, char **papszOptions );
};

/************************************************************************/
/*                        OGRGPSBabelDriver                             */
/************************************************************************/

class OGRGPSBabelDriver : public OGRSFDriver
{
  public:
                ~OGRGPSBabelDriver();

    virtual const char    *GetName();
    virtual OGRDataSource *Open( const char *, int );
    virtual OGRDataSource *CreateDataSource( const char * pszName,
                                             char **papszOptions );
    virtual OGRErr         DeleteDataSource( const char *pszFilename );

    virtual int            TestCapability( const char * );
};

#endif /* ndef _OGR_GPSBABEL_H_INCLUDED */

