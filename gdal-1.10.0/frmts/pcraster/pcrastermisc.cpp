/******************************************************************************
 * $Id: pcrastermisc.cpp 13149 2007-11-29 15:08:00Z warmerdam $
 *
 * Project:  PCRaster Integration
 * Purpose:  PCRaster driver support functions.
 * Author:   Kor de Jong, k.dejong at geog.uu.nl
 *
 ******************************************************************************
 * Copyright (c) 2004, Kor de Jong
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

// Library headers.
#ifndef INCLUDED_IOSTREAM
#include <iostream>
#define INCLUDED_IOSTREAM
#endif

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

// PCRaster library headers.

// Module headers.
#include "gdal_pam.h"

#ifndef INCLUDED_PCRASTERDATASET
#include "pcrasterdataset.h"
#define INCLUDED_PCRASTERDATASET
#endif



CPL_C_START
void GDALRegister_PCRaster(void);
CPL_C_END



void GDALRegister_PCRaster()
{
    if (! GDAL_CHECK_VERSION("PCRaster driver"))
        return;

if(!GDALGetDriverByName("PCRaster")) {

    GDALDriver* driver = new GDALDriver();

    driver->SetDescription("PCRaster");
    driver->SetMetadataItem(GDAL_DMD_LONGNAME, "PCRaster Raster File");
    driver->SetMetadataItem(GDAL_DMD_CREATIONDATATYPES, "Byte Int32 Float32");
    driver->SetMetadataItem(GDAL_DMD_HELPTOPIC, "frmt_various.html#PCRaster");
    driver->SetMetadataItem( GDAL_DMD_EXTENSION, "map" );

    driver->pfnOpen = PCRasterDataset::open;
    driver->pfnCreateCopy = PCRasterDataset::createCopy;

    GetGDALDriverManager()->RegisterDriver(driver);
}
}
