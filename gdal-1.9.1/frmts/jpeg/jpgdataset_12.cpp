/******************************************************************************
 * $Id: jpgdataset_12.cpp 24159 2012-03-24 14:51:28Z rouault $
 *
 * Project:  JPEG JFIF Driver
 * Purpose:  Implement GDAL JPEG Support based on IJG libjpeg.
 * Author:   Even Rouault, even dot rouault at mines dash paris dot org
 *
 ******************************************************************************
 * Copyright (c) 2009, Even Rouault
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

#if defined(JPEG_DUAL_MODE_8_12)
#define LIBJPEG_12_PATH   "libjpeg12/jpeglib.h" 
#define JPGDataset        JPGDataset12
#define JPGRasterBand     JPGRasterBand12
#define JPGMaskBand       JPGMaskBand12
#include "jpgdataset.cpp"

GDALDataset* JPEGDataset12Open(GDALOpenInfo* poOpenInfo)
{
    return JPGDataset12::Open(poOpenInfo);
}

GDALDataset* JPEGDataset12CreateCopy( const char * pszFilename,
                                    GDALDataset *poSrcDS,
                                    int bStrict, char ** papszOptions,
                                    GDALProgressFunc pfnProgress,
                                    void * pProgressData )
{
    return JPGDataset12::CreateCopy(pszFilename, poSrcDS,
                                    bStrict, papszOptions,
                                    pfnProgress, pProgressData);
}

#endif /* defined(JPEG_DUAL_MODE_8_12) */
