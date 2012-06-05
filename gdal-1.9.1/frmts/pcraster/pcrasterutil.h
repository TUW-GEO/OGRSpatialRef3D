/******************************************************************************
 * $Id: pcrasterutil.h 15451 2008-10-03 10:59:42Z kdejong $
 *
 * Project:  PCRaster Integration
 * Purpose:  PCRaster driver support declarations.
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
#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

// PCRaster library headers.
#ifndef INCLUDED_CSF
#include "csf.h"
#define INCLUDED_CSF
#endif

#ifndef INCLUDED_PCRTYPES
#include "pcrtypes.h"
#define INCLUDED_PCRTYPES
#endif

// Module headers.
#ifndef INCLUDED_GDAL_PRIV
#include "gdal_priv.h"
#define INCLUDED_GDAL_PRIV
#endif


GDALDataType       cellRepresentation2GDALType(CSF_CR cellRepresentation);

CSF_VS             string2ValueScale   (const std::string& string);

std::string        valueScale2String   (CSF_VS valueScale);

std::string        cellRepresentation2String(CSF_CR cellRepresentation);

CSF_VS             GDALType2ValueScale (GDALDataType type);

/*
CSF_CR             string2PCRasterCellRepresentation(
                                        const std::string& string);
                                        */

CSF_CR             GDALType2CellRepresentation(
                                        GDALDataType type,
                                        bool exact);

void*              createBuffer        (size_t size,
                                        CSF_CR type);

void               deleteBuffer        (void* buffer,
                                        CSF_CR type);

bool               isContinuous        (CSF_VS valueScale);

double             missingValue        (CSF_CR type);

void               alterFromStdMV      (void* buffer,
                                        size_t size,
                                        CSF_CR cellRepresentation,
                                        double missingValue);

void               alterToStdMV        (void* buffer,
                                        size_t size,
                                        CSF_CR cellRepresentation,
                                        double missingValue);

MAP*               mapOpen             (std::string const& filename,
                                        MOPEN_PERM mode);

CSF_VS             fitValueScale       (CSF_VS valueScale,
                                        CSF_CR cellRepresentation);

void               castValuesToBooleanRange(
                                        void* buffer,
                                        size_t size,
                                        CSF_CR cellRepresentation);

template<typename T>
struct CastToBooleanRange
{
  void operator()(T& value) {
    if(!pcr::isMV(value)) {
      if(value != 0) {
        value = T(value > T(0));
      }
      else {
        pcr::setMV(value);
      }
    }
  }
};



template<>
struct CastToBooleanRange<UINT1>
{
  void operator()(UINT1& value) {
    if(!pcr::isMV(value)) {
      value = UINT1(value > UINT1(0));
    }
  }
};



template<>
struct CastToBooleanRange<UINT2>
{
  void operator()(UINT2& value) {
    if(!pcr::isMV(value)) {
      value = UINT2(value > UINT2(0));
    }
  }
};



template<>
struct CastToBooleanRange<UINT4>
{
  void operator()(UINT4& value) {
    if(!pcr::isMV(value)) {
      value = UINT4(value > UINT4(0));
    }
  }
};
