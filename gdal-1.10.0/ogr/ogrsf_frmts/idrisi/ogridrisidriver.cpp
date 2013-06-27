/******************************************************************************
 * $I$
 *
 * Project:  Idrisi Translator
 * Purpose:  Implements OGRIdrisiDriver.
 * Author:   Even Rouault, even dot rouault at mines dash paris dot org
 *
 ******************************************************************************
 * Copyright (c) 2011, Even Rouault <even dot rouault at mines dash paris dot org>
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

#include "ogr_idrisi.h"
#include "cpl_conv.h"

CPL_CVSID("$Id: ogridrisidriver.cpp 23040 2011-09-04 14:42:12Z rouault $");

// g++ ogr/ogrsf_frmts/idrisi/*.cpp -Wall -g -fPIC -shared -o ogr_Idrisi.so -Iport -Igcore -Iogr -Iogr/ogrsf_frmts/idrisi -Iogr/ogrsf_frmts -Ifrmts/idrisi

extern "C" void RegisterOGRIdrisi();

/************************************************************************/
/*                       ~OGRIdrisiDriver()                         */
/************************************************************************/

OGRIdrisiDriver::~OGRIdrisiDriver()

{
}

/************************************************************************/
/*                              GetName()                               */
/************************************************************************/

const char *OGRIdrisiDriver::GetName()

{
    return "Idrisi";
}

/************************************************************************/
/*                                Open()                                */
/************************************************************************/

OGRDataSource *OGRIdrisiDriver::Open( const char * pszFilename, int bUpdate )

{
    OGRIdrisiDataSource   *poDS = new OGRIdrisiDataSource();

    if( !poDS->Open( pszFilename, bUpdate ) )
    {
        delete poDS;
        poDS = NULL;
    }

    return poDS;
}

/************************************************************************/
/*                           TestCapability()                           */
/************************************************************************/

int OGRIdrisiDriver::TestCapability( const char * pszCap )

{
    return FALSE;
}

/************************************************************************/
/*                        RegisterOGRIdrisi()                       */
/************************************************************************/

void RegisterOGRIdrisi()

{
    OGRSFDriverRegistrar::GetRegistrar()->RegisterDriver( new OGRIdrisiDriver );
}

