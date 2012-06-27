/******************************************************************************
 *
 * Project:  OpenGIS Simple Features Reference Implementation
 * Purpose:  Classes for manipulating spatial reference systems with 
 *           vetical datum suppurt in a platform non-specific manner.
 * Authors:  Gottfried Mandlburger, Johannes Otepka, Bhargav patel
 *
 ******************************************************************************
 * Copyright (c) 2012,  I.P.F., TU Vieanna.
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

#include "ogr_spatialref3D.h"

/************************************************************************/
/*                         OGRSpatialReference3D                        */
/************************************************************************/

/**
 * This class respresents a OpenGIS 3D Spatial Reference System. 
 * It is derived from OGRSpatialReference class and adds elements to 
 * define the Reference System's Vertical Datum. The class provides:
 *
 *  I  ) a (quasi) geoid raster model for transformation between 
 *       ellipsoidal and (quasi) orthometric heights
 *  II ) an additional height correction raster model to compensate the
 *       errors of historic heigth system realizations (e.g. MGI, DHHN...)
 *  III) a constant offset to deal with local height origins, like tidal 
 *       height systems, local height definitions (e.g. Wiener Null) ...
 *  IV ) a scaling factor (e.g. foot-meter-conversion)
 * 
/************************************************************************/


OGRSpatialReference3D::OGRSpatialReference3D()
{
}

OGRSpatialReference3D::OGRSpatialReference3D(const OGRSpatialReference&)
{
}

OGRSpatialReference3D::OGRSpatialReference3D(const char * pszWKT,
                          const char * pszGeoidModel,
                          const char * pszVCorrModel,
                          double dfVOffset,
                          double dfVScale)
{
}

OGRCoordinateTransformation3D::OGRCoordinateTransformation3D()
{
}