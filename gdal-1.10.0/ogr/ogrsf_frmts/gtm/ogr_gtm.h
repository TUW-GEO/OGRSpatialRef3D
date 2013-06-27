/******************************************************************************
 * $Id: ogr_gtm.h 20996 2010-10-28 18:38:15Z rouault $
 *
 * Project:  GTM Driver
 * Purpose:  Declarations for OGR wrapper classes for GTM, and OGR->GTM
 *           translation of geometry.
 * Author:   Leonardo de Paula Rosa Piga; http://lampiao.lsc.ic.unicamp.br/~piga
 *
 ******************************************************************************
 * Copyright (c) 2009, Leonardo de Paula Rosa Piga
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
#ifndef OGR_GTM_H_INCLUDED
#define OGR_GTM_H_INCLUDED

#include "ogrsf_frmts.h"
#include "cpl_conv.h"
#include "cpl_string.h"

class OGRGTMDataSource;

typedef enum
{
    GTM_NONE,
    GTM_WPT,
    GTM_TRACK
} GTMGeometryType;

#ifndef FILE_OFFSETS
#define FILE_OFFSETS
#define NWPTSTYLES_OFFSET 27
#define NWPTSTYLES_SIZE 4

#define NWPTS_OFFSET 35
#define NWPTS_SIZE 4

#define NTRCKS_OFFSET 39
#define NTRCKS_SIZE 4

#define NMAPS_OFFSET 63
#define NMAPS_SIZE 4

#define NTK_OFFSET 67
#define NTK_SIZE 4

#define BOUNDS_OFFSET 47

#define DATUM_SIZE 58

/* GTM_EPOCH is defined as the unix time for the 31 dec 1989 00:00:00 */
#define GTM_EPOCH 631065600

#endif

#ifndef MAX
#  define MIN(a,b)      ((a<b) ? a : b)
#  define MAX(a,b)      ((a>b) ? a : b)
#endif


#include "gtm.h"


/************************************************************************/
/*                           OGRGTMLayer                                */
/************************************************************************/
class OGRGTMLayer : public OGRLayer
{
public:
    OGRGTMLayer();
    ~OGRGTMLayer();
    //
    // OGRLayer Interface
    //
    OGRFeatureDefn* GetLayerDefn();
    OGRSpatialReference* GetSpatialRef() {return poSRS;}
    virtual void ResetReading() = 0;
    virtual OGRFeature* GetNextFeature() = 0;
    virtual int GetFeatureCount(int bForce = TRUE) = 0;
    virtual OGRErr CreateFeature (OGRFeature *poFeature) = 0;

    int TestCapability( const char* pszCap );
    
    OGRErr CreateField( OGRFieldDefn *poField, int bApproxOK );

protected:
    OGRGTMDataSource* poDS;
    OGRSpatialReference* poSRS;
    OGRCoordinateTransformation *poCT;
    char* pszName;

    OGRFeatureDefn* poFeatureDefn;
    int nNextFID;
    int nTotalFCount;
    
    int bError;
    
    static OGRErr CheckAndFixCoordinatesValidity( double& pdfLatitude, double& pdfLongitude );

};


/************************************************************************/
/*                           GTMWaypointLayer                           */
/************************************************************************/
class GTMWaypointLayer : public OGRGTMLayer
{
public:
    GTMWaypointLayer( const char* pszName,
                      OGRSpatialReference* poSRSIn,
                      int bWriterIn,
                      OGRGTMDataSource* poDSIn );
    ~GTMWaypointLayer();
    OGRErr CreateFeature (OGRFeature *poFeature);
    void ResetReading();
    OGRFeature* GetNextFeature();
    int GetFeatureCount(int bForce = TRUE);

    enum WaypointFields{NAME, COMMENT, ICON, DATE};
private:
    void WriteFeatureAttributes( OGRFeature *poFeature, float altitude );
};

/************************************************************************/
/*                           GTMTrackLayer                              */
/************************************************************************/
class GTMTrackLayer : public OGRGTMLayer
{
public:
    GTMTrackLayer( const char* pszName,
                   OGRSpatialReference* poSRSIn,
                   int bWriterIn,
                   OGRGTMDataSource* poDSIn );
    ~GTMTrackLayer();
    OGRErr CreateFeature (OGRFeature *poFeature);
    void ResetReading();
    OGRFeature* GetNextFeature();
    int GetFeatureCount(int bForce = TRUE);
    enum TrackFields{NAME, TYPE, COLOR};

private:
    void WriteFeatureAttributes( OGRFeature *poFeature );
    void WriteTrackpoint( double lat, double lon, float altitude, bool start );

};


/************************************************************************/
/*                           OGRGTMDataSource                           */
/************************************************************************/
class OGRGTMDataSource : public OGRDataSource
{
public:

    // OGRDataSource Interface
    OGRGTMDataSource();
    ~OGRGTMDataSource();

    int Open( const char *pszFilename, int bUpdate );
    int Create( const char *pszFilename, char **papszOptions );

    const char* GetName() { return pszName; }
    int GetLayerCount() { return nLayers; }

    OGRLayer* GetLayer( int );

    OGRLayer* CreateLayer (const char *pszName, 
                           OGRSpatialReference *poSpatialRef=NULL, 
                           OGRwkbGeometryType eGType=wkbUnknown, 
                           char **papszOptions=NULL);
    int TestCapability( const char * );

    // OGRGTMDataSource Methods
    VSILFILE* getOutputFP() { return fpOutput; }
    VSILFILE* getTmpTrackpointsFP() { return fpTmpTrackpoints; }
    VSILFILE* getTmpTracksFP() { return fpTmpTracks; }
    bool isFirstCTError() { return !bIssuedCTError; }
    void issuedFirstCTError() { bIssuedCTError = true; }

    /* Functions to handle with waypoints */
    int getNWpts();
    bool hasNextWaypoint();
    Waypoint* fetchNextWaypoint();
    void rewindWaypoint();
  
    /* Functions to handle with tracks */
    int getNTracks();
    bool hasNextTrack();
    Track* fetchNextTrack();
    void rewindTrack();


    /* Functions for writing ne files */
    float getMinLat() { return minlat; }
    float getMaxLat() { return maxlat; }
    float getMinLon() { return minlon; }
    float getMaxLon() { return maxlon; }

    void checkBounds(float newLat,
                     float newLon);
    int getNumWaypoints() { return numWaypoints; } 
    int getNumTrackpoints() { return numTrackpoints; }
    int getTracks() { return numTracks; };

    int incNumWaypoints() { return ++numWaypoints; } 
    int incNumTrackpoints() { return ++numTrackpoints; }
    int incNumTracks() { return ++numTracks; };
private:
    VSILFILE* fpOutput;
  
    /* GTM is not a contiguous file. We need two temporary files because
       trackpoints and tracks are stored separated and we don't know in
       advance how many trackpoints and tracks the new file will
       have. So, we create temporary file and append the at the end of
       the gtm file whe everything is done, that is, in the
       destructor. */
    VSILFILE* fpTmpTrackpoints;
    char* pszTmpTrackpoints;

    VSILFILE* fpTmpTracks;
    char* pszTmpTracks;

    GTM* poGTMFile;
    char* pszName;
  
    OGRGTMLayer **papoLayers;
    int nLayers;
  
    bool bIssuedCTError;

    /* Used for creating a new file */
    float minlat;
    float maxlat;
    float minlon;
    float maxlon;

    int numWaypoints;
    int numTracks;
    int numTrackpoints;

    void AppendTemporaryFiles();
    void WriteWaypointStyles();
};

/************************************************************************/
/*                             OGRGTMDriver                             */
/************************************************************************/

class OGRGTMDriver : public OGRSFDriver
{
public:
    ~OGRGTMDriver();

    //
    // OGRSFDriver Interface
    //
    const char* GetName();
    OGRDataSource* Open( const char * pszName_, int bUpdate );
    OGRDataSource* CreateDataSource( const char *pszName_, char** papszOptions );

    int TestCapability( const char* pszCap );
};

#endif //OGR_GTM_H_INCLUDED
