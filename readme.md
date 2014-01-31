## Rigorous Vertical Datum Support in OGR Spatial Reference ##
GDAL/OGR is widely used in open source and closed source applications. The OGRSpatialReference in conjunction with the PROJ.4 library is quasi standard for 2D coordinate transformation. Although, 3D transformations are - to some extend - supported within the library, the implemented solution is unsatisfactory from a geodetic point of view. The objective is to establish a rigorous 3D transformation chain, to support Vertical Datum definitions in a generic way.

### What's in this package ###
 * [PROJ](http://trac.osgeo.org/proj/) version 4.8.0
 * [GDAL/OGR](http://www.gdal.org/) version 1.10.0
 * [SpatialRef3D](https://github.com/ottointhesky/OGRSpatialRef3D)
 * Test Utilities : coordinate transformation test (main.cpp), transformation performance measurements (perfmain.cpp), transformation correctness validation (validate.cpp)
 
---
 
What is SpatialRef3D ?
====================
SpatialRef3D is a library that extends some classes in OGR for transforming one spatial reference to another spatial reference. 
The extension implemented in this library support for additional information about height in the coordinate (spatial reference) transformation.
With this extension, the third component in coordinate vector may represent not only ellipsoidal height, but also orthometric height and (possibly inhomogeneous) local height.

### Height Model ###
Height model are additional file(s) containing vertical information in raster format and scalar value. SpatialRef3D support 2 (two) types of information which are Geoid undulation file and vertical correction model which used to correcting inhomogeneous height anomalies. further explanation of the model used can be read in [`doc`]() directory.

### Adding Height Model into the Spatial Reference System ###
Adding height model to the extended Spatial Reference System can be done in two ways : 

 1. using custom WKT format supported by the extension
 
    *Please be aware that the custom format specification is subject to change. Therefore user should be cautious when using this format in conjunction with other library*
   
    The height model are implemented as custom node inside the `GEOGCS` node: 
        
     * a `GEOID` node, pass a string filename denoting path to GDAL compatible raster
        
            GEOID["GeoidFileIdentifier",["BEV\GEOID_2008_GRS80.bil"]],
        
     * a `VCORR` node, pass a string filename denoting path to GDAL compatible raster
        
            VCORR["VertCorrectionIdentifier",["BEV\GV_Hoehengrid_V1.bil"]],
        
     * a `VSHIFT` node, pass a floating-point number inside the node
        
            VSHIFT[-156.68],
        
    When passing the information in custom WKT format, please use `importFromWkt3D` instead of standard `importFromWkt`. `importFromWkt3D` works like `importFromWkt` with additional processing implemented to interpret custom nodes described above.
 
 
 2. programmatically while instantiating `OGRSpatialReference3D` class
     
     To use height model programmatically, use these methods : 
      
      * `SetGeoidModel( const char * pszGeoidModel )`
      * `SetVCorrModel( const char * pszGeoidModel )`
      * `SetVOffset( double dVOffset )`
      
How to build Documentation
====================
* Download and install [Doxygen](http://www.doxygen.org)
* from the command line in the project directory (`OGRSpatialRef3D`) run `doxygen doc\Doxyfile`