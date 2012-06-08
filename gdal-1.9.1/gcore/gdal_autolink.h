#ifndef GDAL_AUTO_LINK_H_INCLUDED
#define GDAL_AUTO_LINK_H_INCLUDED

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

// Warning:
// GDAL 1.7.3 is generally build with nmake, in which case the lib 
// is called "gdal_i.lib", regardless of the configuration (debug/release), 
// and obj-files might not be found  by the debugger. 

#ifndef GDAL_EXPORTS

// Enables autolink of GDAL library
  #pragma comment(lib, "gdal_i.lib")
  #ifdef IPF_LIB_DIAGNOSTIC
  #     pragma message ("Linking to lib file: gdal_i.lib")  
  #endif                                                              

#endif

#endif //GDAL_AUTO_LINK_H_INCLUDED
