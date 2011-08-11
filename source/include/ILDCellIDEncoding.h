#ifndef INCLUDE_ILDCellIDEncoding
#define INCLUDE_ILDCellIDEncoding 1

#include <string>

namespace ILDCellIDEncoding{

  namespace Fields{
    static const std::string subdet   = "subdet" ;
    static const std::string subdet_nbits = "5" ;
    static const std::string side     = "side" ;
    static const std::string side_nbits   = "1" ;
    static const std::string layer    = "layer" ;
    static const std::string layer_nbits  = "10" ;
    static const std::string module   = "module" ;
    static const std::string module_nbits = "8" ;
  }

  static const std::string encoder_string = 
    Fields::subdet + ":" + Fields::subdet_nbits 
    + "," + Fields::side   + ":" + Fields::side_nbits 
    + "," + Fields::layer  + ":" + Fields::layer_nbits  
    + "," + Fields::module + ":" + Fields::module_nbits ; 

  
  namespace DetID{

    static const int VXD =  1 ;
    static const int SIT =  2 ;
    static const int TPC =  3 ;
    static const int SET =  4 ;
    static const int ETD =  5 ;
    static const int FTD =  6 ; 
  }

}

#endif
