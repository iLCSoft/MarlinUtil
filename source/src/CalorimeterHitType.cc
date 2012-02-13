#include "CalorimeterHitType.h"

#include <algorithm>


/** detailed string for calo type */
std::ostream& operator<<(std::ostream& os, const CHT& cht){
  os << " calo hit type: " ;
  
  switch ( cht.caloType() ) {
  case CHT::em :    os << " em,   " ;    break ;
  case CHT::had :   os << " had,  " ;    break ;
  case CHT::muon :  os << " muon, " ;    break ;
  default:  os << "  -  ," ;
  }
  switch ( cht.caloID() ){
  case CHT::ecal :    os << "ecal,  " ;    break ;
  case CHT::hcal :    os << "hcal,  " ;    break ;
  case CHT::yoke :    os << "yoke,  " ;    break ;
  case CHT::lcal :    os << "lcal,  " ;    break ;
  case CHT::lhcal :   os << "lhcal, " ;    break ;
  case CHT::bcal :    os << "bcal,  " ;    break ;
  default:  os << "  -  ," ;
  }
  switch ( cht.layout() ) {
  case CHT::any    :  os << "any,    " ;    break ;
  case CHT::ring   :  os << "ring,   " ;    break ;
  case CHT::endcap :  os << "endcap, " ;    break ;
  case CHT::barrel :  os << "barrel, " ;    break ;
  case CHT::plug   :  os << "plug,   " ;    break ;
  default:  os << "  -  ," ;
  }
  os << " layer: " << cht.layer() ;
  
  return os ;
}

/** Helper functions that should go to Marlinutil/CalorimeterHitTypes.hh */

CHT::Layout layoutFromString(const std::string& name){

  std::string str( name ) ;
  std::transform( str.begin() , str.end() , str.begin(), ::tolower ) ;

  if( str.find("ring" )   != std::string::npos )  return CHT::ring ;
  if( str.find("plug" )   != std::string::npos )  return CHT::plug ;
  if( str.find("endcap" ) != std::string::npos )  return CHT::endcap ;
  if( str.find("barrel" ) != std::string::npos )  return CHT::barrel ;

  std::cout << " not found :" << str << " in " << name << std::endl;
  return CHT::any ;
}

CHT::CaloID caloIDFromString(const std::string& name){

  std::string str( name ) ;
  std::transform( str.begin() , str.end() , str.begin(), ::tolower ) ;

  if( str.find("ecal" ) != std::string::npos )  return CHT::ecal ;
  if( str.find("hcal" ) != std::string::npos )  return CHT::hcal ;
  if( str.find("yoke" ) != std::string::npos )  return CHT::yoke ;
  if( str.find("lcal" ) != std::string::npos )  return CHT::lcal ;
  if( str.find("lhcal") != std::string::npos )  return CHT::lhcal ;
  if( str.find("bcal" ) != std::string::npos )  return CHT::bcal ;

  return CHT::unknown ;
}

