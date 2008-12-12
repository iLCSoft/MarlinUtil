#include "CalorimeterHitType.h"


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
  case CHT::endcap :  os << "endcap, " ;    break ;
  case CHT::barrel :  os << "barrel, " ;    break ;
  case CHT::plug   :  os << "plug,   " ;    break ;
  case CHT::ring   :  os << "ring,   " ;    break ;
  default:  os << "  -  ," ;
  }
  os << " layer: " << cht.layer() ;
  
  return os ;
}
