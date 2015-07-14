#ifndef CalorimeterHitType_h 
#define CalorimeterHitType_h 1

#include <iostream>

/** Helper class for decoding/encoding lcio::CalorimeterHit types for the ILD
 *  detector. The encoding is: caloType + 10 * caloID + 1000 * layout + 10000 * layerNum <br>
 *  (see enums: CaloType, CaloID and Layout for possible values).<br>
 *  Example usage: <br>
 *  <pre>
 *     lcio::CalorimeterHit* cHit = .... ;
 * 
 *     // set the type (e.g. in digitization ) 
 *     cHit->setType( CHT( CHT::em ,  CHT::ecal , CHT::plug , 12 ) ) ;
 *
 *     ...
 *
 *     CHT cht = cHit->getType() ;
 * 
 *     //   sum energies for electromagentic, hadronic and tailcatcher:
 *     if( cht.is( CHT::em ) ) 
 *          e_em +=  cHit->getEnergy() ;
 *     else 
 *       if ( cht.is(CHT::had ) ) 
 *          e_had += cHit->getEnergy() ;
 *       else
 *          e_muon += cHit->getEnergy() ;
 * 
 *     // use only EcalPlug hits:
 *     if( cht.is( CHT::ecal) && cht.is( CHT::plug) ) 
 *      
 *     // get the layer number (e.g. for calibration or clustering) 
 *     unsigned l = cht.layer() ; 
 *     // or directly :
 *     unsigned l = CHT(  cHit->getType() ).layer()  ;
 *
 *     // detailed print:
 *     std::cout <<  CHT(  cHit->getType() ) << std::endl ;
 *
 *  </pre>
 * 
 *  F.Gaede, DESY, 12/2008
 */

class CHT{
  
public:
  
  /** calorimeter types */
  typedef enum {
    em = 0 ,
    had = 1,
    muon = 2 
  } CaloType ;
  
  /** calo ids - specific to ILD */
  typedef enum {
    unknown = 0 ,
    ecal =  1 ,
    hcal =  2 ,
    yoke =  3 ,
    lcal =  4 ,
    lhcal = 5 ,
    bcal =  6 
  } CaloID ;
  
  /** calo layout / subdetector  */ 
  typedef enum {
    any = 0 ,
    barrel = 1 ,
    endcap = 2 ,
    plug = 3 ,
    ring = 4
  } Layout ;
  
  
  /** C'tor for initialization from CalorimeterHit::getType()  */
  CHT(int type) :_type(type) {
  }
  
  /** C'tor  for encoding the calo type inforamtion  */
  CHT(CaloType c, CaloID n, Layout l , unsigned layer) {
    _type = c * fCaloType  + n * fCaloID  + l * fLayout + layer * fLayer ; 
  }
  
  /** calorimeter type: CHT::em , CHT::had, CHT::muon */
  CaloType caloType() const {
    return (CaloType) (_type % fCaloID ) ;
  }

  /** calo ID - see enum CaloID for allowed values */
  CaloID caloID() const {
    return (CaloID) ( (_type % fLayout )  / fCaloID ) ;
  }
  
  /** calo layout - see enum layout for allowed values */
  Layout layout() const {
    return (Layout) ( (_type % fLayer )  / fLayout ) ;
  }
  
  /** calo layer of hit  */
  unsigned layer()  const {
    return unsigned ( _type ) / fLayer ;    
  }

  bool is(CaloType t)  const {
    return caloType() == t ; 
  }
  
  bool is(CaloID n) const  {
    return caloID() == n ; 
  }

  bool is(Layout l) const  {
    return layout() == l ; 
  }
  
  /** automatic conversion to int */
  operator int()  const { return _type ;}

  /** explicit conversion to int */
  int toInt() const  { return _type ;}
  
protected:
  
  int _type ;
  
  static const int fCaloType =     1 ;
  static const int fCaloID   =    10 ;
  static const int fLayout   =  1000 ;
  static const int fLayer    = 10000 ;
  
} ;

/** detailed string for calo type */
std::ostream& operator<<(std::ostream& os, const CHT& cht) ; 


/** Return Layout based on the collection name, e.g. if name contains tolower("endcap") CHT::endcap is returned. In case no known layout
    is found, CHT::any is returned.*/
CHT::Layout layoutFromString(const std::string& name) ;

/** Return caloID based on the collection name, e.g. if name contains tolower("HCal") CHT::hcal is returned. In case no known type
    is found, CHT::unknown is returned.*/
CHT::CaloID caloIDFromString(const std::string& name) ;

/** Return caloType from string, e.g. if name contains tolower("Had") CHT::had is returned. In case no known type
    is found, CHT::em is returned.*/
CHT::CaloType caloTypeFromString(const std::string& name) ;



#endif


// int main(){
  
//   CHT ct( CHT::em ,  CHT::ecal , CHT::plug , 12 )  ;

//   std::cout <<  ct << std::endl ;

//   std::cout << " caloType : " << ct.toInt() << std::endl ;
  
//   std::cout << " caloID " <<    ct.caloID() 
// 	    << " is bcal " << ct.is( CHT::bcal )  
// 	    << " is muon " << ct.is( CHT::muon )  
// 	    << " layer " << ct.layer()   
// 	    << std::endl ;
  
//   int i = CHT( CHT::muon ,  CHT::ecal , CHT::endcap , 123 ) ;
  
//   std::cout << " type i : " << i << std::endl ;

  
//   std::cout <<  " layer from int : " << CHT(i).layer() << std::endl ;

//   CHT it(i) ;

//   std::cout << " caloID " <<    it.caloID() 
// 	    << " is bcal " << it.is( CHT::bcal )  
// 	    << " is muon " << it.is( CHT::muon )  
// 	    << " is endcap " << it.is( CHT::endcap )  
// 	    << " layer " << it.layer()   
// 	    << " layout " << it.layout()   
// 	    << std::endl ;


//   std::cout << it << std::endl ;

//   CHT ct2 = 12345678 ;

//   std::cout << ct2 << std::endl ;

//   std::cout << " sizeof( CHT )  : "  << sizeof( it ) << std::endl ;

// }
