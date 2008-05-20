#include "CutOnCosThetaQQbar.h"

#include <cmath>
//SJA:FIXED:added to make gcc4.3 compliant
#include <cstdlib>


using namespace lcio ;
using namespace marlin ;

CutOnCosThetaQQbar aCutOnCosThetaQQbar ;



CutOnCosThetaQQbar::CutOnCosThetaQQbar() : Processor("CutOnCosThetaQQbar") {

  _description = "The output condition of this processor is false as long as the |cos(theta)| > cutCosTheta and true otherwise.";


  registerProcessorParameter( "cutCosTheta",
			      "cut on |cos(theta)|",
			      _cutCosTheta,
			      (double)1.0);

}


void CutOnCosThetaQQbar::init() {
  
  // usually a good idea to 
  // printParameters();
  
  _nRun = 0 ;
  _nEvt = 0 ;


}


void CutOnCosThetaQQbar::processRunHeader( LCRunHeader* run) {

  ++_nRun;

}


void CutOnCosThetaQQbar::processEvent( LCEvent * evt ) {

  double cosThOfQuarkSystem = getCosThOfQuarkSystem(evt);

  // debug
  std::cout << "cos(theta) of the qqbar system = " << cosThOfQuarkSystem << std::endl;
  
  if ( fabs(cosThOfQuarkSystem) > _cutCosTheta ) {
    
    // debug
    std::cout << "|cos(theta)| > " << _cutCosTheta << "  " << " =>  event discarded, proceed with next event" << std::endl << std::endl;

    setReturnValue(false);

  }
  else setReturnValue(true);

  ++_nEvt;

}


void CutOnCosThetaQQbar::check( LCEvent * evt ) {
 
}


void CutOnCosThetaQQbar::end() {

}


double CutOnCosThetaQQbar::getCosThOfQuarkSystem(const LCEvent* evt) {
  
  try {
    
    std::vector< std::string >::const_iterator iter;
    const std::vector< std::string >* ColNames = evt->getCollectionNames();
    
    for( iter = ColNames->begin() ; iter != ColNames->end() ; iter++) {

      LCCollection* col = evt->getCollection(*iter);

      if ( (col->getTypeName() == LCIO::MCPARTICLE) && (*iter == "MCParticle") ) {
	
	int NMCParticles = col->getNumberOfElements();
	
	for(int j=0; j<NMCParticles; ++j){
	  
	  MCParticle* mcP = dynamic_cast<MCParticle*>(col->getElementAt(j));

	  int pdgMCP = mcP->getPDG();


	  // debug
	  // std::cout << pdgMCP << std::endl;
	  

	  if (pdgMCP==23) {

	    MCParticle* quark1 = dynamic_cast<MCParticle*>(col->getElementAt(j+1));
	    MCParticle* quark2 = dynamic_cast<MCParticle*>(col->getElementAt(j+2));
	    int absPDGQuark1 = abs(quark1->getPDG());
	    int absPDGQuark2 = abs(quark2->getPDG());

	    if ( ( (absPDGQuark1==1) && (absPDGQuark2==1) ) || ( (absPDGQuark1==2) && (absPDGQuark2==2) ) || ( (absPDGQuark1==3) && (absPDGQuark2==3) ) ) {
	      
	      // refer to quark1
	      
  	      const double* p = quark1->getMomentum();
	      const double absP = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
	      double cosTh = p[2]/absP;

	      // debug
	      // std::cout << "p = " << "( " << p[0] << "," << p[1] << "," << p[2] << " )" << "   " << "cosTh = " << cosTh << std::endl;

	      return cosTh;

	    }	    
	    
	  }
	  
	}
	
      }
      
    }

  }
  catch(DataNotAvailableException &e){std::cout << "no valid MC collection in event " << _nEvt << std::endl; };
   
  return 2.0;

}
