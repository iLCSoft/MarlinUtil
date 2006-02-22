#include "SelectEvents.h"
#include <iostream>

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>


using namespace lcio ;
using namespace marlin ;


SelectEvents aSelectEvents ;


SelectEvents::SelectEvents() : Processor("SelectEvents") {
  
  // modify processor description
  _description = "SelectEvent Processor selects certain events from input files" ;
  

  registerProcessorParameter( "FirstEvent" , 
			      "First Event"  ,
			      _firstEvent ,
			      int(0) ) ;

  registerProcessorParameter( "LastEvent" , 
			      "Last Event"  ,
			      _lastEvent ,
			      int(0) ) ;




}


void SelectEvents::init() { 

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void SelectEvents::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void SelectEvents::processEvent( LCEvent * evt ) { 

  if ( (_nEvt >= _firstEvent) && (_nEvt <= _lastEvent) ) setReturnValue(true);
  else setReturnValue(false);
    
  _nEvt ++ ;
}



void SelectEvents::check( LCEvent * evt ) { 

}


void SelectEvents::end(){ 

}

