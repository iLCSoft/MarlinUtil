#include "SelectEvents.h"
#include <iostream>

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

void SelectEvents::processRunHeader( LCRunHeader* ) {

  _nRun++ ;
} 

void SelectEvents::processEvent( LCEvent* ) {

  if ( (_nEvt >= _firstEvent) && (_nEvt <= _lastEvent) ) setReturnValue(true);
  else setReturnValue(false);
    
  _nEvt ++ ;
}



void SelectEvents::check( LCEvent* ) {

}


void SelectEvents::end(){ 

}

