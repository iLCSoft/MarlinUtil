#include "SkipNEvents.h"

using namespace lcio ;
using namespace marlin ;

SkipNEvents aSkipNEvents ;



SkipNEvents::SkipNEvents() : Processor("SkipNEvents") {

  _description = "The output condition of this processor is false for the first n LCEvents. Afterwards it is set to true.";


  registerProcessorParameter( "nSkip",
			      "number of LCEvents to skip",
			      _nSkip,
			      (int)0);

}


void SkipNEvents::init() {
  
  // usually a good idea to 
  // printParameters();
  
  _nRun = 0 ;
  _nEvt = 0 ;


}


void SkipNEvents::processRunHeader( LCRunHeader* ) {

  ++_nRun;

}


void SkipNEvents::processEvent( LCEvent * ) {

  if ( _nEvt < _nSkip ) setReturnValue(false);
  else setReturnValue(true);

  ++_nEvt;

}


void SkipNEvents::check( LCEvent * ) {
 
}


void SkipNEvents::end() {

}
