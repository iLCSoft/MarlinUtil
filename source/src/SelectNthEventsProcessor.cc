#include "SelectNthEventsProcessor.h"
#include <iostream>

using namespace lcio ;
using namespace marlin ;


SelectNthEventsProcessor aSelectNthEventsProcessor ;

SelectNthEventsProcessor::SelectNthEventsProcessor() : Processor("SelectNthEventsProcessor"){
  
  _description = "SelectNthEventsProcessor Processor sets ReturnValue true for every nth event. The selection can be offset and/or inverted.";
  

  registerProcessorParameter( "nEventToSelect", 
			      "Integer number n, every nth event's ReturnValue is set to true; default: 1",
			      _nEventToSelect,
			      int(1) );

  registerProcessorParameter( "SelectionOffset",
			      "Sets an offset m, so the Return Value of the mth event (and every nth before and after it) is set to true; default: 0",
			      _selectionOffset,
			      int(0) );

  registerProcessorParameter( "InvertSelection",
			      "Inverts the ReturnValue; if true, every nth event's ReturnValue is false, all others are true; default: false",
			      _invertSelection,
			      false );
}

void SelectNthEventsProcessor::init(){
  // usually a good idea to
  printParameters() ;

  _nRun = 0;
  _nEvt = 0;
  _nEvtProc = 0;
  _nEvtSkip = 0;
  
}

void SelectNthEventsProcessor::processEvent( LCEvent* ){

  bool skip ( (_nEvt-_selectionOffset)%_nEventToSelect!=0 );
  if( _invertSelection ) skip = !skip;

  setReturnValue(!skip);

  ++_nEvt;
  if (skip) ++_nEvtSkip;
  else ++_nEvtProc;
}

void SelectNthEventsProcessor::end(){ 
  streamlog_out(MESSAGE) << "Total Events: " << _nEvt << ",  Events processed: " << _nEvtProc << ",  Events skipped: " << _nEvtSkip << std::endl;
}
