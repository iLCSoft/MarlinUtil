#ifndef SIMPLETIMER_H
#define SIMPLETIMER_H 1

#include "marlin/Processor.h"
#include "lcio.h"

#include <UTIL/LCTime.h>

#include <iostream>
#include <string>

#include <unistd.h>


using namespace lcio ;
using namespace marlin ;






class SimpleTimer : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new SimpleTimer ; }
  
  SimpleTimer() ;
  
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent* evt ) ; 
  virtual void check( LCEvent* evt ) ; 
  virtual void end() ;
  
  
 protected:

  int _nRun;
  int _nEvt;

  int _time;
  LCTime* _timer;
  int _mode;
  int _secondsToWait;

  void wait(int sleepTime);

};

#endif
