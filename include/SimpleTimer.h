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





/**
 *    Simple timimg processor, which offers a delay after an event has been processed.
 *
 *    @param SecondsToWait : time to wait (in seconds)

 *    @param Mode : toggle for two different timing modes <br>
 *                  0 - no wait, just display time (in sec.) needed per event and summarise time of processing all events
 *                  1 - just wait the time which is specified in the parameter
 *                      SecondsToWait <br>
 *                  2 - calculate the elapsed time since the last call of this processor
 *                      (e.g. in the last event) and wait for the remaining time 
 *                      (SecondsToWait minus elapsed time). If already more time elapsed 
 *                      than SecondsToWait, carry on without any delay.
 *    @author O. Wendt (DESY)
 *    @version $Id: SimpleTimer.h,v 1.4 2006-02-09 18:01:55 owendt Exp $
 *
 */
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

  int _startTime;
  int _time;
  LCTime* _startTimer;
  LCTime* _currentTimer;
  int _mode;
  int _secondsToWait;
  int _secondsPerEvent;
  int _secondsOfJob;
  int _minutesOfJob;

  void wait(int sleepTime);

};

#endif
