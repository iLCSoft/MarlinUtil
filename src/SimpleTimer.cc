#include "SimpleTimer.h"




using namespace lcio ;
using namespace marlin ;
using namespace UTIL;







SimpleTimer aSimpleTimer; // create instance


SimpleTimer::SimpleTimer() : Processor("SimpleTimer") {

  _description = "MARLIN Processor 'SimpleTimer', offers simple timer utilities" ;

  registerProcessorParameter("Mode","Mode",_mode,(int)0);
  registerProcessorParameter("SecondsToWait","Seconds to Wait",_secondsToWait,(int)0);

} 



void SimpleTimer::init() {
  
  _nRun = -1;
  _nEvt = 0;

  _timer = new LCTime();
  _time = _timer->unixTime();

}


void SimpleTimer::processRunHeader(LCRunHeader* run) { 
  _nRun++ ;
  _nEvt = 0;
} 


void SimpleTimer::processEvent(LCEvent* evt) {

  if (_mode == 0) wait(_secondsToWait);

  else if (_mode == 1) {

    if ( (_timer->unixTime() - _time) < _secondsToWait ) { 
      wait(_secondsToWait - (_timer->unixTime() - _time));
    }

  }
  else std::cout << "Wrong mode specified in processor 'SimpleTimer'" << std::endl;

  _time = _timer->unixTime();

}


void SimpleTimer::check(LCEvent* evt) { }

  
void SimpleTimer::end() {

  _time = 0;
  delete _timer;
  _timer = 0;

} 


void SimpleTimer::wait(int sleepTime) {
  sleep(sleepTime);
}
