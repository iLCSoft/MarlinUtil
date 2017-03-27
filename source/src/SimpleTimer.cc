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

  _startTimer = new LCTime();
  _startTime = _startTimer->unixTime();
  _time = _startTime;

}


void SimpleTimer::processRunHeader(LCRunHeader*) {
  _nRun++ ;
  _nEvt = 0;
} 


void SimpleTimer::processEvent(LCEvent*) {

  _currentTimer = new LCTime();

  if (_mode == 0) {
    
    _secondsPerEvent = _currentTimer->unixTime() - _time;
    
    std::cout << "Processing time for event " << _nEvt << ":  " << _secondsPerEvent << " sec" << std::endl << std::endl;
  
    _time = _currentTimer->unixTime();
  
  }

  else if (_mode == 1) {

    wait(_secondsToWait);
    _time = _currentTimer->unixTime();

  }

  else if (_mode == 2) {

    if ( (_currentTimer->unixTime() - _time) < _secondsToWait ) { 
      wait(_secondsToWait - (_currentTimer->unixTime() - _time));
    }

    _time = _currentTimer->unixTime();

  }
  else std::cout << "Wrong mode specified in processor 'SimpleTimer'" << std::endl;

  delete _currentTimer;
  _currentTimer = 0;

  ++_nEvt;

}


void SimpleTimer::check(LCEvent*) { }

  
void SimpleTimer::end() {

  _currentTimer = new LCTime();

  _secondsOfJob = _currentTimer->unixTime() - _startTime;
  int remainder = _secondsOfJob%60;
  _minutesOfJob = _secondsOfJob/60;

  std::cout << "Processing time for " << _nEvt << " events in " << _nRun+1 << " run(s):  " << _minutesOfJob << ":" << remainder << " min (" << _secondsOfJob << " sec)" 
	    << std::endl << std::endl;

  _startTime = 0;;
  _time = 0;
  _secondsPerEvent = 0;
  _secondsOfJob = 0;
  _minutesOfJob = 0;

  delete _currentTimer;
  _currentTimer = 0;

  delete _startTimer;
  _startTimer = 0;

} 


void SimpleTimer::wait(int sleepTime) {
  sleep(sleepTime);
}
