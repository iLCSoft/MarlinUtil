#ifndef SkipNEvents_h
#define SkipNEvents_h 1

#include <iostream>
#include <marlin/Processor.h>


using namespace lcio ;
using namespace marlin ;

class SkipNEvents : public Processor {

 public:

  virtual Processor* newProcessor() { return new SkipNEvents ; }

  SkipNEvents();

  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ;   
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;

  
 private:

  
 protected:

  int _nSkip;

  int _nRun;
  int _nEvt;

} ;

#endif
