#ifndef CutOnCosThetaQQbar_h
#define CutOnCosThetaQQbar_h 1

#include <iostream>

#include <marlin/Processor.h>
#include <lcio.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>


using namespace lcio ;
using namespace marlin ;

class CutOnCosThetaQQbar : public Processor {

 public:

  virtual Processor* newProcessor() { return new CutOnCosThetaQQbar ; }

  CutOnCosThetaQQbar();

  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ;   
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;

  
 private:

  
 protected:

  double _cutCosTheta;

  int _nRun;
  int _nEvt;


  double getCosThOfQuarkSystem(const LCEvent* evt);


} ;

#endif
