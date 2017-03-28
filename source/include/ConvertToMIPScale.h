#ifndef CONVERTTOMIPSCALE_H
#define CONVERTTOMIPSCALE_H 1

#include <iostream>
#include <string>
#include <vector>
#include "marlin/Processor.h"
#include "lcio.h"
#include <EVENT/LCCollection.h>
#include <EVENT/CalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>


using namespace lcio;
using namespace marlin;



class ConvertToMIPScale : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new ConvertToMIPScale; }
  
  
  ConvertToMIPScale();
  
  virtual void init();
  
  virtual void processRunHeader( LCRunHeader* run );
  
  virtual void processEvent( LCEvent * evt ); 
  
  
  virtual void check( LCEvent * evt ); 
  
  
  virtual void end();
  
  
 protected:

  int _nRun=-1;
  int _nEvt=-1;
  
  std::string _inputEcalCollection{};
  std::string _inputHcalCollection{};

  std::string _outputEcalCollection{};
  std::string _outputHcalCollection{};


  float _cutEcal=0.0;
  float _cutHcal=0.0;

  std::vector<float> _mipCoeffEcal{};
  std::vector<float> _mipCoeffHcal{};

};

#endif
