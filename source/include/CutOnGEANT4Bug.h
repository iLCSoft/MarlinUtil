#ifndef CutOnGEANT4Bug_h
#define CutOnGEANT4Bug_h 1

#include <iostream>

#include <vector>

#include <marlin/Processor.h>
#include <lcio.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/CalorimeterHit.h>
#include <UTIL/LCRelationNavigator.h>

#include <MarlinUtil.h>


using namespace lcio ;
using namespace marlin ;

class CutOnGEANT4Bug : public Processor {

 public:

  virtual Processor* newProcessor() { return new CutOnGEANT4Bug ; }

  CutOnGEANT4Bug();

  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ;   
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;

  
 private:

  
 protected:

  double _eMin=5.0;
  double _k=1.75;

  std::string _colNameTracks="Tracks";
  std::string _colNameRelationTrackToMCP="TrueTrackToMCP";
  std::string _colNameRelationCaloHitToSimCaloHit="RelationCaloHit";

  std::vector<float> _calibrCoeffECAL{};
  std::vector<float> _calibrCoeffHCAL{};

  int _nRun=-1;
  int _nEvt=-1;

} ;

#endif
