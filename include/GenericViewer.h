#ifndef VIEWERAR_H
#define VIEWERAR_H 1

#include "marlin/Processor.h"
#include "EVENT/MCParticle.h"
#include "lcio.h"
#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;


/** Viewer Processor <br>
 *  Author : A.Raspereza <br>
 *  This processor displays reconstructed particles. <br>
 *  User has to specify the name of collection <br> 
 *  of reconstructed particles with <br>
 *  processor parameter ParticleCollection <br>
 */
class GenericViewer : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new GenericViewer ; }
  
  
  GenericViewer() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
  
  virtual void check( LCEvent * evt ) ; 
  
  virtual void end() ;
  
  
 protected:

  int _nRun ;
  int _nEvt ;
  
  std::vector<std::string> _caloHitCollections;
  std::vector<std::string> _simCaloHitCollections;
  std::vector<std::string> _trackerHitCollections;
  std::vector<std::string> _simTrackerHitCollections;
  std::string _trueClustersCollection;
  std::string _trueTracksCollection;
  std::string _clustersCollection;
  std::string _tracksCollection;
  std::string _particleCollection;

  int _layerCaloHit;
  int _layerSimCaloHit;
  int _layerTrackerHit;
  int _layerSimTrackerHit;
  int _layerTrueClusters;
  int _layerTrueTracks;
  int _layerClusters;
  int _layerTracks;
  int _layerMCP;
  int _layerBosons;
  int _layerReco;

  int _detModel;

  std::map<MCParticle *, int > _mcpList;
  
  int returnColor(int counter);

  float _bField;


} ;

#endif



