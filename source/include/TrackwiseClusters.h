#ifndef TRACKWISECLUSTERS_H
#define TRACKWISECLUSTERS_H 1

#include <iostream>

#include <string>
#include <vector>
#include <math.h>

#include "EVENT/LCIO.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/Track.h"
#include "IMPL/CalorimeterHitImpl.h"
#include "IMPL/ClusterImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "CaloHitExtended.h"
#include "ClusterExtended.h"
#include "TrackExtended.h"

#include "ClusterShapes.h"
#include "CalorimeterHitWithAttributes.h"

// GEAR include files
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/CalorimeterParameters.h>
#include "random.h"

// for debugging only
#include <MarlinCED.h>


struct TrackwiseClustersParameters {

  std::vector<float> distanceTrackBack;
  std::vector<float> stepTrackBack;
  std::vector<float> resolutionParameter;
  std::vector<float> distanceMergeForward;
  float distanceToTrackSeed;
  float distanceToDefineDirection;
  float resolutionToMerge;
  int nhit_merge_forward;
  int nhit_minimal;
  int typeOfGenericDistance;

  int doMerging;
  int doMergingForward;
  int displayClusters;
  int NDefineSP;
  int nScanToMergeForward;

};


struct TrackwiseClustersGeometryParameters {

  // z position of ECAL Endcap front face
  float zofendcap;
  // radius of ECAL Barrel
  float rofbarrel;
  // offset in Phi angle of Barrel 
  // (Phi = 0 for canonical LC detector)
  float phiofbarrel;
  // Factor defining N_fold symmetry
  // N = 8 for canonical LC detector
  int nsymmetry;
  // Theta of ENDCAP = atan(rofbarrel/zofendcap)
  float thetaofendcap;

  float weightForReso;
  float weightForDist;

  float bField;

};


using namespace lcio;

 
class TrackwiseClusters { 

 public:

  TrackwiseClusters(const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes,const std::vector<float> startPoint,
		    const float pathLengthOnHelixOfStartPoint, const float distanceToHelixOfStartPoint, const std::vector<float> startDirection,
		    const TrackwiseClustersParameters* trackwiseClustersParameters, const TrackwiseClustersGeometryParameters* trackwiseClustersGeometryParameters);
  TrackwiseClusters(const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes,const float* startPoint,
		    const float pathLengthOnHelixOfStartPoint, const float distanceToHelixOfStartPoint, const float* startDirection,
		    const TrackwiseClustersParameters* trackwiseClustersParameters, const TrackwiseClustersGeometryParameters* trackwiseClustersGeometryParameters);

  ~TrackwiseClusters();


  std::vector<ClusterImpl*> doClustering();
  




  

 private:

  std::vector<float> _distanceTrackBack{};
  std::vector<float> _stepTrackBack{};
  std::vector<float> _resolutionParameter{};
  std::vector<float> _distanceMergeForward{};
  float _distanceToTrackSeed=0.0;
  float _distanceToDefineDirection=0.0;
  float _resolutionToMerge=0.0;
  int _nhit_merge_forward=0;
  int _use_tracks=0;
  int _nhit_minimal=0;
  int _typeOfGenericDistance=0;

  int _doMerging=0;
  int _doMergingForward=0;
  int _displayClusters=0;
  int _NDefineSP=0;
  int _nScanToMergeForward=0;

  ClusterExtendedVec _allSuperClusters{};
  ClusterExtendedVec _allClusters{};
  CaloHitExtendedVec _allHits{};

  /** Parameters specifying generic geometry of 
   *  calorimeter system
   */

  // z position of ECAL Endcap front face
  float _zofendcap=0.0;
  // radius of ECAL Barrel
  float _rofbarrel=0.0;
  // offset in Phi angle of Barrel 
  // (Phi = 0 for canonical LC detector)
  float _phiofbarrel=0.0;
  // Factor defining N_fold symmetry
  // N = 8 for canonical LC detector
  int _nsymmetry=0;
  // Theta of ENDCAP = atan(_rofbarrel/_zofendcap)
  float _thetaofendcap=0.0;

  float _weightForReso=0.0;
  float _weightForDist=0.0;

  float _const_pi=M_PI;
  float _const_2pi=2.0*M_PI;
  float _const_pi_n=0.0;
  float _const_2pi_n=0.0;

  float _xmax_in_distance=0.0;
  float _xmin_in_distance=0.0;

  float _bField=0.0;

  int _debugLevel=0;


  std::vector<CalorimeterHitWithAttributes*> _calorimeterHitsWithAttributes{};
  std::vector<float> _startPoint{};
  float _pathLengthOnHelixOfStartPoint=0.0;
  float _distanceToHelixOfStartPoint=0.0;
  std::vector<float> _startDirection{};
  
  void initialiseCollections();
  float findResolutionParameter(CaloHitExtended* fromHit, CaloHitExtended* toHit);
  void CalculateGenericDistance(CaloHitExtended* calohit, float* dist); 
  void BubbleSort(CaloHitExtendedVec& input);  
  float DistanceBetweenPoints(float* x1, float* x2);
  void DisplayClusters(ClusterExtendedVec clusterVec);
  void GlobalSorting();
  void GlobalClustering();
  std::vector<ClusterImpl*> CreateClusterCollection(ClusterExtendedVec clusterVec);
  void mergeForward();
  void mergeLowMultiplicity();
  void calculateProperties(ClusterExtended* Cl);
  void propertiesForAll();
  void CleanUp();
  CalorimeterHitWithAttributes* getCalorimeterHitWithAttributes(CaloHitExtended* calorimeterHitExtended);

} ;


#endif

