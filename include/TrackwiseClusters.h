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

  TrackwiseClusters(const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes,const std::vector<float> startPoint,const std::vector<float> startDirection,
		    const TrackwiseClustersParameters* trackwiseClustersParameters, const TrackwiseClustersGeometryParameters* trackwiseClustersGeometryParameters);
  TrackwiseClusters(const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes,const float* startPoint,const float* startDirection,
		    const TrackwiseClustersParameters* trackwiseClustersParameters, const TrackwiseClustersGeometryParameters* trackwiseClustersGeometryParameters);

  ~TrackwiseClusters();


  std::vector<ClusterImpl*> doClustering();
  




  

 private:

  std::vector<float> _distanceTrackBack;
  std::vector<float> _stepTrackBack;
  std::vector<float> _resolutionParameter;
  std::vector<float> _distanceMergeForward;
  float _distanceToTrackSeed;
  float _distanceToDefineDirection;
  float _resolutionToMerge;
  int _nhit_merge_forward;
  int _use_tracks;
  int _nhit_minimal;
  int _typeOfGenericDistance;

  int _doMerging;
  int _doMergingForward;
  int _displayClusters;
  int _NDefineSP;
  int _nScanToMergeForward;

  ClusterExtendedVec _allSuperClusters;
  ClusterExtendedVec _allClusters;
  CaloHitExtendedVec _allHits;

  /** Parameters specifying generic geometry of 
   *  calorimeter system
   */

  // z position of ECAL Endcap front face
  float _zofendcap;
  // radius of ECAL Barrel
  float _rofbarrel;
  // offset in Phi angle of Barrel 
  // (Phi = 0 for canonical LC detector)
  float _phiofbarrel;
  // Factor defining N_fold symmetry
  // N = 8 for canonical LC detector
  int _nsymmetry;
  // Theta of ENDCAP = atan(_rofbarrel/_zofendcap)
  float _thetaofendcap;

  float _weightForReso;
  float _weightForDist;
  
  float _const_pi  ;
  float _const_2pi ;
  float _const_pi_n  ;
  float _const_2pi_n ;

  float _xmax_in_distance;
  float _xmin_in_distance;

  float _bField;


  std::vector<CalorimeterHitWithAttributes*> _calorimeterHitsWithAttributes;
  std::vector<float> _startPoint;
  std::vector<float> _startDirection;

  
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

} ;


#endif

