#include "TrackwiseClusters.h"

#include <cfloat>

// make gcc > 4.7 compliant
#include <unistd.h>


using namespace lcio;



// ____________________________________________________________________________________________________________________________________________________________________________



TrackwiseClusters::TrackwiseClusters(const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes, std::vector<float> const startPoint, 
				     const float pathLengthOnHelixOfStartPoint, const float distanceToHelixOfStartPoint, 
				     const std::vector<float> startDirection, const TrackwiseClustersParameters* trackwiseClustersParameters, 
				     const TrackwiseClustersGeometryParameters* trackwiseClustersGeometryParameters) {
  

  _calorimeterHitsWithAttributes = calorimeterHitsWithAttributes;
  _startPoint = startPoint;
  _pathLengthOnHelixOfStartPoint = pathLengthOnHelixOfStartPoint;
  _distanceToHelixOfStartPoint = distanceToHelixOfStartPoint;
  _startDirection = startDirection;


  _distanceTrackBack         = trackwiseClustersParameters->distanceTrackBack;
  _stepTrackBack             = trackwiseClustersParameters->stepTrackBack;
  _resolutionParameter       = trackwiseClustersParameters->resolutionParameter;
  _distanceMergeForward      = trackwiseClustersParameters->distanceMergeForward;
  _distanceToTrackSeed       = trackwiseClustersParameters->distanceToTrackSeed;
  _distanceToDefineDirection = trackwiseClustersParameters->distanceToDefineDirection;
  _resolutionToMerge         = trackwiseClustersParameters->resolutionToMerge;
  _nhit_merge_forward        = trackwiseClustersParameters->nhit_merge_forward;

  _nhit_minimal              = trackwiseClustersParameters->nhit_minimal;
  _typeOfGenericDistance     = trackwiseClustersParameters->typeOfGenericDistance;
  _doMerging                 = trackwiseClustersParameters->doMerging;
  _doMergingForward          = trackwiseClustersParameters->doMergingForward;
  _displayClusters           = trackwiseClustersParameters->displayClusters;
  _NDefineSP                 = trackwiseClustersParameters->NDefineSP;
  _nScanToMergeForward       = trackwiseClustersParameters->nScanToMergeForward;

  
  _zofendcap                 = trackwiseClustersGeometryParameters->zofendcap;
  _rofbarrel                 = trackwiseClustersGeometryParameters->rofbarrel;
  _phiofbarrel               = trackwiseClustersGeometryParameters->phiofbarrel;
  _nsymmetry                 = trackwiseClustersGeometryParameters->nsymmetry;
  _thetaofendcap             = trackwiseClustersGeometryParameters->thetaofendcap;

  _weightForReso             = trackwiseClustersGeometryParameters->weightForReso;
  _weightForDist             = trackwiseClustersGeometryParameters->weightForDist;

  _bField                    = trackwiseClustersGeometryParameters->bField;



  _allHits.clear();
  _allSuperClusters.clear();
  _allClusters.clear();
  

  _const_2pi = 2.0*M_PI;
  _const_pi_n  = _const_pi/float(_nsymmetry);
  _const_2pi_n = 2.0*_const_pi/float(_nsymmetry);
  _thetaofendcap = (float)atan((double)(_rofbarrel/_zofendcap));

   
  _xmin_in_distance = 1.0e+10;
  _xmax_in_distance = -1.0e+10;


  _debugLevel = 0;

}



TrackwiseClusters::TrackwiseClusters(const std::vector<CalorimeterHitWithAttributes*> calorimeterHitsWithAttributes,const float* startPoint,
				     const float pathLengthOnHelixOfStartPoint, const float distanceToHelixOfStartPoint,
				     const float* startDirection, const TrackwiseClustersParameters* trackwiseClustersParameters, 
				     const TrackwiseClustersGeometryParameters* trackwiseClustersGeometryParameters) {



  _calorimeterHitsWithAttributes = calorimeterHitsWithAttributes;  
  
  _pathLengthOnHelixOfStartPoint = pathLengthOnHelixOfStartPoint;
  _distanceToHelixOfStartPoint = distanceToHelixOfStartPoint;
  
  for (int i = 0; i < 3; ++i) {
    
    _startPoint.push_back(startPoint[i]);
    _startDirection.push_back(startDirection[i]);

  }


  _distanceTrackBack         = trackwiseClustersParameters->distanceTrackBack;
  _stepTrackBack             = trackwiseClustersParameters->stepTrackBack;
  _resolutionParameter       = trackwiseClustersParameters->resolutionParameter;
  _distanceMergeForward      = trackwiseClustersParameters->distanceMergeForward;
  _distanceToTrackSeed       = trackwiseClustersParameters->distanceToTrackSeed;
  _distanceToDefineDirection = trackwiseClustersParameters->distanceToDefineDirection;
  _resolutionToMerge         = trackwiseClustersParameters->resolutionToMerge;
  _nhit_merge_forward        = trackwiseClustersParameters->nhit_merge_forward;

  _nhit_minimal              = trackwiseClustersParameters->nhit_minimal;
  _typeOfGenericDistance     = trackwiseClustersParameters->typeOfGenericDistance;
  _doMerging                 = trackwiseClustersParameters->doMerging;
  _doMergingForward          = trackwiseClustersParameters->doMergingForward;
  _displayClusters           = trackwiseClustersParameters->displayClusters;
  _NDefineSP                 = trackwiseClustersParameters->NDefineSP;
  _nScanToMergeForward       = trackwiseClustersParameters->nScanToMergeForward;

  
  _zofendcap                 = trackwiseClustersGeometryParameters->zofendcap;
  _rofbarrel                 = trackwiseClustersGeometryParameters->rofbarrel;
  _phiofbarrel               = trackwiseClustersGeometryParameters->phiofbarrel;
  _nsymmetry                 = trackwiseClustersGeometryParameters->nsymmetry;
  _thetaofendcap             = trackwiseClustersGeometryParameters->thetaofendcap;

  _weightForReso             = trackwiseClustersGeometryParameters->weightForReso;
  _weightForDist             = trackwiseClustersGeometryParameters->weightForDist;

  _bField                    = trackwiseClustersGeometryParameters->bField;



  _allHits.clear();
  _allSuperClusters.clear();
  _allClusters.clear();
  

  _const_2pi = 2.0*M_PI;
  _const_pi_n  = _const_pi/float(_nsymmetry);
  _const_2pi_n = 2.0*_const_pi/float(_nsymmetry);
  _thetaofendcap = (float)atan((double)(_rofbarrel/_zofendcap));

   
  _xmin_in_distance = 1.0e+10;
  _xmax_in_distance = -1.0e+10;


  _debugLevel = 0;

}



TrackwiseClusters::~TrackwiseClusters() {
  
  _allHits.clear();
  _allSuperClusters.clear();
  _allClusters.clear();
  
}



// ____________________________________________________________________________________________________________________________________________________________________________



std::vector<ClusterImpl*> TrackwiseClusters::doClustering() {
  
  std::vector<ClusterImpl*> resultingClusters;

  initialiseCollections();
  GlobalSorting();
  GlobalClustering();
  // propertiesForAll();

  if (_doMergingForward == 1) mergeForward();

  if (_displayClusters == 1) DisplayClusters(_allClusters);

  resultingClusters = CreateClusterCollection(_allClusters);

  CleanUp();

  return resultingClusters;

}




void TrackwiseClusters::initialiseCollections() {


  int nHits = 0;
  std::cout << "Trackwise Cluster!!!!" << std::endl;
  sleep(200);

  nHits = _calorimeterHitsWithAttributes.size();

  // debug
  if ( _debugLevel > 5 ) { 
    std::cout << "Start Point: " << "(" << _startPoint.at(0) << "," << _startPoint.at(1) << "," << _startPoint.at(2) << ")" << "  " 
	      << "s on helix: " << _pathLengthOnHelixOfStartPoint << "  " << "dist to helix: " << _distanceToHelixOfStartPoint << std::endl;
    ced_hit ( _startPoint.at(0),_startPoint.at(1),_startPoint.at(2), 2 | 5 << CED_LAYER_SHIFT, 5, 0xff22c8 );
    ced_send_event();
  }
  
  
  for (int j =0; j < nHits; ++j) {
    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(_calorimeterHitsWithAttributes.at(j)->getCalorimeterHit());

    int type = hit->getType();
    
    // convert to be compatible to the TrackwiseClustering Processor
    if ( (type == 0) || (type==1) ) type = 0;
    else type = 1;

    CaloHitExtended *calohit = new CaloHitExtended(hit,type);
    float dist[2];
    CalculateGenericDistance(calohit, dist);		    
    
    if (_typeOfGenericDistance == 0) calohit->setGenericDistance(dist[0]);
    else if (_typeOfGenericDistance == 1) calohit->setGenericDistance(dist[1]);
    else { 
     
      float pathLengthOnHelixOfCaloHit = _calorimeterHitsWithAttributes.at(j)->getPathLengthOnHelix();
      float distanceToHelixOfCaloHit = _calorimeterHitsWithAttributes.at(j)->getDistanceToHelix();

      float distanceInHelixCoordinates = sqrt( pow( (pathLengthOnHelixOfCaloHit - _pathLengthOnHelixOfStartPoint),2) + 
					       pow( (distanceToHelixOfCaloHit - _distanceToHelixOfStartPoint),2) );


      calohit->setGenericDistance(distanceInHelixCoordinates);

      // debug 
      if ( _debugLevel > 5 ) { 
	std::cout << "CaloHit: " << hit << "  " << "_typeOfGenericDistance = " << _typeOfGenericDistance << "  " 
		  << "s on helix: " << _calorimeterHitsWithAttributes.at(j)->getPathLengthOnHelix() << "  " 
		  << "dist to helix: " << _calorimeterHitsWithAttributes.at(j)->getDistanceToHelix() << "  "
		  << "dist to start point: " << distanceInHelixCoordinates << std::endl;
      }

    }
  
    calohit->setDistanceToCalo(dist[1]); // never used !!!!
    float distance = calohit->getGenericDistance();
    if (distance < _xmin_in_distance) 
      _xmin_in_distance = distance;
    if (distance > _xmax_in_distance) 
      _xmax_in_distance = distance;		
    _allHits.push_back(calohit);

  }

}






void TrackwiseClusters::CalculateGenericDistance(CaloHitExtended * calohit, float * dist) {

  float xDistance =0.0;
  float rDistance =0.0;
  
  for (int i(0); i < 3; ++i) {
    float x = calohit->getCalorimeterHit()->getPosition()[i];
    rDistance += x*x; 	
  }
  rDistance = sqrt(rDistance);
  
  float x = calohit->getCalorimeterHit()->getPosition()[0];
  float y = calohit->getCalorimeterHit()->getPosition()[1];
  float z = calohit->getCalorimeterHit()->getPosition()[2];
  float phi = atan2(y,x) - _phiofbarrel + _const_pi_n;
  int nZone = (int)(phi/_const_2pi_n);
  if (phi < 0.)
    phi = phi + _const_2pi;
  phi = phi - nZone * _const_2pi_n - _const_pi_n;
  float radius = sqrt(x*x + y*y);
  float rdist = radius * cos(phi) - _rofbarrel;
  float zdist = fabs(z) - _zofendcap;
  if (rdist > 0 && zdist < 0) {
    xDistance = rdist;
  }
  else if (rdist < 0 && zdist > 0 ) {
    xDistance = zdist;
  }
  else {
    float theta = (float)atan((float)(rdist/zdist));
    if (theta > _thetaofendcap) {
      xDistance = rdist;
    }
    else {
      xDistance = zdist;
    }
    
  }
  xDistance = xDistance + 1.0e-10*rDistance ;
  dist[0] = rDistance;
  dist[1] = xDistance;
  
} 



void TrackwiseClusters::GlobalSorting() {


  int _NGLAYERS = 1000;
  
  float inverse = ((float)_NGLAYERS)/(_xmax_in_distance - _xmin_in_distance);
  
  _allSuperClusters.clear();
  
  for (int i=0; i < _NGLAYERS; ++i) {
    ClusterExtended * cluster = new ClusterExtended();
    _allSuperClusters.push_back(cluster);
  }
  
  for (unsigned int i = 0; i < _allHits.size(); ++i) {
    CaloHitExtended * calohit = _allHits[i];
    float distToIP = calohit->getGenericDistance();
    int index = (int)((distToIP - _xmin_in_distance) * inverse) ;
    if (index >= _NGLAYERS)
      index = _NGLAYERS - 1;
    if (index < 0) 
      index = 0;
    
    _allSuperClusters[index]->addCaloHitExtended(calohit);
  }
  
  _allHits.clear();

  int counter = 0;

  for (int i(0); i < _NGLAYERS; ++i) {
    ClusterExtended * cluster = _allSuperClusters[i];
    CaloHitExtendedVec calohitvec = cluster->getCaloHitExtendedVec();
    int nhits = calohitvec.size();
    if (nhits > 0) {
      BubbleSort(calohitvec);
      for (int ihit(0); ihit < nhits; ++ihit) {
	CaloHitExtended * calohit = calohitvec[ihit];
	calohit->setIndex(counter);
	_allHits.push_back(calohit);
	counter++;
      }
    }
    
  }
  
}



void TrackwiseClusters::BubbleSort(CaloHitExtendedVec& input) {

  unsigned int sizeOfVector = input.size();
  CaloHitExtended *one,*two,*Temp;
  
  for (unsigned int i = 0 ; i < sizeOfVector-1; i++)
    for (unsigned int j = 0; j < sizeOfVector-i-1; j++)
      {
	one = input.at(j);
	two = input.at(j+1);
	if(one->getGenericDistance() > two->getGenericDistance())
	  {
	    Temp = input[j];
	    input[j] = input[j+1];
	    input[j+1] = Temp;
	  }
      }
  
}



void TrackwiseClusters::GlobalClustering() {
  
  
  _allClusters.clear();

  // insert a pseudo hit with the position and direction _startPoint and _startDirection;
  /*    
  CalorimeterHitImpl* pseudoHitImpl = new CalorimeterHitImpl();
  pseudoHitImpl->setEnergy(0.0);
  float position[3];
  for (unsigned int i = 0; i < 3; ++i) position[i] = _startPoint.at(i);
  pseudoHitImpl->setPosition(position);
  CalorimeterHit* pseudoHit = static_cast<CalorimeterHit*>(pseudoHitImpl);  
  CaloHitExtended* pseudoHitExtended = new CaloHitExtended(pseudoHit,0); // change type!!!!!!!!!!!!!!!!!!!!!!!!!!!! needs to be given by MIP stub finder
  pseudoHitExtended->setGenericDistance(0.0);  // change generic distance !!!!!!!!!!!!!!!!!!!!!!! needs to be given by MIP stub finder

  _allHits.insert(_allHits.begin(),pseudoHitExtended);
  */


  // debug
  if ( _debugLevel > 5 ) {

    std::cout << "n of hits in Trackwise Clustering = " << _allHits.size() << std::endl << std::endl;
    for (unsigned int iHit = 0; iHit < _allHits.size(); ++iHit) {

      CaloHitExtended* caloHit = _allHits[iHit];
      CalorimeterHitWithAttributes* calorimeterHitWithAttributesCorrespondingToCaloHitTo = getCalorimeterHitWithAttributes(caloHit);

      std::cout << "CaloHit: " << caloHit << "  " << "_typeOfGenericDistance = " << _typeOfGenericDistance << "  " 
		<< "type of CaloHit: " << caloHit->getType() << "  " 
		<< "s on helix: " << ((calorimeterHitWithAttributesCorrespondingToCaloHitTo->getPathLengthOnHelix()) - _pathLengthOnHelixOfStartPoint) << "  " 
		<< "dist to helix: " << ((calorimeterHitWithAttributesCorrespondingToCaloHitTo->getDistanceToHelix()) - _distanceToHelixOfStartPoint) << "  "
		<< "dist to start point: " << caloHit->getGenericDistance() << std::endl;

    }

  }


  for (unsigned int ihitTo(0); ihitTo < _allHits.size(); ++ihitTo) {

    CaloHitExtended * CaloHitTo = _allHits[ihitTo];

    int ihitFrom = ihitTo - 1;
    int ifound = 0;
    int idTo = CaloHitTo->getType();

    float r_step = _stepTrackBack[idTo];
    float r_min  = r_step;
    float r_dist = _distanceTrackBack[idTo]; 
    
    float YResMin = 1.0e+10;
    float YResCut = _resolutionParameter[idTo];
    float YDistMin = 1.0e+10;

    
    // debug 
    if ( _debugLevel > 5 ) { 

      CalorimeterHitWithAttributes* calorimeterHitWithAttributesCorrespondingToCaloHitTo = getCalorimeterHitWithAttributes(CaloHitTo);

      std::cout << "before while loop: " << std::endl
		<< "ihitTo: " << ihitTo << "  " << "Hit : " << CaloHitTo->getCalorimeterHit() << "  " << "type: " << CaloHitTo->getType() << "  " 
		<< "assigned generic distance: " << CaloHitTo->getGenericDistance() << "  " 
		<< "s: " << ((calorimeterHitWithAttributesCorrespondingToCaloHitTo->getPathLengthOnHelix()) - _pathLengthOnHelixOfStartPoint) << "  "
		<< "d: " << ((calorimeterHitWithAttributesCorrespondingToCaloHitTo->getDistanceToHelix()) - _distanceToHelixOfStartPoint) << std::endl
		<< "ihitFrom: " << ihitFrom << "  " << "ifound: " << ifound << "  " << "r_step: " << r_step << "  " << "r_min: " << r_min << "  " << "r_dist: " << r_dist 
		<< std::endl
		<< "YResMin: " << YResMin << "  " << "YResCut: " << YResCut << "  " << "YDistMin: " << YDistMin << std::endl;

      ced_hit ( CaloHitTo->getCalorimeterHit()->getPosition()[0],CaloHitTo->getCalorimeterHit()->getPosition()[1],CaloHitTo->getCalorimeterHit()->getPosition()[2], 
		2 | 1 << CED_LAYER_SHIFT, 6, 0xff0000 );
      ced_send_event();
      getchar();
    }

    // debug
    int ihitFrominitialValue = ihitFrom;
    
    while (ihitFrom >=0) {
      
      CaloHitExtended * CaloHitFrom = _allHits[ihitFrom];
      float dist_in_generic = DBL_MAX;

      float XDist = 0.0;
      float YRes  = 0.0;

      if ( (_typeOfGenericDistance == 0) ||  (_typeOfGenericDistance == 1) ) { 

	dist_in_generic = CaloHitTo->getGenericDistance()-CaloHitFrom->getGenericDistance();

	float pos1[3];
	float pos2[3];
	for (int iposi=0; iposi<3; ++iposi) {
	  pos1[iposi] = 
	    (float)CaloHitTo->getCalorimeterHit()->getPosition()[iposi];
	  pos2[iposi] = 
	    (float)CaloHitFrom->getCalorimeterHit()->getPosition()[iposi];
	}
      
	XDist = DistanceBetweenPoints(pos1,pos2);

	YRes = findResolutionParameter(CaloHitFrom, CaloHitTo);

      }      
      else { 
     
	CalorimeterHitWithAttributes* calorimeterHitWithAttributesCorrespondingToCaloHitTo = getCalorimeterHitWithAttributes(CaloHitTo);
	CalorimeterHitWithAttributes* calorimeterHitWithAttributesCorrespondingToCaloHitFrom = getCalorimeterHitWithAttributes(CaloHitFrom);
	
	float pathLengthOnHelixOfCaloHitTo = calorimeterHitWithAttributesCorrespondingToCaloHitTo->getPathLengthOnHelix();
	float pathLengthOnHelixOfCaloHitFrom = calorimeterHitWithAttributesCorrespondingToCaloHitFrom->getPathLengthOnHelix();

	float distanceToHelixOfCaloHitTo = calorimeterHitWithAttributesCorrespondingToCaloHitTo->getDistanceToHelix();
	float distanceToHelixOfCaloHitFrom = calorimeterHitWithAttributesCorrespondingToCaloHitFrom->getDistanceToHelix();

	dist_in_generic = CaloHitTo->getGenericDistance()-CaloHitFrom->getGenericDistance();

	XDist = sqrt( pow( (pathLengthOnHelixOfCaloHitTo - pathLengthOnHelixOfCaloHitFrom),2) + 
		      pow( (distanceToHelixOfCaloHitTo - distanceToHelixOfCaloHitFrom),2) );


	YRes = findResolutionParameter(CaloHitFrom, CaloHitTo);

      }
      



      if ( _debugLevel > 5 ) { 

	CalorimeterHitWithAttributes* calorimeterHitWithAttributesCorrespondingToCaloHitTo = getCalorimeterHitWithAttributes(CaloHitTo);
	CalorimeterHitWithAttributes* calorimeterHitWithAttributesCorrespondingToCaloHitFrom = getCalorimeterHitWithAttributes(CaloHitFrom);

	std::cout << "in while loop: " << std::endl
		  << "ihitTo: " << ihitTo << "  " << "Hit : " << CaloHitTo->getCalorimeterHit() << "  " << "type: " << CaloHitTo->getType() << "  " 
		  << "assigned generic distance: " << CaloHitTo->getGenericDistance()  << "  " 
		  << "s: " << calorimeterHitWithAttributesCorrespondingToCaloHitTo->getPathLengthOnHelix() - _pathLengthOnHelixOfStartPoint << "  "
		  << "d: " << calorimeterHitWithAttributesCorrespondingToCaloHitTo->getDistanceToHelix() - _distanceToHelixOfStartPoint << std::endl
		  << "ihitFrom: " << ihitFrom << "  " << "Hit: " << CaloHitFrom->getCalorimeterHit() << "  " << "type: " << CaloHitFrom->getType() << "  " 
		  << "assigned generic distance: " << CaloHitFrom->getGenericDistance()  << "  " 
		  << "s: " << calorimeterHitWithAttributesCorrespondingToCaloHitFrom->getPathLengthOnHelix() - _pathLengthOnHelixOfStartPoint << "  "
		  << "d: " << calorimeterHitWithAttributesCorrespondingToCaloHitFrom->getDistanceToHelix() - _distanceToHelixOfStartPoint << std::endl
		  << "ifound: " << ifound << "  " << "r_step: " << r_step << "  " << "r_min: " << r_min << "  " << "r_dist: " << r_dist << std::endl
		  << "YResMin: " << YResMin << "  " << "YResCut: " << YResCut << "  " << "YDistMin: " << YDistMin << std::endl;
	std::cout << "dist_in_generic = " << dist_in_generic << std::endl;
	
	if ( ihitFrom == ihitFrominitialValue ) {
	  ced_hit ( CaloHitFrom->getCalorimeterHit()->getPosition()[0],CaloHitFrom->getCalorimeterHit()->getPosition()[1],CaloHitFrom->getCalorimeterHit()->getPosition()[2], 
		    0 | 1 << CED_LAYER_SHIFT, 2, 0xf2ff00 );
	  ced_send_event();
	  getchar();
	}

      }


	    
    


      // relaxing cut on 'distance'
      if (dist_in_generic > r_min) {
	if (ifound ==1) {

	  if ( _debugLevel > 5 ) { 	  
	    std::cout << "dist_in_generic: " << dist_in_generic << "  " << "r_min: " << r_min << "  " << "dist_in_generic > r_min: " << (dist_in_generic > r_min) << "  "
		      << "and ifound: " << ifound << " => BREAK while loop" << std::endl;
	  }

	  break;

	}

	r_min +=  r_step; 

      }

      // break on max distance allowed      
      if (dist_in_generic > r_dist) {

	if ( _debugLevel > 5 ) { 
	  std::cout << "dist_in_generic: " << dist_in_generic << "  " << "r_dist: " << r_dist << "  " << "dist_in_generic > r_dist: " << (dist_in_generic > r_dist) << "  "
		    << " => BREAK while loop" << std::endl;
	}

	break;

      }

      if (YRes < 0.)
	std::cout << "Resolution parameter < 0" << std::endl; 
      
      
      float YDist = 1 + _weightForReso*YRes + _weightForDist*XDist;
      
      bool proxCriterion = YDist < YDistMin;
         
      if ( _debugLevel > 5 ) { 
	std::cout << "dist_in_generic: " << dist_in_generic << "  " << "r_min: " << r_min << "  " << "(dist_in_generic > r_min): " << (dist_in_generic > r_min) << "  " 
		  << "r_dist: " << r_dist << "  " << "(dist_in_generic > r_dist): " << (dist_in_generic > r_dist) << std::endl
		  << "XDist: " << XDist << "  " << "YRes: " << YRes << "  " << "YResCut: " << YResCut << "  " << "YDist: " << YDist << "  " << "YDistMin: " << YDistMin << "  " 
		  << "YDist < YDistMin (proxCriterion): " << proxCriterion << std::endl;
      }

      if ( proxCriterion ) {
	YResMin = YRes;
	YDistMin = YDist;
	CaloHitTo->setCaloHitFrom(CaloHitFrom);
	CaloHitTo->setYresFrom(YRes);
	
	  if ( _debugLevel > 5 ) { 
	    std::cout << "proxCriterion fullfilled (" << proxCriterion << "): " << "YResMin = YRes = " << YResMin << "  " << "YDistMin = YDist = " << YDistMin << std::endl;
	  }


      }
      
      if ( proxCriterion && YRes < YResCut ) {
	YResMin = YRes;
	YDistMin = YDist;
	CaloHitTo->setCaloHitFrom(CaloHitFrom);
	CaloHitTo->setYresFrom(YRes);		

	if ( _debugLevel > 5 ) {
	  std::cout << "proxCriterion and YRes < YResCut fullfilled (" << (proxCriterion && YRes < YResCut) << "): " << "YResMin = YRes = " << YResMin << "  " 
		    << "YDistMin = YDist = " << YDistMin << std::endl;
	}

      }
      
      
      if ( YRes < YResCut ) {

	ifound = 1;

	if ( _debugLevel > 5 ) std::cout << "YRes < YResCut fullfilled (" << (YRes < YResCut) << "): " << "ifound = " << ifound << std::endl;
	
      }


      if ( _debugLevel > 5 ) { 
	std::cout << "END of while loop, reducing ihitFrom by one: " << ihitFrom << " -> " << (ihitFrom-1) << "  " << "(if < 0 -> while loop will end)" 
		  << std::endl << std::endl;
      }

      
      ihitFrom --;		
    }
    
    
    if (ifound == 1) { // Attach to already existing cluster
      CaloHitExtended * calohit_AttachTo = CaloHitTo->getCaloHitFrom();
      ClusterExtended * cluster = calohit_AttachTo->getClusterExtended();
      CaloHitExtendedVec calohitvec = cluster->getCaloHitExtendedVec();
      CaloHitTo->setClusterExtended(cluster);
      cluster->addCaloHitExtended(CaloHitTo);

      // debug
      if ( _debugLevel > 5 ) { 
	std::cout << "assign hit to existing cluster ... " << std::endl;
	int color = static_cast<int>(16*reinterpret_cast<long int>(cluster)+0x009988);
	ced_hit ( CaloHitTo->getCalorimeterHit()->getPosition()[0],CaloHitTo->getCalorimeterHit()->getPosition()[1],CaloHitTo->getCalorimeterHit()->getPosition()[2], 
		  0 | 1 << CED_LAYER_SHIFT, 8, color );
	ced_send_event();
	getchar();
      }

      float distanceToHit = 0.0;
      for (int ii=0;ii<3;++ii) {
	float xx =  CaloHitTo->getCalorimeterHit()->getPosition()[ii]
	  -  calohit_AttachTo->getCalorimeterHit()->getPosition()[ii];
	distanceToHit += xx*xx;
      }
      distanceToHit = sqrt(distanceToHit);
      CaloHitTo->setDistanceToNearestHit(distanceToHit);
      float xDir[3];
      bool redefineSP ;
      float dif_in_dist = CaloHitTo->getGenericDistance() - calohit_AttachTo->getGenericDistance();	    
      if (_typeOfGenericDistance == 0) {
	
	redefineSP = (int)calohitvec.size() < _NDefineSP;
	
      }
      
      else {
	
	redefineSP = dif_in_dist < _distanceToDefineDirection;
	
      }
      
      if (redefineSP ) {
	float xx = 0.;
	float yy = 0.;
	float zz = 0.;
	float ee = 0.;
	for (unsigned int i(0); i < calohitvec.size(); ++i) {
	  CaloHitExtended * chit = calohitvec[i];
	  float ene = chit->getCalorimeterHit()->getEnergy();
	  xx += chit->getCalorimeterHit()->getPosition()[0]*ene;
	  yy += chit->getCalorimeterHit()->getPosition()[1]*ene;
	  zz += chit->getCalorimeterHit()->getPosition()[2]*ene;
	  ee += ene;
	}		
	float xSP[3];
	xSP[0] = xx/ee;
	xSP[1] = yy/ee;
	xSP[2] = zz/ee;
	cluster->setStartingPoint(xSP);
      }
      
      float dist_to_SP(0.);
      if (_typeOfGenericDistance == 0) {
	//dist_to_SP = CaloHitTo->getDistanceToCalo();
	for (int i(0); i < 3; ++i) {
	  float xx = CaloHitTo->getCalorimeterHit()->getPosition()[i]-cluster->getStartingPoint()[i];
	  dist_to_SP += xx*xx;
	}
	dist_to_SP = sqrt(dist_to_SP);
      }
      else {
	dist_to_SP = dif_in_dist;
      }
      
      
      if (dist_to_SP < _distanceToDefineDirection ) {
	for (int i(0); i < 3; ++i)
	  xDir[i] =  CaloHitTo->getCalorimeterHit()->getPosition()[i];		
      }
      else {
	for (int i(0); i < 3; ++i) 
	  xDir[i] = CaloHitTo->getCalorimeterHit()->getPosition()[i] - cluster->getStartingPoint()[i];
      }
      
      
      
      CaloHitTo->setDirVec(xDir);
    }
    else { // Create new cluster
      
      if ( ihitTo==0 ) {
	
	ClusterExtended * cluster = new ClusterExtended(CaloHitTo);

	// debug
	if ( _debugLevel > 5 ) { 
	  std::cout << "create new cluster ... " << std::endl;
	  int color = static_cast<int>(16*reinterpret_cast<long int>(cluster)+0x009988);
	  ced_hit ( CaloHitTo->getCalorimeterHit()->getPosition()[0],CaloHitTo->getCalorimeterHit()->getPosition()[1],CaloHitTo->getCalorimeterHit()->getPosition()[2], 
		    0 | 1 << CED_LAYER_SHIFT, 8, color );
	  ced_send_event();
	  getchar();
	}
	

	CaloHitTo->setClusterExtended(cluster);
	CaloHitTo->setDistanceToNearestHit(0.0);
	float xDir[3];

	if ( (_typeOfGenericDistance == 0) ||  (_typeOfGenericDistance == 1) ) {

	  for (int i(0); i < 3; ++i) xDir[i] = CaloHitTo->getCalorimeterHit()->getPosition()[i];

	}
	else {

	  CalorimeterHitWithAttributes* calorimeterHitWithAttributesCorrespondingToCaloHitTo = getCalorimeterHitWithAttributes(CaloHitTo);
	  
	  float pathLengthOnHelixOfCaloHitTo = calorimeterHitWithAttributesCorrespondingToCaloHitTo->getPathLengthOnHelix();
	  float distanceToHelixOfCaloHitTo = calorimeterHitWithAttributesCorrespondingToCaloHitTo->getDistanceToHelix();

	  float distanceInHelixCoordinates = sqrt( pow( (pathLengthOnHelixOfCaloHitTo - _pathLengthOnHelixOfStartPoint),2) + 
						   pow( (distanceToHelixOfCaloHitTo - _distanceToHelixOfStartPoint),2) );

	  for (int i(0); i < 3; ++i) xDir[i] = distanceInHelixCoordinates*_startDirection.at(i);

	}
	
	CaloHitTo->setDirVec(xDir);	    
	_allClusters.push_back( cluster );


      }
      else {

	ClusterExtended * cluster = new ClusterExtended(CaloHitTo);

	// debug
	if ( _debugLevel > 5 ) { 
	  std::cout << "create new cluster ... " << std::endl;
	  int color = static_cast<int>(1024*reinterpret_cast<long int>(cluster)+0xFF9988);
	  ced_hit ( CaloHitTo->getCalorimeterHit()->getPosition()[0],CaloHitTo->getCalorimeterHit()->getPosition()[1],CaloHitTo->getCalorimeterHit()->getPosition()[2], 
		    0 | 1 << CED_LAYER_SHIFT, 8, color );
	  ced_send_event();
	  getchar();
	}



	CaloHitTo->setClusterExtended(cluster);
	CaloHitTo->setDistanceToNearestHit(0.0);
	float xDir[3];

	if ( (_typeOfGenericDistance == 0) ||  (_typeOfGenericDistance == 1) ) {

	  for (int i(0); i < 3; ++i) xDir[i] = CaloHitTo->getCalorimeterHit()->getPosition()[i];

	}
	else {

	  CalorimeterHitWithAttributes* calorimeterHitWithAttributesCorrespondingToCaloHitTo = getCalorimeterHitWithAttributes(CaloHitTo);
	  
	  float pathLengthOnHelixOfCaloHitTo = calorimeterHitWithAttributesCorrespondingToCaloHitTo->getPathLengthOnHelix();
	  float distanceToHelixOfCaloHitTo = calorimeterHitWithAttributesCorrespondingToCaloHitTo->getDistanceToHelix();

	  float distanceInHelixCoordinates = sqrt( pow( (pathLengthOnHelixOfCaloHitTo - _pathLengthOnHelixOfStartPoint),2) + 
						   pow( (distanceToHelixOfCaloHitTo - _distanceToHelixOfStartPoint),2) );

	  for (int i(0); i < 3; ++i) xDir[i] = distanceInHelixCoordinates*(CaloHitTo->getCalorimeterHit()->getPosition()[i]);

	}

	CaloHitTo->setDirVec(xDir);	    
	_allClusters.push_back( cluster );
	
      }

    }
        
  }

}



float TrackwiseClusters::DistanceBetweenPoints(float* x1, float* x2) {

  float xDistance(0.);
  for (int i(0); i < 3; i++) {
    float xx = x1[i] - x2[i];
    xDistance += xx*xx ;
  }
  xDistance = sqrt(xDistance);
  return xDistance;

}



void TrackwiseClusters::propertiesForAll() 
{
  int nclusters = (int)_allClusters.size();
  for (int i=0;i<nclusters;++i) {
    ClusterExtended * Cl = _allClusters[i];
    calculateProperties(Cl);
  }
  
}



void TrackwiseClusters::calculateProperties(ClusterExtended* Cl) {


  CaloHitExtendedVec calohitvec = Cl->getCaloHitExtendedVec();
  int nhcl = (int)calohitvec.size();
  if (nhcl > 0) {
    float * xhit = new float[nhcl];
    float * yhit = new float[nhcl];
    float * zhit = new float[nhcl];
    float * ahit = new float[nhcl];
    float * exhit = new float[nhcl];
    float * eyhit = new float[nhcl];
    float * ezhit = new float[nhcl];    
    float totene = 0.0;
    float totecal = 0.0;
    float tothcal = 0.0;
    RandomNumberGenerator random;
    float zmin = 1.0e+20;
    float zmax = -1.0e+20;
    int jhit = 0;
    for (int ihit(0); ihit < nhcl; ++ihit) {
      CalorimeterHit * calhit = 
	calohitvec[ihit]->getCalorimeterHit();
      if (calohitvec[ihit]->getDistanceToNearestHit() < 100.) {
	xhit[jhit] = calhit->getPosition()[0] + random.EqualDistribution(1.0)[0];
	yhit[jhit] = calhit->getPosition()[1] + random.EqualDistribution(1.0)[0];
	zhit[jhit] = calhit->getPosition()[2];
	ahit[jhit] = calhit->getEnergy();
	exhit[jhit] = 4.0;
	eyhit[jhit] = 4.0;
	ezhit[jhit] = 4.0;
	totene += ahit[jhit];
	if (calohitvec[jhit]->getType() == 0) {
	  totecal += ahit[jhit];
	}
	else {
	  tothcal += ahit[jhit];
	}	
	if (zhit[jhit]<zmin )
	  zmin = zhit[jhit];
	if (zhit[jhit]>zmax)
	  zmax = zhit[jhit];	
	jhit++;
      }	      
    }
    float zBeg;
    float signPz = 1.0;
    if (fabs(zmin)>fabs(zmax)) {
      signPz = -1.0;
      zBeg = zmax;
    }
    else {
      zBeg = zmin;
    }


    ClusterShapes * shapes 
      = new ClusterShapes(jhit,ahit,xhit,yhit,zhit);	    
    shapes->setErrors(exhit,eyhit,ezhit);
    float par[5];
    float dpar[5];
    float chi2 = 1.0e+10;
    float distmax = 1.0e+20;
    float x0 = 1;
    float y0 = 1;
    float r0 = 1;
    float bz = 1;
    float phi0 = 1;
    if (jhit > 3) {
      shapes->FitHelix(500, 0, 1, par, dpar, chi2, distmax);
      x0 = par[0];
      y0 = par[1];
      r0 = par[2];
      bz = par[3];
      phi0 = par[4];
    }
    HelixClass helix;
    helix.Initialize_BZ(x0, y0, r0, 
			bz, phi0, _bField, signPz,
			zBeg);
    
    Cl->setHelix(helix);
    Cl->setHelixChi2R(chi2);
    Cl->setHelixChi2Z(distmax);

    float axis[3];
    float pos[3];
    float axisMod = 0.0;
    //float low[3];
    //float up[3];
    for (int i=0; i<3; ++i) {
      pos[i]  = shapes->getCentreOfGravity()[i];
      axis[i] = shapes->getEigenVecInertia()[i];
      axisMod += axis[i]*axis[i];
      //low[i] = Cl->getLowEdge()[i];
      //up[i]  = Cl->getUpEdge()[i];
    }
    axisMod = sqrt(axisMod);
    Cl->setAxis(axis);
    Cl->setPosition(pos);	    
    Cl->setEccentricity(shapes->getElipsoid_eccentricity());

//    Cl->setLowEdge(low);
//    Cl->setUpEdge(up);
    if (nhcl > 40  ) {      
      //debug
      /*
	std::cout << nhcl << " " << shapes->getElipsoid_eccentricity() << std::endl;
	std::cout << "Chi2 : " << chi2 << " DistMax = " << distmax << std::endl;
	std::cout << "Pos : " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	std::cout << "Low : " << low[0] << " " << low[1] << " " << low[2] << std::endl;
	std::cout << "Up : " << up[0] << " " << up[1] << " " << up[2] << std::endl;     
	std::cout << "Par : " << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] << std::endl;
      */

       //float xx[3];
       //int nn = nhcl / 2;
       // xx[0] = xhit[nn];
       // xx[1] = yhit[nn];
       // xx[2] = zhit[nn];
       // float dd[3];
       // float time = helix.getDistanceToPoint(xx,dd);
       // float dist1 = dd[2];
    }
    
    delete shapes;
    delete[] exhit;
    delete[] eyhit;
    delete[] ezhit;
    delete[] xhit;
    delete[] yhit;
    delete[] zhit;
    delete[] ahit;
    
  }
  
}



void TrackwiseClusters::mergeForward() {
  

  int nClusters = (int)_allClusters.size();
  int nTotHits = (int)_allHits.size();
  
  int iCluster = 0;

  while (iCluster < nClusters) {
    ClusterExtended * clusterAR = _allClusters[iCluster];
    CaloHitExtendedVec hitvec = clusterAR->getCaloHitExtendedVec();
    int nHits = (int)hitvec.size();
    int iforw(0);
    int iback(0);
    CaloHitExtended * calohitAttachTo = 0  ;
    if (nHits > _nhit_minimal && nHits < _nhit_merge_forward) {
      //      std::cout << "attempt to merge forward" << std::endl;
      //      for (int i=1; i<nHits; ++i) {
      //      	float vec[3];
      //	for (int j=0;j<3;++j) 
      //	  vec[j] = hitvec[i]->getCalorimeterHit()->getPosition()[j]-
      //	    hitvec[0]->getCalorimeterHit()->getPosition()[j];
      //	hitvec[i]->setDirVec(vec);
      //      }
      int LowerBound = std::min(_nScanToMergeForward,nHits);
      for (int iCounterHit=0; iCounterHit < LowerBound; ++iCounterHit) {
	CaloHitExtended * calohit = hitvec[nHits-iCounterHit-1];
	int index = calohit->getIndex() + 1;
	int type = calohit->getType();
	float distance = 0.0;
	//	std::cout << " " << index << " " << nTotHits << std::endl; 
	while (distance < _distanceMergeForward[type] && index < nTotHits) {
	  CaloHitExtended * calohitTo = _allHits[index];	
	  distance = calohitTo->getGenericDistance() - calohit->getGenericDistance();
	  ClusterExtended * cluster_dummy = calohitTo->getClusterExtended();
	  CaloHitExtendedVec dummy = cluster_dummy->getCaloHitExtendedVec();
	  float yres = findResolutionParameter(calohit, calohitTo);
	  bool considerHit = yres < 2.0*_resolutionParameter[type]; 
	  //      int ndummy = (int)dummy.size();
	  //	  considerHit = considerHit && (ndummy > _nhit_minimal);
	  considerHit = considerHit && (cluster_dummy != clusterAR);
	  if (considerHit) {
	    calohitAttachTo = calohitTo;
	    calohit->setCaloHitTo(calohitTo);
	    calohit->setYresTo(yres);
	    iforw = 1;
	  }
	  if (iforw == 1)
	    break;
	  index++;
	}
	if (iforw == 1)
	  break;
      }
      if (iforw == 0) {
	for (int iCounterHit=0; iCounterHit < LowerBound; ++iCounterHit) {
	  CaloHitExtended * calohit = hitvec[nHits-iCounterHit-1];
	  int index = calohit->getIndex() - 1;
	  int type = calohit->getType();
	  float distance = 0.0;
	  //	std::cout << " " << index << " " << nTotHits << std::endl; 
	  while (distance < _stepTrackBack[type] && index >= 0) {
	    CaloHitExtended * calohitTo = _allHits[index];	
	    distance = - calohitTo->getGenericDistance() + calohit->getGenericDistance();
	    ClusterExtended * cluster_dummy = calohitTo->getClusterExtended();
	    CaloHitExtendedVec dummy = cluster_dummy->getCaloHitExtendedVec();
	    float yres = findResolutionParameter(calohitTo, calohit);
	    bool considerHit = yres < 2.0*_resolutionParameter[type]; 
	    //      int ndummy = (int)dummy.size();
	    //	  considerHit = considerHit && (ndummy > _nhit_minimal);
	    considerHit = considerHit && (cluster_dummy != clusterAR);
	    if (considerHit) {
	      calohitAttachTo = calohitTo;
	      calohit->setCaloHitTo(calohitTo);
	      calohit->setYresTo(yres);
	      iback = 1;
	    }
	    if (iback == 1)
	      break;
	    index--;
	  }
	  if (iback == 1)
	    break;
	}
      }
      if (iforw == 1) {
	//	std::cout << "Merging forward " << std::endl;
	ClusterExtended * clusterTo = calohitAttachTo->getClusterExtended();
	CaloHitExtendedVec hitvecAttached = clusterTo->getCaloHitExtendedVec();
	int nHitsAttached = (int)hitvecAttached.size();
	for (int jHit = 0; jHit < nHitsAttached; ++jHit) {
	  CaloHitExtended * hitToAttach = hitvecAttached[jHit];
	  clusterAR->addCaloHitExtended( hitToAttach );
	  hitToAttach->setClusterExtended( clusterAR );
	}
      clusterTo->Clear();
      }
      if (iback == 1) {
	//     std::cout << "Merging backward " << std::endl;
	ClusterExtended * clusterTo = calohitAttachTo->getClusterExtended();
	CaloHitExtendedVec hitvecAttached = clusterTo->getCaloHitExtendedVec();
	int nHitsAttached = (int)hitvecAttached.size();
	clusterAR->Clear();
	for (int jHit = 0; jHit < nHitsAttached; ++jHit) {
	  CaloHitExtended * hitToAttach = hitvecAttached[jHit];
	  clusterAR->addCaloHitExtended( hitToAttach );
	  hitToAttach->setClusterExtended( clusterAR );
	}
	for (int jHit =0; jHit < nHits; ++jHit) {
	  CaloHitExtended * hitToAttach = hitvec[jHit];
	  clusterAR->addCaloHitExtended( hitToAttach );
	  hitToAttach->setClusterExtended( clusterAR );
	}
      clusterTo->Clear();
      }
      if (iforw == 0 && iback == 0) {
	iCluster++;
      }
    }
    else {
      iCluster++;
    }
  }

}



void TrackwiseClusters::DisplayClusters(ClusterExtendedVec clusterVec) {


  std::cout << " " << std::endl;
  std::cout << "Debug Output" << std::endl;
  int nhits  = _allHits.size();
  std::cout << "Total number of Hits : " << nhits << std::endl;
  int nclust = clusterVec.size();
  std::cout << "Number of clusters : " << nclust << std::endl;
  int ntot(0);
  int ninbig(0);
  
  

  float totenergy = 0.0;
  float energyinbig = 0.0;
  for (int iclust(0); iclust < nclust; ++iclust) {
    ClusterExtended * Cl = clusterVec[iclust];
    CaloHitExtendedVec calohitvec = Cl->getCaloHitExtendedVec();
    int nhcl = calohitvec.size();
    ntot += nhcl;
    float ene=0.0;
    for (int i=0; i<nhcl; ++i) {
      ene += calohitvec[i]->getCalorimeterHit()->getEnergy();
    }
    totenergy += ene;
    if (nhcl > _nhit_minimal) {
      std::cout << "Cluster  " << iclust << " Number of hits = " << nhcl << std::endl;
      ninbig += nhcl;
      energyinbig += ene;
    }
    
  }
  
  float fraction = energyinbig/totenergy;
  std::cout << std::endl;
  std::cout << "Hits in big clusters     : " << ninbig << std::endl;
  std::cout << "Total energy : " << totenergy << std::endl;
  std::cout << "Fraction in big clusters : " << fraction << std::endl;
  std::cout << "Sumcheck: number of hits : " <<  ntot << std::endl << std::endl;
  
}


float TrackwiseClusters::findResolutionParameter(CaloHitExtended* fromHit, CaloHitExtended* toHit) {
  
  
  float xdistvec[3];
  float dirvec[3];
  float xdist(0.);
  float product(0.);
  float dir(0.);

  if ( (_typeOfGenericDistance == 0) ||  (_typeOfGenericDistance == 1) ) { 

    for (int i(0); i < 3; i++) {
      xdistvec[i] = toHit->getCalorimeterHit()->getPosition()[i] - fromHit->getCalorimeterHit()->getPosition()[i];
      xdist += xdistvec[i]*xdistvec[i];
      dirvec[i] = fromHit->getDirVec()[i];
      dir += dirvec[i]*dirvec[i];
      product += xdistvec[i]*dirvec[i]; 
    }
  }      
  else { 

    /*     
    CalorimeterHitWithAttributes* calorimeterHitWithAttributesCorrespondingToCaloHitTo = getCalorimeterHitWithAttributes(toHit);
    CalorimeterHitWithAttributes* calorimeterHitWithAttributesCorrespondingToCaloHitFrom = getCalorimeterHitWithAttributes(fromHit);
    
    float pathLengthOnHelixOfCaloHitTo = calorimeterHitWithAttributesCorrespondingToCaloHitTo->getPathLengthOnHelix();
    float pathLengthOnHelixOfCaloHitFrom = calorimeterHitWithAttributesCorrespondingToCaloHitFrom->getPathLengthOnHelix();
    
    float distanceToHelixOfCaloHitTo = calorimeterHitWithAttributesCorrespondingToCaloHitTo->getDistanceToHelix();
    float distanceToHelixOfCaloHitFrom = calorimeterHitWithAttributesCorrespondingToCaloHitFrom->getDistanceToHelix();


    xdist = pow( (pathLengthOnHelixOfCaloHitTo - pathLengthOnHelixOfCaloHitFrom),2) + pow( (distanceToHelixOfCaloHitTo - distanceToHelixOfCaloHitFrom),2);
    */

    for (int i(0); i < 3; i++) {
      xdistvec[i] = toHit->getCalorimeterHit()->getPosition()[i] - fromHit->getCalorimeterHit()->getPosition()[i];     
      xdist += xdistvec[i]*xdistvec[i];
      dirvec[i] = fromHit->getDirVec()[i];
      dir += dirvec[i]*dirvec[i];
      product += xdistvec[i]*dirvec[i]; 
    }



  }
  
  xdist = sqrt(xdist);
  dir = sqrt(dir);
  product=product/fmax(1.0e-6,(xdist*dir));
  
  if (product > 1.) {
    product = 0.999999;
  }
  if (product < -1.) {
    product = -0.999999;
  }
  
  float angle = acos(product);
  
  return xdist*angle;

}



std::vector<ClusterImpl*> TrackwiseClusters::CreateClusterCollection(ClusterExtendedVec clusterVec) {
  

  std::vector<ClusterImpl*> resultingClusters;


  int nclust = (int)clusterVec.size();
  
  for (int iclust(0); iclust < nclust; ++iclust) {
    ClusterExtended * Cl = clusterVec[iclust];
    CaloHitExtendedVec calohitvec = Cl->getCaloHitExtendedVec();
    int nhcl = (int)calohitvec.size();
    if (nhcl > _nhit_minimal) {
      ClusterImpl * cluster = new ClusterImpl();
      float * xhit = new float[nhcl];
      float * yhit = new float[nhcl];
      float * zhit = new float[nhcl];
      float * ahit = new float[nhcl];
      float totene = 0.0;
      float totecal = 0.0;
      float tothcal = 0.0;
      for (int ihit(0); ihit < nhcl; ++ihit) {
	CalorimeterHit * calhit = 
	  calohitvec[ihit]->getCalorimeterHit();
	cluster->addHit(calhit,(float)1.0);
	xhit[ihit] = calhit->getPosition()[0];
	yhit[ihit] = calhit->getPosition()[1];
	zhit[ihit] = calhit->getPosition()[2];
	ahit[ihit] = calhit->getEnergy();
	totene += ahit[ihit];
	if (calohitvec[ihit]->getType() == 0) {
	  totecal += ahit[ihit];
	}
	else {
	  tothcal += ahit[ihit];
	}		
      }
      
      ClusterShapes * shape 
	= new ClusterShapes(nhcl,ahit,xhit,yhit,zhit);	    
      cluster->setEnergy(totene);
      cluster->subdetectorEnergies().resize(2);
      cluster->subdetectorEnergies()[0] = totecal;
      cluster->subdetectorEnergies()[1] = tothcal;	    
      cluster->setPosition(shape->getCentreOfGravity());
      float PhiCluster = atan2(shape->getEigenVecInertia()[1],shape->getEigenVecInertia()[0]);
      float ThetaCluster = acos(shape->getEigenVecInertia()[2]);
      cluster->setIPhi(PhiCluster);
      cluster->setITheta(ThetaCluster);	    
      resultingClusters.push_back(cluster);

      delete shape;
      delete[] xhit;
      delete[] yhit;
      delete[] zhit;
      delete[] ahit;

    }

  }

  return resultingClusters;

}


void TrackwiseClusters::CleanUp() {

  
  for (unsigned int i(0); i < _allHits.size(); ++i) {
    CaloHitExtended * calohit = _allHits.at(i);
    delete calohit;
  }
    
  _allHits.clear();

  
  for (unsigned int i(0); i < _allSuperClusters.size(); ++i) {
    ClusterExtended * cluster = _allSuperClusters.at(i);
    delete cluster;
  }
  
  _allSuperClusters.clear();

  for (unsigned int i(0); i < _allClusters.size(); ++i) {
    ClusterExtended * cluster = _allClusters.at(i);
    delete cluster;
  }
  
  _allClusters.clear();

}



CalorimeterHitWithAttributes* TrackwiseClusters::getCalorimeterHitWithAttributes(CaloHitExtended* calorimeterHitExtended) {


  CalorimeterHitWithAttributes* correspondingCalorimeterHitWithAttributes = 0;

  for(std::vector<CalorimeterHitWithAttributes*>::const_iterator i = _calorimeterHitsWithAttributes.begin(); i != _calorimeterHitsWithAttributes.end(); ++i) {

    if ( ((*i)->getCalorimeterHit()) == (calorimeterHitExtended->getCalorimeterHit()) ) {

      correspondingCalorimeterHitWithAttributes = (*i);
      return correspondingCalorimeterHitWithAttributes;

    }

  }

  std::cout << "WARNING: Corresponding CalorimeterHitWithAttributes not found in TrackwiseClusters::getCalorimeterHitWithAttributes(), returning zero pointer" << std::endl;
  return correspondingCalorimeterHitWithAttributes;

}
