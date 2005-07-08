#include "ClusterExtended.h"
#include <math.h>

ClusterExtended::ClusterExtended() {
    _hitVector.clear();
    _trackVector.clear();

    _direction[0] = 0.;
    _direction[1] = 0.;
    _direction[2] = 0.;

    _startingPoint[0] = 0.;
    _startingPoint[1] = 0.;
    _startingPoint[2] = 0.;	

}

ClusterExtended::ClusterExtended( Cluster * cluster ) {
    _hitVector.clear();
    _trackVector.clear();

    _direction[0] = 0.;
    _direction[1] = 0.;
    _direction[2] = 0.;

    _startingPoint[0] = 0.;
    _startingPoint[1] = 0.;
    _startingPoint[2] = 0.;	

    _cluster = cluster;

}




ClusterExtended::ClusterExtended(CaloHitExtended * calohit) {

    _hitVector.clear();
    _hitVector.push_back(calohit);
    _trackVector.clear();

    float rad(0);

    _startingPoint[0] = calohit->getCalorimeterHit()->getPosition()[0];
    _startingPoint[1] = calohit->getCalorimeterHit()->getPosition()[1];
    _startingPoint[2] = calohit->getCalorimeterHit()->getPosition()[2];
    
    for (int i(0); i < 3; ++i) {
	rad += _startingPoint[i]*_startingPoint[i];
    }

    rad = sqrt(rad);

    _direction[0] = _startingPoint[0]/rad;
    _direction[1] = _startingPoint[1]/rad;
    _direction[2] = _startingPoint[2]/rad;
    
}

ClusterExtended::ClusterExtended(TrackExtended * track) {
    _hitVector.clear();
    _trackVector.clear();
    _trackVector.push_back(track);
    for (int i(0); i < 3; ++i) {
	_startingPoint[i] = track->getSeedPosition()[i];
	_direction[i] = track->getSeedDirection()[i];
    }
}


ClusterExtended::~ClusterExtended() {
   _hitVector.clear(); 
   _trackVector.clear();
}

CaloHitExtendedVec& ClusterExtended::getCaloHitExtendedVec() {
    return _hitVector;
}

TrackExtendedVec& ClusterExtended::getTrackExtendedVec() {
    return _trackVector;
}

const float* ClusterExtended::getStartingPoint() {
    return _startingPoint;
}

const float* ClusterExtended::getDirection() {
    return _direction;
}

void ClusterExtended::setStartingPoint(float* sPoint) {
    _startingPoint[0] = sPoint[0];
    _startingPoint[1] = sPoint[1];
    _startingPoint[2] = sPoint[2];
    
}

void ClusterExtended::setDirection(float* direct) {
    _direction[0] = direct[0];
    _direction[1] = direct[1];
    _direction[2] = direct[2];
}

void ClusterExtended::addCaloHitExtended(CaloHitExtended* calohit) {
    _hitVector.push_back(calohit);
}

void ClusterExtended::addTrackExtended(TrackExtended * track) {
    _trackVector.push_back(track);
}

void ClusterExtended::Clear() {
    _hitVector.clear();
    _trackVector.clear();

}

void ClusterExtended::setType( int type ) {
  _type = type;
}

int ClusterExtended::getType() {
  return _type;
}

void ClusterExtended::setCluster(Cluster * cluster) {
  _cluster = cluster;
}

Cluster * ClusterExtended::getCluster() {
  return _cluster;
}
