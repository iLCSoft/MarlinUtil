#include "TrackExtended.h"
#include <math.h>

TrackExtended::TrackExtended( ) {
    _track = NULL;
    _cluster = NULL;
    _trackerHitVector.clear();
}

TrackExtended::TrackExtended( Track * track) {
    _track = track;
    _cluster = NULL;
    _trackerHitVector.clear();
}

TrackExtended::TrackExtended( TrackerHitExtended * trackerhit) {
    _trackerHitVector.clear();
    _trackerHitVector.push_back(trackerhit);
    _track = NULL;
    _cluster = NULL;
}


TrackExtended::~TrackExtended() {}

Track * TrackExtended::getTrack() {
    return _track;
}

const float * TrackExtended::getSeedPosition() {
    return _seedPosition;
}

const float * TrackExtended::getSeedDirection() {
    return _seedDirection;
}

ClusterExtended * TrackExtended::getCluster() {
    return _cluster;
}

ClusterExtended * TrackExtended::getSuperCluster() {
    return _superCluster;
}

TrackerHitExtendedVec & TrackExtended::getTrackerHitExtendedVec() {

    return _trackerHitVector;
}

void TrackExtended::setCluster(ClusterExtended * cluster) {
    _cluster = cluster;
}

void TrackExtended::setSuperCluster(ClusterExtended * superCluster) {
    _superCluster = superCluster;
}

void TrackExtended::setSeedDirection( float * seedDirection ) {
    _seedDirection[0] = seedDirection[0];
    _seedDirection[1] = seedDirection[1];
    _seedDirection[2] = seedDirection[2];
}

void TrackExtended::setSeedPosition( float * seedPosition) {
    _seedPosition[0] = seedPosition[0];
    _seedPosition[1] = seedPosition[1];
    _seedPosition[2] = seedPosition[2];
}

void TrackExtended::addTrackerHitExtended( TrackerHitExtended * trackerhit) {
  _trackerHitVector.push_back(trackerhit);

}

void TrackExtended::ClearTrackerHitExtendedVec() {
  _trackerHitVector.clear();
}

void TrackExtended::setX0( float x0 ) {
  _x0 = x0;
}

void TrackExtended::setY0( float y0 ) {
  _y0 = y0;
}

void TrackExtended::setR0( float r0 ) {
  _r0 = r0;
}

void TrackExtended::setBz( float bz ) {
  _bz = bz;
}

void TrackExtended::setPhi0( float phi0 ) {
  _phi0 = phi0;
}

float TrackExtended::getX0() {
  return _x0;
}

float TrackExtended::getY0() {
  return _y0;
}

float TrackExtended::getR0() {
  return _r0;
}

float TrackExtended::getBz() {
  return _bz;
}

float TrackExtended::getPhi0() {
  return _phi0;
}

void TrackExtended::setStart(float * xStart) {
  _xStart[0] = xStart[0];
  _xStart[1] = xStart[1];
  _xStart[2] = xStart[2];
}

void TrackExtended::setEnd(float * xEnd) {
  _xEnd[0] = xEnd[0];
  _xEnd[1] = xEnd[1];
  _xEnd[2] = xEnd[2];
}

float * TrackExtended::getStart() {
  return _xStart;
}

float * TrackExtended::getEnd() {
  return _xEnd;
}

void TrackExtended::setGroupTracks( GroupTracks * group ) {
  _group = group;
}

GroupTracks * TrackExtended::getGroupTracks() {
  return _group;
}

