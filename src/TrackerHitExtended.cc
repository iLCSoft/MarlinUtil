#include "TrackerHitExtended.h"

TrackerHitExtended::TrackerHitExtended(TrackerHit * trackerhit) {
    _trackerHit = trackerhit;
}

TrackerHitExtended::~TrackerHitExtended() {

}

void TrackerHitExtended::setTrackExtended(TrackExtended * trackAR) {
    _trackAR = trackAR;
} 

void TrackerHitExtended::setTrackerHitTo(TrackerHitExtended * hitTo) {
    _hitTo = hitTo;
}

void TrackerHitExtended::setTrackerHitFrom(TrackerHitExtended * hitFrom) {
    _hitFrom = hitFrom;

} 

void TrackerHitExtended::setGenericDistance(float genericDistance) {
    _genericDistance = genericDistance; 
}

void TrackerHitExtended::setTrackerHit(TrackerHit * hit) {
    _trackerHit = hit;
}

void TrackerHitExtended::setYresTo(float yresTo) {
    _yresTo = yresTo;
}

void TrackerHitExtended::setYresFrom(float yresFrom) {
    _yresFrom = yresFrom;
}

void TrackerHitExtended::setDirVec(float * dirVec) {
    _dirVec[0] = dirVec[0];
    _dirVec[1] = dirVec[1];
    _dirVec[2] = dirVec[2];
}

TrackerHit * TrackerHitExtended::getTrackerHit() {
    return _trackerHit;
}

TrackExtended * TrackerHitExtended::getTrackExtended() {
    return _trackAR;
}

TrackerHitExtended * TrackerHitExtended::getTrackerHitFrom() {
    return _hitFrom;
}
TrackerHitExtended * TrackerHitExtended::getTrackerHitTo() {
    return _hitTo;
}
float TrackerHitExtended::getGenericDistance() {
    return _genericDistance;
}
float TrackerHitExtended::getYresTo() {
    return _yresTo;
}
float TrackerHitExtended::getYresFrom() {
    return _yresFrom;
}

float * TrackerHitExtended::getDirVec() {
    return _dirVec;
}
