#ifndef TRACKERHITAR_H
#define TRACKERHITAR_H 1

#include "lcio.h"
#include "EVENT/LCIO.h"
#include "EVENT/TrackerHit.h"
#include "TrackExtended.h"

using namespace lcio;

class TrackExtended;

class TrackerHitExtended {

 public:

    TrackerHitExtended(TrackerHit * trackerhit);
    ~TrackerHitExtended();
    void setTrackExtended(TrackExtended * trackAR);
    void setTrackerHitTo(TrackerHitExtended * hitTo);
    void setTrackerHitFrom(TrackerHitExtended * hitFrom);
    void setGenericDistance(float genericDistance);
    void setTrackerHit(TrackerHit * hit);
    void setYresTo(float yresTo);
    void setYresFrom(float yresFrom);
    void setDirVec(float * dirVec);

    TrackerHit * getTrackerHit();
    TrackExtended * getTrackExtended();
    TrackerHitExtended * getTrackerHitFrom();
    TrackerHitExtended * getTrackerHitTo();
    float getGenericDistance();
    float getYresTo();
    float getYresFrom();
    float * getDirVec();

 private:

    TrackerHit * _trackerHit;
    TrackExtended * _trackAR;
    TrackerHitExtended * _hitTo;
    TrackerHitExtended * _hitFrom;
    float _yresTo;
    float _yresFrom;
    float _genericDistance;
    float _dirVec[3];
	
};



#endif
