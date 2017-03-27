#ifndef TRACKHITPAIR_H
#define TRACKHITPAIR_H 1

#include "TrackExtended.h"
#include "TrackerHitExtended.h"
#include <vector>

class TrackHitPair;

typedef std::vector<TrackHitPair*> TrackHitPairVec;
/**
 * Class implementing association of TrackExtended and TrackerHitExtended objects. <br>
 * @author A. Raspereza (MPI-Munich)<br>
 */

class TrackHitPair {

 public:

    TrackHitPair(TrackExtended * trkExt, TrackerHitExtended * hitExt, float distance);
    ~TrackHitPair();

    //objects do not own pointers, we can use defaults, or the do own and things
    //leak like mad because the dtor is not correctly implemented
    TrackHitPair(TrackHitPair& trackhitpair) = default;
    TrackHitPair& operator=(TrackHitPair& trackhitpair) = default;

    void setTrackExtended(TrackExtended * trkExt);
    void setTrackerHitExtended(TrackerHitExtended * hitExt);
    void setDistance(float distance);
    TrackExtended * getTrackExtended();
    TrackerHitExtended * getTrackerHitExtended();
    float getDistance();
    

 private:
    TrackExtended * _trackExtended=NULL;
    TrackerHitExtended * _trackerHitExtended=NULL;
    float _distance=0.0;
    

};

#endif
