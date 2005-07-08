#ifndef TRACKAR_H 
#define TRACKAR_H 1

#include "lcio.h"
#include "EVENT/LCIO.h"
#include "EVENT/Track.h"
#include <vector>
#include "ClusterExtended.h"
#include "TrackerHitExtended.h"
#include "GroupTracks.h"

using namespace lcio;

class ClusterExtended;
class GroupTracks;
class TrackerHitExtended;
typedef std::vector<TrackerHitExtended*> TrackerHitExtendedVec;

class TrackExtended {

 public:


    TrackExtended( );
    TrackExtended( TrackerHitExtended * trackerhit );
    TrackExtended( Track * track );
    ~TrackExtended();
    
    Track * getTrack();
    const float * getSeedDirection();
    const float * getSeedPosition();
    ClusterExtended * getCluster();
    ClusterExtended * getSuperCluster();
    TrackerHitExtendedVec & getTrackerHitExtendedVec();
    void setCluster(ClusterExtended * cluster);
    void setSuperCluster(ClusterExtended * superCluster);
    void setSeedDirection( float * seedDirection );
    void setSeedPosition( float * seedPosition);
    void addTrackerHitExtended( TrackerHitExtended * trackerhit);
    void ClearTrackerHitExtendedVec();

    void setX0(float x0);
    void setY0(float y0);
    void setR0(float r0);
    void setBz(float bz);
    void setPhi0(float phi0);

    void setStart(float * xStart);
    void setEnd(float * xEnd);


    float getX0();
    float getY0();
    float getR0();
    float getBz();
    float getPhi0();
    
    float * getStart();
    float * getEnd();

    void setGroupTracks( GroupTracks * group );
    GroupTracks * getGroupTracks();


 private:

    ClusterExtended *_superCluster;
    ClusterExtended *_cluster;
    GroupTracks * _group;
    Track * _track;
    float _seedDirection[3];
    float _seedPosition[3];
    TrackerHitExtendedVec _trackerHitVector;    

    float _x0;
    float _y0;
    float _r0;
    float _bz;
    float _phi0;

    float _xStart[3];
    float _xEnd[3];

};

#endif
