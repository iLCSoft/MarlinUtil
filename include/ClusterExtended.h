#ifndef ClusterExtended_H
#define ClusterExtended_H 1

#include "CaloHitExtended.h"
#include "TrackExtended.h"
#include "EVENT/Cluster.h"


using namespace lcio;

class TrackExtended;
typedef std::vector<TrackExtended*> TrackExtendedVec;

class CaloHitExtended;
typedef std::vector<CaloHitExtended*> CaloHitExtendedVec;

class ClusterExtended;
typedef std::vector<ClusterExtended*> ClusterExtendedVec;


class ClusterExtended {

 public:

    ClusterExtended();
    ClusterExtended( Cluster * cluster );
    ClusterExtended(CaloHitExtended * calohit);
    ClusterExtended(TrackExtended * track);

    ~ClusterExtended();
    
    CaloHitExtendedVec & getCaloHitExtendedVec();
    TrackExtendedVec & getTrackExtendedVec();
    const float* getStartingPoint();
    const float* getDirection();
    void setStartingPoint(float* sPoint);
    void setDirection(float* direct);
    void addCaloHitExtended(CaloHitExtended * calohit);
    void addTrackExtended(TrackExtended * track);
    void setType( int type );
    int getType();

    void Clear();

    void setCluster(Cluster * cluster);
    Cluster * getCluster();

 private:

    TrackExtendedVec _trackVector;
    CaloHitExtendedVec _hitVector;
    float _startingPoint[3];
    float _direction[3];

    int _type;
    Cluster * _cluster;


};

#endif
