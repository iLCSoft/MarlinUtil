#ifndef CaloHitExtended_H
#define CaloHitExtended_H 1

#include "lcio.h"
#include "EVENT/LCIO.h"
#include "EVENT/CalorimeterHit.h"
#include "ClusterExtended.h"

using namespace lcio ;

class ClusterExtended;

/**
 * Class extending native LCIO class CalorimeterHit. <br>
 * Class CaloHitExtended is used in TrackwiseClustering <br>
 * and Wolf processors. <br>
 * @author A. Raspereza (DESY)<br>
 * @version $ld: $<br>
 */
class CaloHitExtended {

    public:

    CaloHitExtended(CalorimeterHit* calhit, int type);
    
    ~CaloHitExtended();

    CalorimeterHit * getCalorimeterHit();
    CaloHitExtended * getCaloHitTo();
    CaloHitExtended * getCaloHitFrom();
    ClusterExtended * getClusterExtended();
    int    getIndex();
    int    getType();
    const float* getDirVec();
    float  getYresTo();
    float  getYresFrom();
    float  getGenericDistance();

    void setCalorimeterHit(CalorimeterHit* calhit);
    void setCaloHitTo(CaloHitExtended* calhitto);
    void setCaloHitFrom(CaloHitExtended* calohitfrom);
    void setClusterExtended(ClusterExtended* cluster);
    void setIndex(int index);
    void setType(int type);
    void setDirVec(float* dirVec);
    void setYresTo(float yresto);
    void setYresFrom(float yresfrom);
    void setGenericDistance(float genericDistance);



    private:
    
    CalorimeterHit * _calohit;
    CaloHitExtended * _calohitTo;
    CaloHitExtended * _calohitFrom;
    ClusterExtended * _clusterAR;
    int _index;
    int _type;
    float _dirVec[3];
    float _yresTo;
    float _yresFrom;
    float _genericDistance;

};

#endif
