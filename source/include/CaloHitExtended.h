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
 * @version $Id: CaloHitExtended.h,v 1.4 2006-02-22 14:53:27 owendt Exp $<br>
 */
class CaloHitExtended {

    public:

    CaloHitExtended(CalorimeterHit* calhit, int type);
    
    ~CaloHitExtended();

    //objects do not own pointers, we can use defaults, or they do own and
    //things leak like mad because the dtor is not correctly implemented
    CaloHitExtended(CaloHitExtended& calohitextended) = default;
    CaloHitExtended& operator=(CaloHitExtended& calohitextended) = default;

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
    void setDistanceToCalo(float distanceToCalo);
    float getDistanceToCalo();
    void setDistanceToNearestHit(float distanceToNearest);
    float getDistanceToNearestHit();
      

    private:
    
    CalorimeterHit * _calohit;
    int _type;
    CaloHitExtended * _calohitTo=NULL;
    CaloHitExtended * _calohitFrom=NULL;
    ClusterExtended * _clusterAR=NULL;
    int _index=0;
    float _dirVec[3];
    float _yresTo=0.0;
    float _yresFrom=0.0;
    float _genericDistance=0.0;
    float _caloDistance=0.0;
    float _distanceToNearestHit=0.0;

};

#endif
