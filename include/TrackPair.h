#ifndef TRACKPAIR_H 
#define TRACKPAIR_H 1
#include "lcio.h"
#include "EVENT/LCIO.h"
#include "EVENT/Track.h"
#include <vector>

using namespace lcio;

class TrackPair;

typedef std::vector<TrackPair*> TrackPairVec;

class TrackPair {

 public:

  TrackPair( );
  ~TrackPair( );

  TrackVec & getTracks();
  Track * getFirstTrack();
  Track * getSecondTrack();
  float getDistance();
  float * getVertex();
  float * getMomentum();
  float getMass();
  int getCode();

  void setFirstTrack( Track * track );
  void setSecondTrack( Track * track );
  void setDistance(float distance);
  void setVertex( float * vertex );
  void setMomentum( float * momentum );
  void setMass( float mass);
  void setCode( int code );

 private:

  TrackVec  _trkVec;
  float _distance;
  float  _vertex[3];
  float  _momentum[3];
  float _mass;
  int _code;

};


#endif
