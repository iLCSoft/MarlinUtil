 
#ifndef CALORIMETERHITWITHATTRIBUTES_H
#define CALORIMETERHITWITHATTRIBUTES_H 1



#include <EVENT/CalorimeterHit.h>


 
class CalorimeterHitWithAttributes { 

 public:

  CalorimeterHitWithAttributes(CalorimeterHit* calorimeterHit, float distanceToHelix, float pathLengthOnHelix) {
  
    _calorimeterHit = calorimeterHit;
    _distanceToHelix = distanceToHelix;
    _pathLengthOnHelix = pathLengthOnHelix;
  
  };
  
  ~CalorimeterHitWithAttributes() {};
  
  CalorimeterHit* getCalorimeterHit() { return _calorimeterHit; };
  float getDistanceToHelix() { return _distanceToHelix; };
  float getPathLengthOnHelix() { return _pathLengthOnHelix; };
  
  void setCalorimeterHit(CalorimeterHit* calorimeterHit) { _calorimeterHit = calorimeterHit; };
  void setDistanceToHelix(float distanceToHelix) { _distanceToHelix = distanceToHelix; };
  void setPathLengthOnHelix(float pathLengthOnHelix) { _pathLengthOnHelix = pathLengthOnHelix; };
  
  
 private:

  CalorimeterHit* _calorimeterHit;
  float _distanceToHelix;
  float _pathLengthOnHelix;

} ;


#endif
