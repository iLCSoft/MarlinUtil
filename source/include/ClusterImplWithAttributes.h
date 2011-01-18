 
#ifndef CLUSTERIMPLWITHATTRIBUTES_H
#define CLUSTERIMPLWITHATTRIBUTES_H 1


#include <vector>
#include <IMPL/ClusterImpl.h>
#include "LCGeometryTypes.h"


// FIXME: put the implementations into a *.cc file 060906 OW

 
class ClusterImplWithAttributes { 

 public:
  
  ClusterImplWithAttributes() { 

    _clusterImpl = new ClusterImpl();

    _pos = new double[3];

    _posEndHit = new double[3];
    _posStartHit = new double[3];

    _directionEndHit = new double[3];
    _directionStartHit = new double[3];
  };
  
  ~ClusterImplWithAttributes() {

    delete _clusterImpl;
    _clusterImpl = 0;

    delete[] _pos;
    _pos = 0;

    delete[] _posEndHit;
    _posEndHit = 0;
    delete[] _posStartHit;
    _posStartHit = 0;

    delete[] _directionEndHit;
    _directionEndHit = 0;
    delete[] _directionStartHit;
    _directionStartHit = 0;
  };
  
  
  void setPositionEndHit(double* posEndHit) { for (unsigned int i = 0; i < 3; ++i) _posEndHit[i] = posEndHit[i]; };
  void setPositionEndHitVec(std::vector<double> posEndHit) { for (unsigned int i = 0; i < 3; ++i) _posEndHit[i] = posEndHit.at(i); };
  void setPositionEndHitLCVec(LCVector3D posEndHit) { for (unsigned int i = 0; i < 3; ++i) _posEndHit[i] = posEndHit[i]; };

  void setPositionStartHit(double* posStartHit) { for (unsigned int i = 0; i < 3; ++i) _posStartHit[i] = posStartHit[i]; };
  void setPositionStartHitVec(std::vector<double> posStartHit) { for (unsigned int i = 0; i < 3; ++i) _posStartHit[i] = posStartHit.at(i); };
  void setPositionStartHitLCVec(LCVector3D posStartHit) { for (unsigned int i = 0; i < 3; ++i) _posStartHit[i] = posStartHit[i]; };
  
  void setPosition(double* position) {    
    float* pos = new float[3];
    for (unsigned int i = 0; i < 3; ++i) pos[i] = (float)(position[i]);
    _clusterImpl->setPosition(pos);
    delete[] pos;
    pos = 0;

  };

  void setPositionVec(std::vector<double> position) {    
    float* pos = new float[3];
    for (unsigned int i = 0; i < 3; ++i) pos[i] = (float)position.at(i);
    _clusterImpl->setPosition(pos);
    delete[] pos;
    pos = 0;
    
  };

  void setPositionLCVec(LCVector3D position) {    
    float* pos = new float[3];
    for (unsigned int i = 0; i < 3; ++i) pos[i] = (float)position[i];
    _clusterImpl->setPosition(pos);
    delete[] pos;
    pos = 0;
    
  };

  void setEnergy(double energy) {_clusterImpl->setEnergy((float)energy); };
  
  void setDirectionEndHit(double* directionEndHit) { for (unsigned int i = 0; i < 3; ++i) _directionEndHit[i] = directionEndHit[i]; };
  void setDirectionEndHitVec(std::vector<double> directionEndHit) { for (unsigned int i = 0; i < 3; ++i) _directionEndHit[i] = directionEndHit.at(i); };
  void setDirectionEndHitLCVec(LCVector3D directionEndHit) { for (unsigned int i = 0; i < 3; ++i) _directionEndHit[i] = directionEndHit[i]; };

  void setDirectionStartHit(double* directionStartHit) { for (unsigned int i = 0; i < 3; ++i) _directionStartHit[i] = directionStartHit[i]; };
  void setDirectionStartHitVec(std::vector<double> directionStartHit) { for (unsigned int i = 0; i < 3; ++i) _directionStartHit[i] = directionStartHit.at(i); };
  void setDirectionStartHitLCVec(LCVector3D directionStartHit) { for (unsigned int i = 0; i < 3; ++i) _directionStartHit[i] = directionStartHit[i]; };

  void setTypeEndHit(int typeEndHit) { _typeEndHit = typeEndHit; };
  void setTypeStartHit(int typeStartHit) { _typeStartHit = typeStartHit; };

  void setIsMIPStub(bool isMIPStub) { _isMIPStub = isMIPStub; };
  void setIsMuon(bool isMuon) { _isMuon = isMuon; };
  
  void setClusterImpl(ClusterImpl* clusterImpl) { _clusterImpl = clusterImpl; };

  
  const double* getPositionEndHit() { return _posEndHit; };
  const std::vector<double> getPositionEndHitVec() { 
    std::vector<double> posEndHit;
    for (unsigned int i = 0; i < 3; ++i) posEndHit.push_back(_posEndHit[i]);
    return posEndHit;
  };
  const LCVector3D getPositionEndHitLCVec() { 
    LCVector3D posEndHit(_posEndHit[0],_posEndHit[1],_posEndHit[2]);
    return posEndHit;
  }; 

  const double* getPositionStartHit() { return _posStartHit; };
  const std::vector<double> getPositionStartHitVec() { 
    std::vector<double> posStartHit;
    for (unsigned int i = 0; i < 3; ++i) posStartHit.push_back(_posStartHit[i]);
    return posStartHit;
  };
  const LCVector3D getPositionStartHitLCVec() { 
    LCVector3D posStartHit(_posStartHit[0],_posStartHit[1],_posStartHit[2]);
    return posStartHit;
  }; 

  const double* getPosition() {     
    for (unsigned int i = 0; i < 3; ++i) _pos[i] = (double)(_clusterImpl->getPosition()[i]);    
    return _pos; 
  };
  const std::vector<double> getPositionVec() {
    std::vector<double> position;
    for (unsigned int i = 0; i < 3; ++i) position.push_back((double)(_clusterImpl->getPosition()[i]));
    return position;
  };
  const LCVector3D getPositionLCVec() {
    LCVector3D position((double)(_clusterImpl->getPosition()[0]),(double)(_clusterImpl->getPosition()[1]),(double)(_clusterImpl->getPosition()[2]));
    return position;
  };

  float getEnergy() {return _clusterImpl->getEnergy();};

  const double* getDirectionEndHit() { return _directionEndHit; };
  const std::vector<double> getDirectionEndHitVec() {
    std::vector<double> directionEndHit;
    for (unsigned int i = 0; i < 3; ++i) directionEndHit.push_back(_directionEndHit[i]);
    return directionEndHit;
  };
  const LCVector3D getDirectionEndHitLCVec() {
    LCVector3D directionEndHit(_directionEndHit[0],_directionEndHit[1],_directionEndHit[2]);
    return directionEndHit;
  };

  const double* getDirectionStartHit() { return _directionStartHit; };
  const std::vector<double> getDirectionStartHitVec() {
    std::vector<double> directionStartHit;
    for (unsigned int i = 0; i < 3; ++i) directionStartHit.push_back(_directionStartHit[i]);    
    return directionStartHit;
  };
  const LCVector3D getDirectionStartHitLCVec() {
    LCVector3D directionStartHit(_directionStartHit[0],_directionStartHit[1],_directionStartHit[2]);
    return directionStartHit;
  };

  int getTypeEndHit() { return _typeEndHit; };
  int getTypeStartHit() { return _typeStartHit; };

  bool isMIPStub() { return _isMIPStub; };
  bool isMuon() { return _isMuon; };
  
  ClusterImpl* getClusterImpl() { return _clusterImpl; };
  void addHit(CalorimeterHit* hit, double contribution) { _clusterImpl->addHit(hit,(float)contribution); };



 private:

  ClusterImpl* _clusterImpl;

  double* _pos;
  double* _posEndHit;
  double* _posStartHit;

  double* _directionEndHit;
  double* _directionStartHit;
  
  int _typeEndHit;
  int _typeStartHit;

  bool _isMIPStub;
  bool _isMuon;

} ;


#endif
