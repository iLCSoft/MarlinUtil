
// FIXME: all methods concerning positions (or more general, concerning return values of LCIO class methods in double precision) need to be overloaded as double precision methods as well (2006/06/19).

#include "HelixClass.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "ced_cli.h"


HelixClass::HelixClass() {
    _const_pi = acos(-1.0);
    _const_2pi = 2.0*_const_pi;
    _const_pi2 = 0.5*_const_pi;
    _FCT = 2.99792458E-4;
}

HelixClass::~HelixClass() {}

void HelixClass::Initialize_VP(float * pos, float * mom, float q, float B) {

  // FIXME: the initialisation is not consistent for all the parametrisations
    _referencePoint[0] = pos[0];
    _referencePoint[1] = pos[1];
    _referencePoint[2] = pos[2];
    _momentum[0] = mom[0];
    _momentum[1] = mom[1];
    _momentum[2] = mom[2];
    _charge = q;
    _bField = B;
    _pxy = sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
    _radius = _pxy / (_FCT*B);
    _omega = q/_radius;
    _tanLambda = mom[2]/_pxy;
    _phiMomRefPoint = atan2(mom[1],mom[0]);
    _xCentre = pos[0] + _radius*cos(_phiMomRefPoint-_const_pi2*q);
    _yCentre = pos[1] + _radius*sin(_phiMomRefPoint-_const_pi2*q);
    _phiRefPoint = atan2(pos[1]-_yCentre,pos[0]-_xCentre);
    _phiAtPCA = atan2(-_yCentre,-_xCentre);
    _phi0 = -_const_pi2*q + _phiAtPCA;
    if (_phi0 < -_const_pi)
	_phi0 += _const_2pi;
    if (_phi0 > _const_pi)
	_phi0 -= _const_2pi;
    _xAtPCA = _xCentre + _radius*cos(_phiAtPCA);
    _yAtPCA = _yCentre + _radius*sin(_phiAtPCA);
    _d0 = -_xAtPCA*sin(_phi0) + _yAtPCA*cos(_phi0);
    _pxAtPCA = _pxy*cos(_phi0);
    _pyAtPCA = _pxy*sin(_phi0);
    float deltaPhi = _phiRefPoint - _phiAtPCA;    
    float xCircles = -pos[2]*q/(_radius*_tanLambda) - deltaPhi;
    xCircles = xCircles/_const_2pi;
    int nCircles;
    int n1,n2;

    if (xCircles >= 0.) {
	n1 = int(xCircles);
	n2 = n1 + 1;
    }
    else {
	n1 = int(xCircles) - 1;
	n2 = n1 + 1;
    }
    
    if (fabs(n1-xCircles) < fabs(n2-xCircles)) {
	nCircles = n1;
    }
    else {
	nCircles = n2;
    }
    _z0 = pos[2] + _radius*_tanLambda*q*(deltaPhi + _const_2pi*nCircles);


    // debug
    /*
    std::cout << "Initialize_VP called : " << std::endl << std::endl;
    printParameters();
    */

}

void HelixClass::Initialize_Canonical(float phi0, float d0, float z0, 
				      float omega, float tanLambda, float B) {
    _omega = omega;
    _d0 = d0;
    _phi0 = phi0;
    _z0 = z0;
    _tanLambda = tanLambda;
    _charge = omega/fabs(omega);
    _radius = 1./fabs(omega);
    _xAtPCA = -_d0*sin(_phi0);
    _yAtPCA = _d0*cos(_phi0);    
    _referencePoint[0] = _xAtPCA;
    _referencePoint[1] = _yAtPCA;
    _referencePoint[2] = _z0;
    _pxy = _FCT*B*_radius;
    _momentum[0] = _pxy*cos(_phi0);
    _momentum[1] = _pxy*sin(_phi0);
    _momentum[2] = _tanLambda * _pxy;    
    _pxAtPCA = _momentum[0];
    _pyAtPCA = _momentum[1];
    _phiMomRefPoint = atan2(_momentum[1],_momentum[0]);
    _xCentre = _referencePoint[0] + 
      _radius*cos(_phi0-_const_pi2*_charge);
    _yCentre = _referencePoint[1] + 
      _radius*sin(_phi0-_const_pi2*_charge);
    _phiAtPCA = atan2(-_yCentre,-_xCentre);
    _phiRefPoint =  _phiAtPCA ;
    _bField = B;

    _bZ = _omega/_tanLambda;


    // debug
    /*
    std::cout << "Initialize_Canonical called: " << std::endl << std::endl;
    printParameters();  
    */

}


void HelixClass::Initialize_BZ(float xCentre, float yCentre, float radius, 
			       float bZ, float phi0, float B, float signPz,
			       float zBegin) {

  // FIXME: _bZ and phi0 are NOT initialised correctely here and within the other two init methods
  _radius = radius;
  _pxy = _FCT*B*_radius;
  _charge = -(bZ*signPz)/fabs(bZ*signPz);
  _momentum[2] = -_charge*_pxy/(bZ*_radius);
  _xCentre = xCentre;
  _yCentre = yCentre;
  _omega = _charge/radius;
  _phiAtPCA = atan2(-_yCentre,-_xCentre);
  _phi0 = -_const_pi2*_charge + _phiAtPCA; 
  if (_phi0 < -_const_pi)
    _phi0 += _const_2pi;
  if (_phi0 > _const_pi)
    _phi0 -= _const_2pi;
  _xAtPCA = _xCentre + _radius*cos(_phiAtPCA);
  _yAtPCA = _yCentre + _radius*sin(_phiAtPCA);
  _d0 = -_xAtPCA*sin(_phi0) + _yAtPCA*cos(_phi0);
  _pxAtPCA = _pxy*cos(_phi0);
  _pyAtPCA = _pxy*sin(_phi0);
  _referencePoint[2] = zBegin;
  _referencePoint[0] = xCentre + radius*cos(bZ*zBegin+phi0);
  _referencePoint[1] = yCentre + radius*sin(bZ*zBegin+phi0);
  _phiRefPoint = atan2(_referencePoint[1]-_yCentre,_referencePoint[0]-_xCentre);
  _phiMomRefPoint =  -_const_pi2*_charge + _phiRefPoint;
  _tanLambda = _momentum[2]/_pxy;
  _momentum[0] = _pxy*cos(_phiMomRefPoint);
  _momentum[1] = _pxy*sin(_phiMomRefPoint);
  
  float deltaPhi = _phiRefPoint - _phiAtPCA;    
  float xCircles = bZ*_referencePoint[2] - deltaPhi;
  xCircles = xCircles/_const_2pi;
  int nCircles;
  int n1,n2;

  if (xCircles >= 0.) {
    n1 = int(xCircles);
    n2 = n1 + 1;
  }
  else {
    n1 = int(xCircles) - 1;
    n2 = n1 + 1;
  }
  
  if (fabs(n1-xCircles) < fabs(n2-xCircles)) {
    nCircles = n1;
  }
  else {
    nCircles = n2;
  }  
  _z0 = _referencePoint[2] - (deltaPhi + _const_2pi*nCircles)/bZ;  
  _bField = B;

  _bZ = bZ;


  // debug
  /*
  std::cout << "Initialize_BZ called: " << std::endl << std::endl;
  printParameters();  
  */
  
}

const float * HelixClass::getMomentum() {
    return _momentum;
}
const float * HelixClass::getReferencePoint() {
    return _referencePoint;
}
float HelixClass::getPhi0() {
    return _phi0;
}
float HelixClass::getD0() {
    return _d0;
}
float HelixClass::getZ0() {
    return _z0;
}
float HelixClass::getOmega() {
    return _omega;
}
float HelixClass::getTanLambda() {
    return _tanLambda;
}
float HelixClass::getPXY() {
    return _pxy;
}
float HelixClass::getXC() {
  return _xCentre;
}

float HelixClass::getYC() {
  return _yCentre;
}

float HelixClass::getRadius() {
  return _radius;
}

float HelixClass::getBz() {
  return _bZ;
}

float HelixClass::getPhiZ() {
  return _phiZ;
}

float HelixClass::getCharge() {
    return _charge;
}

float HelixClass::getPointInXY(float x0, float y0, float ax, float ay, 
			      float * ref , float * point) {

  float time;

  float AA = sqrt(ax*ax+ay*ay);


  if (AA <= 0) {
    time = 1.0e+10; 
    return time;
  }


  float BB = ax*(x0-_xCentre) + ay*(y0-_yCentre);
  BB = BB / AA;

  float CC = (x0-_xCentre)*(x0-_xCentre) 
    + (y0-_yCentre)*(y0-_yCentre) - _radius*_radius;

  CC = CC / AA;

  float DET = BB*BB - CC;
  float tt1 = 0.;
  float tt2 = 0.;
  float xx1,xx2,yy1,yy2; 


  if (DET < 0 ) {
    time = 1.0e+10;
    point[0]=0.0;
    point[1]=0.0;
    point[2]=0.0;
    return time;
  }
  
  
  tt1 = - BB + sqrt(DET);
  tt2 = - BB - sqrt(DET);

  xx1 = x0 + tt1*ax;
  yy1 = y0 + tt1*ay;
  xx2 = x0 + tt2*ax;
  yy2 = y0 + tt2*ay;
  
  float phi1 = atan2(yy1-_yCentre,xx1-_xCentre);
  float phi2 = atan2(yy2-_yCentre,xx2-_xCentre);
  float phi0 = atan2(ref[1]-_yCentre,ref[0]-_xCentre);

  float dphi1 = phi1 - phi0;
  float dphi2 = phi2 - phi0;

  if (dphi1 < 0 && _charge < 0) {
    dphi1 = dphi1 + _const_2pi;
  }
  else if (dphi1 > 0 && _charge > 0) { 
    dphi1 = dphi1 - _const_2pi;
  }

  if (dphi2 < 0 && _charge < 0) {
    dphi2 = dphi2 + _const_2pi;
  }
  else if (dphi2 > 0 && _charge > 0) { 
    dphi2 = dphi2 - _const_2pi;
  }

  // Times
  tt1 = -_charge*dphi1*_radius/_pxy;
  tt2 = -_charge*dphi2*_radius/_pxy;

  if (tt1 < 0. )
    std::cout << "WARNING " << tt1 << std::endl;
  if (tt2 < 0. )
    std::cout << "WARNING " << tt2 << std::endl;
  

  if (tt1 < tt2) {
    point[0] = xx1;
    point[1] = yy1;
    time = tt1;
  }
  else {
    point[0] = xx2;
    point[1] = yy2;
    time = tt2;
  }

  point[2] = ref[2]+time*_momentum[2];

  

  return time;

}


float HelixClass::getPointOnCircle(float Radius, float * ref, float * point) {

  float distCenterToIP = sqrt(_xCentre*_xCentre + _yCentre*_yCentre);

  point[0] = 0.0;
  point[1] = 0.0;
  point[2] = 0.0;

  if ((distCenterToIP+_radius)<Radius) {
    float xx = 1.0e+15;
    return xx;
  }

  if ((_radius+Radius)<distCenterToIP) {
    float xx = 1.0e+15;
    return xx;
  }

  float phiCentre = atan2(_yCentre,_xCentre);
  float phiStar   = Radius*Radius + distCenterToIP*distCenterToIP 
                                    - _radius*_radius;

  phiStar = 0.5*phiStar/fmax(1.0e-20,Radius*distCenterToIP);
  
  if (phiStar > 1.0) 
    phiStar = 0.9999999;
  
  if (phiStar < -1.0)
    phiStar = -0.9999999;
  
  phiStar = acos(phiStar);

  float tt1,tt2,time;

  float xx1 = Radius*cos(phiCentre+phiStar);
  float yy1 = Radius*sin(phiCentre+phiStar);

  float xx2 = Radius*cos(phiCentre-phiStar);
  float yy2 = Radius*sin(phiCentre-phiStar);


  float phi1 = atan2(yy1-_yCentre,xx1-_xCentre);
  float phi2 = atan2(yy2-_yCentre,xx2-_xCentre);
  float phi0 = atan2(ref[1]-_yCentre,ref[0]-_xCentre);

  float dphi1 = phi1 - phi0;
  float dphi2 = phi2 - phi0;

  if (dphi1 < 0 && _charge < 0) {
    dphi1 = dphi1 + _const_2pi;
  }
  else if (dphi1 > 0 && _charge > 0) { 
    dphi1 = dphi1 - _const_2pi;
  }

  if (dphi2 < 0 && _charge < 0) {
    dphi2 = dphi2 + _const_2pi;
  }
  else if (dphi2 > 0 && _charge > 0) { 
    dphi2 = dphi2 - _const_2pi;
  }

  // Times
  tt1 = -_charge*dphi1*_radius/_pxy;
  tt2 = -_charge*dphi2*_radius/_pxy;

  if (tt1 < 0. )
    std::cout << "WARNING " << tt1 << std::endl;
  if (tt2 < 0. )
    std::cout << "WARNING " << tt2 << std::endl;
  

  float time2;
  if (tt1 < tt2) {
    point[0] = xx1;
    point[1] = yy1;
    point[3] = xx2;
    point[4] = yy2;
    time = tt1;
    time2 = tt2;
  }
  else {
    point[0] = xx2;
    point[1] = yy2;
    point[3] = xx1;
    point[4] = yy1;
    time = tt2;
    time2 = tt1;
  }

  point[2] = ref[2]+time*_momentum[2];
  point[5] = ref[2]+time2*_momentum[2];
  

  return time;

}


float HelixClass::getPointInZ(float zLine, float * ref, float * point) {

  float time = zLine - ref[2];

  if (_momentum[2] == 0.) {
    time = 1.0e+10;
    point[0] = 0.;
    point[1] = 0.;
    point[2] = 0.;
    return time;
  }

  time = time/_momentum[2];

  float phi0 = atan2(ref[1] - _yCentre , ref[0] - _xCentre);
  float phi = phi0 - _charge*_pxy*time/_radius;
  float xx = _xCentre + _radius*cos(phi);
  float yy = _yCentre + _radius*sin(phi);

  point[0] = xx;
  point[1] = yy;
  point[2] = zLine;

  return time;


}



float HelixClass::getDistanceToPoint(float* xPoint, float* Distance) {

  float pointOnHelix[3] = {0.0,0.0,0.0};
  
  return getDistanceToPoint(xPoint,Distance,pointOnHelix);

}


float HelixClass::getDistanceToPoint(const float* xPoint, float* Distance) {

  float pointOnHelix[3] = {0.0,0.0,0.0};

  float x[3] = {0.0,0.0,0.0};   
  for (unsigned int i = 0; i < 3; ++i) x[i] = xPoint[i];

  return getDistanceToPoint(x,Distance,pointOnHelix);

}


float HelixClass::getDistanceToPoint(const float* xPoint, float* Distance, float* pointOnHelix) {

  float x[3] = {0.0,0.0,0.0};   
  for (unsigned int i = 0; i < 3; ++i) x[i] = xPoint[i];

  return getDistanceToPoint(x,Distance,pointOnHelix);

}


float HelixClass::getDistanceToPoint(float* xPoint, float* Distance, float* pointOnHelix) {

  float xOnHelix, yOnHelix, zOnHelix;
  float phi = atan2(xPoint[1]-_yCentre,xPoint[0]-_xCentre);
  xOnHelix = _xCentre + _radius*cos(phi);
  yOnHelix = _yCentre + _radius*sin(phi);
  float phi0 = atan2(_referencePoint[1]-_yCentre,_referencePoint[0]-_xCentre);
  float DistXY = (_xCentre-xPoint[0])*(_xCentre-xPoint[0]) + (_yCentre-xPoint[1])*(_yCentre-xPoint[1]);
  DistXY = sqrt(DistXY);
  DistXY = fabs(DistXY - _radius);
  
  int nCircles = 0;

  if (fabs(_tanLambda*_radius)>1.0e-20) {
    float xCircles = phi0 - phi -_charge*(xPoint[2]-_referencePoint[2])/(_tanLambda*_radius);
    xCircles = xCircles/_const_2pi;
    int n1,n2;

    if (xCircles >= 0.) {
	n1 = int(xCircles);
	n2 = n1 + 1;
    }
    else {
	n1 = int(xCircles) - 1;
	n2 = n1 + 1;
    }
    
    if (fabs(n1-xCircles) < fabs(n2-xCircles)) {
	nCircles = n1;
    }
    else {
	nCircles = n2;
    }

  }
  
  float DPhi = _const_2pi*((float)nCircles) + phi - phi0;

  zOnHelix = _referencePoint[2] - _charge*_radius*_tanLambda*DPhi;

  float DistZ = fabs(zOnHelix - xPoint[2]);
  float Time;

  if (fabs(_momentum[2]) > 1.0e-20) {
    Time = (zOnHelix - _referencePoint[2])/_momentum[2];
  }
  else {
    Time = _charge*_radius*DPhi/_pxy;
  }
  
  Distance[0] = DistXY;
  Distance[1] = DistZ;
  Distance[2] = sqrt(DistXY*DistXY+DistZ*DistZ);
  
  pointOnHelix[0] = xOnHelix;
  pointOnHelix[1] = yOnHelix;
  pointOnHelix[2] = zOnHelix;

  return Time;


}

void HelixClass::setHelixEdges(float * xStart, float * xEnd) {

  // FIXME: This method is not used at all. 2006/06/21 OW
  for (int i=0; i<3; ++i) {
    _xStart[i] = xStart[i];
    _xEnd[i] = xEnd[i];
  }

}

float * HelixClass::getDistanceToHelix(HelixClass * helix, float * pos, float * mom) {

  float x01 = _xCentre;
  float y01 = _yCentre;
  
  float x02 = helix->getXC();
  float y02 = helix->getYC();
  
  float rad1 = _radius;
  float rad2 = helix->getRadius();

  float distance = sqrt((x01-x02)*(x01-x02)+(y01-y02)*(y01-y02));

  // FIXME: change return value of this method
  return pos;

}

float HelixClass::getDistanceToLine(LineClass * line) {

  return 1;
}

void HelixClass::getExtrapolatedMomentum(float * pos, float * momentum) {

  float phi = atan2(pos[1]-_yCentre,pos[0]-_xCentre);
  if (phi <0.) phi += _const_2pi;
  phi = phi - _phiAtPCA + _phi0;
  momentum[0] = _pxy*cos(phi);
  momentum[1] = _pxy*sin(phi);
  momentum[2] = _momentum[2];


}


float HelixClass::getPathLength(float* point1, float* point2) {

  // FIXME: check mechanism whether the points are on the helix or not is needed
  
  float pathLength = sqrt( (_bZ*_radius)*(_bZ*_radius) + 1 ) * (point2[2] - point1[2]);

  // for integration in negativ z direction change the sign:
  if ( pathLength < 0.0 ) pathLength *= -1.0;

  return pathLength;

}


bool HelixClass::isOnHelix(const float* point) {

  // still missing

  return false;

}


void HelixClass::printParameters() {

  
  std::cout << "Printout of the different helix parameters: " << std::endl << std::endl;
  std::cout << "0. Parametrisation, i.e. x_vec, p_vec, q and B ( done with 'Initialize_VP(float * pos, float * mom, float q, float B)' ): " << std::endl;
  std::cout << "pos = " << "( " << _referencePoint[0] << ", " << _referencePoint[1] << ", " << _referencePoint[2] << " )" << " [ pos = referencePoint ]" << std::endl 
	    << "p @ pos = " << "( " << _momentum[0] << ", " << _momentum[1] << ", " << _momentum[2] << " )" << std::endl 
	    << "q = " << _charge << std::endl
	    << "B = " << _bField << std::endl << std::endl;

  std::cout << "1. Parametrisation, i.e. xc, yc, r, bZ, phi0, B, sin(pz) and z0 ( done with 'Initialize_BZ(float xCentre, float yCentre, float radius, float bZ, float phi0, float B, float signPz, float zBegin)' ): " << std::endl;
  std::cout << "xc = " << _xCentre << std::endl
	    << "yc = " << _yCentre << std::endl
	    << "r  = " << _radius  << std::endl
	    << "bZ = " << _bZ      << std::endl
	    << "phi0 = " << _phi0  << std::endl
	    << "B = " << _bField   << std::endl
            << "signPz = " << _momentum[2]/fabs(_momentum[2]) << " ( ??? )" << std::endl
	    << "z0 = " << _z0 << " ( ??? )" << std::endl << std::endl;

  std::cout << "2. Parametrisation, i.e. xc, yc, zc, R, omega, phi0, v0, is still missing." << std::endl << std::endl;

  std::cout << "3. Parametrisation ( L3 like ), i.e. phi0, d0, z0, Omega, tanlam and  B ( done with 'Initialize_Canonical(float phi0, float d0, float z0, float omega, float tanlambda, float B)' ): " << std::endl;
  std::cout << "phi0 = " << _phi0  << std::endl
	    << "d0 = " << _d0 << std::endl
	    << "z0 = " << _z0 << std::endl
	    << "Omega = " << _omega << std::endl
	    << "tanlambda = " << _tanLambda << std::endl
	    << "B = " <<  _bField << std::endl
	    << "refPoint = " << "( " << _referencePoint[0] << ", " << _referencePoint[1] << ", " << _referencePoint[2] << " )" << std::endl << std::endl;

  std::cout << "Other member variables:" << std::endl;
  std::cout << "pxy = " << _pxy << std::endl
	    << "phiRefPoint = " << _phiRefPoint << std::endl
	    << "phiAtPCA = " << _phiAtPCA << std::endl
	    << "xAtPCA = " << _xAtPCA << std::endl
	    << "yAtPCA = " << _yAtPCA << std::endl
	    << "pxAtPCA = " << _pxAtPCA << std::endl
	    << "pyAtPCA = " << _pyAtPCA << std::endl
	    << "phiMomRefPoint = " << _phiMomRefPoint << std::endl
	    << "xStart = " << "( " << _xStart[0] << ", " << _xStart[1] << ", " << _xStart[2] << " )" << std::endl 
	    << "xEnd = " << "( " << _xEnd[0] << ", " << _xEnd[1] << ", " << _xEnd[2] << " )" << std::endl << std::endl;
    
}
