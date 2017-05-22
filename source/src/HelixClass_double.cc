#include "HelixClass_double.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "ced_cli.h"

HelixClass_double::HelixClass_double() {}

HelixClass_double::~HelixClass_double() {}

void HelixClass_double::Initialize_VP(double * pos, double * mom, double q, double B) {
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
    while (_phi0<0) _phi0+=_const_2pi;
    while (_phi0>=_const_2pi) _phi0-=_const_2pi;
    _xAtPCA = _xCentre + _radius*cos(_phiAtPCA);
    _yAtPCA = _yCentre + _radius*sin(_phiAtPCA);
    //    _d0 = -_xAtPCA*sin(_phi0) + _yAtPCA*cos(_phi0);
    double pxy = double(_pxy);
    double radius = pxy/double(_FCT*B);
    double xCentre = double(pos[0]) + radius*double(cos(_phiMomRefPoint-_const_pi2*q));
    double yCentre = double(pos[1]) + radius*double(sin(_phiMomRefPoint-_const_pi2*q));
    
    double d0;

    if (q>0) {
      d0 = double(q)*radius - double(sqrt(xCentre*xCentre+yCentre*yCentre));
    }
    else {
      d0 = double(q)*radius + double(sqrt(xCentre*xCentre+yCentre*yCentre));
    }

    _d0 = double(d0);

//     if (fabs(_d0)>0.001 ) {
//       std::cout << "New helix : " << std::endl;
//       std::cout << " Position : " << pos[0] 
// 		<< " " << pos[1]
// 		<< " " << pos[2] << std::endl;
//       std::cout << " Radius = " << _radius << std::endl;
//       std::cout << " RC = " << sqrt(_xCentre*_xCentre+_yCentre*_yCentre) << std::endl;  
//       std::cout << " D0 = " << _d0 << std::endl;
//     }

    _pxAtPCA = _pxy*cos(_phi0);
    _pyAtPCA = _pxy*sin(_phi0);
    double deltaPhi = _phiRefPoint - _phiAtPCA;    
    double xCircles = -pos[2]*q/(_radius*_tanLambda) - deltaPhi;
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

}

void HelixClass_double::Initialize_Canonical(double phi0, double d0, double z0, 
				      double omega, double tanLambda, double B) {
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
}


void HelixClass_double::Initialize_BZ(double xCentre, double yCentre, double radius, 
			       double bZ, double phi0, double B, double signPz,
			       double zBegin) {

  _radius = radius;
  _pxy = _FCT*B*_radius;
  _charge = -(bZ*signPz)/fabs(bZ*signPz);
  _momentum[2] = -_charge*_pxy/(bZ*_radius);
  _xCentre = xCentre;
  _yCentre = yCentre;
  _omega = _charge/radius;
  _phiAtPCA = atan2(-_yCentre,-_xCentre);
  _phi0 = -_const_pi2*_charge + _phiAtPCA;
  while (_phi0<0) _phi0+=_const_2pi;
  while (_phi0>=_const_2pi) _phi0-=_const_2pi;
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
  
  double deltaPhi = _phiRefPoint - _phiAtPCA;    
  double xCircles = bZ*_referencePoint[2] - deltaPhi;
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

}

const double * HelixClass_double::getMomentum() {
    return _momentum;
}
const double * HelixClass_double::getReferencePoint() {
    return _referencePoint;
}
double HelixClass_double::getPhi0() {
  if (_phi0<0.0)
    _phi0 += 2*M_PI;
  return _phi0;
}
double HelixClass_double::getD0() {
    return _d0;
}
double HelixClass_double::getZ0() {
    return _z0;
}
double HelixClass_double::getOmega() {
    return _omega;
}
double HelixClass_double::getTanLambda() {
    return _tanLambda;
}
double HelixClass_double::getPXY() {
    return _pxy;
}
double HelixClass_double::getXC() {
  return _xCentre;
}

double HelixClass_double::getYC() {
  return _yCentre;
}

double HelixClass_double::getRadius() {
  return _radius;
}

double HelixClass_double::getBz() {
  return _bZ;
}

double HelixClass_double::getPhiZ() {
  return _phiZ;
}

double HelixClass_double::getCharge() {
    return _charge;
}

double HelixClass_double::getPointInXY(double x0, double y0, double ax, double ay, 
			      double * ref , double * point) {

  double time;

  double AA = sqrt(ax*ax+ay*ay);


  if (AA <= 0) {
    time = -1.0e+20; 
    return time;
  }


  double BB = ax*(x0-_xCentre) + ay*(y0-_yCentre);
  BB = BB / AA;

  double CC = (x0-_xCentre)*(x0-_xCentre) 
    + (y0-_yCentre)*(y0-_yCentre) - _radius*_radius;

  CC = CC / AA;

  double DET = BB*BB - CC;
  double tt1 = 0.;
  double tt2 = 0.;
  double xx1,xx2,yy1,yy2; 


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
  
  double phi1 = atan2(yy1-_yCentre,xx1-_xCentre);
  double phi2 = atan2(yy2-_yCentre,xx2-_xCentre);
  double phi0 = atan2(ref[1]-_yCentre,ref[0]-_xCentre);

  double dphi1 = phi1 - phi0;
  double dphi2 = phi2 - phi0;

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


double HelixClass_double::getPointOnCircle(double Radius, double * ref, double * point) {

  double distCenterToIP = sqrt(_xCentre*_xCentre + _yCentre*_yCentre);

  point[0] = 0.0;
  point[1] = 0.0;
  point[2] = 0.0;

  if ((distCenterToIP+_radius)<Radius) {
    double xx = -1.0e+20;
    return xx;
  }

  if ((_radius+Radius)<distCenterToIP) {
    double xx = -1.0e+20;
    return xx;
  }

  double phiCentre = atan2(_yCentre,_xCentre);
  double phiStar   = Radius*Radius + distCenterToIP*distCenterToIP 
                                    - _radius*_radius;

  phiStar = 0.5*phiStar/fmax(1.0e-20,Radius*distCenterToIP);
  
  if (phiStar > 1.0) 
    phiStar = 0.9999999;
  
  if (phiStar < -1.0)
    phiStar = -0.9999999;
  
  phiStar = acos(phiStar);

  double tt1,tt2,time;

  double xx1 = Radius*cos(phiCentre+phiStar);
  double yy1 = Radius*sin(phiCentre+phiStar);

  double xx2 = Radius*cos(phiCentre-phiStar);
  double yy2 = Radius*sin(phiCentre-phiStar);


  double phi1 = atan2(yy1-_yCentre,xx1-_xCentre);
  double phi2 = atan2(yy2-_yCentre,xx2-_xCentre);
  double phi0 = atan2(ref[1]-_yCentre,ref[0]-_xCentre);

  double dphi1 = phi1 - phi0;
  double dphi2 = phi2 - phi0;

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
  

  double time2;
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


double HelixClass_double::getPointInZ(double zLine, double * ref, double * point) {

  double time = zLine - ref[2];

  if (_momentum[2] == 0.) {
    time = -1.0e+20;
    point[0] = 0.;
    point[1] = 0.;
    point[2] = 0.;
    return time;
  }

  time = time/_momentum[2];

  double phi0 = atan2(ref[1] - _yCentre , ref[0] - _xCentre);
  double phi = phi0 - _charge*_pxy*time/_radius;
  double xx = _xCentre + _radius*cos(phi);
  double yy = _yCentre + _radius*sin(phi);

  point[0] = xx;
  point[1] = yy;
  point[2] = zLine;

  return time;


}

double HelixClass_double::getDistanceToPoint(double const* xPoint, double * Distance) {

  double zOnHelix;
  double phi = atan2(xPoint[1]-_yCentre,xPoint[0]-_xCentre);
  double phi0 = atan2(_referencePoint[1]-_yCentre,_referencePoint[0]-_xCentre);
  //calculate distance to XYprojected centre of Helix, comparing this with distance to radius around centre gives DistXY
  double DistXY = (_xCentre-xPoint[0])*(_xCentre-xPoint[0]) + (_yCentre-xPoint[1])*(_yCentre-xPoint[1]);
  DistXY = sqrt(DistXY);
  DistXY = fabs(DistXY - _radius);
  
  int nCircles = 0;

  if (fabs(_tanLambda*_radius)>1.0e-20) {
    double xCircles = phi0 - phi -_charge*(xPoint[2]-_referencePoint[2])/(_tanLambda*_radius);
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
  
  double DPhi = _const_2pi*((double)nCircles) + phi - phi0;

  zOnHelix = _referencePoint[2] - _charge*_radius*_tanLambda*DPhi;

  double DistZ = fabs(zOnHelix - xPoint[2]);
  double Time;

  if (fabs(_momentum[2]) > 1.0e-20) {
    Time = (zOnHelix - _referencePoint[2])/_momentum[2];
  }
  else {
    Time = _charge*_radius*DPhi/_pxy;
  }

  Distance[0] = DistXY;
  Distance[1] = DistZ;
  Distance[2] = sqrt(DistXY*DistXY+DistZ*DistZ);

  return Time;


}

//When we are not interested in the exact distance, we can check if we are
//already far enough away in XY, before we start calculating in Z as the
//distance will only increase
double HelixClass_double::getDistanceToPoint(const std::vector<double>& xPoint, double distCut) {
  //calculate distance to XYprojected centre of Helix, comparing this with distance to radius around centre gives DistXY
  double tempx = xPoint[0]-_xCentre;
  double tempy = xPoint[1]-_yCentre;
  double tempsq = sqrt(tempx*tempx + tempy*tempy);
  double tempdf = tempsq - _radius;
  double DistXY = fabs( tempdf );
  //If this is bigger than distCut, we dont have to know how much bigger this is
  if( DistXY > distCut) {
    return DistXY;
  }

  int nCircles = 0;
  double phi = atan2(tempy,tempx);
  double phi0 = atan2(_referencePoint[1]-_yCentre,_referencePoint[0]-_xCentre);
  double phidiff = phi0-phi;
  double  tempz = xPoint[2] - _referencePoint[2];//Yes referencePoint
  double tanradius = _tanLambda*_radius;
  if (fabs(tanradius)>1.0e-20) {
    double xCircles = (phidiff -_charge*tempz/tanradius)/_const_2pi;
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
  double DistZ = - tempz - _charge*tanradius*(_const_2pi*((double)nCircles) - phidiff);
  return sqrt(DistXY*DistXY+DistZ*DistZ);
}//getDistanceToPoint(vector,double)

double HelixClass_double::getDistanceToPoint(const double* xPoint, double distCut) {
  std::vector<double> xPosition(xPoint, xPoint + 3 );//We are expecting three coordinates, must be +3, last element is excluded!
  return getDistanceToPoint(xPosition, distCut);
}//getDistanceToPoint(double*,double)



void HelixClass_double::setHelixEdges(double * xStart, double * xEnd) {
  for (int i=0; i<3; ++i) {
    _xStart[i] = xStart[i];
    _xEnd[i] = xEnd[i];
  }

}

double HelixClass_double::getDistanceToHelix(HelixClass_double * helix, double * pos, double * mom) {

  double x01 = getXC();
  double y01 = getYC();
  
  double x02 = helix->getXC();
  double y02 = helix->getYC();
  
  double rad1 = getRadius();
  double rad2 = helix->getRadius();

  double distance = sqrt((x01-x02)*(x01-x02)+(y01-y02)*(y01-y02));

  bool singlePoint = true;

  double phi1 = 0;
  double phi2 = 0;

  if (rad1+rad2<distance) {
    phi1 = atan2(y02-y01,x02-x01);
    phi2 = atan2(y01-y02,x01-x02);
  }
  else if (distance+rad2<rad1) {
    phi1 = atan2(y02-y01,x02-x01);
    phi2 = phi1;
  }
  else if (distance+rad1<rad2) {
    phi1 = atan2(y01-y02,x01-x02);
    phi2 = phi1;
  }
  else {
    singlePoint = false;
    double cosAlpha = 0.5*(distance*distance+rad2*rad2-rad1*rad1)/(distance*rad2);
    double alpha = acos(cosAlpha);
    double phi0 = atan2(y01-y02,x01-x02);
    phi1 = phi0 + alpha;
    phi2 = phi0 - alpha;
  }
  

  double ref1[3];
  double ref2[3];
  for (int i=0;i<3;++i) {
    ref1[i]=_referencePoint[i];
    ref2[i]=helix->getReferencePoint()[i];
  }
  
  double pos1[6]; //last three accessed in getPointOnCircle
  double pos2[6]; //last three accessed in getPointOnCircle
  double mom1[3];
  double mom2[3];


  if (singlePoint ) {

    double xSect1 = x01 + rad1*cos(phi1);
    double ySect1 = y01 + rad1*sin(phi1);
    double xSect2 = x02 + rad2*cos(phi2);
    double ySect2 = y02 + rad2*sin(phi2);
    double R1 = sqrt(xSect1*xSect1+ySect1*ySect1);
    double R2 = sqrt(xSect2*xSect2+ySect2*ySect2);

    getPointOnCircle(R1,ref1,pos1);
    helix->getPointOnCircle(R2,ref2,pos2);
    
  }
  else {    

    double xSect1 = x02 + rad2*cos(phi1);
    double ySect1 = y02 + rad2*sin(phi1);
    double xSect2 = x02 + rad2*cos(phi2);
    double ySect2 = y02 + rad2*sin(phi2);

//     std::cout << "(xSect1,ySect1)=(" << xSect1 << "," << ySect1 << ")" << std::endl;
//     std::cout << "(xSect2,ySect2)=(" << xSect2 << "," << ySect2 << ")" << std::endl;

    double temp21[3];
    double temp22[3];

    double phiI2  = atan2(ref2[1]-y02,ref2[0]-x02); 
    double phiF21 = atan2(ySect1-y02,xSect1-x02);
    double phiF22 = atan2(ySect2-y02,xSect2-x02);
    double deltaPhi21 = phiF21 - phiI2;
    double deltaPhi22 = phiF22 - phiI2;
    double charge2 = helix->getCharge();
    double pxy2 = helix->getPXY();
    double pz2   = helix->getMomentum()[2];
    if (deltaPhi21 < 0 && charge2 < 0) {
      deltaPhi21 += _const_2pi;
    }
    else if (deltaPhi21 > 0 && charge2 > 0) { 
      deltaPhi21 -= _const_2pi;
    }

    if (deltaPhi22 < 0 && charge2 < 0) {
      deltaPhi22 += _const_2pi;
    }
    else if (deltaPhi22 > 0 && charge2 > 0) { 
      deltaPhi22 -= _const_2pi;
    }

    double time21 = -charge2*deltaPhi21*rad2/pxy2;
    double time22 = -charge2*deltaPhi22*rad2/pxy2;

    double Z21 = ref2[2] + time21*pz2;
    double Z22 = ref2[2] + time22*pz2;

    temp21[0] = xSect1; temp21[1] = ySect1; temp21[2] = Z21;
    temp22[0] = xSect2; temp22[1] = ySect2; temp22[2] = Z22;
    

//     std::cout << "temp21 = " << temp21[0] << " " << temp21[1] << " " << temp21[2] << std::endl;
//     std::cout << "temp22 = " << temp22[0] << " " << temp22[1] << " " << temp22[2] << std::endl;


    double temp11[3];
    double temp12[3];

    double phiI1  = atan2(ref1[1]-y01,ref1[0]-x01); 
    double phiF11 = atan2(ySect1-y01,xSect1-x01);
    double phiF12 = atan2(ySect2-y01,xSect2-x01);
    double deltaPhi11 = phiF11 - phiI1;
    double deltaPhi12 = phiF12 - phiI1;
    double charge1 = _charge;
    double pxy1 = _pxy;
    double pz1   = _momentum[2];
    if (deltaPhi11 < 0 && charge1 < 0) {
      deltaPhi11 += _const_2pi;
    }
    else if (deltaPhi11 > 0 && charge1 > 0) { 
      deltaPhi11 -= _const_2pi;
    }

    if (deltaPhi12 < 0 && charge1 < 0) {
      deltaPhi12 += _const_2pi;
    }
    else if (deltaPhi12 > 0 && charge1 > 0) { 
      deltaPhi12 -= _const_2pi;
    }

    double time11 = -charge1*deltaPhi11*rad1/pxy1;
    double time12 = -charge1*deltaPhi12*rad1/pxy1;

    double Z11 = ref1[2] + time11*pz1;
    double Z12 = ref1[2] + time12*pz1;

    temp11[0] = xSect1; temp11[1] = ySect1; temp11[2] = Z11;
    temp12[0] = xSect2; temp12[1] = ySect2; temp12[2] = Z12;
    
//     std::cout << "temp11 = " << temp11[0] << " " << temp11[1] << " " << temp11[2] << std::endl;
//     std::cout << "temp12 = " << temp12[0] << " " << temp12[1] << " " << temp12[2] << std::endl;

    double Dist1 = 0;
    double Dist2 = 0;

    for (int j=0;j<3;++j) {
      Dist1 += (temp11[j]-temp21[j])*(temp11[j]-temp21[j]);
      Dist2 += (temp12[j]-temp22[j])*(temp12[j]-temp22[j]);
    }

    if (Dist1<Dist2) {
      for (int l=0;l<3;++l) {
	pos1[l] = temp11[l];
	pos2[l] = temp21[l];
      }
    }
    else {
       for (int l=0;l<3;++l) {
	pos1[l] = temp12[l];
	pos2[l] = temp22[l];
      }
    }

  }

  getExtrapolatedMomentum(pos1,mom1);
  helix->getExtrapolatedMomentum(pos2,mom2);

  double helixDistance = 0;

  for (int i=0;i<3;++i) {
    helixDistance += (pos1[i] - pos2[i])*(pos1[i] - pos2[i]);
    pos[i] = 0.5*(pos1[i]+pos2[i]);
    mom[i] = mom1[i]+mom2[i];
  }
  helixDistance = sqrt(helixDistance);

  return helixDistance;

}

void HelixClass_double::getExtrapolatedMomentum(double * pos, double * momentum) {

  double phi = atan2(pos[1]-_yCentre,pos[0]-_xCentre);
  if (phi <0.) phi += _const_2pi;
  phi = phi - _phiAtPCA + _phi0;
  momentum[0] = _pxy*cos(phi);
  momentum[1] = _pxy*sin(phi);
  momentum[2] = _momentum[2];


}
