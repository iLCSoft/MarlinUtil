#ifndef __B_UTIL_H
#define __B_UTIL_H

#include <math.h>

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//   Common independent classes
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

class Point2D {
//-----------------------------------------------------------------------
 public:
  double x;  // Cartesian coordinate 
  double y;  // Cartesian coordinate

  double rz;
  double rz2;
  
  void calc() { rz2=x*x+y*y; rz=sqrt(rz2); }
  void calc(double _x,double _y) { x=_x; y=_y; rz2=x*x+y*y; rz=sqrt(rz2); }
  void calc(Point2D const &a) { x=a.x; y=a.y; rz2=a.rz2; rz=a.rz; }
  double dist2(Point2D &a) { return (x-a.x)*(x-a.x)+(y-a.y)*(y-a.y); }
  double dist(Point2D &a) { return sqrt(dist2(a)); }

  Point2D() : x(0.),y(0.),rz(0.),rz2(0.) { }
  Point2D(double _x,double _y) { calc(_x,_y); }
  Point2D(Point2D const & a) { calc(a); }
  Point2D& operator=(const Point2D &a) { calc(a); return *this; }
};

class Vector3D;

//-----------------------------------------------------------------------
class Point3D : public Point2D {
//-----------------------------------------------------------------------
 public:
  double z;  // Cartesian coordinate
  double r;  // Radius-vector length

  void calc() { Point2D::calc(); r=sqrt(rz2+z*z); }
  void calc(double _x,double _y,double _z) { Point2D::calc(_x,_y); z=_z; r=sqrt(rz2+z*z); }
  void calc(Point3D const &_a) { Point2D::calc(_a); z=_a.z; r=_a.r; }
  
  double dist2(Point3D &a) { return Point2D::dist2(a)+(z-a.z)*(z-a.z); }
  double dist(Point3D &a) { return sqrt(dist2(a)); }

  //  Point3D operator=(Point3D &a) { return a; }
  
  Point3D operator+(Point3D &a){
    Point3D t(x+a.x,y+a.y,z+a.z);
    return t;
  }
  
  Point3D() : Point2D(),z(0.),r(0.) {}
  Point3D(double _x,double _y,double _z) : Point2D(_x,_y) { z=_z; r=sqrt(rz2+z*z); }
  Point3D(Point3D const &a) : Point2D(a) { z=a.z; r=a.r; }

  double rdist(Point3D &a) { return fabs(a.r-r); }

  bool IsBetween(Point3D &p0,Point3D &p1);
  bool IsNear(float px,float py,float pz,float tolerance=5){
    if(fabs(px-x)+fabs(py-y)+fabs(pz-z)<tolerance)
      return true;
    return false;
  }
  Point3D& operator=(const Point3D &_a) { calc(_a); return *this; }
};

//-----------------------------------------------------------------------
class Vector3D : public Point3D {
//-----------------------------------------------------------------------
 public:
      
  Vector3D(Point3D &a,Point3D &b) : Point3D(b.x-a.x,b.y-a.y,b.z-a.z) { }
  Vector3D(Point3D const &a) : Point3D(a) { }
  Vector3D() : Point3D() { }
  Vector3D(double _x,double _y,double _z) : Point3D(_x,_y,_z) { }
  double dot(Vector3D &a) { return x*a.x+y*a.y+z*a.z; }
  double module() { return sqrt(x*x+y*y+z*z); }
  void divide(double a) { x/=a; y/=a; z/=a; }
  void mult(double a) { x*=a; y*=a; z*=a; }
  
  double length2() const { return x*x+y*y+z*z; }
  
  Point3D linear(double a,Point3D &b) {
    Point3D p(a*x+b.x,a*y+b.y,a*z+b.z);
    return p;
  }

  Point3D add(const Point3D &_a){
    Point3D p(x+_a.x,y+_a.y,z+_a.z);
    return p;
  }

  Point3D rotate(const Point3D &v1,const Point3D &v2,const Point3D &v3){
    Point3D p(x*v1.x+y*v2.x+z*v3.x,
	      x*v1.y+y*v2.y+z*v3.y,
	      x*v1.z+y*v2.z+z*v3.z);
    return p;
  }
};

Point3D Projection(Point3D &p,Point3D &p0,Point3D &p1);

//-----------------------------------------------------------------------
class Sphere3D {
//-----------------------------------------------------------------------
 public:
  static const double RADDEG;// = 57.2957795130823209;

  double r;  // Radius
  double t;  // Theta
  double p;  // Phi
  Sphere3D() : r(0.),t(0.),p(0.) {}
  void calc(Point3D &a) {
//    Fill Spherical coordinate system
    double rsq = hypot(a.x,a.y);
    r = hypot(rsq,a.z);
    t = atan2(rsq,a.z);
    p = atan2(a.y,a.x);
    if (p < 0.0) 
      p = 2*M_PI + p;
  }
};
//-----------------------------------------------------------------------
class SimpleLinearFit {
//-----------------------------------------------------------------------
public:
  double sum_x,sum_y,sum_x2,sum_y2,sum_xy;
  unsigned n;
  
  double a;
  double b;
  double md; // mean displacement

  void init(){
    sum_x=sum_y=sum_x2=sum_y2=sum_xy=0.;
    a=b=md=0.;
    n=0;
  }
  void add(double x,double y){
    sum_x+=x;
    sum_x2+=x*x;
    sum_y+=y;
    sum_xy+=x*y;
    n++;
  }
  void calc(){
    if(n>1){
      a=(n*sum_xy-sum_x*sum_y)/(n*sum_x2-sum_x*sum_x);
      b=(sum_y-a*sum_x)/n;
      md=(sum_y2+a*a*sum_x2+n*b*b+2*a*b*sum_x-2*a*sum_xy-2*b*sum_y)/n;
    }
  }
};

//-----------------------------------------------------------------------
template <class _Tp>
//-----------------------------------------------------------------------
class MinMax {
//-----------------------------------------------------------------------
public:
    _Tp min;
    _Tp max;
    MinMax() : min(100000), max(-100000) { }
    void init() { min=100000; max=-100000; }
    void set(_Tp v){
	if(v<min)
	    min=v;
	if(v>max)
	    max=v;
    }
};

#define ARSIZE(x) (sizeof(x)/sizeof(x[0]))

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//   Some useful functions
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
double angle_dist(double t1,double p1,double t2,double p2);

#endif /* __B_UTIL_H */
