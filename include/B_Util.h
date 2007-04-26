#ifndef __B_UTIL_H
#define __B_UTIL_H

#include <math.h>

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

class Point3D : public Point2D {
//-----------------------------------------------------------------------
 public:
  double z;  // Cartesian coordinate
  double r;

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
#endif /* __B_UTIL_H */
