#include <LCLine3D.h>
#include <LCPlane3D.h>
#include <iostream>

#include <cmath>
#include <float.h>
#include <exception>

LCLine3D::LCLine3D()
{
  _point.set(0.,0.,0.);
  _direction.set(1.,0.,0.);
}

LCLine3D::LCLine3D(const LCVector3D point, const LCVector3D direction) 
{
  _direction = direction.unit();
  _point = point - ( point * _direction ) * _direction ;
}

LCLine3D::LCLine3D(const LCLine3D & line) 
{
  _point     = line._point;
  _direction = line._direction;
}

LCLine3D & LCLine3D::operator=(const LCLine3D & rhs) 
{
  _point     = rhs._point;
  _direction = rhs._direction;

  return *this;
}

LCVector3D LCLine3D::position(const double s) const 
{
  return (_point + s*_direction) ;
}

LCVector3D LCLine3D::direction() const 
{
  return _direction;
}

double LCLine3D::distance(const LCVector3D & point) const 
{
  return ( point - position( projectPoint( point ) ) ).mag() ;
}

double LCLine3D::projectPoint(const LCVector3D & point) const 
{
  // the last therm : (...) / _direction.mag2() is not there becaus 
  // the _direction vector is normalised.
  //  return ( 2*point*_direction - _point*_direction ) / _direction.mag2() ;
  return ( point*_direction - _point*_direction ) / _direction.mag2() ;
}

bool LCLine3D::operator==(const LCLine3D & rhs) const 
{
  return (_point == rhs._point &&
	  _direction == rhs._direction);
}

bool LCLine3D::operator!=(const LCLine3D & rhs) const 
{
  return (_point != rhs._point ||
	  _direction != rhs._direction) ;
}

double LCLine3D::intersectionWithPlane(const LCPlane3D plane, bool& pointExists) const 
{
  double c = direction() * plane.normal() ;

  if (c == 0)
    { // no interaction 
      pointExists = false;
      return DBL_MAX;
    }

  pointExists = true;
  return - ( position() * plane.normal() + plane.d() ) / c ;
}

std::ostream & operator << (std::ostream &os, const LCLine3D &l)
{
  return os << l.position() << "+s*" << l.direction() ;
}
