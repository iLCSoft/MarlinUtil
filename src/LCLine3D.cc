#include <LCLine3D.h>

#include <iostream>

#include <cmath>
#include <float.h>
#include <exception>

LCLine3D::LCLine3D(LCVector3D point, LCVector3D direction) 
{
  _point = point;
  _direction = direction.unit();

  LCVector3D origin(0.,0.,0.);
  _point = point( projectPoint(origin) );
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

LCVector3D LCLine3D::point(double s) const 
{
  return (_point + s*_direction) ;
}

LCVector3D LCLine3D::direction() const 
{
  return _direction;
}

double LCLine3D::distance(const LCVector3D & point) const 
{
  return (point-_point)*_direction ;
}

double LCLine3D::projectPoint(const LCVector3D & point) const 
{
  // the last therm : (...) / _direction.mag2() is not there becaus 
  // the _direction vector is normalised.
  return ( 2*point*_direction - _point*_direction ) / _direction.mag2() ;
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

