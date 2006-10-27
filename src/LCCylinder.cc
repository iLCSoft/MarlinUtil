#include <LCCylinder.h>
#include <LCPlane3D.h>
#include <LCLine3D.h>

#include <iostream>

#include <cmath>
#include <float.h>
#include <exception>

LCCylinder::LCCylinder(LCVector3D point1, 
		       LCVector3D point2, 
		       double radius,
		       bool endPlane) 
{
  _axisSstartPoint = point1;
  _axisEndPoint = point2;

  _radius = fabs(radius);
  _endPlane = endPlane;
}

LCCylinder::LCCylinder(double radius, 
		       LCVector3D point, 
		       LCVector3D axis, 
		       bool endPlane) 
{
  _axisSstartPoint = point - axis ;
  _axisEndPoint    = point + axis ;

  _radius = fabs(radius);
  _endPlane = endPlane;
}

LCCylinder::LCCylinder(const LCCylinder & cylinder) 
{
  _axisSstartPoint = cylinder._axisSstartPoint ;
  _axisEndPoint    = cylinder._axisEndPoint ;

  _radius = cylinder._radius ;
  _endPlane = cylinder._endPlane;
}

LCCylinder & LCCylinder::operator=(const LCCylinder & rhs) 
{
  _axisSstartPoint = rhs._axisSstartPoint ;
  _axisEndPoint    = rhs._axisEndPoint ;

  _radius = rhs._radius ;
  _endPlane = rhs._endPlane;

  return *this;
}

LCVector3D LCCylinder::startPoint() const 
{
  return _axisSstartPoint; 
}

LCVector3D LCCylinder::endPoint() const 
{
  return _axisEndPoint  ;
}

LCVector3D LCCylinder::axisDirection() const 
{
  return (_axisEndPoint - _axisSstartPoint).unit();
}

double LCCylinder::length() const 
{
  return (_axisEndPoint - _axisSstartPoint).mag();
}

double LCCylinder::radius() const 
{
  return _radius;
}

double LCCylinder::distance(const LCVector3D & point) const 
{
  int dummy ;
  std::cout << "punkt1 " <<  point << " punkt2 " << projectPoint( point, dummy ) << " code " << dummy << std::endl ;
  return (point - projectPoint( point, dummy ) ).mag() ;
}

LCVector3D LCCylinder::projectPoint(const LCVector3D & point, int & code) const 
{
  LCLine3D a( _axisSstartPoint , axisDirection() );
  std::cout << " a: " << a.position() << " d: " << a.direction() << std::endl;
  double s = a.projectPoint( _axisSstartPoint ) ;
  double e = a.projectPoint( _axisEndPoint ) ;
  double p = a.projectPoint( point ) ;
  double d = a.distance( point ) ;

  std::cout << "s " << s << " e " << e << " p " << p << " d " << d << "point: " << point <<  std::endl;

  double drp = fabs( d - radius() ) ;
  double dsp = fabs( s - p ) ;
  double dep = fabs( e - p ) ;

  std::cout << "drp: " << drp << " dsp: " << dsp << " dep: " <<  dep << std::endl;

  LCVector3D projection;

  // classify in which region the point is located: 
  // point between the two planes at both ends 
  if ( p >= s && p <= e ) 
    { std::cout << "punkt innherhalb des zylinders." << std::endl;
    std::cout << "endplane " << _endPlane << " d " << d << " radius " << radius() << std::endl;
      if (_endPlane && (d <= radius()) )
	{
	  if ( (drp <= dsp) && (drp <= dep) )
	    { std::cout << "mantel " << drp << " " << dsp << " " << dep << std::endl;
	      projection = ( point - a.position(p) ).unit() ;
	      if (projection.mag() < 0.00001) projection = axisDirection().orthogonal().unit() ;
	      projection *= radius() ;
	      projection = a.position(p) + projection;
	      code = 3;
	      return projection;
	    }
	  else if (dsp <= dep)
	    { std::cout << "vorn " << drp << " " << dsp << " " << dep << std::endl;
	      LCPlane3D sPlane(-axisDirection(),_axisSstartPoint);
	      projection = sPlane.projectPoint( point );
	      code = 1;
	      return projection;
	    }
	  else // if ( dep < dsp )
	    { std::cout << "hinten " << drp << " " << dsp << " " << dep << std::endl;
	      LCPlane3D ePlane(axisDirection(),_axisEndPoint);
	      projection = ePlane.projectPoint( point );
	      code = 2;
	      return projection;
	    }
	}
      else
	{
	  projection = ( point - a.position(p) ).unit() ;
	  if (projection.mag() < 0.00001) projection = axisDirection().orthogonal().unit() ;
	  std::cout << "projection1: " << projection << std::endl;
	  projection *= radius() ;
	  std::cout << "projection2: " << projection << std::endl;
	  projection = a.position(p) + projection;
	  std::cout << "projection3: " << projection << std::endl;
	  code = 3;
	  std::cout << "hier return!" << std::endl;
	  return projection;
	}
    }
  else // outside the two planes at the end
    { std::cout << "punkt auserhalb des zylinders." << std::endl;
      if (_endPlane && (d <= radius()) )
	{
	  if ( p < s)
	    {
	      LCPlane3D sPlaneo(-axisDirection(),_axisSstartPoint);
	      projection = sPlaneo.projectPoint( point );
	      code = 1;
	      return projection;
	    }
	  else // if ( p > e )
	    {
	      LCPlane3D ePlaneo(axisDirection(),_axisEndPoint);
	      projection = ePlaneo.projectPoint( point );
	      code = 2;
	      return projection;
	    }
	}
      else
	{
	  if ( p < s)
	    {
	      projection = ( point - a.position(p) ).unit() ;
	      if (projection.mag() < 0.00001) projection = axisDirection().orthogonal().unit() ;
	      projection *= radius() ;
	      projection = a.position(s) + projection;
	      code = 0;
	      return projection;
	    }
	  else // if ( p > e )
	    {
	      projection = ( point - a.position(p) ).unit() ;
	      if (projection.mag() < 0.00001) projection = axisDirection().orthogonal().unit() ;
	      projection *= radius() ;
	      projection = a.position(e) + projection;
	      code = 0;
	      return projection;
	    }
	}
    }

  std::cout << "Classification faild!!!" << std::endl;

  return projection;
}

bool LCCylinder::isInside(const LCVector3D & point) const 
{
  LCLine3D a( _axisSstartPoint , axisDirection() );

  if ( radius() < a.distance( point ) ) return false ;

  double s = a.projectPoint( _axisSstartPoint ) ;
  double e = a.projectPoint( _axisEndPoint ) ;
  double p = a.projectPoint( point ) ;
  if (p < s || p > e) return false;

  return true ;
}

bool LCCylinder::operator==(const LCCylinder & rhs) const 
{
  return (_radius == rhs._radius &&
	  _axisSstartPoint == rhs._axisSstartPoint &&
	  _axisEndPoint == rhs._axisEndPoint &&
	  _endPlane == rhs._endPlane ) ;
}

bool LCCylinder::operator!=(const LCCylinder & rhs) const 
{
  return (_radius != rhs._radius ||
	  _axisSstartPoint != rhs._axisSstartPoint ||
	  _axisEndPoint != rhs._axisEndPoint || 
	  _endPlane != rhs._endPlane) ;
}

