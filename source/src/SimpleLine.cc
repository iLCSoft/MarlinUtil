#include <SimpleLine.h>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/Vector.h>
#include <float.h>

SimpleLine::SimpleLine( LCVector3D ref , LCVector3D direction ) 
{
  _line = LCLine3D(ref,direction);  
}

LCVector3D SimpleLine::getPosition(double s, LCErrorMatrix* errors) const 
{
  if( errors != 0 ) 
    *errors = LCErrorMatrix( 3 , 0 ) ;

  return _line.position(s) ;
}

LCVector3D SimpleLine::getDirection(double s,  LCErrorMatrix* errors) const
{
  if( errors != 0 ) 
    *errors = LCErrorMatrix( 3 , 0 ) ;

  return _line.direction() ;
}

LCErrorMatrix SimpleLine::getCovarianceMatrix( double s) const  
{
  return LCErrorMatrix( 6 , 0 ) ;
}

double SimpleLine::getPathAt(const LCVector3D position ) const  
{
  return _line.projectPoint(position) ;
}

double SimpleLine::getIntersectionWithPlane( LCPlane3D p, bool& pointExists) const
{
  double s = _line.intersectionWithPlane(p,pointExists) ;
  return s ;
}

double SimpleLine::getIntersectionWithCylinder(const LCCylinder & cylinder,
					       bool & pointExists) const
{
  LCVector3D middlePoint = 
    (cylinder.startPoint() + cylinder.endPoint())/2.;
  double minDistance = sqrt( 0.25*cylinder.length()*cylinder.length() 
			     + cylinder.radius()*cylinder.radius() ) ;
  double sProject = _line.projectPoint(middlePoint);
  double pointLineDistance = ( middlePoint - _line.position( sProject ) ).mag();
  if (pointLineDistance > minDistance)
    { // Line does not hit the cylinder 
      pointExists = false ;
      return 0 ;
    }

  double sStart = sProject - minDistance ; 
  double sEnd   = sProject + minDistance ;

  if (sStart > sEnd)
    {
      double temp = sStart;
      sStart = sEnd;
      sEnd = temp;
    }

  if ( (sStart < 0.) && (sEnd < 0.) )
    { // Intersection is in backwards direction
      pointExists = false ;
      return 0;
    }
  else if ( (sStart < 0.) && (sEnd > 0.) )
    { // intersection region starts in backwards direction
      sStart = 0;
    }

  double s = sStart, sOld = -DBL_MAX;
  double epsilon = 0.0000001;

  while ( (s-sOld) > epsilon )
    {
      double d = fabs( cylinder.distance( _line.position(s) ) ) ;
      sOld = s;
      if (s > sEnd)
        { // Problem! no intersection !
          pointExists = false ;
          return 0;
        }
      s += d;
    }

  pointExists = true ;
  return s ;
}
