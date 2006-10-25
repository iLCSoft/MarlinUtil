#include "SimpleLine.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"

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
  //FIXME: needs to be implemented
  pointExists = false ;
  return 0 ;
}
