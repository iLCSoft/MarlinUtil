#include "SimpleLine.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"

SimpleLine::SimpleLine( LCVector3D ref , LCVector3D direction ) {
  
  // set reference point so that s==0 is closest to origin
  _a = direction.unit() ;
  _r = ref - (ref * _a ) * _a ;
}

LCVector3D SimpleLine::getPosition(double s, LCErrorMatrix* errors) const {
  
  if( errors != 0 ) 
    *errors = LCErrorMatrix( 3 , 0 ) ;

  return _r + s * _a  ;
}

LCVector3D SimpleLine::getDirection(double s,  LCErrorMatrix* errors) const{
  
  if( errors != 0 ) 
    *errors = LCErrorMatrix( 3 , 0 ) ;

  return  _a  ;
}


LCErrorMatrix SimpleLine::getCovarianceMatrix( double s) const  {
  return LCErrorMatrix( 6 , 0 ) ;
}


double SimpleLine::getPathAt(const LCVector3D position ) const  {

  return ( position - _r ) * _a  ;
}


double SimpleLine::getIntersectionWithPlane( LCPlane3D p, bool& pointExists) const{

  // ---- Implemented solving the Matrix equation------------

  CLHEP::HepMatrix m( 3, 3 )  ;
  CLHEP::HepVector c( 3 ) ;
  
  m(1,1) =  p.a() ;
  m(1,2) =  p.b() ;
  m(1,3) =  p.c() ;
  m(2,1) =  1. ;
  m(2,2) =  0. ;
  m(2,3) =  - _a.x() / _a.z() ;
  m(3,1) =  0. ;
  m(3,2) =  1. ;
  m(3,3) =  - _a.y() / _a.z();
  
  c(1) =   - p.d() ;
  c(2) =  _r.x() - _r.z() * _a.x() / _a.z() ;
  c(3) =  _r.y() - _r.z() * _a.y() / _a.z() ;
  
  int invertError(0) ;
  m.invert( invertError ) ;

  pointExists = !invertError ; 
  
  CLHEP::HepVector xV = m * c ;

//   std::cout << " point computed: " << xV << std::endl ;

  double x = ( xV(3) - _r(2) ) / _a(2) ;  
  // this is sick: z is HepVector(3) and Hep3Vector(2) !!!!!
  
  return x ;
}

// double SimpleLine::getIntersectionWithPlane( LCPlane3D p, bool& pointExists) const{



//  // ---- Implemented using vector algebra ------------

//   LCVector3D pp0 = p.point( _r ) ;       // projection of s==0
//   LCVector3D pp1 = p.point( _r + _a ) ;  // projection of s==1

//   LCVector3D ap = pp1 - pp0  ;

//   // ---------------------------------------------
//   //    pp0 + x * ap = _r + x * _a
//   // --------------------------------------------

//   pointExists = ( ap != _a )  ;

//   double x = 1e99 ;

//   if( pointExists )

//   // ----  (_r - pp0 )  and ( ap - _a ) are parallel if solution exists
//     x =  LCVector3D( _r - pp0 ).mag() / LCVector3D( ap - _a ).mag() ;

//   return x ;
// }




double SimpleLine::getIntersectionWithCylinder(LCVector3D center, 
					       LCVector3D axis, 
					       double radius,
					       bool & pointExists) const  {

  //FIXME: needs to be implemented
  pointExists = false ;
  return 0 ;
}


  
