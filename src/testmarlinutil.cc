
#include <iostream>
#include "SimpleLine.h"
#include "CLHEP/Geometry/Transform3D.h"
#include <cstdlib>



double drand() ;

bool isequal(  double x , double y , double epsilon=1.e-9) {

  return ( 2.*fabs(x-y)/fabs(x+y) < epsilon ) ;  
}

void printLine( const SimpleLine& l , double s=0) {

  int nStep = 10 ;
  double s0 =  s + 0. ;
  double s1 = s0 + 1.  ;

  for(int i=0 ; i < nStep ; ++i ){

//     double s = s0 + double(i) * ( s1 -s0 ) / double( nStep ) ;

    double s = s1 - double(i) * ( s1 -s0 ) / double( nStep ) ;

    LCPoint3D p = l.getPosition( s ) ;
    std::cout << " s: " << s 
	      <<  " p : "  <<   p
	      <<  "  - s (p) " << l.getPathAt( p ) 
	      << std::endl  ;
  }
}

int main(){
  
  std::cout  << std::endl
	     << " ------- MarlinUtil test program ---------- " 
	     << std::endl
	     << std::endl ;

//   LCPoint3D ref( 0. , 2. , 3. ) ;
//   LCVector3D dir( 0. , 1. , 1. ) ;
//   SimpleLine l( ref , - dir  ) ;
//   printLine( l ) ;
//   std::cout << " --------------- " << std::endl ;
//   HepGeom::Translate3D   trans( ref ) ;
//   LCPoint3D ref2 = trans * ref  ;
//   SimpleLine l2( ref2 , dir  ) ;
//   printLine( l2 ) ;
//   std::cout << " --------------- " << std::endl ;
//   LCErrorMatrix err ;
//   LCVector3D p2 = l2.getPosition( 0. , &err ) ;
//   std::cout << "l2.getPosition( 0. , err ) ; " <<  p2 << " +/- " << err  << std::endl;


//---------------------------------------------------------------------------------------
  std::cout << "  - SimpleLine  :  " ;

  int nError = 0 ;

  for(int j = 0 ; j < 10000 ; j ++  ) {

    // random line l :
    LCPoint3D ref( drand() , drand() ,  drand() ) ;
    LCVector3D dir( drand() , drand() , drand() ) ;
    double s = drand() ;
    SimpleLine l( ref , dir  ) ;


    //---- check that closest point to random point on line is the point itself:
    LCPoint3D p = l.getPosition(s) ;
    double t = l.getPathAt( p ) ;

    if( ! isequal( t, s ) ) nError++ ;


    //---- check that random point in orthogonal plane through point p is p: 
    // plane through p orthogonal to line
    LCPlane3D plane( dir, p ) ;  
    LCPoint3D rp0( drand() , drand() ,  drand() ) ;

    // a random point on that plane:
    LCPoint3D rp1 = plane.point( rp0 ) ;

    double u = l.getPathAt( rp1 ) ;
  
    if( ! isequal( u, s ) ) nError++ ;
    

    //------  check intersection with random plane
    double sp = drand() ;
    // random plane through random point on the line:
    LCPoint3D np(  drand() , drand() ,  drand() ) ;
    LCPoint3D rp =  l.getPosition( sp ) ;
    LCPlane3D rplane( np , rp ) ;

//     std::cout << " random point : " << rp << std::endl ;

    bool ok = false ;
    double spc = l.getIntersectionWithPlane( rplane , ok ) ;

    if( ! isequal( sp, spc ) ) nError++ ;
    
    if( ! ok ) {
      nError ++ ;
      std::cout << " no intersection found : " <<  l.getPosition( sp ) 
		<< " --- " <<  l.getPosition( spc )  << std::endl ;
    }
  }
  
  if( nError != 0 ) 
    std::cout << " nError : " << nError << std::endl ;
  else 
    std::cout << "  Ok " << std::endl ;

  std::cout << std::endl ;

}
double drand() { return  1000. * double( rand() ) / RAND_MAX  ; } 
