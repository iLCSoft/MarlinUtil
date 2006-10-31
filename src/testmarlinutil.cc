
#include <iostream>
#include "SimpleLine.h"
#include "SimpleHelix.h"
#include <LCCylinder.h>
#include <LCLine3D.h>
#include <LCPlane3D.h>
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

    LCVector3D p = l.getPosition( s ) ;
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

//   LCVector3D ref( 0. , 2. , 3. ) ;
//   LCVector3D dir( 0. , 1. , 1. ) ;
//   SimpleLine l( ref , - dir  ) ;
//   printLine( l ) ;
//   std::cout << " --------------- " << std::endl ;
//   HepGeom::Translate3D   trans( ref ) ;
//   LCVector3D ref2 = trans * ref  ;
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
    LCVector3D ref( drand() , drand() ,  drand() ) ;
    LCVector3D dir( drand() , drand() , drand() ) ;
    double s = drand() ;
    SimpleLine l( ref , dir  ) ;


    //---- check that closest point to random point on line is the point itself:
    LCVector3D p = l.getPosition(s) ;
    double t = l.getPathAt( p ) ;

    if( ! isequal( t, s ) ) nError++ ;


    //---- check that random point in orthogonal plane through point p is p: 
    // plane through p orthogonal to line
    LCPlane3D plane( dir, p ) ;  
    LCVector3D rp0( drand() , drand() ,  drand() ) ;

    // a random point on that plane:
    LCVector3D rp1 = plane.projectPoint( rp0 ) ;

    double u = l.getPathAt( rp1 ) ;
  
    if( ! isequal( u, s ) ) nError++ ;
    

    //------  check intersection with random plane
    double sp = drand() ;
    // random plane through random point on the line:
    LCVector3D np(  drand() , drand() ,  drand() ) ;
    LCVector3D rp =  l.getPosition( sp ) ;
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


// ----- Testing the LCCylinder class

  LCVector3D cs(1,0,0), ce(3,0,0),cp(0.9,0.1,0);
  double cr = 1;
  LCCylinder cylinder(cs,ce,cr,false);

    int code;
    std::cout << "Projection: " << cylinder.projectPoint(cp,code) 
  	    << " " << code << std::endl;

    LCPlane3D p1(1,1,1,-1);
    LCVector3D po1(0,0,0), ax1(1,1,1);
    LCLine3D l1(po1,ax1);
    bool yesno;

    std::cout << "s inter: " << l1.intersectionWithPlane(p1,yesno) << " " << yesno << std::endl;
    std::cout << "point: " << l1.position(l1.intersectionWithPlane(p1,yesno)) << std::endl;
    std::cout << "project: " << p1.projectPoint(po1) << std::endl;

// ----- Testing the Helix implementation of the Trajectory:

  LCVector3D ref( 0. , 0. , 0. );

  double d0 = 10;
  double phi0 = 0;
  double omega = 0.1;
  double z0 = 0;
  double tanLambda = 1;

  SimpleHelix helix(d0,phi0,omega,z0,tanLambda,ref);
  helix.printProperties();
  LCVector3D nWall(1.,0.,0.);
  LCPlane3D myWall( nWall, 0. );
  bool pWallExists;
  double sWall = helix.getIntersectionWithPlane(myWall,pWallExists);
  std::cout << " s: " << sWall << " p: " << helix.getPosition(sWall) 
	    << " exists: " << pWallExists << std::endl;
  LCVector3D nWall2(1.,-1.,0.);
  LCVector3D mypoint(0.,0.,0.);
  LCPlane3D myWall2( nWall2, 1. );
  std::cout << "Ebene2: " << myWall2 << " Abstand " << myWall2.distance(mypoint) 
	    << " punkt: " << myWall2.projectPoint(mypoint) << std::endl;

  LCVector3D myp1(0.,0.,0.);
  LCVector3D myn1(1.,1.,1.);
  SimpleLine myl1(myp1,myn1);

  LCVector3D mypz1(3.,3.,3.);
  LCVector3D mypz2(9.,9.,9.);
  LCCylinder myz1(mypz1,mypz2,2.,true);

  bool bmyl1;
  double sl1 = myl1.getIntersectionWithCylinder(myz1,bmyl1);
  std::cout << "sl1: " << sl1 << " point: " << myl1.getPosition(sl1) << " exists: " << bmyl1 << std::endl;

// ----------- Helix intersection with cylinder

  double xd0 = 0;
  double xphi0 = 0;
  double xomega = 0.009;
  double xz0 = 0;
  double xtanLambda = 1;

  SimpleHelix helix2(xd0,xphi0,xomega,xz0,xtanLambda,ref);

  LCVector3D myxpz1(0.,0.,-400.);
  LCVector3D myxpz2(0.,0.,400.);
  LCCylinder myxz1(myxpz1,myxpz2,200.);

  bool bmyh1;
  double sh1 = helix2.getIntersectionWithCylinder(myxz1,bmyh1);
  std::cout << "sh1: " << sh1 << " point: " << helix2.getPosition(sh1) << " exists: " << bmyh1 << std::endl;

}

double drand() { return  1000. * double( rand() ) / RAND_MAX  ; } 
