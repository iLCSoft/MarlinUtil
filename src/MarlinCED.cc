#include "MarlinCED.h"



MarlinCED* MarlinCED::_me = 0 ;

MarlinCED* MarlinCED::instance() {
  
  if( _me == 0 )
    _me = new MarlinCED ;
  return _me ;
}

void MarlinCED::init( Processor* proc ) {
  
  if( instance()->_first == 0 ){
    
    instance()->_first = proc ;
    
    ced_client_init("localhost",7286);
    ced_register_elements();
  }

  instance()->_last = proc ;
}

void MarlinCED::newEvent( Processor* proc ) {
  if( proc == instance()->_first ) {

    ced_new_event(); 


    // FIXME: for now have hard codes Tesla geometry ....
    // here we would read some geometry description ....
    
    // from V.Morgunov : The Simplest geometry 
    static CED_GeoCylinder geoCylinders[] = {       // for TESLA Detector Geometry
      {    50.0,  3,  0.0, 5658.5, -5658.5, 0xff     }, // beam tube
      {   380.0,  4,  0.0, 2658.5, -2658.5, 0xff     }, // inner TPC
      {  1840.0,  8, 22.5, 2700.0, -2700.0, 0x7f7f1f }, // inner ECAL
      {  2045.7,  8, 22.5, 2700.0, -2700.0, 0x7f7f1f }, // outer ECAL
      {  2045.7,  8, 22.5, 101.00,  2820.0, 0x7f7f1f }, // endcap ECAL
      {  2045.7,  8, 22.5, 101.00, -3022.0, 0x7f7f1f }, // endcap ECAL
      {  3000.0, 16,  0.0, 2658.5, -2658.5, 0xcf00   }, // outer HCAL
      {  3000.0,  8, 22.5, 702.25,  2826.0, 0xcf00   }, // endcap HCAL
      {  3000.0,  8, 22.5, 702.25, -4230.5, 0xcf00   }, // endcap HCAL
    }; 
    /*
      static CED_GeoCylinder geoCylinders[] = {    // for Prototype
      {    180.0,  4,  45.0, 110.0, 0.0, 0xff }, // beam tube
      {    500.0,  4,  45.0, 250.0, 220., 0xff } // inner TPC
      };
    */
    ced_geocylinders(sizeof(geoCylinders)/sizeof(CED_GeoCylinder),geoCylinders);


  }
}

void MarlinCED::draw( Processor* proc ) {
  
  if( proc == instance()->_last ) {

    ced_draw_event();

    std::cout << "        [ Press return for next event ] " << std::endl ; 

    getchar();
  }
//   else
//     ced_send_event();





}
void MarlinCED::drawHelix(float b, float charge, float x, float y, float z,
			  float px, float py, float pz, int marker, int size, int col,
			  float rmin, float rmax, float zmax ) {
  
  double cFactor = 2.9979251e-4;
  
  double pt = hypot(px,py);  //sqrt(px*px + py*py);

  if( pt < 0.0000001 ) return ;

//   double p  = hypot(pt,pz);  //sqrt(px*px + py*py + pz*pz);
  
  double r =  pt / ( cFactor * b * std::abs( charge )  ) ;
//   double phi = std::atan2( x , y ) ; 

  double sign =  charge > 0 ? 1 : -1 ;
  sign = - sign  ; // FIXME: need to check the convention - but this works !?
 
  double phi = std::atan2( py , px ) + ( 2. + sign ) * M_PI / 2. ;

//   std::cout << "  atan2( py , px ) : " 
// 	    << std::atan2( py , px ) 
// 	    << " px,py,pz: " << px << ", " << py << ", " << pz 
// 	    << " charge:  " << charge 
// 	    << std::endl ;

  //center of helix
  double cx = x - ( sign * py * r / pt ) ;
  double cy = y + ( sign * px * r / pt ) ;
  double cz = z ; 
  
  double x1 =  x ; 
  double y1 =  y ;
  double z1 =  z ;
  double step = 0.005 ; 

  int nSteps  = 50 + int( 50. / pt ) ;

  for (int j = 0; j < nSteps ; j++) {

    double alpha = step*j ;  
  
    double x2 = cx + r * cos( phi + sign * alpha ) ;
    double y2 = cy + r * sin( phi + sign * alpha ) ;
    double z2 = cz + r * alpha * pz / pt ;
    
    double r_current  = hypot( x2, y2 ) ;



    if( std::abs(z2) > zmax || r_current > rmax  ) 
      break ;

    if( r_current > rmin ) 
      ced_line( x1, y1, z1, x2, y2, z2 , marker , size, col);	 
    
    x1 = x2;
    y1 = y2;
    z1 = z2;
  }
}

