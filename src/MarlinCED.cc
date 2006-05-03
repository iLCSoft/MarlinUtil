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


void MarlinCED::newEvent( Processor* proc , int modelID ) {
  if( proc == instance()->_first ) {

    ced_new_event(); 
    drawDetector(modelID);

  }
}


void MarlinCED::draw( Processor* proc , int waitForKeyboard ) {
  
  if( proc == instance()->_last ) {

    ced_draw_event();
   
    if ( waitForKeyboard == 1 ) {
      
      std::cout << "        [ Press return for next event ] " << std::endl ; 
      
      getchar();

    }
    
  }
  //   else
  //   ced_send_event();

}


void MarlinCED::drawHelix(float b, float charge, float x, float y, float z,
			  float px, float py, float pz, int marker, int size, int col,
			  float rmin, float rmax, float zmax ) {
  
  double cFactor = 2.9979251e-4;
  
  double pt = sqrt(px*px + py*py); // hypot(px,py)

  if( pt < 0.0000001 ) return ;

//   double p  = sqrt(px*px + py*py + pz*pz); // hypot(pt,pz);
  
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
  double step = 0.05 ; 

  int nSteps  = 50 + int( 150. / pt ) ;

  for (int j = 0; j < nSteps ; j++) {

    double alpha = step*j ;  
  
    double x2 = cx + r * cos( phi + sign * alpha ) ;
    double y2 = cy + r * sin( phi + sign * alpha ) ;
    double z2 = cz + r * alpha * pz / pt ;
    
    double r_current  = sqrt( x2*x2 + y2*y2); // hypot( x2, y2 ) 

    if( std::abs(z2) > zmax || r_current > rmax  ) 
      break ;

    if( r_current > rmin ) 
      ced_line( x1, y1, z1, x2, y2, z2 , marker , size, col);	 
    
    x1 = x2;
    y1 = y2;
    z1 = z2;
  }
}



void MarlinCED::drawSpike( float x0, float y0, float z0,float x1, float y1, float z1,unsigned int color, unsigned int layer ) {

  //    const float s0 = 0.;
  const float s1 = .92;
  const float s2 = .94;
  const float s3 = .96;
  const float s4 = .98;
  //    const float s5 = 1. ;
  
  float p0[3]= { x0, y0, z0 };
  float p1[3]= { (1-s1)*x0 + s1*x1 , (1-s1)*y0 + s1*y1 , (1-s1)*z0 + s1*z1 };
  float p2[3]= { (1-s2)*x0 + s2*x1 , (1-s2)*y0 + s2*y1 , (1-s2)*z0 + s2*z1 };
  float p3[3]= { (1-s3)*x0 + s3*x1 , (1-s3)*y0 + s3*y1 , (1-s3)*z0 + s3*z1 };
  float p4[3]= { (1-s4)*x0 + s4*x1 , (1-s4)*y0 + s4*y1 , (1-s4)*z0 + s4*z1 };
  float p5[3]= { x1, y1, z1 };
  
  unsigned int layty = layer << CED_LAYER_SHIFT ;
  
  ced_line( p0[0],p0[1],p0[2], p1[0],p1[1],p1[2], layty ,6 , color );
  ced_line( p1[0],p1[1],p1[2], p2[0],p2[1],p2[2], layty ,5 , color );
  ced_line( p2[0],p2[1],p2[2], p3[0],p3[1],p3[2], layty ,4 , color );
  ced_line( p3[0],p3[1],p3[2], p4[0],p4[1],p4[2], layty ,3 , color );
  ced_line( p4[0],p4[1],p4[2], p5[0],p5[1],p5[2], layty ,2 , color );
  
  ced_hit ( p0[0],p0[1],p0[2], CED_HIT_POINT | layer << CED_LAYER_SHIFT, 0, color );
  ced_hit ( p5[0],p5[1],p5[2], CED_HIT_POINT | layer << CED_LAYER_SHIFT, 0, color );
  
}


void MarlinCED::drawMCParticle(MCParticle* MCP, bool drawSimHits, LCEvent* event, int marker, int size, int color, int layer, double BField) {


  if ( MCP == 0 ) return;

  // FIXME scales the lenght of the symbolised momentum line
  double scale = 500.0;
 
  float x1 = (float)MCP->getVertex()[0];
  float y1 = (float)MCP->getVertex()[1];
  float z1 = (float)MCP->getVertex()[2];
  
  float x2 = (float)MCP->getEndpoint()[0];
  float y2 = (float)MCP->getEndpoint()[1];
  float z2 = (float)MCP->getEndpoint()[2];
	
  float p1 = (float)MCP->getMomentum()[0];
  float p2 = (float)MCP->getMomentum()[1];
  float p3 = (float)MCP->getMomentum()[2];
  float p  = sqrt(p1*p1 + p2*p2 + p3*p3);
  
	
  ced_hit(x1,y1,z1,0|layer<<CED_LAYER_SHIFT,4,0x00ff00);
  ced_hit(x2,y2,z2,1|layer<<CED_LAYER_SHIFT,10,0xff0000);
  ced_line(x1,y1,z1,x2,y2,z2,layer<<CED_LAYER_SHIFT,1,0xee9a49);
  ced_line(x1,y1,z1,scale*(p1/p),scale*(p2/p),scale*(p3/p),layer<<CED_LAYER_SHIFT,1,0xeed8ae);

  float charge = MCP->getCharge();

  if ( charge != 0 ) drawHelix(BField,charge,x1,y1,z1,p1,p2,p3,0,0,0xee9a49);

  if (drawSimHits) drawHitCollectionsByMCContribution(event,MCP,marker,size,color,layer);

}



void MarlinCED::drawSimTrackerHits(LCEvent* event, int marker, int size, int color, int layer) {

  drawHitCollectionsByType(event,LCIO::SIMTRACKERHIT,marker,size,color,layer);

}



void MarlinCED::drawSimCalorimeterHits(LCEvent* event, int marker, int size, int color, int layer) {

  drawHitCollectionsByType(event,LCIO::SIMCALORIMETERHIT,marker,size,color,layer);

}



void MarlinCED::drawSimHits(LCEvent* event, int marker, int size, int color, int layer) {

  drawHitCollectionsByType(event,LCIO::SIMTRACKERHIT,marker,size,color,layer);
  drawHitCollectionsByType(event,LCIO::SIMCALORIMETERHIT,marker,size,color,layer);

}



void MarlinCED::drawTrackerHits(LCEvent* event, int marker, int size, int color, int layer) {

  drawHitCollectionsByType(event,LCIO::TRACKERHIT,marker,size,color,layer);

}



void MarlinCED::drawCalorimeterHits(LCEvent* event, int marker, int size, int color, int layer) {

  drawHitCollectionsByType(event,LCIO::CALORIMETERHIT,marker,size,color,layer);

}


void MarlinCED::drawHits(LCEvent* event, int marker, int size, int color, int layer) {

  drawHitCollectionsByType(event,LCIO::TRACKERHIT,marker,size,color,layer);
  drawHitCollectionsByType(event,LCIO::CALORIMETERHIT,marker,size,color,layer);

}


void MarlinCED::drawTrack(Track* track, int marker, int size, int color, int layer) {

  const TrackerHitVec trackerHits = track->getTrackerHits();
  
  drawObjectsWithPosition(trackerHits.begin(),trackerHits.end(),marker,size,color,layer);

}


void MarlinCED::drawCluster(Cluster* cluster, int marker, int size, int color, int layer) {

  const CalorimeterHitVec clusterHits = cluster->getCalorimeterHits();
  
  drawObjectsWithPosition(clusterHits.begin(),clusterHits.end(),marker,size,color,layer);

}


void MarlinCED::drawRecoParticle(ReconstructedParticle* reco, int marker, int size, int color, int layer) {


  if ( reco == 0 ) return;

  unsigned int NofTracks   = reco->getTracks().size();
  unsigned int NofClusters = reco->getClusters().size();


  // FIXME: A track might be composed of several other tracks => insert a second loop
  for (unsigned int i = 0; i < NofTracks; ++i) {
    
    Track* track = reco->getTracks()[i];

    drawTrack(track,marker,size,color,layer);
    
  }

  for (unsigned int i = 0; i < NofClusters; ++i) {
    
    Cluster* cluster = reco->getClusters()[i];
    drawCluster(cluster,marker,size,color,layer);

  }

}
