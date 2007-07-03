
#ifdef USE_CLHEP  // only if CLHEP is available !
#include "CLHEP/HepPDT/ParticleID.hh"
#endif

#include "MarlinCED.h"
#include "gear/GearMgr.h" 
#include <gear/TPCParameters.h>
#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>
#include <gear/PadRowLayout2D.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>
#include <LCGeometryTypes.h>
#include "ced_cli.h"

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
    //drawDetector(modelID);
    if ( modelID == 99999 ) drawGEARTelescope();
    else drawGEARDetector();

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
			  float px, float py, float pz, int marker, int size, unsigned int col,
			  float rmin, float rmax, float zmax) {

  // FIXME : check for zmin as well, i.e. cylindrical coordinates

  
  double cFactor = 2.9979251e-4;
    
  double pt = sqrt(px*px + py*py); // hypot(px,py)
  double pl = pz;
  double absp = sqrt(px*px + py*py + pz*pz);

  // debug
  // std::cout << "|p| = " << absp << "  " << "pt = " << pt << "  " << "pl = " << pl << std::endl;



  // FIXME: use a parameter for this cut or better this should be a function of the B field, charge and momentum 2006/07/04 OW    
  if ( (pt >= 0.001) && (pt <= 40.0) ) {
  
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
    double step = 0.01;  // initial 0.05
    
    int nSteps  = 50 + int( 150. / pt ) ;
    
    for (int j = 0; j < nSteps ; j++) {
      
      double alpha = step*j ;  
      
      double x2 = cx + r * cos( phi + sign * alpha ) ;
      double y2 = cy + r * sin( phi + sign * alpha ) ;
      double z2 = cz + r * alpha * pz / pt ;
      
      double r_current  = sqrt( x2*x2 + y2*y2); // hypot( x2, y2 ) 

      
      // debug
      /*    
      std::cout << "step = " << step << "  " << "nSteps = " << nSteps << "  " << "alpha = " << alpha << "  " << "|z2| = " << std::abs(z2) << "  " << "zmax = " << zmax 
		<< "  " << "r_current = " << r_current << "  " << "rmax = " << rmax << "  " << "rmin = " << rmin << std::endl;
      */


      if( std::abs(z2) > zmax || r_current > rmax  ) 
	break ;
    
      if( r_current > rmin ) 
	ced_line( x1, y1, z1, x2, y2, z2 , marker , size, col);	 
      
      x1 = x2;
      y1 = y2;
      z1 = z2;
    }

  }
  else { //if pt < 0.001, draw a straight line from start point to the intersection point of the momentum extrapolation with the outer cylinder shell (rmax,zmax)
      

    float absP =sqrt(px*px + py*py + pz*pz);
    float k = 0.0;
    float kr = 0.0;
    float kz = 0.0;
    float summand = 0.0;
    float radicant = 0.0;

    // find intersection with rmax
    summand = (-1)*( absP*(px*x + py*y)/(pow(px,2) + pow(py,2)) );
    radicant = summand*summand - ( (pow(absP,2)*(pow(x,2)+pow(y,2)-pow(rmax,2)))/(pow(px,2) + pow(py,2)) );
    
    if (radicant < 0) {

      std::cout << "Error in 'MarlinCED::drawHelix()': Startpoint beyond (rmax,zmax)" << std::endl;
      return;

    }
    
    kr = summand + sqrt(radicant);
    kz = ((zmax-z)*absP)/pz;

    if ( kr >= kz ) k = kr;
    else k = kz;

    if (k < 0.0) {
      
      std::cout << "Error in 'MarlinCED::drawHelix()': No intersection point with the outer cylinder shell (rmax,zmax) found" << std::endl;
      return;

    }

    float xEnd = x + (k*px)/absP;
    float yEnd = y + (k*py)/absP;
    float zEnd = z + (k*pz)/absP;


    // debug
    std::cout << "start point = " << "( " << x << ", " << y << ", " << z << " )" << "  " 
	      << "end point = " << "( " << xEnd << ", " << yEnd << ", " << zEnd << " )" << std::endl;
    

    ced_line( x, y, z, xEnd, yEnd, zEnd , marker , size, col);	 
    
  }

} 

void MarlinCED::drawTrajectory(const Trajectory* t, const int marker,
			       const int size, const unsigned int col,
			       const float rmin, const float rmax,
			       const float zmax) 
{
  if (rmax <= rmin || zmax == 0 ) return;
  double stepSize = 5.0; // initial 0.05
  double nStepKill = int(2*3.14*rmax/stepSize);

  double rmax2 = rmax*rmax, rmin2 = rmin*rmin,xmagxy2;
  double s = - stepSize;
  int nz = 0;
  LCVector3D x,xold;
  x = t->getPosition(s);
  for (;;)
    {
      s += stepSize;
      xold = x ;
      x = t->getPosition(s);
      if (x.z() == xold.z()) 
	{
	  ++nz;
	  if (nz > nStepKill) break;
	}
      xmagxy2 = x.x()*x.x() + x.y()*x.y() ;
      if (fabs(x.z()) > zmax) break;
      if (xmagxy2 < rmin2) continue;
      if (xmagxy2 > rmax2) break;
      ced_line( xold.x(), xold.y(), xold.z(), x.x(), x.y(), x.z(), 
		marker, size, col);
    }
}

void MarlinCED::drawSpike( float x0, float y0, float z0,float x1, float y1, float z1, unsigned int color, unsigned int layer ) {

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


void MarlinCED::drawMCParticle(MCParticle* MCP, bool drawSimHits, LCEvent* event, int marker, int size, unsigned int color, unsigned int layer, double bField,
			       double rmin, double zmin, double rmax, double zmax, bool drawOnDifferentLayers) {


  if ( MCP == 0 ) return;
 
  double x1 = MCP->getVertex()[0];
  double y1 = MCP->getVertex()[1];
  double z1 = MCP->getVertex()[2];
  
  double x2 = MCP->getEndpoint()[0];
  double y2 = MCP->getEndpoint()[1];
  double z2 = MCP->getEndpoint()[2];
	
  double p1 = MCP->getMomentum()[0];
  double p2 = MCP->getMomentum()[1];
  double p3 = MCP->getMomentum()[2];
  
  float charge = 0.0;

  unsigned int l = 0;
  

  #ifdef USE_CLHEP  // only if CLHEP is available !
  charge =   HepPDT::ParticleID( MCP->getPDG() ).threeCharge() / 3.00;
  #else
  charge = MCP->getCharge();
  #endif

  // debug
  // std::cout << bField << "  " << charge << "  " << x1 << "  " << y1 << "  " << z1 << "  " << p1 << "  " << p2 << "  " << p3 << "  " << color << std::endl;

  bool isCharged       = charge != 0.0;    
  bool isNeutrino      = (abs(MCP->getPDG())==12) || (abs(MCP->getPDG())==14) || (abs(MCP->getPDG())==16);
  bool isBackscattered = MCP->isBackscatter();
    

  // charged MC Particles are displayed on layer and their SimHits optionally on (layer + 10)
  if (isCharged && !isNeutrino && !isBackscattered) {

    drawHelix(bField,charge,x1,y1,z1,p1,p2,p3, marker | ( layer << CED_LAYER_SHIFT ), size, color, (float)rmin, (float)rmax, (float)zmax);
    if (drawSimHits) drawHitCollectionsByMCContribution(event,MCP,marker,size+2,color,layer+10);

  }

  // neutral MC Particles are displayed on (layer+1) and their SimHits optionally on (layer + 11)
  else if (!isCharged && !isNeutrino && !isBackscattered) {

    if (drawOnDifferentLayers) l = layer+1;
    else l = layer;

    ced_line(x1,y1,z1,x2,y2,z2, marker | ( l << CED_LAYER_SHIFT ), size, color);

    if (drawSimHits) {
      if (drawOnDifferentLayers) l = layer+11;
      else l = layer+10;
      drawHitCollectionsByMCContribution(event,MCP,marker,size+2,color,l);
    }

  }
  // backscatterd charged particles and neutrinos are displayed on (layer+2) and their SimHits optionally on (layer + 12)
  else if (isCharged && !isNeutrino && isBackscattered) {

    if (drawOnDifferentLayers) l = layer+2;
    else l = layer;

    drawHelix(bField,charge,x1,y1,z1,p1,p2,p3, marker | ( l << CED_LAYER_SHIFT ), size, color, (float)rmin, (float)rmax, (float)zmax);

    if (drawSimHits) {
      if (drawOnDifferentLayers) l = layer+12;
      else l = layer+10;
      drawHitCollectionsByMCContribution(event,MCP,marker,size+2,color,l);
    }

  }
  // backscatterd charged particles and neutrinos are displayed on (layer+2) and their SimHits optionally on (layer + 12)
  else if (!isCharged && isNeutrino || isBackscattered) {
  
    if (drawOnDifferentLayers) l = layer+2;
    else l = layer;

    ced_line(x1,y1,z1,x2,y2,z2, marker | ( l << CED_LAYER_SHIFT ), size, color);

    if (drawSimHits) {
      if (drawOnDifferentLayers) l = layer+12;
      else l = layer+10;
      drawHitCollectionsByMCContribution(event,MCP,marker,size+2,color,l);
    }

  }
  else std::cout << "This MC Particle has not been displayed : id = " << MCP->id() << "  " << "PDG Code : " << MCP->getPDG() << std::endl;
  
}



void MarlinCED::drawMCParticleTree(LCEvent* event, std::string colNameMC, double energyCut,  double bField, double rIn, double zIn, double rOut, double zOut) {

  
  try {
    
    std::vector< std::string >::const_iterator i;
    const std::vector< std::string >* ColNames = event->getCollectionNames();
    
    for( i = ColNames->begin() ; i != ColNames->end() ; i++) {
      
      LCCollection* col = event->getCollection( *i ) ;
      
      if ( (col->getTypeName() == LCIO::MCPARTICLE) && (*i == colNameMC) ) {
	
	int nMCP = col->getNumberOfElements();

	for(int j=0; j<nMCP; ++j){

	  MCParticle* mcP = dynamic_cast<MCParticle*> ( col->getElementAt( j ) ) ;

	  double energy = mcP->getEnergy();

	  if ( energy >= energyCut ) {

	    const double* rStart = mcP->getVertex();
	    const double* rEnd   = mcP->getEndpoint();

	    double r2Start_rp = rStart[0]*rStart[0] + rStart[1]*rStart[1];
	    double zStart     = rStart[2];
	    
	    double r2End_rp   = rEnd[0]*rEnd[0] + rEnd[1]*rEnd[1];
	    double zEnd       = rEnd[2];
	    
	    // completely inside innermost detector
	    bool withinInnerDetector = (r2Start_rp < rIn*rIn) && (r2End_rp < rIn*rIn) && (zStart < zIn) && (zEnd < zIn);
	    // completely beyond outermost detector not taking into account the calorimeters (because of showering)
	    bool beyondOuterDetector = (r2Start_rp > rOut*rOut) && (r2End_rp > rOut*rOut) && (zStart > zOut) && (zEnd > zOut);
	  
	    if (!withinInnerDetector && !beyondOuterDetector) {

	      unsigned int color = MarlinDrawUtil::getColor(mcP->getPDG());   
	  
	      MarlinCED::drawMCParticle(mcP,true,event,2,1,color,1,bField,rIn,zIn,rOut,zOut, true);
	     
	    }

	  }

	}

      }

    }

  }
  catch(DataNotAvailableException &e){

    std::cout << "no valid MC collection in event." << std::endl ;

  };

}



void MarlinCED::drawSimTrackerHits(LCEvent* event, int marker, int size, unsigned int color, unsigned int layer) {

  drawHitCollectionsByType(event,LCIO::SIMTRACKERHIT,marker,size,color,layer);

}



void MarlinCED::drawSimCalorimeterHits(LCEvent* event, int marker, int size, unsigned int color, unsigned int layer) {

  drawHitCollectionsByType(event,LCIO::SIMCALORIMETERHIT,marker,size,color,layer);

}



void MarlinCED::drawSimHits(LCEvent* event, int marker, int size, unsigned int color, unsigned int layer) {

  drawHitCollectionsByType(event,LCIO::SIMTRACKERHIT,marker,size,color,layer);
  drawHitCollectionsByType(event,LCIO::SIMCALORIMETERHIT,marker,size,color,layer);

}



void MarlinCED::drawTrackerHits(LCEvent* event, int marker, int size, unsigned int color, unsigned int layer) {

  drawHitCollectionsByType(event,LCIO::TRACKERHIT,marker,size,color,layer);

}



void MarlinCED::drawCalorimeterHits(LCEvent* event, int marker, int size, unsigned int color, unsigned int layer) {

  drawHitCollectionsByType(event,LCIO::CALORIMETERHIT,marker,size,color,layer);

}


void MarlinCED::drawHits(LCEvent* event, int marker, int size, unsigned int color, unsigned int layer) {

  drawHitCollectionsByType(event,LCIO::TRACKERHIT,marker,size,color,layer);
  drawHitCollectionsByType(event,LCIO::CALORIMETERHIT,marker,size,color,layer);

}


void MarlinCED::drawTrack(Track* track, int marker, int size, unsigned int color, unsigned int layer) {

  const TrackerHitVec trackerHits = track->getTrackerHits();
  
  drawObjectsWithPosition(trackerHits.begin(),trackerHits.end(),marker,size,color,layer);

}


void MarlinCED::drawCluster(Cluster* cluster, int marker, int size, unsigned int color, unsigned int layer) {

  const CalorimeterHitVec clusterHits = cluster->getCalorimeterHits();
  
  drawObjectsWithPosition(clusterHits.begin(),clusterHits.end(),marker,size,color,layer);

}


void MarlinCED::drawClusterImpl(const ClusterImpl* cluster, int marker, int size, unsigned int color, unsigned int layer) {

  const CalorimeterHitVec clusterHits = cluster->getCalorimeterHits();
  
  drawObjectsWithPosition(clusterHits.begin(),clusterHits.end(),marker,size,color,layer);

}


void MarlinCED::drawRecoParticle(ReconstructedParticle* reco, int marker, int size, unsigned int color, unsigned int layer) {


  if ( reco == 0 ) return;

  unsigned int NofTracks   = reco->getTracks().size();
  unsigned int NofClusters = reco->getClusters().size();


  // FIXME: A track might be composed of several other tracks => insert a second, third, ... loop
  for (unsigned int i = 0; i < NofTracks; ++i) {
    
    Track* track = reco->getTracks()[i];
    drawTrack(track,marker,size,color,layer);
    
  }

  for (unsigned int i = 0; i < NofClusters; ++i) {

    Cluster* cluster = reco->getClusters()[i];
    drawCluster(cluster,marker,size,color,layer);

  }

}


void MarlinCED::drawGEARTelescope() {
  
  const gear::SiPlanesParameters&  siPlanesParameters  = Global::GEAR->getSiPlanesParameters();
  const gear::SiPlanesLayerLayout& siPlanesLayerLayout = siPlanesParameters.getSiPlanesLayerLayout();
  
  double * sizes  = new double[3];
  double * center = new double[3];
  unsigned int color = 0xFFFFFF;

  for ( int iLayer = 0 ; iLayer < siPlanesLayerLayout.getNLayers() ; iLayer++ ) {
    center[0] = siPlanesLayerLayout.getSensitivePositionX(iLayer);
    center[1] = siPlanesLayerLayout.getSensitivePositionY(iLayer);
    center[2] = siPlanesLayerLayout.getSensitivePositionZ(iLayer);
    sizes[0]  = siPlanesLayerLayout.getSensitiveSizeX(iLayer);
    sizes[1]  = siPlanesLayerLayout.getSensitiveSizeY(iLayer);
    sizes[2]  = siPlanesLayerLayout.getSensitiveThickness(iLayer) ;
    ced_geobox( sizes, center, color );
  }
  delete [] center;
  delete [] sizes;
}



void MarlinCED::drawGEARDetector() {
  //
  // code from V.Morgunov, DESY
  //

  const gear::TPCParameters&  pTPC      = Global::GEAR->getTPCParameters();
  const gear::PadRowLayout2D& padLayout = pTPC.getPadLayout();
  const gear::DoubleVec&      planeExt  = padLayout.getPlaneExtent();
  float r_min_tpc = planeExt[0];
  float r_max_tpc = planeExt[1];
  float z_max_tpc = pTPC.getMaxDriftLength();
  // ===================================================================
  //               Get gear file parameters  for calorimeters
  // ===================================================================

  const gear::CalorimeterParameters& pECAL_B = 
    Global::GEAR->getEcalBarrelParameters();

  float r_min_ecal_bar = pECAL_B.getExtent()[0];
  float r_max_ecal_bar = pECAL_B.getExtent()[1]; // wrong, takes the ECAL_1 structure only
  float z_min_ecal_bar = pECAL_B.getExtent()[2];
  float z_max_ecal_bar = pECAL_B.getExtent()[3];
  
  
  const gear::CalorimeterParameters& pECAL_E = 
    Global::GEAR->getEcalEndcapParameters();
  float r_min_ecal_ecap = pECAL_E.getExtent()[0];
  float r_max_ecal_ecap = pECAL_E.getExtent()[1];
  float z_min_ecal_ecap = pECAL_E.getExtent()[2];
  float z_max_ecal_ecap = pECAL_E.getExtent()[3]; // wrong, takes the ECAL_1 structure only
  
  
  const gear::CalorimeterParameters& pHCAL_B = Global::GEAR->getHcalBarrelParameters();
  //  _innerHcalRadius = float(pHcalBarrel.getExtent()[0]);
  float r_min_hcal_bar = pHCAL_B.getExtent()[0];
  float r_max_hcal_bar = pHCAL_B.getExtent()[1];
  float z_min_hcal_bar = pHCAL_B.getExtent()[2];
  float  z_max_hcal_bar = pHCAL_B.getExtent()[3];
  

  const gear::CalorimeterParameters& pHCAL_E = Global::GEAR->getHcalEndcapParameters();
  float r_min_hcal_ecap = pHCAL_E.getExtent()[0];
  float r_max_hcal_ecap = pHCAL_E.getExtent()[1];
  float z_min_hcal_ecap = pHCAL_E.getExtent()[2];
  float z_max_hcal_ecap = pHCAL_E.getExtent()[3];
  
  /*
  std::cout<<"  ++++++++++++ Parameters will be used for GeoCylinders to draw geometry   +++++++++++++++ "<<std::endl;
  std::cout<<" TPC       "<< r_min_tpc      <<"; "<<r_max_tpc      <<"; "<<z_max_tpc<<std::endl;
  std::cout<<" Ecal Bar  "<< r_min_ecal_bar <<"; "<<r_max_ecal_bar <<"; "<<z_min_ecal_bar <<"; "<<z_max_ecal_bar<<std::endl;
  std::cout<<" Hcal Bar  "<< r_min_hcal_bar <<"; "<<r_max_hcal_bar <<"; "<<z_min_hcal_bar <<"; "<<z_max_hcal_bar<<std::endl;
  std::cout<<" Ecal Ecap "<< r_min_ecal_ecap<<"; "<<r_max_ecal_ecap<<"; "<<z_min_ecal_ecap<<"; "<<z_max_ecal_ecap<<std::endl;
  std::cout<<" Hcal Ecap "<< r_min_hcal_ecap<<"; "<<r_max_hcal_ecap<<"; "<<z_min_hcal_ecap<<"; "<<z_max_hcal_ecap<<std::endl;
  */  

  //                   It was printed for LDC00Sc +0 +0
  //   ++++++++++++ Parameters will be used for GeoCylinders to draw geometry   +++++++++++++++
  //  TPC        386;   1626;           2500
  //  Ecal Bar  1700;   1840;     0;    2730
  //  Hcal Bar  1910;   2890;     0;    2800
  //  Ecal Ecap  300;   1884;  2829.5;  2969.5
  //  Hcal Ecap  300;   2920;  3018;    3998
  
  //                   It was printed for LDC00Sc +20 +20
  //  ++++++++++++ Parameters will be used for GeoCylinders to draw geometry   +++++++++++++++
  // TPC         386;   1826;           2700
  // Ecal Bar   1900;   2040;     0;    2930
  // Hcal Bar   2110;   3090;     0;    3000
  // Ecal Ecap   300;   2084;  3029.5;  3169.5
  // Hcal Ecap   300;   3120;  3218;    4198
  // ===============================================================================
  // To convert inner radius of polygone to its outer radius to draw correct geocylinders
  float Cos4  = cos(M_PI/4.0);
  float Cos8  = cos(M_PI/8.0);
  float Cos16 = cos(M_PI/16.);
  
  float r_inn_ecal_bar     = r_min_ecal_bar/Cos8;   // convertion of  inner radius of polygone to its outer radius
  float r_out_ecal_bar     = (r_max_ecal_bar+56.0)/Cos8 ; // + second ECAL structure thickness = 56 mm
  float r_inn_ecal_ecap    = r_min_ecal_ecap/Cos4;
  float r_out_ecal_ecap    = r_max_ecal_ecap/Cos8;
  float thick_ecal_ecap    = 0.5*(z_max_ecal_ecap + 56.0 - z_min_ecal_ecap); // + second ECAL structure thickness = 56 mm
  float shift_ecal_z_plus  = z_min_ecal_ecap;
  float shift_ecal_z_minus = z_min_ecal_ecap + 2.0*thick_ecal_ecap;
  
  float r_inn_hcal_bar     = r_min_hcal_bar/Cos8;
  float r_out_hcal_bar     = r_max_hcal_bar/Cos16;
  
  float r_inn_hcal_ecap    = r_min_hcal_ecap/Cos4;
  float r_out_hcal_ecap    = r_max_hcal_ecap/Cos8;
  float thick_hcal_ecap    = 0.5*(z_max_hcal_ecap - z_min_hcal_ecap + 20.0); // +20 by hand to see hits inside
  float shift_hcal_z_plus  = z_min_hcal_ecap;
  float shift_hcal_z_minus = z_min_hcal_ecap + 2.0*thick_hcal_ecap;
  // ===============================================================================
  /**
     GeoCylinder
     typedef struct {
     float d;           // radius (outer) of polygone
     unsigned  sides;   // poligon order
     float rotate;      // angle degree
     float z;           // 1/2 length
     float shift;       // in z
     unsigned color;
     } CED_GeoCylinder;
  */
  static CED_GeoCylinder geoCylindersANY[] = {       // for ANY Detector Geometry
    { r_min_tpc,        4, 45.0, z_max_tpc, -z_max_tpc, 0xff      }, // inner TPC
    { r_inn_ecal_bar ,  8, 22.5, z_max_ecal_bar,  -z_max_ecal_bar,     0x7f7f1f  }, // inner ECAL Barrel
    { r_out_ecal_bar ,  8, 22.5, z_max_ecal_bar,  -z_max_ecal_bar,     0x7f7f1f  }, // outer ECAL Barrel
    //      { r_inn_ecal_ecap,  4,  0.0, thick_ecal_ecap,  shift_ecal_z_plus,  0x7f7f1f  }, // inner endcap ECAL +Z
    //      { r_inn_ecal_ecap,  4,  0.0, thick_ecal_ecap, -shift_ecal_z_minus, 0x7f7f1f  }, // inner endcap ECAL -Z
    { r_out_ecal_ecap,  8, 22.5, thick_ecal_ecap,  shift_ecal_z_plus,  0x7f7f1f  }, // outer endcap ECAL +Z
    { r_out_ecal_ecap,  8, 22.5, thick_ecal_ecap, -shift_ecal_z_minus, 0x7f7f1f  }, // outer endcap ECAL -Z
    { r_inn_hcal_bar ,  8, 22.5, z_max_hcal_bar,  -z_max_hcal_bar,     0x00cf00  }, // inner HCAL Barrel
    { r_out_hcal_bar , 16,  0.0, z_max_hcal_bar,  -z_max_hcal_bar,     0x00cf00  }, // outer HCAL Barrel
    //      { r_inn_hcal_ecap,  4,  0.0, thick_hcal_ecap,  shift_hcal_z_plus,  0x00cf00  }, // inner endcap HCAL +Z
    //      { r_inn_hcal_ecap,  4,  0.0, thick_hcal_ecap, -shift_hcal_z_minus, 0x00cf00  }, // inner endcap HCAL -Z
    { r_out_hcal_ecap,  8,  0.0, thick_hcal_ecap,  shift_hcal_z_plus,  0x00cf00  }, // outer endcap HCAL +Z
    { r_out_hcal_ecap,  8,  0.0, thick_hcal_ecap, -shift_hcal_z_minus, 0x00cf00  }  // outer endcap HCAL -Z
  };
  
  ced_geocylinders(sizeof(geoCylindersANY)/sizeof(CED_GeoCylinder),geoCylindersANY);
  
  
  
  // ===============================================================================
  // ===============================================================================
  
  
}
