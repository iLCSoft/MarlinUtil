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

#include <gear/VXDLayerLayout.h>
#include <gear/VXDParameters.h>
#include <gear/GEAR.h>
#include <gear/BField.h>
#include <gearimpl/Vector3D.h>

//SJA:FIXED:added to make gcc4.3 compliant
#include <cstdlib>

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


void MarlinCED::newEvent( Processor* proc , int modelID, LCEvent* evt) {
  if( proc == instance()->_first ) {
    ced_new_event(); 
    if(evt!=0) 
        instance()->_currEvent = evt;
//     //drawDetector(modelID);
//     if ( modelID == 99999 ) drawGEARTelescope();
//     else drawGEARDetector();

    switch(modelID) {

    case 0:
      drawGEARDetector();
      break ;

    case 99999:
      drawGEARTelescope();
      break;

    default:
      // don't draw anything
      break;
    }
    
    
  }
}

MCParticle* MarlinCED::getMCParticleFromID(int partID, LCEvent* evt) {
    //SM-H This could probably be improved...
    LCCollection * mcpcol ;
    if(partID>0) {
        if(evt!=0) {
            try {
                mcpcol = evt->getCollection("MCParticle");
                int nelem = mcpcol->getNumberOfElements();
                for (int ielem(0); ielem < nelem; ++ielem) {
                    MCParticle * mcp = 
                    dynamic_cast<MCParticle*>(mcpcol->getElementAt(ielem));
                    if(mcp->id()==partID) {
                     //Found particle with appropriate ID
                         return mcp;
                    }
                }             
            }
            catch(DataNotAvailableException &e) {
                 std::cout << "Data not available" << std::endl;
            }
        }
    }
    //If particle is not in collection, return a null pointer
    return NULL;
}

void MarlinCED::printMCParticle(MCParticle* part, int daughterIndent, int motherIndent) {
    //Define width of numbers, and decimal precision 
    const int width = 10;
    const int prec = 2;
    std::string daughterIndent_str;
    std::string motherIndent_str;
    for (int i = 0; i < daughterIndent; i++) {
        daughterIndent_str.append("->");
    }
    for (int i = 0; i < motherIndent; i++) {
        motherIndent_str.append("<-");
    }
    streamlog_out(MESSAGE) << std::endl
        << motherIndent_str
        << daughterIndent_str
        <<  "[   id   ]PDG    |  px      ,  py      ,  pz      |  energy  |gen|[simstat]|  vertex x,      y   ,    z     |  endpoint x,      y  ,   z       |    mass  |  charge  | [parents] - [daughters] |"    
        << std::endl;

    std::setiosflags(std::ios::fixed);
    streamlog_out(MESSAGE) << motherIndent_str;
    streamlog_out(MESSAGE) << daughterIndent_str;
    streamlog_out(MESSAGE) << std::setfill(' ');
    streamlog_out(MESSAGE) << "[" << std::setw(8) << part->id() << "]"; 
    streamlog_out(MESSAGE) << std::setw(7) << part->getPDG() << "|";
    streamlog_out(MESSAGE) << std::scientific << std::setw(width) << std::setprecision(prec)
        << std::setw(width) << part->getMomentum()[0] << ","
        << std::setw(width) << part->getMomentum()[1] << ","
        << std::setw(width) << part->getMomentum()[2] << "|";
    streamlog_out(MESSAGE) << std::scientific << std::setw(width) << std::setprecision(prec)
        << part->getEnergy() << "|";
    streamlog_out(MESSAGE) << part->getGeneratorStatus() << "  |";
    streamlog_out(MESSAGE) << LCTOOLS::getSimulatorStatusString( part ).c_str() << "|";
    streamlog_out(MESSAGE) << std::scientific << std::setw(width) << std::setprecision(prec)
        << std::setw(width) << part->getVertex()[0] << ","
        << std::setw(width) << part->getVertex()[1] << ","
        << std::setw(width) << part->getVertex()[2] << "|";
    streamlog_out(MESSAGE) << std::scientific << std::setw(width) << std::setprecision(prec)
        << std::setw(width) << part->getEndpoint()[0] << ", "
        << std::setw(width) << part->getEndpoint()[1] << ", "
        << std::setw(width) << part->getEndpoint()[2] << "|";
    streamlog_out(MESSAGE) << std::scientific << std::setw(width) << std::setprecision(prec)
        << part->getMass() << "|";
    streamlog_out(MESSAGE) << std::scientific << std::setw(width) << std::setprecision(prec)
        << part->getCharge() << "|";
    streamlog_out(MESSAGE) << "    [" ;

    streamlog_out(MESSAGE) << part->getParents().size();
    streamlog_out(MESSAGE) << "]    -    [" ;
    streamlog_out(MESSAGE) << part->getDaughters().size();
    streamlog_out(MESSAGE) << "] " << std::endl ;
    streamlog_out(MESSAGE) << std::endl 
        << "-------------------------------------------------------------------------------- " 
        << std::endl;
}

//SM-H: Loop recursively through each level of the family hierarchy, up to some specified limit
void MarlinCED::printMCFamily(MCParticle* part, unsigned int daughterBranches, unsigned int motherBranches,
                               unsigned int daughterIndent, unsigned int motherIndent) {

    printMCParticle(part, daughterIndent, motherIndent);
    for (unsigned int i = 0; i < daughterBranches; i++) {
        for (unsigned int j = 0; j < part->getDaughters().size(); j++) {
            MCParticle* part_daughter = part->getDaughters()[j];
            printMCFamily(part_daughter, daughterBranches-1, 0, daughterIndent+1, 0);
        }
    }
    for (unsigned int i = 0; i < motherBranches; i++) {
        for (unsigned int j = 0; j < part->getParents().size(); j++) {
            MCParticle* part_mother = part->getParents()[j];
            printMCFamily(part_mother, 0, motherBranches-1, 0, motherIndent+1);
        }
    }
    return;
}

//SM-H: Loop recursively through each level of the family hierarchy, up to some specified limit
//Draw each particle in the hierarchy
void MarlinCED::printAndDrawMCFamily(MCParticle* part, LCEvent* evt, unsigned int daughterBranches, 
                                      unsigned int motherBranches, unsigned int daughterIndent, unsigned int motherIndent) {

    double bField = Global::GEAR->getBField().at(gear::Vector3D(0,0,0)).z() ;
    const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;

    //int colour = 0xff00ff;
    int colour = abs(0xff00ff-abs((daughterIndent-motherIndent)*128));
    double endpoint_r = gearTPC.getPadLayout().getPlaneExtent()[1];
    double endpoint_z = gearTPC.getMaxDriftLength();
    if(gearTPC.getPadLayout().getPlaneExtent()[1] > sqrt(part->getEndpoint()[0]*part->getEndpoint()[0] + part->getEndpoint()[1]*part->getEndpoint()[1]))
        endpoint_r = sqrt(part->getEndpoint()[0]*part->getEndpoint()[0] + part->getEndpoint()[1]*part->getEndpoint()[1]);
    if(gearTPC.getMaxDriftLength() > part->getEndpoint()[2])
        endpoint_z = part->getEndpoint()[2];
    if(part->getPDG() < 81 || part->getPDG() > 100) {
        //Ignore internal MC particles and neutrals

        //MarlinCED::newEvent( this, _detModel ) ;
        int layer = daughterIndent + 1;
        MarlinCED::drawMCParticle(part, false, evt, 1, 1, colour, layer, bField, 0, 0,endpoint_r, 
                endpoint_z, false);
        printMCParticle(part, daughterIndent, motherIndent);
        //MarlinCED::draw( this, _waitForKeyboard ) ;

    }
    for (unsigned int i = 0; i < daughterBranches; i++) {
        for (unsigned int j = 0; j < part->getDaughters().size(); j++) {
            std::cout << "Daughter of " << part->id() << std::endl;
            MCParticle* part_daughter = part->getDaughters()[j];
            printAndDrawMCFamily(part_daughter, evt, daughterBranches-1, 0, daughterIndent+1, 0);
        }
    }
    for (unsigned int i = 0; i < motherBranches; i++) {
        for (unsigned int j = 0; j < part->getParents().size(); j++) {
            std::cout << "Mother of " <<  part->id()<< std::endl;
            MCParticle* part_mother = part->getParents()[j];
            printAndDrawMCFamily(part_mother, evt, 0, motherBranches-1, 0, motherIndent+1);
        }
    }

    return;
}

void MarlinCED::draw( Processor* proc , int waitForKeyboard ) {
  
  if( proc == instance()->_last ) {

    ced_draw_event();
   
    if ( waitForKeyboard == 1 ) {
      std::cout << "[ Press any key for next event ]" << std::endl ; 
      //SM-H: TODO: Experimental picking code 
      //std::cout << "[ Press 'p' to enter picking mode, or any key for next event ]" << std::endl ; 
      char c = getchar();
      //if(c=='p' || c=='P') {
      //  std::cout << "Entered picking mode" << std::endl;
      //  int old_id = 0; //Check whether this was the previous particle selected: if yes, don't reprint it
      //  for(;;) {
      //    int id = ced_selected_id();
      //    //std::cout << id << std::endl;
      //    if(id>0 && id!=old_id) {
      //      //Need to get MCParticle from ID -> require knowledge of event
      //      old_id = id;
      //      if(instance()->_currEvent!=0) {
      //        MCParticle* mcp = MarlinCED::getMCParticleFromID(id, instance()->_currEvent);
      //        if(mcp!=0) {
      //            printMCParticle(mcp);
      //        }
      //      }
      //    }
      //    else if(id == 0) { 
      //        std::cout << "No particle selected" << std::endl;
      //    }
      //    //else {
      //    //  break;
      //    //}
      //  }
      //}
    }
  }
  //   else
  //   ced_send_event();
}


/**
 * Improved drawHelix() method. Draws straight lines as well.
 */
//SM-H: Added id to drawHelix (default zero), which allows for implementation of picking
void MarlinCED::drawHelix(float b, float charge, float x, float y, float z,
			  float px, float py, float pz, int marker, int size, unsigned int col,
			  float rmin, float rmax, float zmax, unsigned int id)  {

  	// FIXME : check for zmin as well, i.e. cylindrical coordinates

	double cFactor = 2.9979251e-4;
    const double high_pt = 100.0;//Transverse momentum high enough for the particle not to curve noticeably 
  	double pt = sqrt(px*px + py*py);

  // FIXME: use a parameter for this cut or better this should be a function of the B field, charge and momentum 2006/07/04 OW  
  
  // SD: FIXME: Adaptive step-number (or get rid of it!) and adaptive draw step!  
  if ( (pt >= 0.01) && (pt <= high_pt && charge!=0) ) {
    double r =  pt / ( cFactor * b * std::abs( charge )  ) ;
    double sign =  charge > 0 ? 1 : -1 ;
    sign = - sign  ; // FIXME: need to check the convention - but this works !?
    double phi = std::atan2( py , px ) + ( 2. + sign ) * M_PI / 2. ;
    //center of helix
    double cx = x - ( sign * py * r / pt ) ;
    double cy = y + ( sign * px * r / pt ) ;
    double cz = z ; 
    
    double x1 =  x ; 
    double y1 =  y ;
    double z1 =  z ;
    double step = 0.05;  // initial 0.05
    
    // FIX ME: do the adaptive step number...
	
	// cheap adaptive algorithms
	if (px>1 || py >1 || px <-1 || py <-1 ){
		step = 0.005;
		if (px>5 || py >5 || px <-5 || py <-5){
			step = 0.001;
			//std::cout << "Above 5 momenta" << std::endl;
		}
	}
	//std::cout << step << std::endl;
    
    int nSteps = 1000000;
    
    for (int j = 0; j < nSteps ; j++) {
      
      double alpha = step*j ;  
      
      double x2 = cx + r * cos( phi + sign * alpha ) ;
      double y2 = cy + r * sin( phi + sign * alpha ) ;
      double z2 = cz + r * alpha * pz / pt ;
      
      double r_current  = sqrt(x2*x2 + y2*y2); // hypot( x2, y2 ) 

	/*
	 *  interpolation and loop break
	 */
      if( std::abs(z2) > zmax || r_current > rmax  ) {
      	
      	double alpha = step*(j+0.5); 
      	
     	 x2 = cx + r * cos( phi + sign * alpha ) ;
      	 y2 = cy + r * sin( phi + sign * alpha ) ;
      	 z2 = cz + r * alpha * pz / pt ;
      	//std::cout 	<< "Number of steps = " << j << std::endl;
		break ;
      }
      
		if( r_current >= (rmin+step)) {
			ced_line_ID( x1, y1, z1, x2, y2, z2 , marker , size, col, id);
      	}
    x1 = x2;
    y1 = y2;
    z1 = z2;

    }

  }
  //For high momentum tracks, just draw straight line 
  else if (pt > high_pt) { 
        std::cout << "pt = " << pt << std::endl;
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
        
        // this has been improved
        
        if (z + (kr*pz)/absP > zmax || z + (kr*pz)/absP < -zmax){
        	k = kz;
        }
        else k = kr;
        
        if (k < 0.0 && k!=kz) {
          std::cout << "Error in 'MarlinCED::drawHelix()': No intersection point with the outer cylinder shell (rmax,zmax) found" << std::endl;
          return;
        }
        
        float xEnd = x + (k*px)/absP;
        float yEnd = y + (k*py)/absP;
        float zEnd = z + (k*pz)/absP;
        
        if (rmin != 0){
        	std::cout << "FIX ME: Inner cylinder not taken into account!" << std::endl;
        	return;
        }
        
        ced_line_ID(x, y, z, xEnd, yEnd, zEnd , marker , size, col, id);
    
  	}
    else {
        std::cout << "Low momentum particle given point instead of helix" << std::endl;
        const double delta = 0.0001;
        ced_line_ID(x, y, z, x+delta, y+delta, z+delta, marker , size, col, id);

        //ced_hit ( x,y,z, marker,size , col);
    }
}
  


void MarlinCED::drawTrajectory(const Trajectory* t, const int marker,
			       const int size, const unsigned int col,
			       const float rmin, const float rmax,
			       const float zmax) 
{
  if (rmax <= rmin || zmax == 0 ) return;
  double stepSize = 5.0; // initial 0.05
  double nStepKill = int(2*M_PI*rmax/stepSize);

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


  //SM-H: Calls drawHelix with MCP->Iid(), which allows for implementation of picking
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
  

  charge = MCP->getCharge();

  // debug
  std::cout << bField << "  " << charge << "  " << x1 << "  " << y1 << "  " << z1 << "  " << x2 << "  " << y2 << "  " << z2 << "  " << color << " " 
            << p1 << "  " << p2 << "  " << p3 << "  " << std::endl;


  bool isCharged       = charge != 0.0;    
  bool isNeutrino      = (abs(MCP->getPDG())==12) || (abs(MCP->getPDG())==14) || (abs(MCP->getPDG())==16);
  bool isBackscattered = MCP->isBackscatter();
    
  //std::cout << "isCharged = " << isCharged << " isNeutrino = " << isNeutrino << " isBackscattered " << isBackscattered << std::endl;

  // charged MC Particles are displayed on layer and their SimHits optionally on (layer + 10)
  if (isCharged && !isNeutrino && !isBackscattered) {

    drawHelix(bField,charge,x1,y1,z1,p1,p2,p3, marker | ( layer << CED_LAYER_SHIFT ), size, color, (float)rmin, (float)rmax, (float)zmax, MCP->id());
    if (drawSimHits) drawHitCollectionsByMCContribution(event,MCP,marker,size+2,color,layer+10);

  }

  // neutral MC Particles are displayed on (layer+1) and their SimHits optionally on (layer + 11)
  else if (!isCharged && !isNeutrino && !isBackscattered) {

    if (drawOnDifferentLayers) l = layer+1;
    else l = layer;

    ced_line_ID(x1,y1,z1,x2,y2,z2, marker | ( l << CED_LAYER_SHIFT ), size, color, MCP->id());

    if (drawSimHits) {
      if (drawOnDifferentLayers) l = layer+11;
      else l = layer+10;
      drawHitCollectionsByMCContribution(event,MCP,marker,size+2,color,l);
    }

  }
  // backscattered charged particles and neutrinos are displayed on (layer+2) and their SimHits optionally on (layer + 12)
  else if (isCharged && !isNeutrino && isBackscattered) {

    if (drawOnDifferentLayers) l = layer+2;
    else l = layer;

    drawHelix(bField,charge,x1,y1,z1,p1,p2,p3, marker | ( l << CED_LAYER_SHIFT ), size, color, (float)rmin, (float)rmax, (float)zmax, MCP->id());

    if (drawSimHits) {
      if (drawOnDifferentLayers) l = layer+12;
      else l = layer+10;
      drawHitCollectionsByMCContribution(event,MCP,marker,size+2,color,l);
    }

  }
  // backscattered charged particles and neutrinos are displayed on (layer+2) and their SimHits optionally on (layer + 12)
  else if (!isCharged && isNeutrino || isBackscattered) {
  
    if (drawOnDifferentLayers) l = layer+2;
    else l = layer;

    ced_line_ID(x1,y1,z1,x2,y2,z2, marker | ( l << CED_LAYER_SHIFT ), size, color, MCP->id());

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




void MarlinCED::drawGEARDetector(){
//
// code from V.Morgunov, MPI
//

   const gear::TPCParameters&  pTPC      = Global::GEAR->getTPCParameters();
   const gear::PadRowLayout2D& padLayout = pTPC.getPadLayout();
   const gear::DoubleVec&      planeExt  = padLayout.getPlaneExtent();
   float r_min_tpc = planeExt[0];
   float r_max_tpc = planeExt[1];
   float z_max_tpc = pTPC.getMaxDriftLength();
   const gear::CalorimeterParameters& pECAL_B = 
            Global::GEAR->getEcalBarrelParameters();
   float r_min_ecal_bar = pECAL_B.getExtent()[0];
   float r_max_ecal_bar = pECAL_B.getExtent()[1];
// float z_min_ecal_bar = pECAL_B.getExtent()[2];
   float z_max_ecal_bar = pECAL_B.getExtent()[3];
   const gear::CalorimeterParameters& pECAL_E = 
            Global::GEAR->getEcalEndcapParameters();
// float r_min_ecal_ecap = pECAL_E.getExtent()[0];
   float r_max_ecal_ecap = pECAL_E.getExtent()[1];
   float z_min_ecal_ecap = pECAL_E.getExtent()[2];
   float z_max_ecal_ecap = pECAL_E.getExtent()[3];
   const gear::CalorimeterParameters& pHCAL_B = 
            Global::GEAR->getHcalBarrelParameters();
   //  _innerHcalRadius = float(pHcalBarrel.getExtent()[0]);
   float r_min_hcal_bar = pHCAL_B.getExtent()[0];
   float r_max_hcal_bar = pHCAL_B.getExtent()[1];
// float z_min_hcal_bar = pHCAL_B.getExtent()[2];
   float z_max_hcal_bar = pHCAL_B.getExtent()[3];
   const gear::CalorimeterParameters& pHCAL_R = 
             Global::GEAR->getHcalRingParameters();
// float r_min_hcal_ring = pHCAL_R.getExtent()[0];
    float r_max_hcal_ring = pHCAL_R.getExtent()[1];
    float z_min_hcal_ring = pHCAL_R.getExtent()[2];
    float z_max_hcal_ring = pHCAL_R.getExtent()[3];
   const gear::CalorimeterParameters& pHCAL_E = 
            Global::GEAR->getHcalEndcapParameters();
// float r_min_hcal_ecap = pHCAL_E.getExtent()[0];
   float r_max_hcal_ecap = pHCAL_E.getExtent()[1];
   float z_min_hcal_ecap = pHCAL_E.getExtent()[2];
   float z_max_hcal_ecap = pHCAL_E.getExtent()[3];
   
     // **************************************** //
  // ** Building VTX Detector ** //
  // **************************************** //

  //--Get GEAR Parameters--
  const gear::VXDParameters& pVXDDetMain = Global::GEAR->getVXDParameters();
  const gear::VXDLayerLayout& pVXDLayerLayout = pVXDDetMain.getVXDLayerLayout();

  int nLayersVTX = pVXDLayerLayout.getNLayers();
  
  float rad2deg = 180.0 / M_PI;

  for (int i=0; i<nLayersVTX; ++i) {

    int nLadders = pVXDLayerLayout.getNLadders(i);

    float _ladder_phi0 = float(pVXDLayerLayout.getPhi0(i));

    float _sensitive_distance = float(pVXDLayerLayout.getSensitiveDistance(i));
    float _sensitive_thickness = float(pVXDLayerLayout.getSensitiveThickness(i));
    float _sensitive_width = float(pVXDLayerLayout.getSensitiveWidth(i));
    float _sensitive_length = float(pVXDLayerLayout.getSensitiveLength(i));
    float _sensitive_offset = float (pVXDLayerLayout.getSensitiveOffset(i));

    float currPhi;
    float angleLadders = 2*M_PI / nLadders;
    float cosphi, sinphi;

    _sensitive_distance +=0.5* _sensitive_thickness;

    for (int j=0; j<nLadders; ++j) {

      currPhi = _ladder_phi0 + (angleLadders * j);
      cosphi = cos(currPhi);
      sinphi = sin(currPhi);

      double * sizes  = new double[3];
      double * center = new double[3];
      unsigned int color = 0xFFFFFF;

      center[0] = (_sensitive_distance*cosphi - _sensitive_offset*sinphi);
      center[1] = (_sensitive_distance*sinphi + _sensitive_offset*cosphi);
      center[2] = 0.0;
      sizes[0]  = _sensitive_thickness;
      sizes[1]  = _sensitive_width;
      sizes[2]  = _sensitive_length;

      unsigned int layer = 11<<CED_LAYER_SHIFT;
      
      double *rotate = new double[3];
      rotate[2] = currPhi*rad2deg;

      ced_geobox_r( sizes, center, rotate, color, layer);

      delete [] center;
      delete [] sizes;



    }

  }


   
   
// =======================================================================
//To convert inner radius of polygone to its outer radius
// float Cos4  = cos(M_PI/4.0);
   float Cos8  = cos(M_PI/8.0);
   float Cos16 = cos(M_PI/16.);
// convertion of  inner radius of polygone to its outer radius
   float r_inn_ecal_bar     = r_min_ecal_bar/Cos8;  
   float r_out_ecal_bar     = (r_max_ecal_bar)/Cos8 ;
//   float r_inn_ecal_ecap    = r_min_ecal_ecap/Cos4;
   float r_out_ecal_ecap    = r_max_ecal_ecap/Cos8;
   float thick_ecal_ecap    = 0.5*(z_max_ecal_ecap - z_min_ecal_ecap);
   float shift_ecal_z_plus  = z_min_ecal_ecap;
   float shift_ecal_z_minus = z_min_ecal_ecap + 2.0*thick_ecal_ecap;
   float r_inn_hcal_bar     = r_min_hcal_bar/Cos8;
   float r_out_hcal_bar     = r_max_hcal_bar/Cos16;
//   float r_inn_hcal_ring    = r_min_hcal_ring/Cos4;
   float r_out_hcal_ring    = r_max_hcal_ring/Cos8;
   float thick_hcal_ring    = 0.5*(z_max_hcal_ring - 
		 z_min_hcal_ring + 20.0); // +20 by hand to see hits inside
   float shift_hcalr_z_plus  = z_min_hcal_ring;
   float shift_hcalr_z_minus = z_min_hcal_ring + 2.0*thick_hcal_ring;
//   float r_inn_hcal_ecap    = r_min_hcal_ecap/Cos4;
   float r_out_hcal_ecap    = r_max_hcal_ecap/Cos8;
   float thick_hcal_ecap    = 0.5*(z_max_hcal_ecap - 
		 z_min_hcal_ecap + 20.0); // +20 by hand to see hits inside
   float shift_hcal_z_plus  = z_min_hcal_ecap;
   float shift_hcal_z_minus = z_min_hcal_ecap + 2.0*thick_hcal_ecap;
// ========================================================================
    static CED_GeoCylinder geoCylindersANY[] = {       // for ANY Detector Geometry
      { r_min_tpc,        40, 0.0, z_max_tpc,       -z_max_tpc,          0xff      }, // inner TPC  40 also temporary
      { r_max_tpc,        40, 0.0, z_max_tpc,       -z_max_tpc,          0xff      }, // outer TPC  temporary
      { r_inn_ecal_bar ,  8, 22.5, z_max_ecal_bar,  -z_max_ecal_bar,     0x7f7f1f  }, // inner ECAL Barrel
      { r_out_ecal_bar ,  8, 22.5, z_max_ecal_bar,  -z_max_ecal_bar,     0x7f7f1f  }, // outer ECAL Barrel
      { r_out_ecal_ecap,  8, 22.5, thick_ecal_ecap,  shift_ecal_z_plus,  0x7f7f1f  }, // outer endcap ECAL +Z
      { r_out_ecal_ecap,  8, 22.5, thick_ecal_ecap, -shift_ecal_z_minus, 0x7f7f1f  }, // outer endcap ECAL -Z
      { r_inn_hcal_bar ,  8, 22.5, z_max_hcal_bar,  -z_max_hcal_bar,     0x00cf00  }, // inner HCAL Barrel
      { r_out_hcal_bar , 16,  0.0, z_max_hcal_bar,  -z_max_hcal_bar,     0x00cf00  }, // outer HCAL Barrel
      { r_out_hcal_ring,  8,  0.0, thick_hcal_ring,  shift_hcalr_z_plus,  0x00cf00  }, // outer ring HCAL +Z
      { r_out_hcal_ring,  8,  0.0, thick_hcal_ring, -shift_hcalr_z_minus, 0x00cf00 } , // outer ring HCAL -Z      
      { r_out_hcal_ecap,  8,  0.0, thick_hcal_ecap,  shift_hcal_z_plus,  0x00cf00  }, // outer endcap HCAL +Z
      { r_out_hcal_ecap,  8,  0.0, thick_hcal_ecap, -shift_hcal_z_minus, 0x00cf00  }  // outer endcap HCAL -Z      
    };
// ========================================================================
      ced_geocylinders(sizeof(geoCylindersANY)/sizeof(CED_GeoCylinder),
		       geoCylindersANY);
} // End of class DrawGeometry

