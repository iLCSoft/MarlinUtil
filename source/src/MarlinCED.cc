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

//hauke
#include "UTIL/LCTOOLS.h"
#include "UTIL/Operators.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/RawCalorimeterHit.h"
#include "EVENT/SimTrackerHit.h"

#include "UTIL/LCTOOLS.h"
#include "UTIL/Operators.h"

#include "EVENT/LCCollection.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/RawCalorimeterHit.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/TPCHit.h"
#include "EVENT/TrackerRawData.h"
#include "EVENT/TrackerData.h"
#include "EVENT/TrackerPulse.h"
#include "EVENT/LCIO.h"
#include "EVENT/MCParticle.h"
#include "EVENT/LCFloatVec.h"
#include "EVENT/LCIntVec.h"
#include "IMPL/LCFlagImpl.h"
#include "EVENT/Track.h"
#include "EVENT/Cluster.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Vertex.h"
#include "EVENT/LCGenericObject.h" 
#include "EVENT/LCRelation.h"
#include "LCIOSTLTypes.h"
#include <signal.h>

#include "UTIL/LCObjectHandle.h"
#include "UTIL/LCTime.h"
#include "UTIL/CellIDDecoder.h"
#include "UTIL/PIDHandler.h"

/*
using namespace std ;
using namespace EVENT ;
using namespace IMPL ;
*/
using namespace UTIL;



//SJA:FIXED:added to make gcc4.3 compliant
#include <cstdlib>
#include <cstdio>

//hauke
#include <ctime>
#include <time.h>
#include <termios.h>
#include <poll.h>

//for kbhit
#include <sys/select.h>
#include <termios.h>

MarlinCED* MarlinCED::_me = 0;

//hauke hoelbe
CEDPickingMap CEDPickingHandler::map;
CEDFunctionMap CEDPickingHandler::funcMap;
CEDPickingHandler *CEDPickingHandler::instance = NULL;

std::vector<std::string> MarlinCED::_descs(MAX_LAYER, ""); //layer descriptions

CEDPickingHandler& CEDPickingHandler::getInstance() {
    if( !instance){
        std::cout<<"new instance"<<std::endl;
        instance = new CEDPickingHandler();
    }
    return *instance;
}

void CEDPickingHandler::update(LCEvent *evt){
    const std::vector<std::string> *collNames = evt->getCollectionNames();
    const LCCollection *coll;
    const LCObject *obj;
    CEDFunctionMap::iterator iter;
    std::map<std::string, void (*)(const LCObject *)> printDefaultMap;
    std::map<std::string, void (*)(const LCObject *)>::iterator printDefaultIter;   
    clock_t start =  clock(); 

    //default print functions
    printDefaultMap.insert(std::make_pair(LCIO::MCPARTICLE,&defaultOutputFunction<MCParticle>));
    printDefaultMap.insert(std::make_pair(LCIO::TRACKERHIT, &defaultOutputFunction<TrackerHit>));
    printDefaultMap.insert(std::make_pair(LCIO::SIMTRACKERHIT, &defaultOutputFunction<SimTrackerHit>));
    printDefaultMap.insert(std::make_pair(LCIO::CALORIMETERHIT, &defaultOutputFunction<CalorimeterHit>));
    printDefaultMap.insert(std::make_pair(LCIO::SIMCALORIMETERHIT, &defaultOutputFunction<SimCalorimeterHit>));
    printDefaultMap.insert(std::make_pair(LCIO::VERTEX, &defaultOutputFunction<Vertex>));
    printDefaultMap.insert(std::make_pair(LCIO::RECONSTRUCTEDPARTICLE, &defaultOutputFunction<ReconstructedParticle>));
    printDefaultMap.insert(std::make_pair(LCIO::TRACK, &defaultOutputFunction<Track>));
    printDefaultMap.insert(std::make_pair(LCIO::CLUSTER, &defaultOutputFunction<Cluster>));
    printDefaultMap.insert(std::make_pair(LCIO::LCRELATION, &defaultOutputFunction<LCRelation>));
    printDefaultMap.insert(std::make_pair(LCIO::LCFLOATVEC, &defaultOutputFunction<LCFloatVec>));

    //add default print functions. Set it only if the user dont have set his own functions
    for( printDefaultIter = printDefaultMap.begin(); printDefaultIter != printDefaultMap.end(); printDefaultIter++ ) {
        iter=funcMap.find(printDefaultIter->first);
        if(iter == funcMap.end()){ //user dont have registered a function for this collection name
            registerFunction(printDefaultIter->first, printDefaultIter->second);
        }
    } 

    std::string typeName;
    std::string collName;
    CEDMapParticleObject particleObj;
    for(unsigned int i=0;i<collNames->size();i++){
        collName=collNames->at(i);
        coll =  evt->getCollection(collName);
        typeName=coll->getTypeName();
        for (int j=0;j<coll->getNumberOfElements();j++){
            obj=coll->getElementAt(j);
            /*
            //debug
            if(true){ //collName == LCIO::MCPARTICLE){
                streamlog_out(DEBUG) << "------------------ TEST -----------------------" << std::endl;
                streamlog_out(DEBUG) << lcio_long(* (EVENT::TrackerHit *) obj,coll);
                streamlog_out(DEBUG) << "-----------------------------------------------" << std::endl;
            }
            //debug end
            */
            particleObj.obj=obj;

            iter =  funcMap.find(collName);
            if(iter != funcMap.end()){ //user have registered a function for this collection name
                particleObj.function=iter->second;
            }else{
                iter=funcMap.find(typeName);
                if(iter != funcMap.end()){ //user have a registered a function for this _typeName_ 
                    particleObj.function=iter->second;
                }else{ 
                    streamlog_out(DEBUG) << "CEDPickingHandler: cant register " << collName << "/" << typeName 
                                         << " (no function given)" << std::endl;
                    continue;
                }
            }
            /*test output at startup (debug)
            print ALL objects while filling the map*/
            //particleObj.function(particleObj.obj);

	    // streamlog_out( DEBUG ) << "  registering object of type : " <<  typeName << "  with ID= " << obj->id() << std::endl ;

	    map.insert(std::pair<const int, CEDMapParticleObject>(obj->id(),particleObj));
        }
   }
   clock_t end = clock() ; 
   streamlog_out(DEBUG) << "CEDPickingHandler::Map size: " << map.size() << " time: " << double( end - start ) / double(CLOCKS_PER_SEC) << "s" << std::endl;
}

void CEDPickingHandler::printID(int id){
    CEDPickingMap::iterator iter;
    CEDMapParticleObject obj;
    void (*printFunction)(const LCObject *); 

    iter = map.find(id);
    if( iter != map.end() ){ 
        obj=iter->second; 

        printFunction=obj.function;
        printFunction(obj.obj); 
    }else{
        streamlog_out(WARNING) << "No print function registered for this collection- or type name!" << std::endl;
    }

}

//--------------------------------------------------------------------------------------------------------

void CEDPickingHandler::registerFunction(std::string type, void (*printFunction)(const LCObject *)){   
    funcMap.insert(std::make_pair(type,printFunction));
}

int CEDPickingHandler::kbhit(void) {
    //http://stackoverflow.com/questions/448944/c-non-blocking-keyboard-input#448982
    struct timeval tv = { 0L, 0L };
    fd_set fds;
    FD_SET(0, &fds);
    return select(1, &fds, NULL, NULL, &tv);
}
/*
struct termios CEDPickingHandler::orig_termios;

void CEDPickingHandler::reset_terminal_mode()
{
    //http://stackoverflow.com/questions/448944/c-non-blocking-keyboard-input#448982
    tcsetattr(0, TCSANOW, &orig_termios);
}

int CEDPickingHandler::getch(void)
{
    int r;
    unsigned char c;
    if ((r = read(0, &c, sizeof(c))) < 0) {
        return r;
    } else {
        return c;
    }
}

void CEDPickingHandler::set_conio_terminal_mode()
{
    //http://stackoverflow.com/questions/448944/c-non-blocking-keyboard-input#448982
    struct termios new_termios;

    // take two copies - one for now, one for later 
    tcgetattr(0, &CEDPickingHandler::orig_termios);
    memcpy(&new_termios, &orig_termios, sizeof(new_termios));

    // register cleanup handler, and set the new terminal mode 
    atexit(CEDPickingHandler::reset_terminal_mode);
    cfmakeraw(&new_termios);
    tcsetattr(0, TCSANOW, &new_termios);
}

*/

void MarlinCED::add_layer_description(const std::string &desc, int layerID){
    std::string tmp;
    //std::cout<<"LAYER: add: " << desc << std::endl;
    //std::cout << "LAYER: search " << desc << " into " << _descs.at(layerID);
    if(layerID > MAX_LAYER || layerID < 0){return;}
    if( _descs.at(layerID).find(desc.c_str()) == std::string::npos){      
        //std::cout << " found " << std::endl;
        tmp=_descs.at(layerID);
        if(! tmp.empty()){
            tmp.append(", ");
        }
        tmp.append(desc);
        _descs.at(layerID)=tmp;
    }else{
        //std::cout << " NOT found " << std::endl;
    }
}

void MarlinCED::set_layer_description(const std::string &desc, int layerID){
    if(layerID > MAX_LAYER || layerID < 0){return;}
    //std::cout<<"LAYER: set: " << desc <<  std::endl;
    _descs.at(layerID)=desc;

}

void MarlinCED::write_layer_description(void){
    //std::cout<<"LAYER: write all layer in ced" << std::endl;
    unsigned int i;
    //for(i=0;i<25;i++){
    for(i=0; i<_descs.size(); i++){
        //std::cout<<"LAYER " << i << ": "<< _descs.at(i) << std::endl;
        ced_describe_layer(_descs.at(i).c_str(), i);
    }
}

//end hauke hoelbe

MarlinCED* MarlinCED::instance() {
  if( _me == 0 )
    _me = new MarlinCED ;
  return _me ;
}


void MarlinCED::init( Processor* proc ) {
  
  if( instance()->_first == 0 ){
    
    instance()->_first = proc ;
    
    char *port, *host;
    port = getenv("CED_PORT");
    host = getenv("CED_HOST");
    if((port == NULL || port[0] == 0) && (host == NULL || host[0] == 0)){
      ced_client_init("localhost",7286);
    }else if(port == NULL || port[0] == 0){
      streamlog_out(MESSAGE)<< "Use user defined host " << host << std::endl;
      ced_client_init(host,7286);
    }else if(host == NULL|| host[0] == 0){
      streamlog_out(MESSAGE)<< "Use user defined port " << port << std::endl;
      ced_client_init("localhost",atoi(port));
    }else{
      streamlog_out(MESSAGE)<< "Use user defined host " << host << ", port " <<  port << std::endl;
      ced_client_init(host,atoi(port));
    }
    
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


//hauke hoelbe modify 08.02.2010
void MarlinCED::draw( Processor* proc , int waitForKeyboard ) {
    //char message[200];
    int i=0;
    CEDPickingHandler &pHandler=CEDPickingHandler::getInstance();


    if( proc == instance()->_last ) {
    //    ced_draw_event();
        MarlinCED::write_layer_description();
//ced_picking_text("test1 test2 test3");

        ced_send_event();
        if ( waitForKeyboard == 1 ) {
            streamlog_out(MESSAGE) << "Double click for picking. Press <ENTER> for the next event." << std::endl;
            //test:

            
            //streamlog_out(MESSAGE)  << "(Please do not resize this terminal window!)" << std::endl;

/*
            struct pollfd pfd[1];
            pfd[0].fd = 1;
            pfd[0].events=POLLIN;

*/
            signal(SIGWINCH,SIG_IGN);
            //streamlog_out(DEBUG) << "sigwinch ign" << std::endl;

            while(!CEDPickingHandler::kbhit()){
//            while(!poll(pfd,1,0)){

                usleep(100000); //micro seconds
                //timeval Timeout;
                //Timeout.tv_sec = 5;
                //Timeout.tv_usec = 0;
                //select( 0, (fd_set *)0, (fd_set *)0, (fd_set *)0, &Timeout );

                int id = ced_selected_id_noblock();
                if(id>=0) {
                    streamlog_out(DEBUG) << "DEBUG: got id: " << id <<std::endl;
                    if(id == 0){
                        streamlog_out(WARNING) << "Picking nothing, or an object with ID 0!" << std::endl;
                    }else{
                        pHandler.printID(id);
                        ced_picking_text("test1 test2 test3",i++);
                        ced_send_event();
                    }
                }
            }
            
            //streamlog_out(DEBUG) << "try to getchar()" <<std::endl;
            signal(SIGWINCH,SIG_IGN);
            char c = getchar();
            //char c = CEDPickingHandler::getch();
            if(c=='q'||c=='Q'||c==3){ //quit if the user pressed q or strg+c (3 = strg+c)
                exit(0);
            }
            streamlog_out(MESSAGE) << "--------- END ---------------\n";
        }
    }
}

/**
 * Improved drawHelix() method. Draws straight lines as well.
 */
//SM-H: Added id to drawHelix (default zero), which allows for implementation of picking
void MarlinCED::drawHelix(float b, float charge, float x, float y, float z,
			  float px, float py, float pz, int marker, int size, unsigned int col,
			  float rmin, float rmax, float zmax, unsigned int id)  {
    //return; //draw nothing

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
    
    //int nSteps = 1000000;
    int nSteps = int(100/step); //hauke
    //streamlog_out(DEBUG) << "draw helix (nsteps: " << nSteps << ") id= " << id<<std::endl; 

    
    int count_lines=0;
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
      	// streamlog_out(DEBUG) 	<< "Number of steps = " << j << std::endl;
		break ;
      }
      
		if( r_current >= (rmin+step)) {
            count_lines++;
			ced_line_ID( x1, y1, z1, x2, y2, z2 , marker , size, col, id);
			//ced_line( x1, y1, z1, x2, y2, z2 , marker , size, col);
      	}
    x1 = x2;
    y1 = y2;
    z1 = z2;

    }
    // streamlog_out(DEBUG)<<"added " <<count_lines <<"ced_line_ID to CED"<<std::endl;

  }
  //For high momentum tracks, just draw straight line 
  else if (pt > high_pt) { 
         streamlog_out(DEBUG) << "pt = " << pt << std::endl;
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
           streamlog_out(ERROR) << "Error in 'MarlinCED::drawHelix()': Startpoint beyond (rmax,zmax)" << std::endl;
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
           streamlog_out(DEBUG) << "Error in 'MarlinCED::drawHelix()': No intersection point with the outer cylinder shell (rmax,zmax) found" << std::endl;
          return;
        }
        
        float xEnd = x + (k*px)/absP;
        float yEnd = y + (k*py)/absP;
        float zEnd = z + (k*pz)/absP;
        
        if (rmin != 0){
        	 streamlog_out(DEBUG) << "FIX ME: Inner cylinder not taken into account!" << std::endl;
        	return;
        }
        
        ced_line_ID(x, y, z, xEnd, yEnd, zEnd , marker , size, col, id);
    
  	}
    else {
         streamlog_out(DEBUG) << "Low momentum particle given point instead of helix" << std::endl;
        const double delta = 0.0001;
        ced_line_ID(x, y, z, x+delta, y+delta, z+delta, marker , size, col, id);

        //ced_hit ( x,y,z, marker,size , col);
    }
}
  


void MarlinCED::drawTrajectory(const Trajectory* t, const int marker,
			       const int size, const unsigned int col,
			       const float rmin, const float rmax,
			       const float zmax, unsigned int id)  //hauke:: addet optional argument id
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
      ced_line_ID( xold.x(), xold.y(), xold.z(), x.x(), x.y(), x.z(), 
		marker, size, col, id); //hauke: added id
    }
}

void MarlinCED::drawSpike( float x0, float y0, float z0,float x1, float y1, float z1, unsigned int color, unsigned int layer, unsigned int id) {

//hauke: added optional argument id
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
  
  ced_line_ID( p0[0],p0[1],p0[2], p1[0],p1[1],p1[2], layty ,6 , color, id ); //hauke: added id
  ced_line_ID( p1[0],p1[1],p1[2], p2[0],p2[1],p2[2], layty ,5 , color, id ); //hauke ...
  ced_line_ID( p2[0],p2[1],p2[2], p3[0],p3[1],p3[2], layty ,4 , color, id );
  ced_line_ID( p3[0],p3[1],p3[2], p4[0],p4[1],p4[2], layty ,3 , color, id );
  ced_line_ID( p4[0],p4[1],p4[2], p5[0],p5[1],p5[2], layty ,2 , color, id );
  
  //ced_hit_ID ( p0[0],p0[1],p0[2], CED_HIT_POINT | layer << CED_LAYER_SHIFT, 0, color, id );
  //ced_hit_ID ( p5[0],p5[1],p5[2], CED_HIT_POINT | layer << CED_LAYER_SHIFT, 0, color, id );
  ced_hit_ID ( p0[0],p0[1],p0[2], CED_HIT_POINT, layer, 0, color, id );
  ced_hit_ID ( p5[0],p5[1],p5[2], CED_HIT_POINT, layer, 0, color, id );

  
}


void MarlinCED::drawMCParticle(MCParticle* MCP, bool drawSimHits, LCEvent* event, int marker, int size, unsigned int color, unsigned int layer, double bField,
			       double rmin, double zmin, double rmax, double zmax, bool drawOnDifferentLayers) {


   streamlog_out(DEBUG)<<"Hauke: draw mcparticle, id="<<MCP->id() << std::endl;


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
   streamlog_out(DEBUG) << bField << "  " << charge << "  " << x1 << "  " << y1 << "  " << z1 << "  " << x2 << "  " << y2 << "  " << z2 << "  " << color << " " 
            << p1 << "  " << p2 << "  " << p3 << "  " << std::endl;


  bool isCharged       = charge != 0.0;    
  bool isNeutrino      = (abs(MCP->getPDG())==12) || (abs(MCP->getPDG())==14) || (abs(MCP->getPDG())==16);
  bool isBackscattered = MCP->isBackscatter();
    
  // streamlog_out(DEBUG) << "isCharged = " << isCharged << " isNeutrino = " << isNeutrino << " isBackscattered " << isBackscattered << std::endl;

  // charged MC Particles are displayed on layer and their SimHits optionally on (layer + 10)
  if (isCharged && !isNeutrino && !isBackscattered) {

    drawHelix(bField,charge,x1,y1,z1,p1,p2,p3, marker | ( layer << CED_LAYER_SHIFT ), size, color, (float)rmin, (float)rmax, (float)zmax, MCP->id());
    //drawHelix(bField,charge,x1,y1,z1,p1,p2,p3, marker | ( layer), size, color, (float)rmin, (float)rmax, (float)zmax, MCP->id());

    if (drawSimHits) drawHitCollectionsByMCContribution(event,MCP,marker,size+2,color,layer+10);

  }

  // neutral MC Particles are displayed on (layer+1) and their SimHits optionally on (layer + 11)
  else if (!isCharged && !isNeutrino && !isBackscattered) {

    if (drawOnDifferentLayers){ 
        l = layer+1;
    } else{ l = layer;}

    ced_line_ID(x1,y1,z1,x2,y2,z2, marker | ( l << CED_LAYER_SHIFT ), size, color, MCP->id());
    //ced_line_ID(x1,y1,z1,x2,y2,z2, marker | ( l ), size, color, MCP->id());


    if (drawSimHits) {
        if (drawOnDifferentLayers){ 
            l = layer+11;
        }
      else{ 
            l = layer+10;
      }

      drawHitCollectionsByMCContribution(event,MCP,marker,size+2,color,l);
    }

  }
  // backscattered charged particles and neutrinos are displayed on (layer+2) and their SimHits optionally on (layer + 12)
  else if (isCharged && !isNeutrino && isBackscattered) {

    if (drawOnDifferentLayers) {
        l = layer+2;
    }
    else l = layer;

    drawHelix(bField,charge,x1,y1,z1,p1,p2,p3, marker | ( l << CED_LAYER_SHIFT ), size, color, (float)rmin, (float)rmax, (float)zmax, MCP->id());
    //drawHelix(bField,charge,x1,y1,z1,p1,p2,p3, marker | ( l), size, color, (float)rmin, (float)rmax, (float)zmax, MCP->id());


    if (drawSimHits) {
      if (drawOnDifferentLayers) {
            l = layer+12;
      }
      else{ 
            l = layer+10;

      }
      drawHitCollectionsByMCContribution(event,MCP,marker,size+2,color,l);
    }

  }
  // backscattered charged particles and neutrinos are displayed on (layer+2) and their SimHits optionally on (layer + 12)
  else if( (!isCharged && isNeutrino) || isBackscattered ){
  
    if (drawOnDifferentLayers){ 
            l = layer+2;
    }
    else l = layer;

    ced_line_ID(x1,y1,z1,x2,y2,z2, marker | ( l << CED_LAYER_SHIFT ), size, color, MCP->id());
    //ced_line_ID(x1,y1,z1,x2,y2,z2, marker | ( l), size, color, MCP->id());


    if (drawSimHits) {
      if (drawOnDifferentLayers){ 
            l = layer+12;
       }
      else l = layer+10;
      drawHitCollectionsByMCContribution(event,MCP,marker,size+2,color,l);
    }

  }
  else  streamlog_out(DEBUG) << "This MC Particle has not been displayed : id = " << MCP->id() << "  " << "PDG Code : " << MCP->getPDG() << std::endl;
  
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

     streamlog_out(WARNING) << "no valid MC collection in event." << std::endl ;

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
   //   float r_min_ecal_ecap = pECAL_E.getExtent()[0];
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
 float r_min_hcal_ring = pHCAL_R.getExtent()[0];
    float r_max_hcal_ring = pHCAL_R.getExtent()[1];
    float z_min_hcal_ring = pHCAL_R.getExtent()[2];
    float z_max_hcal_ring = pHCAL_R.getExtent()[3];
   const gear::CalorimeterParameters& pHCAL_E = 
            Global::GEAR->getHcalEndcapParameters();
   float r_min_hcal_ecap = pHCAL_E.getExtent()[0];
   float r_max_hcal_ecap = pHCAL_E.getExtent()[1];
   float z_min_hcal_ecap = pHCAL_E.getExtent()[2];
   float z_max_hcal_ecap = pHCAL_E.getExtent()[3];
   
   const gear::CalorimeterParameters& pLHCal = 
            Global::GEAR->getLHcalParameters();
   float r_min_lhcal = pLHCal.getExtent()[0];
   float r_max_lhcal = pLHCal.getExtent()[1];
   float z_min_lhcal = pLHCal.getExtent()[2];
   float z_max_lhcal = pLHCal.getExtent()[3];

   const gear::CalorimeterParameters& pLCal = 
            Global::GEAR->getLcalParameters();
   float r_min_lcal = pLCal.getExtent()[0];
   float r_max_lcal = pLCal.getExtent()[1];
   float z_min_lcal = pLCal.getExtent()[2];
   float z_max_lcal = pLCal.getExtent()[3];

   const gear::CalorimeterParameters& pBeamcal = 
            Global::GEAR->getBeamCalParameters();
   float r_min_beamcal = pBeamcal.getExtent()[0];
   float r_max_beamcal = pBeamcal.getExtent()[1];
   float z_min_beamcal = pBeamcal.getExtent()[2];
   float z_max_beamcal = pBeamcal.getExtent()[3];


   const gear::CalorimeterParameters& pYOKE_B = 
     Global::GEAR->getYokeBarrelParameters();
   //  _innerYokeRadius = float(pYokeBarrel.getExtent()[0]);
   float r_min_yoke_bar = pYOKE_B.getExtent()[0];
   float r_max_yoke_bar = pYOKE_B.getExtent()[1];
   //float z_min_yoke_bar = pYOKE_B.getExtent()[2];
   float z_max_yoke_bar = pYOKE_B.getExtent()[3];
   const gear::CalorimeterParameters& pYOKE_R = 
             Global::GEAR->getYokePlugParameters();
   float r_min_yoke_plug = pYOKE_R.getExtent()[0];
   float r_max_yoke_plug = pYOKE_R.getExtent()[1];
   float z_min_yoke_plug = pYOKE_R.getExtent()[2];
   float z_max_yoke_plug = pYOKE_R.getExtent()[3];
   const gear::CalorimeterParameters& pYOKE_E = 
     Global::GEAR->getYokeEndcapParameters();
   float r_min_yoke_ecap = pYOKE_E.getExtent()[0];
   float r_max_yoke_ecap = pYOKE_E.getExtent()[1];
   float z_min_yoke_ecap = pYOKE_E.getExtent()[2];
   float z_max_yoke_ecap = pYOKE_E.getExtent()[3];
   
 
   const gear::GearParameters&  pCoil      = Global::GEAR->getGearParameters("CoilParameters");
   float coil_half_z         =  pCoil.getDoubleVal("Coil_cryostat_half_z" ) ;
   float coil_inner_radius  =   pCoil.getDoubleVal("Coil_cryostat_inner_radius" ) ;
   float coil_outer_radius  =   pCoil.getDoubleVal("Coil_cryostat_outer_radius" ) ;

   const gear::GearParameters&  pFTD      = Global::GEAR->getGearParameters("FTD");
   //   DoubleVec& ftd_thick  =  pFTD.getDoubleVal("FTDDiskSiThickness" ) ;
   const DoubleVec& ftd_d   =  pFTD.getDoubleVals("FTDDiskSupportThickness" )  ;
   const DoubleVec& ftd_ri  =  pFTD.getDoubleVals("FTDInnerRadius" )  ;
   const DoubleVec& ftd_ro  =  pFTD.getDoubleVals("FTDOuterRadius" )  ;
   const DoubleVec& ftd_z   =  pFTD.getDoubleVals("FTDZCoordinate" )  ;
   

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
      //unsigned int layer = 11;

      
      double *rotate = new double[3];

      rotate[0] = 0.0;
      rotate[1] = 0.0;
      rotate[2] = currPhi*rad2deg;


      ced_geobox_r( sizes, center, rotate, color, layer);

      delete [] center;
      delete [] sizes;
      delete [] rotate;


    }

  }


   
   
// =======================================================================
//To convert inner radius of polygone to its outer radius
//   float Cos4  = cos(M_PI/4.0);
   float Cos8  = cos(M_PI/8.0);
   float Cos12  = cos(M_PI/12.0);
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

   float r_inn_hcal_ring    = r_min_hcal_ring/Cos8; //fg: was cos4
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
   
   float thick_lhcal         = 0.5 * ( z_max_lhcal -  z_min_lhcal ) ; 
   float shift_lhcal_z_plus  = z_min_lhcal;
   float shift_lhcal_z_minus = z_min_lhcal +  2.0 * thick_lhcal ;

   float thick_lcal         = 0.5 * ( z_max_lcal -  z_min_lcal ) ; 
   float shift_lcal_z_plus  = z_min_lcal;
   float shift_lcal_z_minus = z_min_lcal +  2.0 * thick_lcal ;

   float thick_beamcal         = 0.5 * ( z_max_beamcal -  z_min_beamcal ) ; 
   float shift_beamcal_z_plus  = z_min_beamcal;
   float shift_beamcal_z_minus = z_min_beamcal +  2.0 * thick_beamcal ;


   float r_inn_yoke_bar     = r_min_yoke_bar/Cos12;
   float r_out_yoke_bar     = r_max_yoke_bar/Cos12;

   float r_inn_yoke_plug    = r_min_yoke_plug/Cos12; //fg: was cos4
   float r_out_yoke_plug    = r_max_yoke_plug/Cos12;
   float thick_yoke_plug    = 0.5*(z_max_yoke_plug - 
		 z_min_yoke_plug + 20.0); // +20 by hand to see hits inside
   float shift_yoker_z_plus  = z_min_yoke_plug;
   float shift_yoker_z_minus = z_min_yoke_plug + 2.0*thick_yoke_plug;
   //   float r_inn_yoke_ecap    = r_min_yoke_ecap/Cos12;
   float r_out_yoke_ecap    = r_max_yoke_ecap/Cos12;
   float thick_yoke_ecap    = 0.5*(z_max_yoke_ecap - 
		 z_min_yoke_ecap + 20.0); // +20 by hand to see hits inside
   float shift_yoke_z_plus  = z_min_yoke_ecap;
   float shift_yoke_z_minus = z_min_yoke_ecap + 2.0*thick_yoke_ecap;


   //SIT
   const gear::GearParameters& pSITDet = Global::GEAR->getGearParameters("SIT");
   
   const DoubleVec& rSIT = pSITDet.getDoubleVals("SITLayerRadius")  ;
   const DoubleVec& lSIT = pSITDet.getDoubleVals("SITLayerHalfLength") ;

// ========================================================================
/*
    static CED_GeoCylinder geoCylindersANY[] = {       // for ANY Detector Geometry
      //
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
      ced_geocylinders(sizeof(geoCylindersANY)/sizeof(CED_GeoCylinder), geoCylindersANY);

*/
   
   // colors used in Mokka ILD_00
   static const unsigned tpcCol  = 0xf5f300 ;
   static const unsigned ecalCol = 0x7bf300 ;
   static const unsigned hcalCol = 0xc4c231 ;
   static const unsigned yokeCol = 0x18c2c4 ;
   static const unsigned coilCol = 0x4949dd ;
   static const unsigned ftdCol  = 0x651c93 ;
   static const unsigned fcalCol = 0xabaaab ;

   static const unsigned sitCol = 0xdddddd ; // light grey
   
   int fl=NUMBER_DATA_LAYER;
   static CED_GeoTube geoTubesANY[] = {       // for ANY Detector Geometry
     //ATTENTION: Please order the items from the inner to the outer

      //   radius,       inner_radius,outer_edges,inner_edges,o_phi0, i_phi0-o_phi0, z_max, z_min,    color, inner_detector_shape visiable in classic view, outer detector shape visiable in classic view

      { ftd_ri[0],          ftd_ro[0],  40,  40,   0.0   , 0, ftd_d[0]  ,   ftd_z[0] ,  ftdCol,fl+0 },  //  FTD    
      { ftd_ri[1],          ftd_ro[1],  40,  40,   0.0   , 0, ftd_d[1]  ,   ftd_z[1] ,  ftdCol,fl+0 },  //  FTD    
      { ftd_ri[2],          ftd_ro[2],  40,  40,   0.0   , 0, ftd_d[2]  ,   ftd_z[2] ,  ftdCol,fl+0 },  //  FTD    
      { ftd_ri[3],          ftd_ro[3],  40,  40,   0.0   , 0, ftd_d[3]  ,   ftd_z[3] ,  ftdCol,fl+0},  //  FTD    
      { ftd_ri[4],          ftd_ro[4],  40,  40,   0.0   , 0, ftd_d[4]  ,   ftd_z[4] ,  ftdCol,fl+0},  //  FTD    
      { ftd_ri[5],          ftd_ro[5],  40,  40,   0.0   , 0, ftd_d[5]  ,   ftd_z[5] ,  ftdCol,fl+0},  //  FTD    
      { ftd_ri[6],          ftd_ro[6],  40,  40,   0.0   , 0, ftd_d[6]  ,   ftd_z[6] ,  ftdCol,fl+0 },  //  FTD    
      { ftd_ri[0],          ftd_ro[0],  40,  40,   0.0   , 0, ftd_d[0]  ,  - ftd_z[0] ,  ftdCol,fl+0},  //  FTD    
      { ftd_ri[1],          ftd_ro[1],  40,  40,   0.0   , 0, ftd_d[1]  ,  - ftd_z[1] ,  ftdCol,fl+0},  //  FTD    
      { ftd_ri[2],          ftd_ro[2],  40,  40,   0.0   , 0, ftd_d[2]  ,  - ftd_z[2] ,  ftdCol,fl+0},  //  FTD    
      { ftd_ri[3],          ftd_ro[3],  40,  40,   0.0   , 0, ftd_d[3]  ,  - ftd_z[3] ,  ftdCol,fl+0},  //  FTD    
      { ftd_ri[4],          ftd_ro[4],  40,  40,   0.0   , 0, ftd_d[4]  ,  - ftd_z[4] ,  ftdCol,fl+0},  //  FTD    
      { ftd_ri[5],          ftd_ro[5],  40,  40,   0.0   , 0, ftd_d[5]  ,  - ftd_z[5] ,  ftdCol,fl+0},  //  FTD    
      { ftd_ri[6],          ftd_ro[6],  40,  40,   0.0   , 0, ftd_d[6]  ,  - ftd_z[6] ,  ftdCol,fl+0},  //  FTD    

      { r_max_beamcal,      r_min_beamcal,              40, 40,    0., 0, thick_beamcal,  shift_beamcal_z_plus,   fcalCol,fl+1}, //    BEAMCAL +Z
      { r_max_beamcal,      r_min_beamcal,              40, 40, 0.,    0, thick_beamcal, -shift_beamcal_z_minus,  fcalCol,fl+1},  //   BEAMCAL -Z      
         
      { r_max_tpc,          r_min_tpc,                  40, 40,  0.0, 0, z_max_tpc,        -z_max_tpc,            tpcCol,fl+2,0,1}, //  TPC

      { r_out_ecal_bar,     r_inn_ecal_bar,              8,  8, 22.5, 0, z_max_ecal_bar,   -z_max_ecal_bar,       ecalCol,fl+3,0,1}, //  ECAL Barrel

      { r_out_ecal_ecap,    0.5*(r_max_lhcal+r_max_lcal),8, 40, 22.5, 0,    thick_ecal_ecap,   shift_ecal_z_plus,    ecalCol,fl+3}, //  endcap ECAL +Z
      { r_out_ecal_ecap,    0.5*(r_max_lhcal+r_max_lcal),8, 40, 22.5, 0,  thick_ecal_ecap,  -shift_ecal_z_minus,   ecalCol,fl+3}, //  endcap ECAL -Z
      //hauke: the ecal above and below overlaps(??)
      //{ r_out_ecal_ecap,    r_min_ecal_ecap,             8, 40, 22.5, 0, thick_ecal_ecap,   shift_ecal_z_plus,    ecalCol  }, //  endcap ECAL +Z
      //{ r_out_ecal_ecap,    r_min_ecal_ecap,             8, 40, 22.5, 0, thick_ecal_ecap,  -shift_ecal_z_minus,   ecalCol  }, //  endcap ECAL -Z

      { r_out_hcal_bar,     r_inn_hcal_bar,             16,  8, 11.25, 11.25, z_max_hcal_bar, - z_max_hcal_bar,      hcalCol,fl+4,0,1}, //  HCAL Barrel
      { r_out_hcal_ring,    r_inn_hcal_ring,             8,  8,  22.5, 0,  thick_hcal_ring,  shift_hcalr_z_plus,  hcalCol,fl+4,0,1}, //  ring HCAL +Z
      { r_out_hcal_ring,    r_inn_hcal_ring,             8,  8,  22.5, 0,  thick_hcal_ring, -shift_hcalr_z_minus, hcalCol,fl+4,0,1} , //  ring HCAL -Z 
      { r_out_hcal_ecap,    r_min_hcal_ecap,             8, 40,  22.5, 0,  thick_hcal_ecap,  shift_hcal_z_plus,   hcalCol,fl+4,0,1}, //  endcap HCAL +Z
      { r_out_hcal_ecap,    r_min_hcal_ecap,             8, 40,  22.5, 0, thick_hcal_ecap, -shift_hcal_z_minus,  hcalCol,fl+4,0,1},  //  endcap HCAL -Z      
      

      //{ coil_outer_radius,  coil_inner_radius,          40, 40,        0.0,   0, coil_half_z, -coil_half_z,  coilCol  },  //  coil     
      { coil_outer_radius,  coil_inner_radius,          200, 200,        0.0,   0, coil_half_z, -coil_half_z,  coilCol,fl+5},  //  coil     

      { r_out_yoke_plug,    r_inn_yoke_plug,            12, 12,        15.0,   0, thick_yoke_plug,  shift_yoker_z_plus,  yokeCol,fl+5,0,0}, //  plug YOKE +Z
      { r_out_yoke_plug,    r_inn_yoke_plug,            12, 12,        15.0,   0, thick_yoke_plug, -shift_yoker_z_minus, yokeCol,fl+5,0,0} , //  plug YOKE -Z 



      { r_max_lhcal,        r_min_lhcal,                40, 40,    0., 0, thick_lhcal,  shift_lhcal_z_plus,   fcalCol,fl+6,0,0}, //    LHCAL +Z
      { r_max_lhcal,        r_min_lhcal,                40, 40,    0., 0, thick_lhcal, -shift_lhcal_z_minus,  fcalCol ,fl+6,0,0},  //   LHCAL -Z      

      { r_max_lcal,         r_min_lcal,                 40, 40,    0., 0, thick_lcal,  shift_lcal_z_plus,   fcalCol,fl+6,0,0}, //    LCAL +Z
      { r_max_lcal,         r_min_lcal,                 40, 40,    0., 0, thick_lcal, -shift_lcal_z_minus,  fcalCol,fl+6,0,0},  //   LCAL -Z      
 
       

      { r_out_yoke_bar,     r_inn_yoke_bar,             12, 12,        15.0,  0, z_max_yoke_bar, - z_max_yoke_bar,      yokeCol,    fl+7,0,0}, //  YOKE Barrel
      { r_out_yoke_ecap,    r_min_yoke_ecap,            12, 12,        15.0,   0, thick_yoke_ecap,  shift_yoke_z_plus,   yokeCol,   fl+7,0,0}, //  endcap YOKE +Z
      { r_out_yoke_ecap,    r_min_yoke_ecap,            12, 12,        15.0,   0, thick_yoke_ecap, -shift_yoke_z_minus,  yokeCol,   fl+7,0,0},  //  endcap YOKE -Z      


      { rSIT[0],          rSIT[0]-0.1 ,                 40, 40,  0.0, 0, lSIT[0],        -lSIT[0],            sitCol,fl+8,0,1}, //  SIT
      { rSIT[1],          rSIT[1]-0.1 ,                 40, 40,  0.0, 0, lSIT[1],        -lSIT[1],            sitCol,fl+8,0,1}  //  SIT
      
    };
   // const DoubleVec& ftd_d   =  pFTD.getDoubleVals("FTDDiskSupportThickness" )  ;
   // const DoubleVec& ftd_ri  =  pFTD.getDoubleVals("FTDInnerRadius" )  ;
   // const DoubleVec& ftd_ro  =  pFTD.getDoubleVals("FTDOuterRadius" )  ;
   // const DoubleVec& ftd_z   =  pFTD.getDoubleVals("FTDZCoordinate" )  ;
// ========================================================================
      ced_geotubes(sizeof(geoTubesANY)/sizeof(CED_GeoTube), geoTubesANY);
      set_layer_description("FTD",fl+0);
      set_layer_description("Beamcal",fl+1);
      set_layer_description("TPC",fl+2);
      set_layer_description("ECAL",fl+3);
      set_layer_description("HCAL",fl+4);
      set_layer_description("Coil",fl+5);
      set_layer_description("LHCAL",fl+6);
      set_layer_description("Yoke",fl+7);
      set_layer_description("SIT", fl+8);
/*
      set_layer_description(" ",fl+9);
      set_layer_description(" ",fl+10);
      set_layer_description(" ",fl+11);
*/

      write_layer_description();
} // End of class DrawGeometry


/*
void ObjectIDHandler::getIDforIndex(col,int index){
}
int  ObjectIDHandler::getObjectforID(int ID){
}
*/

