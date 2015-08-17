#include "DDMarlinCED.h"

#include <LCGeometryTypes.h>
#include "ced_cli.h"
 
//hauke//
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
//#include "UTIL/CellIDDecoder.h"
#include "UTIL/PIDHandler.h"

using namespace UTIL;

// make gcc > 4.7 compliant
#include <unistd.h>//

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



//for detector drawing (Thorben Quast)
#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h" 
#include "DDRec/DetectorData.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"
#include "TColor.h"
using namespace DD4hep::Geometry ;
using namespace DD4hep;
using namespace DD4hep::DDRec ;


DDMarlinCED* DDMarlinCED::_me = 0;

//hauke hoelbe
DDCEDPickingMap DDCEDPickingHandler::map;
DDCEDFunctionMap DDCEDPickingHandler::funcMap;
DDCEDPickingHandler *DDCEDPickingHandler::instance = NULL;

std::vector<std::string> DDMarlinCED::_descs(CED_MAX_LAYER, ""); //layer descriptions

DDCEDPickingHandler& DDCEDPickingHandler::getInstance() {
    if( !instance){
        std::cout<<"new instance"<<std::endl;
        instance = new DDCEDPickingHandler();
    }
    return *instance;
}

void DDCEDPickingHandler::update(LCEvent *evt){
    const std::vector<std::string> *collNames = evt->getCollectionNames();
    const LCCollection *coll;
    const LCObject *obj;
    DDCEDFunctionMap::iterator iter;
    std::map<std::string, void (*)(const LCObject *)> printDefaultMap;
    std::map<std::string, void (*)(const LCObject *)>::iterator printDefaultIter;
    clock_t start =  clock();
    
    //default print functions
    printDefaultMap.insert(std::make_pair(LCIO::MCPARTICLE,&defaultOutputFunction<MCParticle>));
    printDefaultMap.insert(std::make_pair(LCIO::TRACKERHIT, &defaultOutputFunction<TrackerHit>));
    printDefaultMap.insert(std::make_pair(LCIO::TRACKERHITPLANE, &defaultOutputFunction<TrackerHitPlane>));
    printDefaultMap.insert(std::make_pair(LCIO::TRACKERHITZCYLINDER, &defaultOutputFunction<TrackerHitZCylinder>));
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
    DDCEDMapParticleObject particleObj;
    for(unsigned int i=0;i<collNames->size();i++){
        collName=collNames->at(i);
        coll =  evt->getCollection(collName);
        typeName=coll->getTypeName();
        for (int j=0;j<coll->getNumberOfElements();j++){
            obj=coll->getElementAt(j);
            particleObj.obj=obj;
            
            iter =  funcMap.find(collName);
            if(iter != funcMap.end()){ //user have registered a function for this collection name
                particleObj.function=iter->second;
            }else{
                iter=funcMap.find(typeName);
                if(iter != funcMap.end()){ //user have a registered a function for this _typeName_
                    particleObj.function=iter->second;
                }else{
                    streamlog_out(DEBUG) << "DDCEDPickingHandler: cant register " << collName << "/" << typeName
                    << " (no function given)" << std::endl;
                    continue;
                }
            }
            map.insert(std::pair<const int, DDCEDMapParticleObject>(obj->id(),particleObj));
        }
    }
    clock_t end = clock() ;
    streamlog_out(DEBUG) << "DDCEDPickingHandler::Map size: " << map.size() << " time: " << double( end - start ) / double(CLOCKS_PER_SEC) << "s" << std::endl;
}

void DDCEDPickingHandler::printID(int id){
    DDCEDPickingMap::iterator iter;
    DDCEDMapParticleObject obj;
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

void DDCEDPickingHandler::registerFunction(std::string type, void (*printFunction)(const LCObject *)){
    funcMap.insert(std::make_pair(type,printFunction));
}

int DDCEDPickingHandler::kbhit(void) {
    //http://stackoverflow.com/questions/448944/c-non-blocking-keyboard-input#448982
    struct timeval tv = { 0L, 0L };
    fd_set fds;
    FD_SET(0, &fds);
    return select(1, &fds, NULL, NULL, &tv);
}


void DDMarlinCED::add_layer_description(const std::string &desc, int layerID){
    std::string tmp;
    if(layerID > CED_MAX_LAYER || layerID < 0){return;}
    if( _descs.at(layerID).find(desc.c_str()) == std::string::npos){
        tmp=_descs.at(layerID);
        if(! tmp.empty()){
            tmp.append(", ");
        }
        tmp.append(desc);
        _descs.at(layerID)=tmp;
    }else{
    }
}

void DDMarlinCED::set_layer_description(const std::string &desc, int layerID){
    if(layerID > CED_MAX_LAYER || layerID < 0){return;}
    _descs.at(layerID)=desc;
    
}

void DDMarlinCED::write_layer_description(void){
    //std::cout<<"LAYER: write all layer in ced" << std::endl;
    unsigned int i;
    //for(i=0;i<25;i++){
    for(i=0; i<_descs.size(); i++){
        ced_describe_layer(_descs.at(i).c_str(), i);
    }
}

//end hauke hoelbe

DDMarlinCED* DDMarlinCED::instance() {
    if( _me == 0 )
        _me = new DDMarlinCED ;
    return _me ;
}


void DDMarlinCED::init( Processor* proc ) {
    
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


void DDMarlinCED::newEvent( Processor* proc , LCEvent* evt) {
    if( proc == instance()->_first ) {
        ced_new_event();
        if(evt!=0)
            instance()->_currEvent = evt;    
    }
}

//hauke hoelbe modify 08.02.2010
void DDMarlinCED::draw( Processor* proc , int waitForKeyboard ) {
    int i=0;
    DDCEDPickingHandler &pHandler=DDCEDPickingHandler::getInstance();
    
    
    if( proc == instance()->_last ) {
        //    ced_draw_event();
        DDMarlinCED::write_layer_description();
        //ced_picking_text("test1 test2 test3");
        
        ced_send_event();
        if ( waitForKeyboard == 1 ) {
            streamlog_out(MESSAGE) << "Double click for picking. Press <ENTER> for the next event." << std::endl;
            //test:
            
            signal(SIGWINCH,SIG_IGN);
            
            while(!DDCEDPickingHandler::kbhit()){
                //            while(!poll(pfd,1,0)){
                
                usleep(100000); //micro seconds
                
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
            
            signal(SIGWINCH,SIG_IGN);
            char c = getchar();
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
void DDMarlinCED::drawHelix(float b, float charge, float x, float y, float z,
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
            }
        }
        
        int nSteps = int(100/step); //hauke
        
        
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
                break ;
            }
            
            if( r_current >= (rmin+step)) {
                count_lines++;
                ced_line_ID( x1, y1, z1, x2, y2, z2 , marker , size, col, id);
            }
            x1 = x2;
            y1 = y2;
            z1 = z2;
            
        }
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
            streamlog_out(ERROR) << "Error in 'DDMarlinCED::drawHelix()': Startpoint beyond (rmax,zmax)" << std::endl;
            return;
        }
        
        kr = summand + sqrt(radicant);
        kz = ((zmax-z)*absP)/pz;
        
        // this has been improved
        
        if (z + (kr*pz)/absP > zmax || z + (kr*pz)/absP < -zmax){
            k = kz;
        }
        else k = kr;
        
        if (k < 0.0 ) {
            streamlog_out(DEBUG2) << "DDMarlinCED::drawHelix(): negative intersection parameter - will revert sign ... "
            << std::endl;
            //fg: k cannot be negativ ( particle is moving along its 3-momentum ....)
            k = -k ;
        }
        
        float xEnd = x + (k*px)/absP;
        float yEnd = y + (k*py)/absP;
        float zEnd = z + (k*pz)/absP;
        
        if (rmin != 0){
            streamlog_out(DEBUG) << "FIX ME: Inner cylinder not taken into account!" << std::endl;
            return;
        }
        
        
        streamlog_out(DEBUG1) << "DDMarlinCED::drawHelix()' - pt : " << pt << " |p| = " << absP
        << ", x " << x
        << ", y " << y
        << ", z " << z
        << ", px " << px
        << ", py " << py
        << ", pz " << pz
        << ", xEnd " << xEnd
        << ", yEnd " << yEnd
        << ", zEnd " << zEnd
        << std::endl ;
        
        
        ced_line_ID(x, y, z, xEnd, yEnd, zEnd , marker , size, col, id);
        
    }
    else {
        streamlog_out(DEBUG) << "Low momentum particle given point instead of helix" << std::endl;
        const double delta = 0.0001;
        ced_line_ID(x, y, z, x+delta, y+delta, z+delta, marker , size, col, id);
    }
}


void DDMarlinCED::drawDD4hepDetector( DD4hep::Geometry::LCDD& lcdd, bool _surfaces, StringVec _detailled){
  typedef std::vector< DD4hep::Geometry::DetElement> DetVec ;
  // get DetElements for the main sub detectors from dd4hep 
  const DetVec& trackers     = lcdd.detectors( "tracker" ) ;
  const DetVec& calorimeters = lcdd.detectors( "calorimeter" ) ;
  const DetVec& passiveDets  = lcdd.detectors( "passive" ) ;
  //allocate reference to the surface manager
  DD4hep::DDRec::SurfaceManager& surfMan = *lcdd.extension<DD4hep::DDRec::SurfaceManager>();
  //some temporary parameters for visualization
  unsigned color; bool visible;
  //temporary objects
  DetElement det; std::string detName;

  std::vector<CEDGeoTube> gTV ; 
  int detLayer = NUMBER_DATA_LAYER; 


  //--- loop over all calorimeters
  for( unsigned i=0,n=calorimeters.size() ; i<n ; ++i ){
    det = calorimeters[i] ;
    detName = det.name() ;
    LayeredCalorimeterData* calo = 0 ;
    streamlog_out( MESSAGE ) << " ......processing " << detName << std::endl;  
    //try to get the appropriate extension
    try{ 
      calo = det.extension<LayeredCalorimeterData>() ; 
    } catch(std::runtime_error& e){
      streamlog_out( MESSAGE ) <<  detName 
             << " has no extension of type LayeredCalorimeterData. "
             <<   std::endl;  
    }

    //get the visAttributes of the detElement's volume or (if not existing) some default values
    getVisAttributes(det, color, visible);
    int layer = detLayer++;
    bool isDrawn = _surfaces && DrawSurfaces(surfMan, detName, color, layer);

    //draw if the object exists and its visibility is set true
    if(visible && !isDrawn) {                    
      if (calo != 0){
        //get the required parameter set for drawing a CEDGeoTube
        CEDGeoTubeParams params;    
        params = CalorimeterParameterConversion(calo);
        //consistently allocated this helper for linking the element correctly to the GUI control via ID
        gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry, params.phi0, params.delta_phi, params.delta_z,  params.z0, color , layer  ,1,1 ) ) ; 
        if (!params.isBarrel){
          //place the second one symmetric to the z-axis. An additional shift by the width of the layer is needed since the appropriate argument represents the left handed start of the geometry
          gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry, params.phi0, params.delta_phi, params.delta_z,  - (params.z0+2.0*params.delta_z), color , layer  ,1,1 ) ) ; 
        }
        isDrawn = true;
      }
      if (!calo)
        isDrawn = DrawSurfaces(surfMan, detName, color, layer);
    }
    if(isDrawn){
      DDMarlinCED::set_layer_description( detName , layer );
      streamlog_out( MESSAGE ) << detName << " has successfully been processed."<< std::endl;
    }
    std::cout<<std::endl;
  }
  //--- draw trackers 
  //mostly repetition of the calorimeter 
  for( unsigned i=0,n=trackers.size() ; i<n ; ++i ){
    det = trackers[i] ;
    detName = trackers[i].name() ;
    ZPlanarData* trkPlanar = 0; ZDiskPetalsData* trkDisk = 0; FixedPadSizeTPCData* trkTPC = 0;
    streamlog_out( MESSAGE ) << " ......processing" <<  detName << std::endl; 
    
    try{ 
      trkPlanar = det.extension<ZPlanarData>();
    } catch(std::runtime_error& e){
      try{
        trkDisk = det.extension<ZDiskPetalsData>();
      }catch(std::runtime_error& e){
        try{
          trkTPC = det.extension<FixedPadSizeTPCData>();
        } catch(std::runtime_error& e){
            streamlog_out( MESSAGE ) <<  detName 
             << " has no extension of type ZPlanarData/ZDiskPetalsData. "
             <<   std::endl;           
        }
      }
    }

    int layer = detLayer++;  
    getVisAttributes(det, color, visible);
    bool isDrawn = _surfaces && DrawSurfaces(surfMan, detName, color, layer);
    //get the visAttributes of the detElement's volume or (if not existing) some default values
    if (!isDrawn && visible){
      //the following if statements are exclusive, i.e. only one may apply
      if(trkPlanar){
        for (std::vector<DDRec::ZPlanarData::LayerLayout>::iterator thisLayer = trkPlanar->layers.begin(); thisLayer != trkPlanar->layers.end(); thisLayer++){
          LayerGeometry Geo;
          Geo = TrackerLayerParameterConversion(thisLayer);
          if (detailledDrawing(_detailled, detName)){
            for( unsigned stave_i=0; stave_i<Geo.staves.size() ; ++stave_i ){
              ced_geobox_r_ID( Geo.staves[stave_i].sizes, Geo.staves[stave_i].center, Geo.staves[stave_i].rotate, color, layer,0);
              ced_geobox_r_solid( Geo.staves[stave_i].sizes, Geo.staves[stave_i].center, Geo.staves[stave_i].rotate, color, layer);
            }
          }
          else{ 
            gTV.push_back( CEDGeoTube( Geo.tube.Rmax, Geo.tube.Rmin, Geo.tube.inner_symmetry, Geo.tube.outer_symmetry, Geo.tube.phi0, Geo.tube.delta_phi, Geo.tube.delta_z,  Geo.tube.z0, color , layer  ,1,1 ) );
          }
        }
        isDrawn = true;
      }

      if(trkDisk){ 
        if (detailledDrawing(_detailled, detName)){
          streamlog_out( MESSAGE )<<detName<<": Not drawn for now (appropriate geometry does not exist)"<<std::endl;
        }
        else{
          for (std::vector<DDRec::ZDiskPetalsData::LayerLayout>::iterator thisLayer = trkDisk->layers.begin(); thisLayer != trkDisk->layers.end(); thisLayer++){
            CEDGeoTubeParams params;    
            params = PetalParameterConversion(thisLayer);
            gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry, params.phi0, params.delta_phi, params.delta_z,  params.z0, color , layer  ,1,1 ) ) ; 
            //place the second one symmetric to the z-axis. An additional shift by the width of the layer is needed since the appropriate argument represents the left handed start of the geometry
            gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry, params.phi0, params.delta_phi, params.delta_z,  -(params.z0+2*params.delta_z), color , layer  ,1,1 ) ) ; 
            DDMarlinCED::set_layer_description( detName , layer );
          }
        }
        isDrawn = true;
      }
      
      if(trkTPC) {
       CEDGeoTubeParams params;    
       params = TPCParameterConversion(trkTPC);
       gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry, params.phi0, params.delta_phi, params.delta_z,  params.z0, color , layer  ,1,1 ) ) ; 
       isDrawn = true;
      }

      if(!trkTPC && !trkDisk && !trkPlanar)//if no simplified geometry is given, try surfaces
        isDrawn = DrawSurfaces(surfMan, detName, color, layer);
      
    }
    if(isDrawn){
      DDMarlinCED::set_layer_description( detName , layer );
      streamlog_out( MESSAGE ) << detName << " has successfully been processed."<< std::endl;
    }
    std::cout<<std::endl;
  }


  //--- draw passive 
  for( unsigned i=0,n=passiveDets.size() ; i<n ; ++i ){
    det = passiveDets[i] ;
    detName = passiveDets[i].name() ;
    ConicalSupportData* passiveConical = 0; LayeredCalorimeterData* passiveCalo = 0;
    streamlog_out( MESSAGE ) << " ......processing " <<  detName << std::endl;  
    try{ 
        passiveConical = det.extension<ConicalSupportData>();
    } catch(std::runtime_error& e){
      try{
        passiveCalo = det.extension<LayeredCalorimeterData>();
      } catch(std::runtime_error& e){
          streamlog_out( MESSAGE ) <<  detName 
             << " has no extension of type ConicalSupportData/LayeredCalorimeterData. "
             <<   std::endl;  
      }
    }
    //get the visAttributes of the detElement's volume or (if not existing) some default values
    getVisAttributes(det, color, visible);
    int layer = detLayer++;
    bool isDrawn = _surfaces && DrawSurfaces(surfMan, detName, color, layer);
    if (!isDrawn && visible){
      if (passiveConical){
        streamlog_out( MESSAGE )<<detName<<" is not drawn for now (not needed)."<<std::endl;
        for (std::vector<DDRec::ConicalSupportData::Section>::iterator thisSection = passiveConical->sections.begin(); thisSection != passiveConical->sections.end(); thisSection++){
        }
        isDrawn = true;
      }
      if (passiveCalo){
        for (std::vector<DDRec::LayeredCalorimeterData::Layer>::iterator thisLayer = passiveCalo->layers.begin(); thisLayer != passiveCalo->layers.end(); thisLayer++){
          CEDGeoTubeParams params;
          params = CalorimeterLayerParameterConversion(thisLayer);
          gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry,  params.phi0, params.delta_phi, params.delta_z,  params.z0, color , layer  ,1,1 ) ) ; 
        }
        isDrawn = true;
      }
      if (!passiveConical && !passiveCalo)
        isDrawn = DrawSurfaces(surfMan, detName, color, layer);
    }
    if (isDrawn){
      DDMarlinCED::set_layer_description( detName , layer );
      streamlog_out( MESSAGE ) << detName << " has successfully been processed."<< std::endl;
    }
    std::cout<<std::endl;
  }
  // ========================================================================
  //Draw the tubes:
  ced_geotubes( gTV.size() ,  (CED_GeoTube*) &gTV[0] );

  // ========================================================================
  
  DDMarlinCED::write_layer_description();
}

void DDdraw_helix( float b, float charge, float x, float y, float z,
                float px, float py, float pz,
                int marker, int size, unsigned int col,
                float rmin, float rmax, float zmax, unsigned int id){
    
    
    DDMarlinCED::drawHelix( b,  charge,  x,  y,  z, px,  py,  pz,  marker,  size, col, rmin,  rmax,  zmax,  id) ;
    
}



/***detector draw helpers***/

bool detailledDrawing(StringVec _detailled, std::string detName){
  if(_detailled.size() == 0) 
    return false;
  unsigned index = 0;
  while(index < _detailled.size()){
    if (detName.compare(_detailled[index])==0)
      return true;
    index++;
  }
  return false;
}

//converts the parameters in LayeredCalorimeterData given by the appropriate drivers
//into those required by the CEDGeoTube
CEDGeoTubeParams CalorimeterParameterConversion (LayeredCalorimeterData *calo){
  //get all the information from the lcdd class
  double rMin = calo->extent[0]/dd4hep::mm ;
  double rMax = calo->extent[1]/dd4hep::mm ;
  //attention! the meaning of these two variables depends on their context 
  //(see e.g. ECalEndcap_o2_v01_geo.cpp lines 78/79 vs. ECalBarrel_o1_v01.cpp lines 74/75)
  double zMin = calo->extent[2]/dd4hep::mm ; 
  double zMax = calo->extent[3]/dd4hep::mm ;
  
  double inner_symmetry = calo->inner_symmetry;
  double outer_symmetry = calo->outer_symmetry;
  //Note that the phi0's are given in degrees rather than rad in this implementation
  double inner_phi0 = calo->inner_phi0*180/M_PI;
  double outer_phi0 = calo->outer_phi0*180/M_PI;
  //new convention to prevent weird overlaps
  int NCircle = 36;
  //correct for the case of circle approximation, the given is conventional
  if (outer_symmetry < 1)  outer_symmetry = NCircle;
  if (inner_symmetry < 1)  inner_symmetry = NCircle;

  CEDGeoTubeParams returnParams;
  //given: distance middle point - center of edge section
  //required: distance middle point - edge intersection (corner) to prevent weird overlaps
  returnParams.Rmax = rMax/cos(M_PI/outer_symmetry); 
  returnParams.Rmin = rMin/cos(M_PI/inner_symmetry);
  
  //nothing to convert
  returnParams.inner_symmetry = inner_symmetry;
  returnParams.outer_symmetry = outer_symmetry;
  
  //by default, CED draws tubes with the corner facing downward
  //what we want: phi0 = 0 <=> straight edge parallel to x-y plane
  //therefore, the tube must be rotated by 360 - (90 + 360./(2*number of sides)) degrees since 
  //respecting the rotation symmetry in 360/n, the minimal phi0 is calculated to prevent interference with phi cuts in the CEDViewer
  //phi0 == 0 implies a normal vector of the first cell parallel to the +x-axis 
  //but CED starts drawing symmetrically from the +y-axis
  returnParams.phi0 = outer_phi0 + 270. - 180./outer_symmetry;
  returnParams.phi0 = returnParams.phi0 - (360./outer_symmetry)*(int (returnParams.phi0/(360./outer_symmetry))); 
  //in ILD and CLIC, both inner and outer shapes agree. For a generic solution, another parameter in LayeredCalorimeterData (e.g. phi_inner) should be introduced
  //i.e. delta_phi that is the rotational angle of the inner angle with respect to the outer layer
  returnParams.delta_phi = -returnParams.phi0 + inner_phi0 + 270. - 180./inner_symmetry;
  returnParams.delta_phi = returnParams.delta_phi - (360./inner_symmetry)*(int (returnParams.delta_phi/(360./inner_symmetry)));

  

  //endcaps and barrels take the same parameters in CED but the interpretation of the z-coordinates differs:
  returnParams.isBarrel = calo->layoutType == LayeredCalorimeterData::BarrelLayout;
  
  //type specific conversions
  if (returnParams.isBarrel){
    //barrels are drawn centered at z=0, i.e. they start at -zMax (z0) and their half length (delta_z) is zMax
    returnParams.delta_z = zMax;
    returnParams.z0 = -zMax;
  }
  else {
    //In the drivers, zMin is given to be the start of one disk. (Symmetry at the z-axis is assumed during placement.)
    //zMax is passed as zMin + thickness such that half the distance is correctly obtained by the line below.
    returnParams.delta_z = 0.5*(zMax - zMin);
    returnParams.z0 = zMin;
  }
  return returnParams;
}
//converts the parameters in ZDiskPetalsData given by the appropriate drivers
//into those required by the CEDGeoTube
CEDGeoTubeParams PetalParameterConversion (std::vector<DDRec::ZDiskPetalsData::LayerLayout>::iterator thisLayer){

  double phi0 = thisLayer->phi0*180/M_PI;
  double distanceSensitive = thisLayer->distanceSensitive/dd4hep::mm;
  //Supposedly, this is the actual width of the ring. The expansion along the z-axis is hard coded to 0.3.
  double lengthSensitive = thisLayer->lengthSensitive/dd4hep::mm;
  double zPosition = thisLayer->zPosition/dd4hep::mm;
  double thicknessSupport = thisLayer->thicknessSupport/dd4hep::mm;
  int petalNumber = thisLayer->petalNumber;

  CEDGeoTubeParams returnParams;
  /*   
  //Old implementation with GEAR
  returnParams.Rmax = distanceSensitive;
  returnParams.Rmin = distanceSensitive-1.0*lengthSensitive;
  */
  returnParams.Rmax = distanceSensitive+1.0*lengthSensitive;
  returnParams.Rmin = distanceSensitive;

  //(see comment in line 151)
  returnParams.Rmax = returnParams.Rmax/cos(M_PI/petalNumber);
  returnParams.Rmin = returnParams.Rmin/cos(M_PI/petalNumber);
  //Number of edges = number of petals
  returnParams.inner_symmetry = petalNumber;
  returnParams.outer_symmetry = petalNumber;
  //(see comment in line 159)
  returnParams.phi0 = phi0 + 270. - 180./petalNumber - (360./petalNumber)*(int ((270.-180./petalNumber)/(360./petalNumber)));
  returnParams.delta_phi = 0.0;
  
  //thicknessSupport is negligibly small
  returnParams.delta_z = thicknessSupport; 
  //Again: z0 is the left handed starting point for drawing
  returnParams.z0 = zPosition-returnParams.delta_z;
  return returnParams;
}

//converts the parameters from a LayeredCalorimeterData layer given by the appropriate drivers
//into those required by the CEDGeoTube
CEDGeoTubeParams CalorimeterLayerParameterConversion(std::vector<DDRec::LayeredCalorimeterData::Layer>::iterator thisLayer){
  double distance = thisLayer->distance/dd4hep::mm;
  double thickness = thisLayer->thickness/dd4hep::mm;
  double cellSize0 = thisLayer->cellSize0/dd4hep::mm;
  double cellSize1 = thisLayer->cellSize1/dd4hep::mm;
  
  int NCircle = 36;   //hard coded number of edges to form a circle

  CEDGeoTubeParams returnParams;
  returnParams.Rmax = distance + 1.0*thickness;
  returnParams.Rmin = distance;
  //(see comment in line 151)
  returnParams.Rmax = returnParams.Rmax/cos(M_PI/NCircle);
  returnParams.Rmin = returnParams.Rmin/cos(M_PI/NCircle);
  //assume round edges
  returnParams.inner_symmetry = NCircle;
  returnParams.outer_symmetry = NCircle;
  returnParams.phi0 = 270. - 180./NCircle - (360./NCircle)*(int ((270.-180./NCircle)/(360./NCircle)));
  returnParams.delta_phi = 0.0;

  //cellSize1 is half the length along the z-axis (see e.g. Solenoid_o1_v01_gep.cpp line 58/59)
  returnParams.delta_z = cellSize1;
  //cellSize0 is defined to be the middle point of the geometry 
  returnParams.z0 = -cellSize1 + cellSize0;

  return returnParams;
}

//converts the parameters from a FixedPadSizeTPCData given by the appropriate drivers
//into those required by the CEDGeoTube
CEDGeoTubeParams TPCParameterConversion(FixedPadSizeTPCData *tpc){
  double zHalf = tpc->zHalf/dd4hep::mm;
  //these radii include the insensitive space!
  double rMin = tpc->rMin/dd4hep::mm;
  double rMax = tpc->rMax/dd4hep::mm;
  
  int NCircle = 36;   //hard coded number of edges to form a circle

  CEDGeoTubeParams returnParams;

  returnParams.Rmax = rMax;
  returnParams.Rmin = rMin;
  //(see comment in line 151)
  returnParams.Rmax = returnParams.Rmax/cos(M_PI/NCircle);
  returnParams.Rmin = returnParams.Rmin/cos(M_PI/NCircle);
  //assume round edges
  returnParams.inner_symmetry = NCircle;
  returnParams.outer_symmetry = NCircle;
  returnParams.phi0 = 270. - 180./NCircle - (360./NCircle)*(int ((270.-180./NCircle)/(360./NCircle)));
  returnParams.delta_phi = 0.0;

  
  returnParams.delta_z = zHalf;
  returnParams.z0 = -zHalf;

  return returnParams;
}

//converts the parameters from a ZPlanarData::LayerLayout layer given by the appropriate drivers
//into those required by the CEDGeoBox (for drawing of staves) or by CEDGeoTube (for approximation of the set of staves into tubes)
LayerGeometry TrackerLayerParameterConversion(std::vector<DDRec::ZPlanarData::LayerLayout>::iterator thisLayer){
  int nLadders = thisLayer->ladderNumber;
  double phi0 = thisLayer->phi0*180/M_PI;
  
  double distance_sensitive = thisLayer->distanceSensitive/dd4hep::mm ;
  double thickness_sensitive = thisLayer->thicknessSensitive/dd4hep::mm ;
  double width_sensitive = thisLayer->widthSensitive/dd4hep::mm ;
  double offset_sensitive = thisLayer->offsetSensitive/dd4hep::mm ;
  double zHalf_sensitive = thisLayer->zHalfSensitive/dd4hep::mm ;  
  
  LayerGeometry Geometry;
  //two possibilites to draw: 
  //1) Draw layer consisting of staves
  double currentPhi; double radius;  double deltaPhi = 360./nLadders;
  for (int i=0; i<nLadders; i++){
    //placement shall begin along -z axis like for all other geometries
    currentPhi = phi0 + i*deltaPhi;
    
    //distance_sensitive is passed to be the minimal radius by the driver;
    //for drawing, it should be converted to the radius of the middle line
    radius = distance_sensitive+0.5* thickness_sensitive;

    CEDGeoBox stave;
    //place the center of the box at the appropriate coordinates
    //offset_sensitive is an additonal shift in placement orthogonal to the original 2D radius vector
    stave.center[0] = (radius*cos(currentPhi*M_PI/180) - offset_sensitive*sin(currentPhi*M_PI/180));
    stave.center[1] = (radius*sin(currentPhi*M_PI/180) + offset_sensitive*cos(currentPhi*M_PI/180));
    stave.center[2] = 0.0;  //placed z=0
    //dimensions are straight forward in an xyz coordinate system
    stave.sizes[0]  = thickness_sensitive;
    stave.sizes[1]  = width_sensitive;
    stave.sizes[2]  = zHalf_sensitive * 2;
    //the individual staves are finally rotated in the coordinate system to their appropriate position
    //herby, the rotation is performed around the z-axis
    stave.rotate[0] = 0.0;
    stave.rotate[1] = 0.0;
    stave.rotate[2] = currentPhi;

    Geometry.staves.push_back(stave);
  }


  //2) Summarize the set of staves into a geotube
  //analogous conversions as for the CalorimeterLayerParameterConversion (see above)
  Geometry.tube.Rmax = (distance_sensitive+thickness_sensitive)/cos(M_PI/nLadders); 
  Geometry.tube.Rmin = distance_sensitive/cos(M_PI/nLadders);
  Geometry.tube.inner_symmetry = nLadders;
  Geometry.tube.outer_symmetry = nLadders;
  Geometry.tube.phi0 = phi0 + 270. - 180./nLadders - (360./nLadders)*(int ((270.-180./nLadders)/(360./nLadders)));
  Geometry.tube.delta_phi = 0.0;
  Geometry.tube.delta_z = zHalf_sensitive;
  Geometry.tube.z0 = - zHalf_sensitive;

  //return both possible sets of parameters describing the geometry. The choice for either is implemented in the main draw routine
  return Geometry;
}

//draws the given surfaces as a set of individual lines from indicated start- to the endpoint
bool DrawSurfaces(DD4hep::DDRec::SurfaceManager &surfMan, std::string detName, unsigned color, int layer){
  typedef DD4hep::DDRec::SurfaceMap SMap;
  const SMap* sMap = surfMan.map(detName);
  int lineCounter = 0;
  if(sMap) {
    for (SMap::const_iterator it = sMap->begin(); it != sMap->end(); ++it){
      DD4hep::DDRec::Surface* surf = dynamic_cast<DD4hep::DDRec::Surface*> (it->second);
      if (!surf) continue;
      if (!(surf->type().isVisible())) continue;
      const std::vector<std::pair<Vector3D,Vector3D> > lines = surf->getLines();
      if (lines.empty()){
        streamlog_out( MESSAGE )<<" **** drawSurfaces(): empty lines vector for surface "<< *surf <<std::endl;
        continue;
      }
      for(unsigned i = 0; i <lines.size(); i++){
        unsigned default_width = 2;
        unsigned default_type = layer;
        ced_line(lines[i].first.x()/dd4hep::mm ,lines[i].first.y()/dd4hep::mm ,lines[i].first.z()/dd4hep::mm ,
                lines[i].second.x()/dd4hep::mm ,lines[i].second.y()/dd4hep::mm ,lines[i].second.z()/dd4hep::mm ,
                default_type,default_width, color);
        ced_line_ID(lines[i].first.x()/dd4hep::mm ,lines[i].first.y()/dd4hep::mm ,lines[i].first.z()/dd4hep::mm ,
                lines[i].second.x()/dd4hep::mm ,lines[i].second.y()/dd4hep::mm ,lines[i].second.z()/dd4hep::mm ,
                default_type,default_width, color, layer);
        lineCounter++; 
      }
    }
  } 
  if (lineCounter > 0)
    streamlog_out( MESSAGE )<<"Surfaces have been used."<<std::endl;
  return lineCounter > 0; //at least one line must have been drawn, otherwise the surface is considered to be undrawn
}



void getVisAttributes(DD4hep::Geometry::DetElement det, unsigned &color, bool &visible) {
    DD4hep::Geometry::VisAttr thisVisAttribute = det.volume().visAttributes();  
    if (thisVisAttribute.isValid()){
      TColor* c = new TColor(thisVisAttribute.color(), 1, 1, 1);
      //convert the given TColor into a hexadecimal whereby R_i, G_i, B_i integers
      //color = 0x00|R_1 R_2|G_1 G_2|B_1 B_2|
      color = ((int(255*c->GetRed())<<16) |  //multiply with 2^16 to get the last two bits
        (int(255*c->GetGreen())<<8)|          //multiply with 2^8 to get the middle bits
        (int(255*c->GetBlue())<<0));          //multiply with 2^0 to get the first two bits
                                              //the | operator is a bitwise addition of the number
      delete c;
    }
    else{
      streamlog_out( MESSAGE )<<"color: pointer does not exist"<<std::endl;
        color  = 8947848; //== 0xff999999
        visible = true;   
    }
    //TODO: Hard coded to make every element visible for now
    visible = true;
}
