/***********************************************************************************************
Minimized copy of MarlinCED specified for the use by DDCEDViewer processor (CEDViewer).
Debug comments and unnecessary functionality (for that processor) as well as the GEAR dependence 
are removed from the original code.
It includes helper functions to draw LCIO and the complete draw procedure for the geometry.

author: Thorben Quast (CERN Summer Student 2015)
date: 12 August 2015
***********************************************************************************************/

#ifndef DDMarlinCED_h
#define DDMarlinCED_h 1

#include <cmath>

#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include "marlin/Processor.h"
#include "marlin/Global.h"

#include "ced_cli.h"

#include <lcio.h>
#include <UTIL/LCTypedVector.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/Operators.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <EVENT/Cluster.h>
#include <EVENT/Vertex.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/CalorimeterHit.h>

#include "IMPL/ClusterImpl.h"

#include <MarlinDrawUtil.h>
#include <Trajectory.h>

#include <ctime> //hauke


using namespace marlin ;

//Includes for detector drawing
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h" 
#include "DDRec/DetectorData.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"


struct DDCEDMapParticleObject{const LCObject *obj;  void (*function)(const LCObject *);};

typedef std::map<const int, DDCEDMapParticleObject> DDCEDPickingMap;
typedef std::map<std::string, void (*)(const LCObject *)> DDCEDFunctionMap;


template <class T>
    void printDefault(const LCObject *raw){const T* obj = (T*) raw; std::cout << obj;}

/** Bring the feature to print a LCIO-Objekt given by his ID. 
 *  To Provide this behavior it is required to register all ids first. 
 *  This is an Singelton Class, use it with getInstance(). 
 *  @author Hauke Hoelbe (DESY)
 *  @version May 2010 
 */
class DDCEDPickingHandler{
    private:
        /** The pointer to the DDCEDPickingHandler Objekt 
         */
        static DDCEDPickingHandler *instance;

        /** In this map the id and the DDCEDMapParticleObject is stored. DDCEDMapParticleObject contains a pointer to the LCIO-Object and a function pointer 
         * to print the LCIO-Object 
         */
        static DDCEDPickingMap map;

        /** This map contains the name of the collection or typeName and the function pointer of the function with should be called for this 
          * collection/type. The Update method looks first for the collection name. This means when the collection name overwrite the typeName.
         */
        static DDCEDFunctionMap funcMap;
        //static struct termios orig_termios;

        //DDCEDPickingHandler(); //is not allow to instance this
        //DDCEDPickingHandler(const DDCEDPickingHandler& cc); //not allows to make copies
        //~DDCEDPickingHandler();


    public:
        //class
        /** Returns a pointer to the DDCEDPickingHandler object 
         */
        static DDCEDPickingHandler& getInstance();

        /** The default print method for all LCIO objects 
         */

        template <class T>
        static void defaultOutputFunction(const LCObject *obj){

          streamlog_out(MESSAGE) << " Picked object " << obj << " with id " << obj->id() << std::endl;
          streamlog_out(MESSAGE) << *(T*)obj;

        }

        /** Returns 1 if a key was been pressed, otherwise 0. 
          * (Should be not part of DDCEDPickingHandler)
        */
        static int kbhit(void);
        /*
        static void set_conio_terminal_mode(void);
        static void reset_terminal_mode(void);
        static int getch(void);
        */


        /** Register all LCIO-objects of the given LCEvent in the DDCEDPickingHandler map.
          * This method iterate over all objects how are part of the LCEvent. <br>
          * For each object, the collection name will searched in the funcMap. If this 
          * collection name not found in the funcMap the typeName will be searched. 
          * Now the object pointer, the object id and the given print function from the funcMap 
          * will be stored in the the DDCEDPickingHandler map.<br> 
          * If either the collection name and the typeName not found in the funcMap, this function is not 
          * able to register the id of the object. The result is that print this object with his id is not possible. 
         */
        void update(LCEvent *);
        
        /** Print the LCIO-object given by his ID.
         */
        void printID(int);
        
        /** This method provides to register a (user defined) print function. 
          * The first argument is the collection- or type-Name of the LCIO Object, the second a pointer of the print function.<br>
          * Example: pHandler.registerFunction(LCIO::TRACKERHIT, &yourPrintFunction);
         */
        void registerFunction(std::string, void (*)(const LCObject *));
};
//end hauke hoelbe


/** Singleton class to manage access to CED from several processors. All processors using CED
 *  have to use the methods init(), newEvent() and draw().
 */

class DDMarlinCED {
  
 public:
  static DDMarlinCED* instance() ;
  
  LCEvent* _currEvent=NULL;
  
  /** To be called by every processor that uses CED in intit(). 
   */
  static void init( Processor* proc ) ;

  static void newEvent( Processor* proc , LCEvent* evt=0  ) ;

  /** To be called by every processor that uses CED in processEvent() after drawing everything. 
   *  Actually draws the event. The flag waitForKeyboard indicates if after an event is drawn 
   *  an input from the keyboard is expected (waitForKeyboard=1) or not (waitForKeyboard=0). 
   */
  static void draw( Processor* proc , int waitForKeyboard=1 ) ;
  
  static void getParticleFromID(int, LCEvent*);

  /** Draw all objects in iterator range [first,last) which have a method getPosition() with the given color marker and 
   *  size in the given layer (default 0). The template takes classes providing a class method 'getPosition()' 
   *  as template argument.
   */
  template <class In>
  static void drawObjectsWithPosition(In first, In last, int marker, int size ,unsigned int color, unsigned int layer=0, const char * /*PickingMessage*/="") {
    while( first != last ) {
      int id = (*first)->id();
      ced_hit_ID( (*first)->getPosition()[0],
      (*first)->getPosition()[1],
      (*first)->getPosition()[2],
        marker,layer, size , color, id ) ;

      ++first ;
    }  
  }


//hauke hoelbe
  template <class In>
  static void drawObjectsWithPositionID(LCCollection* /*col*/,In first, In last, int marker, int size ,unsigned int color, unsigned int layer=0) {
    int i=0;
    while( first != last ) {
      int id = (*first)->id(); 
      std::cout << "test!!! " << std::endl;
      ced_hit_ID( (*first)->getPosition()[0],
	       (*first)->getPosition()[1],
	       (*first)->getPosition()[2],
            marker, layer, size , color, id ) ;

      ++first ;
      i++;
    }  
  }

 
  /** Draws a helix from the given point(x,y,z) for momentum(px,py,pz) in a B-field b (in Tesla) 
   */
  static void drawHelix(float b, float charge, float x, float y, float z,
			float px, float py, float pz, int marker, int size, 
			unsigned int col,
			float rmin=10.0, float rmax=3000.0, float zmax=4500.0, unsigned int id = 0);

/*************************/
  static int getIDfromIndex(LCCollection* col, int index);

    static void set_layer_description(const std::string& desc, int layerID);
    static void add_layer_description(const std::string& desc, int layerID);
    static void write_layer_description(void);

    /* Draws the detector geometry for CLIC and ILD. 
     * features:
     * - improved, i.e. more exact, placements 
     * - generic
     * - no GEAR dependence
     * - surface (optionally) drawn as set of lines
     * 
     * author: Thorben Quast, CERN Summer Student 2015
     * date: 31/07/2015
    */
     //TODO: rename
    static void drawDD4hepDetector( dd4hep::Detector& theDetector, bool _surfaces, StringVec _detailled);

private:
    static int _int_count;
    static std::vector<std::string> _descs;


protected:

  //hauke hoelbe: 08.02.2010
  DDMarlinCED() :  _currEvent(NULL), _first(0) , _last(0) {}
  
  static DDMarlinCED* _me ;
  
  Processor* _first=NULL;
  Processor* _last=NULL;
  // helper method to draw hit collections by type
  static void drawHitCollectionsByType(LCEvent* event, const char* type, int marker, int size, 
				       unsigned int color, unsigned int layer=0) {
    
    try {
      
      std::vector< std::string >::const_iterator iter;
      const std::vector< std::string >* ColNames = event->getCollectionNames();
      
      for( iter = ColNames->begin() ; iter != ColNames->end() ; iter++) {
	
	LCCollection* col = event->getCollection( *iter ) ;
	
	if ( col->getTypeName() == type ) {
	  
	  if ( type == LCIO::SIMTRACKERHIT ) {
	    
	    LCTypedVector<SimTrackerHit> v(col);
	    drawObjectsWithPosition(v.begin(),v.end(),marker,size,color,layer);
	    
	  }
	  
	  if ( type == LCIO::SIMCALORIMETERHIT ) {
	    
	    LCTypedVector<SimCalorimeterHit> v(col);
	    drawObjectsWithPosition(v.begin(),v.end(),marker,size,color,layer);
	    
	  }
	  
	  if ( type == LCIO::TRACKERHIT ) {
	    
	    LCTypedVector<TrackerHit> v(col);
	    drawObjectsWithPosition(v.begin(),v.end(),marker,size,color,layer);
	    
	  }
	  
	  if ( type == LCIO::CALORIMETERHIT ) {
	    
	    LCTypedVector<CalorimeterHit> v(col);
	    drawObjectsWithPosition(v.begin(),v.end(),marker,size,color,layer);
	    
	  }
	}
      }
    }    
    catch(DataNotAvailableException &e){}

  }

  // FIXME: Not so elegant, refine! Use iterators, templates etc. See drawHitCollectionsByType(...).
  // helper method to draw hit collections by MC Contribution
  static void drawHitCollectionsByMCContribution(LCEvent* event, MCParticle* MCP, int marker, int size, 
             unsigned int color, unsigned int layer=0) {
    
    std::vector< std::string >::const_iterator iter;
    const std::vector< std::string >* ColNames = event->getCollectionNames();
    
    for( iter = ColNames->begin() ; iter != ColNames->end() ; iter++) {
      
      LCCollection* col = event->getCollection( *iter ) ;
      
      if ( col->getTypeName() == LCIO::SIMTRACKERHIT ) {

      int n = col->getNumberOfElements();
  
      for (int i = 0; i < n; ++i) {

        SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(col->getElementAt(i));
    
        if (hit->getMCParticle() == MCP) {
          double x = hit->getPosition()[0];
          double y = hit->getPosition()[1];
          double z = hit->getPosition()[2];
            ced_hit_ID(x,y,z,marker, layer,size,color, MCP->id());


        }
      }
      }

      if ( col->getTypeName() == LCIO::SIMCALORIMETERHIT ) {

      int n = col->getNumberOfElements();
  
      for (int i = 0; i < n; ++i) {
    
        SimCalorimeterHit* hit = dynamic_cast<SimCalorimeterHit*>(col->getElementAt(i));

        int nMC = hit->getNMCContributions();

        bool found = false;
        for (int j = 0; j < nMC; ++j) {
          if (hit->getParticleCont(j) == MCP) {
            found = true; 
            break;
          }
        }  
        
        if (found) {
          double x = hit->getPosition()[0];
          double y = hit->getPosition()[1];
          double z = hit->getPosition()[2];
            ced_hit_ID(x,y,z,marker,layer,size,color, MCP->id());
        }
      }
      }
    }
  }
} ;

extern "C"
void DDdraw_helix( float b, float charge, float x, float y, float z,
		 float px, float py, float pz, 
		 int marker, int size, unsigned int col, 
		 float rmin=10.0, float rmax=3000.0, float zmax=4500.0, unsigned int id = 0) ;


/******* HELPERS ********/

//read out of the "_detailled" parameter
bool detailledDrawing(StringVec _detailled, std::string detName);

//Set of geometric parameters for initialization of a CEDGeoBox class object
struct CEDGeoBox {
  double  sizes[3] ;
  double  center[3] ;
  double rotate[3];
};
//Set of geometric parameters for initialization of a CEDGeoTube class object
struct CEDGeoTubeParams {
  double Rmax; double Rmin; double inner_symmetry; double outer_symmetry; double phi0; double delta_phi; double delta_z; double z0; 
  //boolean that decides if the GeoTube is drawn twice at two different zPositions
  bool isBarrel;
};

//Convenient summary of both parameter sets above as (tracker) layers may be drawn as one tube or as a sequence of staves (-->GeoBox)
struct LayerGeometry {
  CEDGeoTubeParams tube{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,false};
  std::vector<CEDGeoBox> staves{};
};



void getVisAttributes(dd4hep::DetElement det, unsigned &color, bool &visible);

/***detector draw helpers***/

//converts the parameters in LayeredCalorimeterData given by the appropriate drivers
//into those required by the CEDGeoTube
CEDGeoTubeParams CalorimeterParameterConversion (dd4hep::rec::LayeredCalorimeterData *calo);

//converts the parameters in ZDiskPetalsData given by the appropriate drivers
//into those required by the CEDGeoTube
CEDGeoTubeParams PetalParameterConversion (std::vector<dd4hep::rec::ZDiskPetalsData::LayerLayout>::iterator thisLayer);

//converts the parameters from a LayeredCalorimeterData layer given by the appropriate drivers
//into those required by the CEDGeoTube
CEDGeoTubeParams CalorimeterLayerParameterConversion(std::vector<dd4hep::rec::LayeredCalorimeterData::Layer>::iterator thisLayer);

//converts the parameters from a FixedPadSizeTPCData given by the appropriate drivers
//into those required by the CEDGeoTube
CEDGeoTubeParams TPCParameterConversion(dd4hep::rec::FixedPadSizeTPCData *tpc);

//converts the parameters from a ZPlanarData::LayerLayout layer given by the appropriate drivers
//into those required by the CEDGeoBox (for drawing of staves) or by CEDGeoTube (for approximation of the set of staves into tubes)
LayerGeometry TrackerLayerParameterConversion(std::vector<dd4hep::rec::ZPlanarData::LayerLayout>::iterator thisLayer);

//draws the given surfaces as a set of individual lines
bool DrawSurfaces(const dd4hep::rec::SurfaceManager &surfMan, std::string detName, unsigned color, int layer);

#endif




