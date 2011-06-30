#ifndef MarlinCED_h
#define MarlinCED_h 1

#include <cmath>

#include <vector>
#include <iostream>
#include <iomanip>
#include "marlin/Processor.h"
#include "marlin/Global.h"

#include "ced_cli.h"

#include <lcio.h>
#include <UTIL/LCTypedVector.h>
#include <UTIL/LCTOOLS.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <EVENT/Cluster.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/SimCalorimeterHit.h>

#include "IMPL/ClusterImpl.h"

#include <MarlinDrawUtil.h>
#include <Trajectory.h>

#include <ctime> //hauke


using namespace marlin ;

struct CEDMapParticleObject{const LCObject *obj;  void (*function)(const LCObject *);};

typedef std::map<const int, CEDMapParticleObject> CEDPickingMap;
typedef std::map<std::string, void (*)(const LCObject *)> CEDFunctionMap;


template <class T>
    void printDefault(const LCObject *raw){const T* obj = (T*) raw; std::cout << obj;}

/** Bring the feature to print a LCIO-Objekt given by his ID. 
 *  To Provide this behavior it is required to register all ids first. 
 *  This is an Singelton Class, use it with getInstance(). 
 *  @author Hauke Hoelbe (DESY)
 *  @version May 2010 
 */
class CEDPickingHandler{
    private:
        /** The pointer to the CEDPickingHandler Objekt 
         */
        static CEDPickingHandler *instance;

        /** In this map the id and the CEDMapParticleObject is stored. CEDMapParticleObject contains a pointer to the LCIO-Object and a function pointer 
         * to print the LCIO-Object 
         */
        static CEDPickingMap map;

        /** This map contains the name of the collection or typeName and the function pointer of the function with should be called for this 
          * collection/type. The Update method looks first for the collection name. This means when the collection name overwrite the typeName.
         */
        static CEDFunctionMap funcMap;
        //static struct termios orig_termios;

        //CEDPickingHandler(); //is not allow to instance this
        //CEDPickingHandler(const CEDPickingHandler& cc); //not allows to make copies
        //~CEDPickingHandler();


    public:
        //class
        /** Returns a pointer to the CEDPickingHandler object 
         */
        static CEDPickingHandler& getInstance();

        /** The default print method for all LCIO objects 
         */

        template <class T>
        static void defaultOutputFunction(const LCObject *obj){
            streamlog_out(MESSAGE) << *(T*)obj;
        }

        /** Returns 1 if a key was been pressed, otherwise 0. 
          * (Should be not part of CEDPickingHandler)
        */
        static int kbhit(void);
        /*
        static void set_conio_terminal_mode(void);
        static void reset_terminal_mode(void);
        static int getch(void);
        */


        /** Register all LCIO-objects of the given LCEvent in the CEDPickingHandler map.
          * This method iterate over all objects how are part of the LCEvent. <br>
          * For each object, the collection name will searched in the funcMap. If this 
          * collection name not found in the funcMap the typeName will be searched. 
          * Now the object pointer, the object id and the given print function from the funcMap 
          * will be stored in the the CEDPickingHandler map.<br> 
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

class MarlinCED {
  
 public:
  static MarlinCED* instance() ;
  
  LCEvent* _currEvent;
  
  /** To be called by every processor that uses CED in intit(). 
   */
  static void init( Processor* proc ) ;

  /** To be called by every processor that uses CED in processEvent() before actually drawing 
   *  something. Draws the detector depending on modelID:
   * <ul>
   * <li>0:       use gear file to draw full detector</li>
   * <li>99999:   draw EUTelescope from gear file</li>
   * <li>other:   don't draw a detector model</li>
   * </ul>
   *  
   */
  static void newEvent( Processor* proc , int modelID=0, LCEvent* evt=0  ) ;

  /** To be called by every processor that uses CED in processEvent() after drawing everything. 
   *  Actually draws the event. The flag waitForKeyboard indicates if after an event is drawn 
   *  an input from the keyboard is expected (waitForKeyboard=1) or not (waitForKeyboard=0). 
   */
  static void draw( Processor* proc , int waitForKeyboard=1 ) ;
  
  static void getParticleFromID(int, LCEvent*);

  static void printMCParticle(MCParticle* mcp, int daughterIndent=0, int motherIndent=0);
  void printMCFamily(MCParticle* part, unsigned int daughterBranches, unsigned int motherBranches, 
                     unsigned int daughterIndent=0, unsigned int motherIndent=0);
  void printAndDrawMCFamily(MCParticle* part, LCEvent * evt, unsigned int daughterBranches, unsigned int motherBranches, 
                            unsigned int daughterIndent=0, unsigned int motherIndent=0);


  /** Draw all objects in iterator range [first,last) which have a method getPosition() with the given color marker and 
   *  size in the given layer (default 0). The template takes classes providing a class method 'getPosition()' 
   *  as template argument.
   */
  template <class In>
  static void drawObjectsWithPosition(In first, In last, int marker, int size ,unsigned int color, unsigned int layer=0, const char * PickingMessage="") {
    //std::cout<<"drawObjectsWithPosition register a map id"<< std::endl;
    //int id;
    while( first != last ) {

         //sprintf(pickingMessage, "Position: %f, %f", (*first)->getPosition()[1], (*first)->getPosition()[2]);

         //id = idMap::add(pickingMessage);

         int id = (*first)->id();
      ced_hit_ID( (*first)->getPosition()[0],
	       (*first)->getPosition()[1],
	       (*first)->getPosition()[2],
	       //marker | ( layer << CED_LAYER_SHIFT ) , size , color, id ) ;
            marker,layer, size , color, id ) ;

      ++first ;
    }  
  }


//hauke hoelbe
  template <class In>
  static void drawObjectsWithPositionID(LCCollection* col,In first, In last, int marker, int size ,unsigned int color, unsigned int layer=0) {
    int i=0;
    while( first != last ) {
      int id = (*first)->id(); 
      std::cout << "test!!! " << std::endl;
      //int id = getIDfromIndex(col,i);
      ced_hit_ID( (*first)->getPosition()[0],
	       (*first)->getPosition()[1],
	       (*first)->getPosition()[2],
	       //marker | ( layer << CED_LAYER_SHIFT ) , size , color, id ) ;
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
  
  
  /** Draws a trajectory in the volume described by rmin, rmax, zmax
   */
  static void drawTrajectory(const Trajectory* t, const int marker, 
			     const int size, const unsigned int col,
			     const float rmin=10.0, const float rmax=3000.0, 
			     const float zmax=4500.0, unsigned int id=0) ;

  /** Draws a 'spike', i.e. a bold arrow, from (x0,y0,z0) to (x1,y1,z1) with color on layer e.g. to display jet thrust axes
   */
  static void drawSpike(float x0, float y0, float z0, float x1, float y1, float z1, unsigned int color, unsigned int layer, unsigned int id=0);

  /** Draws the detector using the geometry parameters from GEAR */
  static void drawGEARDetector() ;

  /** Draws the telescope using the geometry parameters from GEAR */
  static void drawGEARTelescope() ;

  /** Draws the detector only. modelID: LDC(0) - default , SID(1), GLD(2), 
   * CaliceTestBeam(3). Needed if you would like to draw the detector 
   * wire frame several times in one event
   * @deprecated
   */
  static void drawDetector(int modelID) {
    
    
    // FIXME: for now have hard coded geometry ....
    // here we would read some geometry description ....
    
    // from V.Morgunov / A.Raspereza: simple outlines of LDC, SID, GLD:
    
    static CED_GeoCylinder geoCylindersLDC[] = {       // for TESLA Detector Geometry
      //      {    50.0,  6,  0.0, 5658.5, -5658.5, 0xff      }, // beam tube
      {   380.0, 24,  0.0, 2658.5, -2658.5, 0xff      }, // inner TPC
      {  1840.0,  8, 22.5, 2700.0, -2700.0, 0x7f7f1f  }, // inner ECAL
      {  2045.7,  8, 22.5, 2700.0, -2700.0, 0x7f7f1f  }, // outer ECAL
      {  2045.7,  8, 22.5, 101.00,  2820.0, 0x7f7f1f  }, // endcap ECAL
      {  2045.7,  8, 22.5, 101.00, -3022.0, 0x7f7f1f  }, // endcap ECAL
      {  3000.0, 16,  0.0, 2658.5, -2658.5, 0xcf00    }, // outer HCAL
      {  3000.0,  8, 22.5, 702.25,  2826.0, 0xcf00    }, // endcap HCAL
      {  3000.0,  8, 22.5, 702.25, -4230.5, 0xcf00    }, // endcap HCAL
	// radius, poligon order, angle degree, 1/2 length, shift in z, color
    }; 
    
    static CED_GeoCylinder geoCylindersSiD[] = {       // for SiD Detector Geometry
      //	{    12.0,  100,  0.0, 5658.5, -5658.5, 0xff      }, // beam tube
	{  186.35,   20,  0.0,  271.0,  -271.0, 0xff      }, // 1st SiD layer
	{  448.85,   20,  0.0,  621.0,  -621.0, 0xff      }, // 2nd SiD layer
	{  711.35,   20,  0.0,  971.0,  -971.0, 0xff      }, // 3d  SiD layer
	{  973.85,   20,  0.0, 1321.0, -1321.0, 0xff      }, // 4th SiD layer
	{ 1236.35,   20,  0.0, 1649.4, -1649.4, 0xff      }, // 5th SiD layer
	{ 1270.00,   30,  0.0, 1682.5, -1682.5, 0x7f7f1f  }, // inner ECAL
	{ 1385.00,   30,  0.0, 1682.5, -1682.5, 0x7f7f1f  }, // outer ECAL
	{ 1385.00,   30,  0.0,  562.5,  1682.5, 0x7f7f1f  }, // endcap ECAL
	{ 1385.00,   30,  0.0,  562.5, -2807.5, 0x7f7f1f  }, // endcap ECAL
	{ 2500.00,   30,  0.0, 1795.0, -1795.0, 0xcf00    }, // outer  HCAL  
	{ 2500.00,   30,  0.0, 487.25 , 1795.0, 0xcf00    }, // endcap HCAL
	{ 2500.00,   30,  0.0, 487.25 ,-2770.0, 0xcf00    }, // endcap HCAL      
    };


    static CED_GeoCylinder geoCylindersGLD[] = {       // for Huge Detector Geometry
      //      {    50.0,  6,  0.0, 5658.5, -5658.5, 0xff      }, // beam tube
      {   400.0, 24,  0.0, 2600.0, -2600.0, 0xff      }, // inner TPC
      {  2100.0, 30, 22.5, 2700.0, -2700.0, 0x7f7f1f  }, // inner CALO
      {  3405.0, 30, 22.5, 2700.0, -2700.0, 0x7f7f1f  }, // outer CAL0
      {  3405.0, 30, 22.5,  737.5,  2700.0, 0x7f7f1f  }, // endcap CAL0
      {  3405.0, 30, 22.5,  737.5, -4175.0, 0x7f7f1f  }, // endcap ECAL
      //      {  3000.0, 16,  0.0, 2658.5, -2658.5, 0xcf00    }, // outer HCAL
      //      {  3000.0,  8, 22.5, 702.25,  2826.0, 0xcf00    }, // endcap HCAL
      //      {  3000.0,  8, 22.5, 702.25, -4230.5, 0xcf00    }, // endcap HCAL
	// radius, poligon order, angle degree, 1/2 length, shift in z, color
    }; 

   
    static CED_GeoCylinder geoCylindersCaliceTestBeam[] = {  // Calice TestBeam Prototype (simple wire frame)
      {  192.25*sqrt(2.0),  4,  45.0,  101.35, -2.0*101.35      , 0x7f7f1f  },    // ECAL (Proto03)
      {  450.00*sqrt(2.0),  4,  45.0,  572.87,        0.00      , 0xcf00    },    // HCAL (TBHcal04)
      {  500.00*sqrt(2.0),  4,  45.0,  457.60,  2.0*572.87+250.0, 0xcf00    },    // TCatcher (TBCatcher04)
	// radius (surrounding (out-circle!) the polygon), poligon order, angle degree, 1/2 length, shift in z, color
    }; 

    if (modelID == 0) {
      ced_geocylinders(sizeof(geoCylindersLDC)/sizeof(CED_GeoCylinder),geoCylindersLDC);
    }
    else if (modelID == 1) {
      ced_geocylinders(sizeof(geoCylindersSiD)/sizeof(CED_GeoCylinder),geoCylindersSiD);
    }
    else if (modelID == 2) {
      ced_geocylinders(sizeof(geoCylindersGLD)/sizeof(CED_GeoCylinder),geoCylindersGLD);
    }
    else if (modelID == 3) {
      ced_geocylinders(sizeof(geoCylindersCaliceTestBeam)/sizeof(CED_GeoCylinder),geoCylindersCaliceTestBeam);
    }


  }
  
  /** Draws a thin line between vertex-point and end-point of a MC particle, another thin line at the vertex-point
   * symbolising the initial momentum vector and all the hits in the SimTrackerHit- and 
   * SimCalorimeterHitCollections which are produced by this MC particle, if toggle drawSimHits is true. 
   * These SimHits are drawn with a marker of size and color. All objects are drawn on the same layer.
   */
  static void drawMCParticle(MCParticle* MCP, bool drawSimHits, LCEvent* event, int marker, int size, 
			     unsigned int color, unsigned int layer=0, double bField = 4.0, 
			     double rmin = 0.0, double zmin = 0.0, double rmax = 3000.0, double zmax = 4500.0, 
			     bool drawOnDifferentLayers = true);
  
  /** Draws the full MC Particle Tree on the layers 1, 2 and 3. On layer 1 the charged particles are displayed and on 
   * layer shift-1 their corresponding SimHits. On layer 2 all the neutral particles are shown, again with their 
   * SimHits on layer shift-2. The layers 3 and shift-3 are used for the backscattered MC Particles. 
   * The variables energyCut, rIn, zIn, rOut and zOut are cut values meant to reduce the number particles displayed; 
   * i.e. only particles with energy larger than 'energyCut' and within the cylinder described 
   * by rIn, zIn, rOut, zOut are shown.
   */
  static void drawMCParticleTree(LCEvent* event, std::string colNameMC, double energyCut,  
				 double bField = 4.0, double rIn = 0.0, double zIn = 0.0, 
				 double rOut = 3000.0 , double zOut = 4500.0);

  /** Draws the hits of all SimTrackerHit Collections of event with a marker of size and color on layer
   */
  static void drawSimTrackerHits(LCEvent* event, int marker, int size, unsigned int color, unsigned int layer=0);

  /** Draws the hits of all SimCalorimeterHit Collections of event with a marker of size and color on layer
   */
  static void drawSimCalorimeterHits(LCEvent* event, int marker, int size, unsigned int color, unsigned int layer=0);

  /** Draws the hits of all SimHit Collections of event with a marker of size and color on layer
   */
  static void drawSimHits(LCEvent* event, int marker, int size, unsigned int color, unsigned int layer=0);

  /** Draws the hits of all TrackerHit Collections of event with a marker of size and color on layer
   */
  static void drawTrackerHits(LCEvent* event, int marker, int size, unsigned int color, unsigned int layer=0);

  /** Draws the hits of all CalorimeterHit Collections of event with a marker of size and color on layer
   */
  static void drawCalorimeterHits(LCEvent* event, int marker, int size, unsigned int color, unsigned int layer=0);

  /** Draws the hits of all Hit Collections of event with a marker of size and color on layer
   */
  static void drawHits(LCEvent* event, int marker, int size, unsigned int color, unsigned int layer=0);

  /** Draws the hits of a track with a marker of size and color on layer
   */
  static void drawTrack(Track* track, int marker, int size, unsigned int color, unsigned int layer=0);

  /** Draws the hits of a cluster with a marker of size and color on layer
   */
  static void drawCluster(Cluster* cluster, int marker, int size, unsigned int color, unsigned int layer=0);

  /** Draws the hits of a cluster implementation with a marker of size and color on layer
   */
  static void drawClusterImpl(const ClusterImpl* cluster, int marker, int size, unsigned int color, unsigned int layer=0);

  /** Draws the hits of a recontructed particle reco with a marker of size and color on layer
   */
  static void drawRecoParticle(ReconstructedParticle* reco, int marker, int size, unsigned int color, unsigned int layer=0);


/*************************/
  static int getIDfromIndex(LCCollection* col, int index);
//{
//      _int_count++;
//      std::cout<<"Registration" << "col.getTypeName(): " << col->getTypeName() << " id: " << index << " test i:" << _int_count << std::endl;
//      return(1);
//}

    static void set_layer_description(const std::string& desc, int layerID);
    static void add_layer_description(const std::string& desc, int layerID);
    static void write_layer_description(void);

private:
    static int _int_count;
    static std::vector<std::string> _descs;


protected:

  //hauke hoelbe: 08.02.2010
  MarlinCED() : _first(0) , _last(0){ _currEvent=0; }
  
  static MarlinCED* _me ;
  
  Processor* _first ;
  Processor* _last ;
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
    	    //ced_hit_ID(x,y,z,marker | ( layer << CED_LAYER_SHIFT ),size,color, MCP->id());
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
    	    //ced_hit_ID(x,y,z,marker | ( layer << CED_LAYER_SHIFT ),size,color, MCP->id());
            ced_hit_ID(x,y,z,marker,layer,size,color, MCP->id());

    
    	  }
    	}
      }
    }
  }


  
} ;



#endif




