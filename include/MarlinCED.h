#ifndef MarlinCED_h
#define MarlinCED_h 1

#include <cmath>

#include <vector>

#include "marlin/Processor.h"
#include "ced.h"
#include "ced_cli.h"

#include <lcio.h>
#include <UTIL/LCTypedVector.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <EVENT/Cluster.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/SimCalorimeterHit.h>



using namespace marlin ;


/** Singleton class to manage access to CED from several processors. All processors using CED
 *  have to use the methods init(), newEvent() and draw().
 */
class MarlinCED {
  
 public:
  static MarlinCED* instance() ;
  
  
  /** To be called by every processor that uses CED in intit(). 
   */
  static void init( Processor* proc ) ;

  /** To be called by every processor that uses CED in processEvent() before actually drawing 
   *  something. Draws the detector with modelID: LDC(0) - default , SID(1), GLD(2).
   *  
   */
  static void newEvent( Processor* proc , int modelID=0  ) ;

  /** To be called by every processor that uses CED in processEvent() after drawing everything. 
   *  Actually draws the event. The flag waitForKeyboard indicates if after an event is drawn an input from the keyboard is expected (waitForKeyboard=1) or not (waitForKeyboard=0). 
   */
  static void draw( Processor* proc , int waitForKeyboard=1 ) ;
  

  /** Draw all objects in iterator range [first,last) which have a method getPosition() with the given color marker and 
   *  size in the given layer (default 0). The template takes classes providing a class method 'getPosition()' as template argument.
   */
  template <class In>
  static void drawObjectsWithPosition(In first, In last, int marker, int size ,int color, int layer=0) {
    while( first != last ) {
      ced_hit( (*first)->getPosition()[0],
	       (*first)->getPosition()[1],
	       (*first)->getPosition()[2],
	       marker | ( layer << CED_LAYER_SHIFT ) , size , color ) ;
      ++first ;
    }  
  }
  

  /** Draws a helix from the given point(x,y,z) for momentum(px,py,pz) in a B-field b (in Tesla) 
   */
  static void drawHelix(float b, float charge, float x, float y, float z,
			float px, float py, float pz, int marker, int size, int col,
			float rmin=10., float rmax=3000.0, float zmax=4500.)  ; 
    

  //  static void drawHelixVM(float b, float charge, float x, float y, float z,
  //	                      float px, float py, float pz, int marker, int size, int col)  ; 


  /** Draws a 'spike', i.e. a bold arrow, from (x0,y0,z0) to (x1,y1,z1) with color on layer e.g. to display jet thrust axes
   */
  static void drawSpike(float x0, float y0, float z0, float x1, float y1, float z1, unsigned int color, unsigned int layer);

  /** Draws the detector only. modelID: LDC(0) - default , SID(1), GLD(2). Needed if you would like to draw the detector wire frame several times in one event
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

    /*
      static CED_GeoCylinder geoCylinders[] = {    // for Prototype
      {    180.0,  4,  45.0, 110.0, 0.0, 0xff },   // beam tube
      {    500.0,  4,  45.0, 250.0, 220., 0xff }   // inner TPC
      };
    */

    if (modelID == 0) {
      ced_geocylinders(sizeof(geoCylindersLDC)/sizeof(CED_GeoCylinder),geoCylindersLDC);
    }
    else if (modelID == 1) {
      ced_geocylinders(sizeof(geoCylindersSiD)/sizeof(CED_GeoCylinder),geoCylindersSiD);
    }
    else if (modelID == 2) {
      ced_geocylinders(sizeof(geoCylindersGLD)/sizeof(CED_GeoCylinder),geoCylindersGLD);
    }

  }
  
  /** Draws a thin line between vertex-point and end-point of a MC particle, another thin line at the vertex-point symbolising the initial momentum vector and all the hits in the SimTrackerHit- and SimCalorimeterHitCollections which are produced by this MC particle, if toggle drawSimHits is true. These SimHits are drawn with a marker of size and color. All objects are drawn on the same layer.
   */
  static void drawMCParticle(MCParticle* MCP, bool drawSimHits, LCEvent* event, int marker, int size, int color, int layer=0, double BField = 4.0);

  /** Draws the hits of all SimTrackerHit Collections of event with a marker of size and color on layer
   */
  static void drawSimTrackerHits(LCEvent* event, int marker, int size, int color, int layer=0);

  /** Draws the hits of all SimCalorimeterHit Collections of event with a marker of size and color on layer
   */
  static void drawSimCalorimeterHits(LCEvent* event, int marker, int size, int color, int layer=0);

  /** Draws the hits of all SimHit Collections of event with a marker of size and color on layer
   */
  static void drawSimHits(LCEvent* event, int marker, int size, int color, int layer=0);

  /** Draws the hits of all TrackerHit Collections of event with a marker of size and color on layer
   */
  static void drawTrackerHits(LCEvent* event, int marker, int size, int color, int layer=0);

  /** Draws the hits of all CalorimeterHit Collections of event with a marker of size and color on layer
   */
  static void drawCalorimeterHits(LCEvent* event, int marker, int size, int color, int layer=0);

  /** Draws the hits of all Hit Collections of event with a marker of size and color on layer
   */
  static void drawHits(LCEvent* event, int marker, int size, int color, int layer=0);

  /** Draws the hits of a track with a marker of size and color on layer
   */
  static void drawTrack(Track* track, int marker, int size, int color, int layer=0);

  /** Draws the hits of a cluster with a marker of size and color on layer
   */
  static void drawCluster(Cluster* cluster, int marker, int size, int color, int layer=0);

  /** Draws the hits of a recontructed particle reco with a marker of size and color on layer
   */
  static void drawRecoParticle(ReconstructedParticle* reco, int marker, int size, int color, int layer=0);

protected:

  MarlinCED() : _first(0) , _last(0) {}
  
  static MarlinCED* _me ;
  
  Processor* _first ;
  Processor* _last ;

  // helper method to draw hit collections by type
  static void drawHitCollectionsByType(LCEvent* event, const char* type, int marker, int size, int color, int layer=0) {

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


  // FIXME: Not so elegant, refine! Use iterators, templates etc. See drawHitCollectionsByType(...).
  // helper method to draw hit collections by MC Contribution
  static void drawHitCollectionsByMCContribution(LCEvent* event, MCParticle* MCP, int marker, int size, int color, int layer=0) {
    
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
	    ced_hit(x,y,z,marker | ( layer << CED_LAYER_SHIFT ),size,color );

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
	    ced_hit(x,y,z,marker | ( layer << CED_LAYER_SHIFT ),size,color );

	  }
	}
      }
    }
    
  }


  
} ;
#endif



