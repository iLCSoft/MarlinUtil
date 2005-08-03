#ifndef MarlinCED_h
#define MarlinCED_h 1

#include <vector>

#include "marlin/Processor.h"
#include "ced.h"
#include "ced_cli.h"

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
   *  something. Draws the detector (currently only the LDC).
   *  
   */
  static void newEvent( Processor* proc ) ;

  /** To be called by every processor that uses CED in processEvent() after drawing everything. 
   *  Actually draws the event.
   */
  static void draw( Processor* proc ) ;
  

  /** Draw all objects in iterator that have a method getPosition() with the given color marker and 
   *  size in the given layer (default 0) 
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
  static void MarlinCED::drawHelix(float b, float charge, float x, float y, float z,
				   float px, float py, float pz, int marker, int size, int col,
				   float rmin=10., float rmax=3000.0, float zmax=4500.)  ; 
    

//   static void MarlinCED::drawHelixVM(float b, float charge, float x, float y, float z,
// 				   float px, float py, float pz, int marker, int size, int col)  ; 
    
protected:

  MarlinCED() : _first(0) , _last(0) {}
  
  static MarlinCED* _me ;
  
  Processor* _first ;
  Processor* _last ;
  
} ;
#endif



