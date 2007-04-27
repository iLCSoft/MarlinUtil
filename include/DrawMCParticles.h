#ifndef DrawMCParticles_h
#define DrawMCParticles_h 1

#include <iostream>
#include <string>

#include "marlin/Processor.h"
#include "lcio.h"

#include <EVENT/LCCollection.h>
#include <MarlinCED.h>

// GEAR include files
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>


using namespace lcio ;
using namespace marlin ;


/**  
 * 
 * @author O. Wendt, DESY
 * @version $Id: DrawMCParticles.h,v 1.1 2007-04-27 14:04:29 owendt Exp $ 
 */

class DrawMCParticles : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new DrawMCParticles ; }
  
  
  DrawMCParticles() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  /** Input collection name.
   */
  std::string _colNameMC ;


  double _rIn;
  double _zIn;
  double _rOut;
  double _zOut;
  double _energyCut;

  double _bField;
  int _waitForKeyboard;

  int _nRun ;
  int _nEvt ;
} ;

#endif



