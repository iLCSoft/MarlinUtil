#ifndef SelectEvents_h
#define SelectEvents_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;


/** Select a range of events by setting the boolean 'SelectEvents' to true 
 *  for the events processed by this processor in the range 
 *  FirstEvent-LastEvent.
 * 
 * @param FirstEvent  first processed event where 'SelectEvents' is true
 * @param LastEvent   last processed event where 'SelectEvents' is true
 * 
 * @author O. Wendt, DESY
 * @version $Id: SelectEvents.h,v 1.2 2008-04-16 15:03:03 gaede Exp $ 
 */

class SelectEvents : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new SelectEvents ; }
  
  
  SelectEvents() ;
  
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

  int _nRun=-1;
  int _nEvt=-1;
  
  int _firstEvent=0;
  int _lastEvent=0;

} ;

#endif



