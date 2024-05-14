#ifndef SelectNthEventsProcessor_h
#define SelectNthEventsProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;


/** Select a number n, every nth event will have the processor ReturnValue set to true, all others to false.
 *  The selection can be inverted and the starting event can be offset from event 0.
 *  Use in steering file in execute section via
 *  <if condition="MySelectNthEventsProcessor">
 *    <processor name="MyAnalysisProcessor">
 *  </if>
 *
 * @param nEventToSelect: Integer number n. Every nth event's ReturnValue is set to true.
 * @param SelectionOffset: Sets an offset m, so the Return Value of the mth event (and every nth before and after it) is set to true.
 * @param InvertSelection: Inverts the ReturnValue. If true, every nth event's ReturnValue is false, all others are true.
 * 
 * @author U Einhaus, DESY
 * @version $Id: SelectNthEventsProcessor.h, v1.0 2024-01-21 $
 */

class SelectNthEventsProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new SelectNthEventsProcessor ; }
  
  
  SelectNthEventsProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* ) {};
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent* );
  
  
  virtual void check( LCEvent* ) {};
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  int _nRun=0;
  int _nEvt=0;
  int _nEvtProc=0;
  int _nEvtSkip=0;
  
  int _nEventToSelect=1;
  int _selectionOffset=0;
  bool _invertSelection=false;

} ;

#endif



