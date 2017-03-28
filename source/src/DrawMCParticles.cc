#include "DrawMCParticles.h"


using namespace lcio ;
using namespace marlin ;


DrawMCParticles aDrawMCParticles ;


DrawMCParticles::DrawMCParticles() : Processor("DrawMCParticles") {
  
  _description = "DrawMCParticles draws the MC Particle prediction w/o calorimeters in the ced display" ;
  


  registerProcessorParameter( "MCCollectionName" , 
			      "Name of the MCParticle collection"  ,
			      _colNameMC ,
			      std::string("MCParticle") ) ;


  // FIXME: take this information from GEAR file
  registerProcessorParameter( "rIn","Radius of the innermost detector component (VTX) in cylindrical coordinates",
			      _rIn,
			      double(15.5) ) ;
  // FIXME: take this information from GEAR file
  registerProcessorParameter( "zIn","z coordinate of the innermost detector component (VTX) in cylindrical coordinates",
			      _zIn,
			      double(50.0) ) ;

  /*
  // FIXME: take this information from GEAR file
  // FIXME: take calo face instead of outer TPC radius
  registerProcessorParameter( "rOut","Radius of the outermost detector component (TPC) inside the calorimeters in cylindrical coordinates",
			      _rOut,
			      double(1690.0) ) ;
  // FIXME: take this information from GEAR file
  registerProcessorParameter( "zOut","z coordinate of the outermost detector component (TPC) inside the calorimeters inside the calorimeters in cylindrical coordinates",
			      _zOut,
			      double(2000.0) ) ;

  */



  registerProcessorParameter( "EnergyCut",
			      "Energy Cut in GeV",
			      _energyCut,
			      double(0.01) ) ;
  
  registerProcessorParameter("WaitForKeyboard","Wait for Keyboard before proceed",
			     _waitForKeyboard,
			     (int)1);


}


void DrawMCParticles::init() { 

  // usually a good idea to
  printParameters() ;

  //FIXME: constant offset
  const double offset = 0.0;

  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;

  _bField = gearTPC.getDoubleVal("BField");
  _rOut = padLayout.getPlaneExtent()[1]+offset;
  _zOut = gearTPC.getMaxDriftLength()+offset;

   

  MarlinCED::init(this);


  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void DrawMCParticles::processRunHeader( LCRunHeader* ) {

  _nRun++ ;
} 

void DrawMCParticles::processEvent( LCEvent * evt ) { 
  
  // Reset drawing buffer and START drawing collection
  MarlinCED::newEvent(this);

  MarlinCED::drawMCParticleTree(evt,_colNameMC,_energyCut,_bField,_rIn,_zIn,_rOut,_zOut);

  MarlinCED::draw(this,_waitForKeyboard);

  _nEvt ++ ;
}



void DrawMCParticles::check( LCEvent* ) {

}


void DrawMCParticles::end(){ 

}
