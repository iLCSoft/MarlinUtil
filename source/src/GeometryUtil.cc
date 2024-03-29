#include "GeometryUtil.h"

#include <streamlog/streamlog.h>

#include <DD4hep/DD4hepUnits.h>

double MarlinUtil::getBzAtOrigin() {

  double bfield(0.0);

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
  if ( not (theDetector.state() == dd4hep::Detector::READY) ) {
    throw std::runtime_error("Detector geometry not initialised, cannot get bfield");
  }
  const double position[3]={0,0,0}; // position to calculate magnetic field at (the origin in this case)
  double magneticFieldVector[3]={0,0,0}; // initialise object to hold magnetic field
  theDetector.field().magneticField(position,magneticFieldVector); // get the magnetic field vector from DD4hep
  bfield = magneticFieldVector[2]/dd4hep::tesla; // z component at (0,0,0)
  return bfield;

}


dd4hep::rec::LayeredCalorimeterData const* MarlinUtil::getLayeredCalorimeterData(unsigned int includeFlag,
                                                                                 unsigned int excludeFlag) {

  dd4hep::rec::LayeredCalorimeterData * theExtension = 0;

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
  std::vector<dd4hep::DetElement> const& theDetectors = dd4hep::DetectorSelector(theDetector).detectors(  includeFlag, excludeFlag ) ;

  if(  theDetectors.size() > 0 ){
    streamlog_out(DEBUG) << " getLayeredCalorimeterData :  includeFlag: " << dd4hep::DetType( includeFlag )
                         << " excludeFlag: " << dd4hep::DetType( excludeFlag )
      //size is > 0 so we can safely use at(0)
                         << "  found : " << theDetectors.size() << "  - first det: " << theDetectors.at(0).name() << std::endl ;
  }

  if( theDetectors.size()  != 1 ){
    std::stringstream es;
    es << " getLayeredCalorimeterData: selection is not unique (or empty) includeFlag: "
       << dd4hep::DetType( includeFlag ) << " excludeFlag: " << dd4hep::DetType( excludeFlag )
       << " --- found detectors : " ;
    for( unsigned i=0, N= theDetectors.size(); i<N ; ++i ){
      es << theDetectors.at(i).name() << ", " ;
    }
    throw std::runtime_error( es.str() ) ;
  }

  //size is 1 or we would have exited before, so we can safely use at(0)
  theExtension = theDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();

  return theExtension;
}
