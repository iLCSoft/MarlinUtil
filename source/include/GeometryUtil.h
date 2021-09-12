#ifndef GEOMETRYUTIL_H
#define GEOMETRYUTIL_H 1

#include <DDRec/DetectorData.h>

namespace MarlinUtil {


  /**
   * Returns the bfield value in Z direction at (0 0 0),
   *
   * Obtains value from DD4hep (lcdd) Throws an exception if geometry is not
   * instantiated correctly
   */
  double getBzAtOrigin();


  /**
   * Returns DDRec detector extension
   *
   * include/excludeFlags must describe unique detector, otherwise exception is thrown
   * @param includeFlag : bitmask describing detector to select
   * @param excludeFlag : bitmask describing detector to exclude (optional)
   *
   */
  dd4hep::rec::LayeredCalorimeterData const* getLayeredCalorimeterData(unsigned int includeFlag,
                                                                       unsigned int excludeFlag=0);


   /**
   * Returns DDRec detector extension ZPlanarData for vertex detector
   * (link to the dd4hep documentation of the class so user knows which parameters he get get from this...)
   */
   dd4hep::rec::ZPlanarData* getVXDData();


   /**
   * Returns DDRec detector extension ZPlanarData for SIT detector
   * (link to the dd4hep documentation of the class so user knows which parameters he get get from this...)
   */
   dd4hep::rec::ZPlanarData* getSITData();


   /**
   * Returns DDRec detector extension ZDiskPetalsStruct for FTD detector
   * (link to the dd4hep documentation of the class so user knows which parameters he get get from this...)
   */
   dd4hep::rec::ZDiskPetalsData* getFTDData();


   /**
   * Returns DDRec detector extension FixedPadSizeTPCData for TPC detector
   * (link to the dd4hep documentation of the class so user knows which parameters he get get from this...)
   */
   dd4hep::rec::FixedPadSizeTPCData* getTPCData();


   /**
   * Returns DDRec detector extension ZPlanarData for SET detector
   * (link to the dd4hep documentation of the class so user knows which parameters he get get from this...)
   */
   dd4hep::rec::ZPlanarData* getSETData();









}//namespace MarlinUtil

#endif // GEOMETRYUTIL_H
