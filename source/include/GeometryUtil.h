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

}//namespace MarlinUtil

#endif // GEOMETRYUTIL_H
