#ifndef ILDGEOMETRYUTIL_H
#define ILDGEOMETRYUTIL_H 1

#include <DDRec/DetectorData.h>


namespace MarlinUtil{
    /**
    * Returns DDRec detector extension ZPlanarData for vertex detector of the ILD
    * (link to the dd4hep documentation of the class so user knows which parameters he get get from this...)
    */
    dd4hep::rec::ZPlanarData* getVXDData();


    /**
    * Returns DDRec detector extension ZPlanarData for SIT detector of the ILD
    * (link to the dd4hep documentation of the class so user knows which parameters he get get from this...)
    */
    dd4hep::rec::ZPlanarData* getSITData();


    /**
    * Returns DDRec detector extension ZDiskPetalsStruct for FTD detector of the ILD
    * (link to the dd4hep documentation of the class so user knows which parameters he get get from this...)
    */
    dd4hep::rec::ZDiskPetalsData* getFTDData();


    /**
    * Returns DDRec detector extension FixedPadSizeTPCData for TPC detector of the ILD
    * (link to the dd4hep documentation of the class so user knows which parameters he get get from this...)
    */
    dd4hep::rec::FixedPadSizeTPCData* getTPCData();


    /**
    * Returns DDRec detector extension ZPlanarData for SET detector of the ILD
    * (link to the dd4hep documentation of the class so user knows which parameters he get get from this...)
    */
    dd4hep::rec::ZPlanarData* getSETData();
}

#endif
