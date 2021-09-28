#ifndef ILDGEOMETRYUTIL_H
#define ILDGEOMETRYUTIL_H 1

#include <DDRec/DetectorData.h>

/** Namespace for the ILD specific geometry functions.

MarlinUtil::ILD contains helper functions which provide easier interface for extracting
geometry information of the ILD detector elements.

They are based on the template function MarlinUtil::getDetData(const std::string& detName)
which can be applied to extract the detector element extension for any generic detector model.
*/
namespace MarlinUtil::ILD{
    /**Get dd4hep::rec::StructExtension of the dd4hep::rec::ZPlanarStruct for the VXD detector.
    */
    dd4hep::rec::ZPlanarData* getVXDData();


    /**Get dd4hep::rec::StructExtension of the dd4hep::rec::ZPlanarStruct for the SIT detector.
    */
    dd4hep::rec::ZPlanarData* getSITData();


    /**Get dd4hep::rec::StructExtension of the dd4hep::rec::ZDiskPetalsStruct for the FTD detector.
    */
    dd4hep::rec::ZDiskPetalsData* getFTDData();


    /**Get dd4hep::rec::StructExtension of the dd4hep::rec::FixedPadSizeTPCStruct for the TPC detector.
    */
    dd4hep::rec::FixedPadSizeTPCData* getTPCData();


    /**Get dd4hep::rec::StructExtension of the dd4hep::rec::ZPlanarStruct for the SET detector.
    */
    dd4hep::rec::ZPlanarData* getSETData();
}

#endif
