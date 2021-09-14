#include "ILDGeometryUtil.h"
#include "GeometryUtil.h"

dd4hep::rec::ZPlanarData* MarlinUtil::getVXDData(){
    return getDetData<dd4hep::rec::ZPlanarData>("VXD");
}


dd4hep::rec::ZPlanarData* MarlinUtil::getSITData(){
    return getDetData<dd4hep::rec::ZPlanarData>("SIT");
}


dd4hep::rec::ZDiskPetalsData* MarlinUtil::getFTDData(){
    return getDetData<dd4hep::rec::ZDiskPetalsData>("FTD");
}


dd4hep::rec::FixedPadSizeTPCData* MarlinUtil::getTPCData(){
    return getDetData<dd4hep::rec::FixedPadSizeTPCData>("TPC");
}


dd4hep::rec::ZPlanarData* MarlinUtil::getSETData(){
    return getDetData<dd4hep::rec::ZPlanarData>("SET");
}
