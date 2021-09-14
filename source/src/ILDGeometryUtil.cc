#include "ILDGeometryUtil.h"
#include "GeometryUtil.h"

dd4hep::rec::ZPlanarData* MarlinUtil::ILD::getVXDData(){
    return getDetData<dd4hep::rec::ZPlanarData>("VXD");
}


dd4hep::rec::ZPlanarData* MarlinUtil::ILD::getSITData(){
    return getDetData<dd4hep::rec::ZPlanarData>("SIT");
}


dd4hep::rec::ZDiskPetalsData* MarlinUtil::ILD::getFTDData(){
    return getDetData<dd4hep::rec::ZDiskPetalsData>("FTD");
}


dd4hep::rec::FixedPadSizeTPCData* MarlinUtil::ILD::getTPCData(){
    return getDetData<dd4hep::rec::FixedPadSizeTPCData>("TPC");
}


dd4hep::rec::ZPlanarData* MarlinUtil::ILD::getSETData(){
    return getDetData<dd4hep::rec::ZPlanarData>("SET");
}
