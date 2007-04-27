#include <iostream>
#include <assert.h>
#include <typeinfo>
#include <math.h>

#include "Phys_Geom_Database.h"

#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gearimpl/TPCParametersImpl.h>
#include <gearimpl/FixedPadSizeDiskLayout.h>
#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>

#include <Elements.h>

using namespace std;
using namespace marlin;

/**  
 *    Physical Geometrical database will be used in the reconstruction
 *    procedure to get fast access to values/variables dependent on
 *    the particular and different calorimeter zones with different
 *    physical properties and geometrical relations between cells. 
 *
 *    Such a database can be created after reading the geometry 
 *    record in LCIO if it will be presented in or another kind of 
 *           outstanding geometrical database like GEAR.
 *
 *    Another way is to create such database during simulation
 *      step as an obligatory for any simulation program, 
 *             during the detector geometry creation.
 *       That is much easy way to form such a database.
 *            and put it as member (or class) of LCIO.
 *
 *
 *
 *
 *       PGdb::ZONE PGdb::get_zone(Point3D &p)
 *
 *     function returns zone "number" (enum type)
 *          for any arbitrary 3-d point
 *
 *
 *    @author V. L. Morgunov, A. Zhelezov  (DESY/ITEP)
 * 
 **/

PGdb PGDB;

PGDBP aPGDB;

//============================================================================
PGDBP::PGDBP() : Processor("PGDBP") {
//============================================================================
  _description = "Physical and Geometrical Database (for Boojum at least)" ;
}
//============================================================================
void PGDBP::init() {             //the begin of the job
//============================================================================
  PGDB.init();
}
//============================================================================
void PGdb::Zone::_init_tpc() {
//============================================================================
  const gear::TPCParameters  &pTPC      = Global::GEAR->getTPCParameters();
  const gear::PadRowLayout2D &padLayout = pTPC.getPadLayout();
  const gear::DoubleVec      &planeExt  = padLayout.getPlaneExtent();
  no=TPC;
  shape = CYLINDER; 
  symmetry  = 0;
  phi0      = 0.;
  r_inner   = planeExt[0]; 
  r_outer   = planeExt[1]; 
  z_inner   = 0.0; 
  z_outer   = pTPC.getMaxDriftLength(); 
  n_sampl   = padLayout.getNRows (); 
  sampling  = (r_outer- r_inner)/n_sampl; 
  min_lay   = 0; 
  max_lay   = n_sampl-1; 
  absorber  = 0.0; 
  abs_mat   = 0; 
  detector  = sampling; 
  det_mat   = 18;       //  Argon has been taken for the moment
  cell_size = sampling; // formaly 
//---------------------------------
  mip_vis   = 30.0e-9 ; // Should be in GEAR
  e_coeff   = 1.0     ; // Should be in GEAR
//---------------------------------
}
//============================================================================
void PGdb::Zone::_init_endpl1(){
//============================================================================
  const gear::TPCParameters  &pTPC      = Global::GEAR->getTPCParameters();
  const gear::PadRowLayout2D &padLayout = pTPC.getPadLayout();
  const gear::DoubleVec      &planeExt  = padLayout.getPlaneExtent();
  no        = ENDPLATE1;
  shape     = CYLINDER;
  symmetry  = 0;
  phi0      = 0.;
  r_inner   = planeExt[0];
  r_outer   = planeExt[1];
  z_inner   =  pTPC.getMaxDriftLength();
  z_outer   =  pTPC.getMaxDriftLength()+160.0;
 //-----------------added KP
  Zeff      =    21.01;
  Aeff      =    45.09;
  Rhoeff    =   0.2643;
  Ieff      = 145.13e-9;
  Rmeff     = 504.54;// mm
  x0eff     = 419.29 ;  // mm
  Eceff     = 17.628 ; // MeV  
  _init_all_common();
}
//============================================================================
void PGdb::Zone::_init_endpl2(){
//============================================================================
  const gear::TPCParameters  &pTPC      = Global::GEAR->getTPCParameters();
  const gear::PadRowLayout2D &padLayout = pTPC.getPadLayout();
  const gear::DoubleVec      &planeExt  = padLayout.getPlaneExtent();
  no        = ENDPLATE2;
  shape     = CYLINDER;
  symmetry  = 0;
  phi0      = 0.;
  r_inner   = planeExt[0];
  r_outer   = planeExt[1];
  z_inner   = -pTPC.getMaxDriftLength();
  z_outer   = -(pTPC.getMaxDriftLength()+160.0);
 //-----------------dodao KP
  Zeff      =    21.01;
  Aeff      =    45.09;
  Rhoeff    =   0.2643;
  Ieff      = 145.13e-9;
  Rmeff     = 504.54;// mm
  x0eff     = 419.29 ;  // mm
  Eceff     = 17.628 ; // MeV
  _init_all_common();
}
//============================================================================
static int _ecal_last_layer(const gear::CalorimeterParameters& pCAL){
//============================================================================
  int i;
  const gear::LayerLayout &lb = pCAL.getLayerLayout() ;
  int nLayerb = lb.getNLayers() ;
  double tlb  = lb.getThickness(0);
  double dlb  = lb.getAbsorberThickness(0);
  for( i=1 ; i < nLayerb ; i++ )
    if(abs(tlb-lb.getThickness(i))>0.0000001 ||
       abs(dlb-lb.getAbsorberThickness(i))>0.0000001)
      return i;
  cerr<<"ERROR: Can't get boundary of CAL"<<endl;
  return i;
}
//============================================================================
void PGdb::Zone::_init_ecal_bar_common(){
//============================================================================
  const gear::CalorimeterParameters& pECAL_B = 
    Global::GEAR->getEcalBarrelParameters();
  shape    = POLYGON;
  symmetry = pECAL_B.getSymmetryOrder();
  phi0     = pECAL_B.getPhi0();
  z_inner  = 0.0;
  z_outer  = pECAL_B.getExtent()[3];
  abs_mat  = 74;  // Tungsten Should be in GEAR
  detector = 0.5; // No chance to extract  Should be in GEAR
  det_mat  = 14;  // Silicon Should be in GEAR
}
//============================================================================
void PGdb::Zone::_init_ecal1_bar(int last_layer){
//============================================================================
  const gear::CalorimeterParameters& pECAL_B = 
    Global::GEAR->getEcalBarrelParameters();
  const gear::LayerLayout &lb = pECAL_B.getLayerLayout() ;
  no        = ECAL1_BAR;
  r_inner   = lb.getDistance(0);
  r_outer   = lb.getDistance(last_layer); 
  n_sampl   = last_layer; 
  sampling  = lb.getThickness(0); 
  min_lay   = 0;
  max_lay   = last_layer-1; 
  absorber  = lb.getAbsorberThickness(0); 
  cell_size = lb.getCellSize0(0); // should be 10 mm in GEAR
//---------------------------------
  mip_vis   = 170.0e-6 ;  // Should be in GEAR
  e_coeff   = 33.02346 ;  // for model LDC00
  // e_coeff   = 46.7703 ;  // for model LDC01  Should be in GEAR
  //-----------------added KP
  Zeff      =    67.41;
  Aeff      =  166.868;
  Rhoeff    =     7.75;
  Ieff      =706.67e-9;
  Rmeff     = 22.0833;  // mm
  x0eff     = 9.4587 ;  // mm
  Eceff     = 9.0804  ; // MeV
  eprime    = 0.72002;

  _init_ecal_bar_common();
}
//============================================================================
void PGdb::Zone::_init_ecal2_bar(int last_ecal1_layer){
//============================================================================
  const gear::CalorimeterParameters& pECAL_B = 
    Global::GEAR->getEcalBarrelParameters();
  const gear::LayerLayout &lb = pECAL_B.getLayerLayout() ;
  int nLayerb = lb.getNLayers() ;
  no        = ECAL2_BAR;
  r_inner   = lb.getDistance(last_ecal1_layer); 
  r_outer   = lb.getDistance(nLayerb-1)+lb.getThickness(nLayerb-1); 
  n_sampl   = nLayerb-last_ecal1_layer; 
  sampling  = lb.getThickness(nLayerb-1); 
  min_lay   = last_ecal1_layer;
  max_lay   = nLayerb-1; 
  absorber  = lb.getAbsorberThickness(nLayerb-1); 
  cell_size = lb.getCellSize0(nLayerb-1); // should be 10 mm
//---------------------------------
  mip_vis   = 170.0e-6 ;  // Should be in GEAR
  e_coeff   = 93.56822 ;  //for model LDC00
  // e_coeff   =  86.3749 ;  //for model LDC01
 //-----------------added KP
  Zeff      =   71.635;
  Aeff      =  177.614;
  Rhoeff    =   12.577;
  Ieff      = 729.09e-9; //GeV
  Rmeff     = 13.8418;// mm
  x0eff     = 8.135953 ;  // mm
  Eceff     = 8.4709  ; // MeV  dou
  eprime    = 0.709672;

  _init_ecal_bar_common();
}
//============================================================================
void PGdb::Zone::_init_ecal_cap_common(){
//============================================================================
  const gear::CalorimeterParameters& pECAL_E = 
    Global::GEAR->getEcalEndcapParameters();
  shape     = POLYGON;
  symmetry  = pECAL_E.getSymmetryOrder();
  phi0      = pECAL_E.getPhi0();
  r_inner   = pECAL_E.getExtent()[0]; 
  r_outer   = pECAL_E.getExtent()[1]; 
  abs_mat   = 74;  // Tungsten
  detector  = 0.5; // No chance to extract 
  det_mat   = 14;  // Silicon
}
//============================================================================
void PGdb::Zone::_init_ecal1_cap(int last_layer){
//============================================================================
  const gear::CalorimeterParameters& pECAL_E = 
    Global::GEAR->getEcalEndcapParameters();
  const gear::LayerLayout &le = pECAL_E.getLayerLayout() ;
  no        = ECAL1_CAP;
  z_inner   = le.getDistance(0);; 
  z_outer   = le.getDistance(last_layer); 
  n_sampl   = last_layer; 
  sampling  = le.getThickness(0); 
  min_lay   = 0;
  max_lay   = last_layer-1; 
  absorber  = le.getAbsorberThickness(0); 
  cell_size = le.getCellSize0(0); // should be 10 mm
//---------------------------------
  mip_vis   = 170.0e-6 ;  // Should be in GEAR
  e_coeff   = 33.02346 ;  // for model LDC00
  //    e_coeff   = 46.7703 ;  // for model LDC01  Should be in parameters
  //-----------------added KP
  Zeff      =    67.41;
  Aeff      =  166.868;
  Rhoeff    =     7.75;
  Ieff      = 706.67e-9;
  Rmeff     = 22.0833;// mm
  x0eff     = 9.4587 ;  // mm
  Eceff     = 9.0804  ; // MeV
  eprime    = 0.72002;

  _init_ecal_cap_common();
}
//============================================================================
void PGdb::Zone::_init_ecal2_cap(int last_ecal1_layer){
//============================================================================
  const gear::CalorimeterParameters& pECAL_E = 
    Global::GEAR->getEcalEndcapParameters();
  const gear::LayerLayout &le = pECAL_E.getLayerLayout() ;
  int nLayere = le.getNLayers() ;
  no        = ECAL2_CAP;
  z_inner   = le.getDistance(last_ecal1_layer); 
  z_outer   = le.getDistance(nLayere-1)+le.getThickness(nLayere-1);; 
  n_sampl   = nLayere-last_ecal1_layer; 
  sampling  = le.getThickness(nLayere-1); 
  min_lay   = last_ecal1_layer;
  max_lay   = nLayere-1; 
  absorber  = le.getAbsorberThickness(nLayere-1); 
  cell_size = le.getCellSize0(nLayere-1); // should be 10 mm
//---------------------------------
  mip_vis   = 170.0e-6 ;  // Should be in GEAR
  e_coeff   = 93.56822 ;  //for model LDC00
//    e_coeff   =  86.3749 ;  //for model LDC01   
 //-----------------added KP
  Zeff      =   71.635;
  Aeff      =  177.614;
  Rhoeff    =   12.577;
  Ieff      =729.09e-9;
  Rmeff     = 13.8418;// mm
  x0eff     = 8.135953 ;  // mm
  Eceff     = 8.4709  ; // MeV
  eprime    = 0.709672;

  _init_ecal_cap_common();
}
//============================================================================
void PGdb::Zone::_init_hcal(PGdb::ZONE zone){
//============================================================================
  const gear::CalorimeterParameters *pHCAL;
  if(zone==HCAL_BAR)
   pHCAL = & Global::GEAR->getHcalBarrelParameters();
  else
   pHCAL = & Global::GEAR->getHcalEndcapParameters();
  const gear::LayerLayout &lhb = pHCAL->getLayerLayout() ;
  no        = zone;
  shape     = POLYGON;
  symmetry  = pHCAL->getSymmetryOrder();
  phi0      = pHCAL->getPhi0();
  r_inner   = pHCAL->getExtent()[0]; 
  r_outer   = pHCAL->getExtent()[1]; 
  z_inner   = pHCAL->getExtent()[2]; 
  z_outer   = pHCAL->getExtent()[3]; 
  n_sampl   = lhb.getNLayers() ;
  sampling  = lhb.getThickness(0); 
  min_lay   = 0;
  max_lay   = n_sampl-1; 
  absorber  = lhb.getAbsorberThickness(0); 
  abs_mat   = 26;  // atomic number = Iron
  detector  = 5.0; // Thickness [mm] -- No chance to extract 
  det_mat   = 6;   //  atomic number = Carbon
  cell_size = lhb.getCellSize0(0); // should be 30 mm
  //---------------------------------
  mip_vis   = 875.0e-6 ;   // Should be in GEAR
  e_coeff   = 21.196262;   // for model LDC00
  //    e_coeff   = 22.01925 ;  //for model LDC01
 //-----------------added KP
  Zeff      =   25.318;
  Aeff      =   54.354;
  Rhoeff    =    6.134;
  Ieff      = 283.34e-9;
}
//============================================================================
void PGdb::Zone::_init_all_common(){
//============================================================================
//   Absolutly voluntary -- never used
  n_sampl  = 1; // Number of layers in calorimeter or TPC rings
  sampling = 1; // Sampling layer thickness
  min_lay  = 1; // Minimal Layer number
  max_lay  = 1; // Maximal Layer number
  absorber = 1; // Absorber layer thickness
  abs_mat  = 1; // Absorber material
  detector = 1; // Detector layer thickness
  det_mat  = 1; // Detector material
  cell_size= 1; // Cell size in X and Y (square shape)
  mip_vis  = 1;
  e_coeff  = 1;
}
//============================================================================
void PGdb::Zone::_init_vtx(){
//============================================================================
  no        = VTX;
  shape     = CYLINDER;
  symmetry  = 0;
  phi0      = 0.;
  r_inner   = 0.; 
  r_outer   = PGDB[ECAL1_CAP].r_inner; 
  z_inner   = 0; 
  z_outer   = PGDB[HCAL_CAP].z_outer; 
  _init_all_common();
}
//============================================================================
void PGdb::Zone::_init_coil(){
//============================================================================
  no        = COIL;
  shape     = CYLINDER;
  symmetry  = 0;
  phi0      = 0.;
  r_inner   = PGDB[HCAL_BAR].r_outer/cos(M_PI/PGDB[HCAL_BAR].symmetry); 
  r_outer   = r_inner+850;  // Should be in GEAR
  z_inner   = 0;
  z_outer   = PGDB[HCAL_CAP].z_outer; 
 //-----------------added KP
  Zeff      =  13.0;
  Aeff      =  27.0;
  Rhoeff    =   2.7;
  Ieff      = 166.0e-9;

  _init_all_common();
}
//============================================================================
void PGdb::Zone::_init_detector(){
//============================================================================
  no        = DETECTOR;
  shape     = CYLINDER;
  symmetry  = 0;
  phi0      = 0.;
  r_inner   = 0.; 
  r_outer   = PGDB[COIL].r_outer; 
  z_inner   = 0; 
  z_outer   = PGDB[HCAL_CAP].z_outer; 
  x0eff     = 100000.0 ;  // mm
  _init_all_common();
}
//============================================================================
void PGdb::Zone::_init_world(){
//============================================================================
  no        = WORLD;
  shape     = CYLINDER;
  symmetry  = 0;
  phi0      = 0.;
  r_inner   = PGDB[COIL].r_outer; 
  r_outer   = 1.e10; 
  z_inner   = PGDB[HCAL_CAP].z_outer; 
  z_outer   = 1.e10; 
  x0eff     = 100000.0 ;  // mm
  _init_all_common();
}
//============================================================================
void PGdb::Zone::_init_final() {
//============================================================================
  //      Predicted distance to neighbor for particular zone
  //               1.41 is ~ sqrt(2)
  r_neighbor = 1.4*((cell_size < sampling) ? sampling : cell_size/2.);
  r_neighbor = r_neighbor * r_neighbor;

  //          Volume including absorber
  cell_vol = cell_size*cell_size*sampling;  

//  Predicted cutoffs  for particular zone around MIP
//             Should be in GEAR
  //   0.5 -- no MIP signal less than this energy
  cut_noise = 0.3*mip_vis;   // in units of RAW data 

  //   MIP peak is within 0.5 MIP and 2.0 MIP
  //    2.0 -- choosen by eyes see page 267 old log book
  cut_mip  = 2.0*mip_vis;   // in units of RAW data 

  //      3.7 -- choosen by eyes see page 267 in old log book
  cut_hadr = 3.7*mip_vis;   // in units of RAW data 

  //        Just only green and red hits -- no blue hits at all
  //      cut_hadr = 2.0*mip_vis;   // in units of RAW data 

  // MIP most probable  physical energy, i.e. including energy lost in absorber
  //       Predicted MIP amplitude and density
  mip_whole = mip_vis*e_coeff;    // in [GeV]
  // MIP energy density in whole cell volume
  mip_dens = mip_whole/cell_vol; // in [GeV]/[mm]^3

  if(shape==POLYGON && symmetry){
    if(symmetry == 2)  
// because GEAR think that endcap has a symmetry = 2 instead of 8
      symmetry = 8;
    a0 = 2*M_PI/symmetry;
    r_max=r_outer/cos(a0/2.);
    r_min=r_inner/cos(a0/2.);
  } else {
    a0 = 2.*M_PI/8;  // default symmetry
    r_max=r_outer;
    r_min=r_inner;
  }
}
//============================================================================
void PGdb::init(){
//============================================================================
  B_Field  = Global::GEAR->getTPCParameters().getDoubleVal("BField");
  
  // Sequence is important !!!
  zone[TPC]._init_tpc();
  zone[ENDPLATE1]._init_endpl1();
  zone[ENDPLATE2]._init_endpl2();

  int last_layer=_ecal_last_layer(Global::GEAR->getEcalBarrelParameters());
  zone[ECAL1_BAR]._init_ecal1_bar(last_layer);
  zone[ECAL2_BAR]._init_ecal2_bar(last_layer);

  last_layer=_ecal_last_layer(Global::GEAR->getEcalEndcapParameters());
  zone[ECAL1_CAP]._init_ecal1_cap(last_layer);
  zone[ECAL2_CAP]._init_ecal2_cap(last_layer);

  zone[HCAL_BAR]._init_hcal(HCAL_BAR);
  zone[HCAL_CAP]._init_hcal(HCAL_CAP);

  zone[VTX]._init_vtx();
  zone[COIL]._init_coil();
  zone[DETECTOR]._init_detector();
  zone[WORLD]._init_world();

  for( int k=0 ; k < ZONE_COUNT ; k++ )
    zone[k]._init_final();

  cout << "    =================================================================" 
       << endl ;
  cout << "          Phys_Geom_Database created and initialized  "
       <<" with "<< ZONE_COUNT <<" zones "<< endl;
  cout << "    =================================================================" 
       << endl ;
}   // End PGdb_init() 
//============================================================================
//      Simple Detector Geometry Zone Finder
//============================================================================
bool PGdb::Zone::in_polygon(double r,Point3D &p){
//============================================================================
  double ph = atan2(p.y,p.x) + phi0;
  if (ph < 0.0) ph = 2.*M_PI + ph;
  return (p.rz > r/cos(ph-(trunc((ph+a0/2.)/a0)*a0) ))?false:true;
}
//============================================================================
bool PGdb::Zone::inside(Point3D &p){
//============================================================================
  if(p.z < z_inner || p.z > z_outer)
    return false;
  if(p.rz > r_max)
    return false;
  if(p.rz > r_outer)
    if(!in_polygon(r_outer,p))
      return false;
  if(p.rz > r_min)
    return true;
  if(p.rz > r_inner)
    return !in_polygon(r_inner,p);
  return false;
}
//============================================================================
PGdb::ZONE PGdb::get_zone(Point3D &p){
//============================================================================
// As well as number of zones is small we can go along all of them
//           can be optimized  -- latter
//  DETECTOR contains of all of them; WORLD is outside region of DETECTOR
//============================================================================
  Point3D ap(abs(p.x),abs(p.y),abs(p.z));
  if( p.x==0.0 && p.y==0.0 && p.z==0.0)// ???????
    return VTX;
  if(!zone[DETECTOR].inside(ap))
    return WORLD;
  for(unsigned i=0;i<DETECTOR;i++)
    if(zone[i].inside(ap))
      return (ZONE)i;
  return DETECTOR;
}
//============================================================================
const char *PGdb::Zone::get_name() const {
//============================================================================
// These names should be the same as in ZONE typedef
  static const char *names[]={ 
    "VTX", "TPC" , "ECAL1_BAR", "ECAL2_BAR", "ECAL1_CAP", "ECAL2_CAP",
    "HCAL_BAR", "HCAL_CAP", "COIL", "DETECTOR", "WORLD"
  };
  if(no>=sizeof(names)/sizeof(names[0]))
    return "Unknown";
  return names[no];
}
//============================================================================
ostream &operator<<(ostream &o,const PGdb &d){
//============================================================================
  o<<"==================================================================="<< endl;
  o<<" ================ Physical Geometrical Database  ================="<< endl;  
  for (unsigned i = 0; i < PGdb::ZONE_COUNT; i++){
    const PGdb::Zone &z=d[static_cast<PGdb::ZONE>(i)];
    o<<"==================================================================="<< endl;
    o<<" \x1b[30;46m                 Zone   "<< i << "  " << z.get_name() 
     <<"                        \x1b[0m\x1b[30m"<< endl;
    o<<"-----------------------------------------------------------------------"<< endl;
    o<<" R inner    R outer    Z inner   Z outer    Shape    N samples"<< endl;
    o<<"   "<<z.r_inner
     <<"       "<<z.r_outer
     <<"       "<<z.z_inner
     <<"       "<<z.z_outer
     <<"       "<<z.shape
     <<"       "<<z.n_sampl<<endl;
    o<<"-----------------------------------------------------------------------"<< endl;
    o<<" Min--Max layers,    Thicknesses,         Cell size,        Volume,       R-neib  "<< endl;
    o<<"    "<< z.min_lay 
     <<"    "<< z.max_lay
     <<"       "<< z.sampling         <<" mm, "
     <<"  "<< z.detector         <<" mm, "
     <<"        "<< z.cell_size        <<" mm, "
     <<"      "<< z.cell_vol         <<" mm^3, "
     <<"      "<< sqrt(z.r_neighbor)  <<" mm"        <<endl; 
    if(i==PGdb::HCAL_CAP  ||
       i==PGdb::HCAL_BAR  ||  
       i==PGdb::ECAL1_BAR ||  
       i==PGdb::ECAL2_BAR ||  
       i==PGdb::ECAL1_CAP ||  
       i==PGdb::ECAL2_CAP ){
      o<<"-----------------------------------------------------------------------"<< endl;
      o<<"      MIP visible       Noise cut,          MIP cut,          HADR cut "<< endl
       <<"      "<< z.mip_vis*1.e6   <<" [keV], "
       <<"      "<< z.cut_noise*1.e6 <<" [keV], "
       <<"      "<< z.cut_mip*1.e6   <<" [keV], "
       <<"      "<< z.cut_hadr*1.e6  <<" [keV]  " << endl;
      o<<"-----------------------------------------------------------------------"<< endl;
      o<<"     E_coeff;                      MIP whole [GeV];               MIP density "<< endl
       <<"   "<< z.e_coeff   <<" Phys./Vis.[GeV],"  
       <<"      "<< z.mip_whole <<" [GeV], "    
       <<"      "<< z.mip_dens  <<" [GeV/mm^3]" << endl <<endl;
    }
  }
  o<<"-----------------------------------------------------------------------"<< endl;
  return o;
} // ------- end Output stream

// +++++++++++ End Phys_Geom_Database definition  +++++++++++++++++++++++

