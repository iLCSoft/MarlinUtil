#ifndef PHYS_GEOM_DATABASE_H
#define PHYS_GEOM_DATABASE_H 1

#include <iostream>
#include "marlin/Processor.h"
#include "lcio.h"
#include <B_Util.h>

/**
 *                                                                                <br>
 *   Phys_Geom_Database class make connection between detector geometry and       <br>
 *             Particle-Flow reconstruction package.                              <br>
 *                                                                                <br>
 *       It might be forms as geometry driver in Mokka.                           <br>
 *       Even more - maybe better create it at the simulation                     <br>
 *       step and put it onto LCIO simulation output file.                        <br>
 *                                                                                <br>
 *   *************************************************************                <br>
 *       It strongly depends on the PARTICULAR detector geometry.                 <br>
 *       GEAR is used to get that particular detector geometry                    <br>
 *         "The Latest" and unique link to geometry before                        <br>
 *                  the PFA reconstruction procedure                              <br>
 *   ****************************************************************             <br>
 *                                                                                <br>
 *   It also may consist of all reconstruction tuning parameters.                 <br>
 *     Some of them depends on zone number, some does not.                        <br>
 *                                                                                <br>
 * Usage:                                                                         <br>
 *                                                                                <br>
 * PGdb::ZONE PGdb::get_zone(Point3D &p);                                         <br>
 *     Function return zone "number" (enum type)                                  <br>
 *      for any arbitrary 3-d point                                               <br>
 *                                                                                <br>
 * if Hit h is a hit than h.cp is Caretsian coordinates of hit;                   <br>
 *                                                                                <br>
 *   ZONE i=PGDB.get_zone(h.cp);                                                  <br>
 *   ...                                                                          <br>
 *   int iel = pgdb[i].det_mat;                                                   <br>
 *   BooElemens[pgdb[i].det_mat];                                                 <br>
 *   double ecal_inn_rad;                                                         <br>
 *   pgdb[PGDB::ECAL1_BAR].r_inner;                                               <br>
 *                                                                                <br>
 *             Created 3 Oct 2004   V.L. Morgunov                                 <br>
 *             Last update 6 August 2006                                          <br>
 *             Last update 26 April 2006                                          <br>
 *                                                                                <br>
 *           Addons from Predrag Krsonos  26-April-2007                           <br>
 *                                                                                <br>
 *                                                                                <br>
 *    @author V. L. Morgunov, A. Zhelezov  (DESY/ITEP)                            <br>
 *                                                                                <br>
 */

//============================================================================
class PGdb {
//============================================================================
 public:
// !! Change get_name() in cc file in case of changes here !!
  typedef enum {
    VTX=0,
    TPC,
    ENDPLATE1, // TPC
    ENDPLATE2, // TPC
    ECAL1_BAR,
    ECAL2_BAR,
    ECAL1_CAP,
    ECAL2_CAP,
    HCAL_BAR,
    HCAL_CAP,
    COIL,
    DETECTOR,
    WORLD,
    ZONE_COUNT
  } ZONE; 
//============================================================================
  class Volume {
//============================================================================
  public:
    typedef enum {
      CYLINDER=0,
      POLYGON
    }SHAPE;
    SHAPE shape;
    unsigned symmetry;
    double   phi0;
    double   r_inner;
    double   r_outer;
    double   z_inner;
    double   z_outer;
  };
//============================================================================
  class Zone : public Volume {
//============================================================================
  private: 
    double r_min;
    double r_max;
    double a0;     // Polygone segment angle = 2Pi/symmetry for speed up
  public:
    unsigned           no; // Zone ID
    //---------------------------------
    //        Detector Parameters
    //---------------------------------
    unsigned      n_sampl; // Number of layers in calorimeter or TPC rings
    double       sampling; // Sampling layer thickness
    unsigned      min_lay; // Minimal Layer number
    unsigned      max_lay; // Maximal Layer number
    double       absorber; // Absorber layer thickness
    int           abs_mat; // Absorber material
    double       detector; // Detector layer thickness
    int           det_mat; // Detector material
    double      cell_size; // Cell size in X and Y (square shape)
    double           Zeff; // effective Z
    double           Aeff; // effective A
    double         Rhoeff; // effective Ro  g/cm^3
    double           Ieff; // effective exc. energy [GeV]
    double          Rmeff ;// Molere radius [mm]
    double          x0eff; // X0 [mm]
    double          Eceff; // Critical energy [MeV]
    double          eprime; // EM shower parameter
//---------------------------------
//     Physical Parameters
//---------------------------------
// MIP most probable detector energy, does not include energy lost in absorber
    double  mip_vis;
//      Coeff  converts visible energy to physical energy 
// i.e. converts energy lost in detector into energy lost in sampling layer  
//                   including absorber 
    double  e_coeff;
    
    double  r_neighbor; // Predicted distance to neighbor for particular zone
    
    double  cell_vol;     // Volume including absorber
//    RGB   Predicted cutoffs  for particular zone around MIP
    double  cut_noise;
    double  cut_mip;
    double  cut_hadr;
// MIP most probable  physical energy, i.e. including energy lost in absorber
    double  mip_whole;
    
    double  mip_dens; // MIP energy density in whole cell volume

    const char *get_name() const;

  private:
    friend class PGdb;
    void _init_final();
    void _init_tpc();
    void _init_endpl1();
    void _init_endpl2();
    void _init_ecal_bar_common();
    void _init_ecal1_bar(int last_layer);
    void _init_ecal2_bar(int last_ecal1_layer);
    void _init_ecal_cap_common();
    void _init_ecal1_cap(int last_layer);
    void _init_ecal2_cap(int last_ecal1_layer);
    void _init_hcal(PGdb::ZONE zone);
    void _init_all_common();
    void _init_vtx();
    void _init_coil();
    void _init_detector();    
    void _init_world();
//    Geometry 
    bool inside(Point3D &p);
    bool in_polygon(double r,Point3D &p);
  };

 private:
  Zone zone[ZONE_COUNT];

 public:
  double B_Field;
  void init();
  ZONE get_zone(Point3D &p);
//  ZONE get_zone(Hit &h); // Hit is not separate class ????? will be less calculations
  const Zone &operator[](ZONE z) const { return zone[z]; }

}; // +++++++++++ End PGdb --  Phys_Geom_Database definition  +++++++++++++++++++++++

extern PGdb PGDB;

std::ostream &operator<<(std::ostream &o,const PGdb &d);

using namespace lcio ;
using namespace marlin ;

//============================================================================
class PGDBP : public Processor {
//============================================================================
 public:
  virtual Processor*  newProcessor() { return new PGDBP ; }
  PGDBP();
  virtual void init() ;
};

#endif
