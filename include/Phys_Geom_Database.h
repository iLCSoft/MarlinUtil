#ifndef PHYS_GEOM_DATABASE_H
#define PHYS_GEOM_DATABASE_H 1

#include <iostream>
#include "marlin/Processor.h"
#include "lcio.h"
#include <B_Util.h>
#include <string>
/**
 *                                                                                <br>
 *       It might be forms as gometry driver in Mokka.                            <br>
 *       Even more - maybe better create it at the simulation                     <br>
 *       step and put it onto LCIO simulation outut file.                         <br>
 *                                                                                <br>
 *   Phys_Geom_Database class make connection between detector geometry and       <br>
 *             Particle-Flow reconstruction package.                              <br>
 *                                                                                <br>
 *   *************************************************************                <br>
 *       It strongly depends on the PARTICULAR detector geometry.                 <br>
 *         "The Lastest" and unique link to geometry before                       <br>
 *                  the PF reconstruction procedure                               <br>
 *   ****************************************************************             <br>
 *                                                                                <br>
 *    NOW it made BY HAND for TESLA detector calorimeter                          <br>
 *                                                                                <br>
 *   It also may conists of all reconstruction tuning parameters.                 <br>
 *     Some of them depends on zone nomber, some does not.                        <br>
 *                                                                                <br>
 *   Each hit of the Rec_Hits class has a zone number                             <br>
 *   that is the reference to this database, which give us                        <br>
 *   physical and geometrical properties of particular zone.                      <br>
 *                                                                                <br>
 *                                                                                <br>
 *                                                                                <br>
 *   PGdb::ZONE PGdb::get_zone(Point3D &p)                                        <br>
 *     Function return zone "number" (enum type)                                  <br>
 *      for any arbitrary 3-d point                                               <br>
 *                                                                                <br>
 * Usage:                                                                         <br>
 *                                                                                <br>
 *                                                                                <br>
 * Hit h;                                                                         <br>
 *                                                                                <br>
 *   ZONE i=PGDB.get_zone(h.cp);                                                  <br>
 *   ZONE i;                                                                      <br>
 *   unsigned j;                                                                  <br>
 *   ...                                                                          <br>
 *   int iel = pgdb[i].det_mat;                                                   <br>
 *   BooElemens[pgdb[i].det_mat];                                                 <br>
 *   double ecal_inn_rad;                                                         <br>
 *   pgdb[PGDB::ECAL1_BAR].r_inner;                                               <br>
 *                                                                                <br>
 *                                                                                <br>
 *    @author V. L. Morgunov, A. Zhelezov  (DESY/ITEP)                            <br>
 *                                                                                <br>
 *    modifications and extensions by P. Krstonosic (DESY)                        <br>
 */

//============================================================================
class PGdb {
//============================================================================
 public:
  typedef enum {
    VTX=0,
    TPC,
    ECAL1_BAR,
    ECAL2_BAR,
    ECAL1_CAP,
    ECAL2_CAP,
    HCAL_BAR,
    HCAL_CAP,
    COIL,
    ENDPLATE1,
    ENDPLATE2,
    DETECTOR,
    WORLD,
    ZONE_COUNT
  } ZONE; //       !! Change get_name() in cc file in case of changes here !!
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
    double a0;
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
    double          Rmeff ;// mm
    double          x0eff; // mm
    double          Eceff; // MeV 
    double          eprime; //


    //---------------------------------
    //     Physical Parameters
    //---------------------------------
    // MIP most probable detector energy, does not include energy lost in absorber
    double  mip_vis;
    //      Coeff  converts visible enery to physical energy 
    // i.e. converts energy lost in detector into energy lost in sampling layer  
    //                   including absorber 
    double  e_coeff;
    
    double  r_neibour; // Predicted distance to neibour for particular zone
    
    double  cell_vol;     // Volume including absorber
    //    RGB   Predicted cutoffs  for particular zone around MIP
    double  cut_noise;
    double  cut_mip;
    double  cut_hadr;
    // MIP most probable  physical energy, i.e. including energy lost in absorber
    double  mip_whole;
    
    double  mip_dens; // MIP energy density in whole cell volume

    const char* get_name() const;

  private:
    friend class PGdb;
    void _init_final();
    void _init_tpc();
    void _init_ecal_bar_common();
    void _init_ecal1_bar(int last_layer);
    void _init_ecal2_bar(int last_ecal1_layer);
    void _init_ecal_cap_common();
    void _init_ecal1_cap(int last_layer);
    void _init_ecal2_cap(int last_ecal1_layer);
    void _init_hcal(PGdb::ZONE zone);
    void _init_endpl1();
    void _init_endpl2();
    void _init_all_common();
    void _init_vtx();
    void _init_coil();
    void _init_detector();    
    void _init_world();

    bool inside(Point3D &p);
    bool in_polygon(double r,Point3D &p);
  
  };

 private:
  Zone zone[ZONE_COUNT];

 public:

  //  Zone zone[ZONE_COUNT];
  double B_Field;
  void init();
  ZONE get_zone(Point3D &p);
  ZONE get_zone2(int i);

  //  ZONE get_zone(Hit &h); // Ht is not separate class ????? will be less calculations
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
