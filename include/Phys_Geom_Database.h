#ifndef PHYS_GEOM_DATABASE_H
#define PHYS_GEOM_DATABASE_H 1

#include "PG_Zone.h"

class Phys_Geom_Database{
/**
 *       It might be forms as gometry driver in Mokka.
 *       Even more - maybe better create it at the simulation 
 *       step and put it onto LCIO simulation outut file. 
 *
 *   Phys_Geom_Database class make connection between detector geometry and
 *             Particle-Flow reconstruction package.
 *
 *   *************************************************************
 *       It strongly depends on the PARTICULAR detector geometry.
 *         "The Lastest" and unique link to geometry before 
 *                  the PF reconstruction procedure
 *   ****************************************************************
 *
 *    NOW it made BY HAND for TESLA detector calorimeter
 *
 *   It also may conists of all reconstruction tuning parameters.
 *     Some of them depends on zone nomber, some does not.
 *
 *   Each hit of the Rec_Hits class has a zone number
 *   that is the reference to this database, which give us
 *   physical and geometrical properties of particular zone.
 *    
 */
  
 public:
    Phys_Geom_Database(int);           // Constructor
   ~Phys_Geom_Database();              // Destructor

inline    PG_Zone& operator[](int idx) { return z[idx]; }
inline    PG_Zone& operator()(int idx) { return z[idx]; }
inline    double p_ecal_z_size(){return ecal_z_size;}
inline    double p_hcal_z_size(){return hcal_z_size;}
inline    double p_tpc_z_size() {return tpc_z_size ;}
inline    double p_tpc_inn_rad(){return tpc_inn_rad;}
inline    double p_tpc_out_rad(){return tpc_out_rad;}

    void print();

 protected:

    PG_Zone *z;
    int ndb;             // Number of zones in the detector
    double ecal_z_size;
    double hcal_z_size;
    double tpc_z_size ;
    double tpc_inn_rad;
    double tpc_out_rad;

} ;    // End description of class Phys_Geom_Database

#endif
