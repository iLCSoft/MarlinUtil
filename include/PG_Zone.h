#ifndef PG_ZONE_H
#define PG_ZONE_H 1

class PG_Zone {

/**  
 *    Physical Geometrical database will be used in the reconstruction
 *    procedure to get fast access to values/varibales dependent on
 *    the particular and different calorimeter zones with different
 *    physical properties and geometrical relations between cells. 
 *
 *    Such a database can be created after reading the geometry 
 *    record in LCIO if it will be presented in or another kind of 
 *           outstanding geometrical database.
 *
 *    Another way is to create such database during simulation
 *      step as an abligatory for any simulation program, 
 *             during the detector geometry creation.
 *       That is much easy way to form such a database.
 *            and put it as member (or class) of LCIO.
 */

public:

    int     min_lay;    // Minimal Layer number
    int     max_lay;    // Minimal Layer number

    double  th_sampl;   // Sampling layer thickness
    double  th_det;     // Detector layer thickness
    double  cell_size;  // Cell size in X and Y (square)
// MIP most probable  detector energy, does not include energy lost in absorber
    double  mip_vis;
//      Coeff  converts visible enery to physical energy 
// i.e. converts energy lost in detector into energy lost in sampling layer  
//                   including absorber 
    double  e_coeff;

    double  r_neibour; // Predicted distance to neibour for particular zone

    double  cell_vol;     // Volume including absorber
//       Predicted cutoffs  for particular zone around MIP
    double  cut_noise;
    double  cut_mip;
    double  cut_hadr;
// MIP most probable  physical energy, i.e. including energy lost in absorber
    double  mip_whole;

    double  mip_dens; // MIP energy density in whole cell volume

};
#endif
