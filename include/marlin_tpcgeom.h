
#ifdef __cplusplus
/**
 * header file to provide common tpc geometry to both fortran and cpp implementations
 * tracking code. In the cpp implementation this is via namespace
 */
#endif

#define NPADROWS 200
#define INNERRAD 38.6
#define OUTERRAD 162.6
#define MAXDRIFT 250.
#define IONPOTEN 0.00000003
#define TPCRPRES 0.016
#define TPCZRES  0.1
#define TPCPIXZ  0.14
#define TPCPIXRP 0.22

#ifdef _LANGUAGE_FORTRAN

      INTEGER LTPDRO, LTTROW, LTSROW, LTWIRE, LTSTYP, LTSLOT, LTCORN,
     +        LTSECT, LTTPAD, LMXPDR, LTTSRW
c(kh) PARAMETER (LTPDRO=118,LTTROW=19,LTSROW=12,LTWIRE=200,LTSTYP=3,
      PARAMETER (LTPDRO=NPADROWS,LTTROW=19,LTSROW=12,LTWIRE=200,LTSTYP=3,
     +           LTSLOT=12,LTCORN=6,LTSECT=LTSLOT*LTSTYP,LTTPAD=4,
     +           LMXPDR=150,LTTSRW=11)
      INTEGER NTSECT, NTPCRN, ITPTYP, ITPSEC, IENDTP
      REAL RTPCMN, DRTPMN, DRTPMX, DZTPMX, TPFRDZ,
     &     TPFRDW, TPAVDZ, TPFOF1, TPFOF2, TPFOF3, TPPROW, TPTROW,
     &     TPCORN, TPPHI0, TPCPH0, TPSPH0
      COMMON /TPGEOM/RTPCMN,DRTPMN,DRTPMX,DZTPMX,
     &               TPFRDZ,TPFRDW,TPAVDZ,TPFOF1,TPFOF2,TPFOF3,
     &               TPPROW(LTPDRO),TPTROW(LTTROW),NTSECT,
     &               NTPCRN(LTSTYP),TPCORN(2,LTCORN,LTSTYP),
     &               TPPHI0(LTSECT),TPCPH0(LTSECT),TPSPH0(LTSECT),
     &               ITPTYP(LTSECT),ITPSEC(LTSECT),IENDTP(LTSECT)
#endif

#ifdef __cplusplus

#ifndef tpcgeom_h
#define tpcgeom_h
namespace tpcgeom
{

#ifdef compile_tpcgeom
  //dimentions should be in mm
  /**number of tpc pad rows */
  extern const int nrtpc = NPADROWS;
  // don't understand why setting these const causes linking error   
  //outer radius of active tpc volume in mm
  extern const float tpcacro = OUTERRAD*10.;
  //inner radius of active tpc volume  in mm
  extern const float tpcacri = INNERRAD*10.;
  //radial pad size in mm
  float tpcpadr = (tpcacro-tpcacri)/(float)(nrtpc);
  //maximum drift length in mm
  float zdrift = MAXDRIFT*10.;
  //radius of given the nth row
  float rrow[nrtpc];
  //rphi parameterised z dependand resolution fixme need function here
  //instead of variable in mm
  float tpc_rphi_res_max = TPCRPRES*10.;
  //z resolution in mm
  float tpc_z_res = TPCZRES*10.;
  //read out pixel size in rphi in mm
  float pix_rp = TPCPIXRP*10.;
  //read out time slice in mm
  float pix_z = TPCPIXZ*10.;
  //number of pixels in phi per row
  int nphi_pix[nrtpc];
  //number of time slices per row
  int nz_pix = (int)(2.*zdrift/pix_z);
  int npix_tpc; 
  // Excitation energy for a single electron (GeV) (10*Z eV; Z=18 for Argon) 
  float ionisation_potential = IONPOTEN;
#endif //compile_tpcgeom

  /*----------------------------------------------------------- */

  extern   const int nrtpc;
  extern   const float tpcacro;
  extern   const float tpcacri ;
  extern   float tpcpadr ;
  extern   float zdrift;
  extern   float rrow[];
  extern   float tpc_rphi_res_max;
  extern   float tpc_z_res;
  extern   float pix_rp;
  extern   float pix_z;
  extern   int nphi_pix[];
  extern   int nz_pix ;
  extern   int npix_tpc; 
  extern   float ionisation_potential;
}


#endif
#endif
