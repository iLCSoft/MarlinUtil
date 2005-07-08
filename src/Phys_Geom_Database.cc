#include <iostream>
#include <assert.h>
#include <typeinfo>
#include <math.h>

#include "Phys_Geom_Database.h"

using namespace std;

/**  
 *    Physical Geometrical database will be used in the reconstruction
 *    procedure to get fast access to values/varibales dependent on
 *    the particular and different calorimeter zones with different
 *    physical properties and geometrical relations between cells. 
 *
 *    Such a database can be created after readind the geometry 
 *    record in LCIO if it will be presented in or another kind of 
 *           outstanding geometrical database.
 *
 *    Another way is to create such database during simulation
 *      step as an abligatory for any simulation program, 
 *             during the detector geometry creation.
 *       That is much easy way to form such a database.
 *            and put it as member (or class) of LCIO.
 */

/*
  ===================== HBARREL geometry ============
  Inner Calorimeter Barrel Radius=  170.cm
  Outer Calorimeter Barrel Radius=  300.cm
  ECAL thickness                    20.7999992cm
  Inner HCAL barrel Radius=         190.800003cm
  Outer HCAL barrel Radius=         296.799988cm
  Z_size of  HCAL barrel=           265.850006cm
  HCAL barrel wall thickness=       0.200000003cm
  HCAL barrel absorber thickness=   2.cm
  HCAL barrel scint thickness=      0.5cm
  HCAL barrel gap thickness=        0.150000006cm
  HCAL bar sampling =               2.6500001cm


----------------------------------------------------------------------
c TPC Dimensions:
C Total inner thickness of TPC material is 3%X_0 and 1X_0 of Al is 8.9 cm
c THickness of material in the Inner wall (Barrel) is 1%X_0
      TPCTHBI=0.01*8.9
c THickness of material in the Outer wall (Barrel) is 3%X_0
      TPCTHBO=0.03*8.9
c Endplate thickness - 30cm.
      TPCTHKE=30.
*     thinner endplates to accomodate the FCH and to not
*     interfere with the ECAL
      TPCTHKE=23.
C Inner Radius
      RTPCINN=32.
C Outer Radius
c MVL      RTPCOUT=170.
      RTPCOUT=169.
C Inner Active radius
      TPCACRI=38.6
C Outer Active radius
      TPCACRO=162.6
c  z-half-length of TPC:
      TPCHLFZ=280.
      TPCHLFZ=273.
C Number of radial readout rings:
c(kh) NRTPC=118
      NRTPC=200
C Maximum drift length:
      ZDRIFT=TPCHLFZ-TPCTHKE
c Radial pad size
      TPCPADR=(TPCACRO-TPCACRI)/REAL(NRTPC)

*/
// --------------- Constructor ----------------------
//-----------------------------------------------------------------------
//     Everything are filled in Constructor now, once per run
//-----------------------------------------------------------------------
Phys_Geom_Database::Phys_Geom_Database(int nn) 
{
    ndb  = nn;             // number of zones in database
    z = new PG_Zone [ndb]; // = malloc(start_size*sizeof(PG_Zone ))

//---------------------------------
//    Initialize ECAL zones
//                 Barrel
//---------------------------------
    z[0].min_lay   =  1       ;  
    z[0].max_lay   = 30       ;
    z[0].th_sampl  =   4.35   ; 
    z[0].th_det    =   0.5    ;
    z[0].cell_size =  10.     ; 
    z[0].mip_vis   = 170.0e-6 ; 
    z[0].e_coeff   =  28.369  ;

    z[1].min_lay   = 31       ;
    z[1].max_lay   = 40       ; 
    z[1].th_sampl  =   7.15   ; 
    z[1].th_det    =   0.5    ; 
    z[1].cell_size =  10.     ; 
    z[1].mip_vis   = 170.0e-6 ; 
    z[1].e_coeff   =  80.378  ;
//---------------------------------
//                  EndCaps
//---------------------------------
    z[2].min_lay   =  1       ;
    z[2].max_lay   = 30       ; 
    z[2].th_sampl  =   4.35   ; 
    z[2].th_det    =   0.5    ;
    z[2].cell_size =  10.     ; 
    z[2].mip_vis   = 170.0e-6 ; 
    z[2].e_coeff   = 28.369   ;

    z[3].min_lay   = 31       ;  
    z[3].max_lay   = 40       ;  
    z[3].th_sampl  =   7.15   ; 
    z[3].th_det    =   0.5    ;
    z[3].cell_size =  10.     ; 
    z[3].mip_vis   = 170.0e-6 ; 
    z[3].e_coeff   =  80.378  ;
    
//---------------------------------
//    Initialize HCAL zones
//                 Barrel
//---------------------------------
    z[4].min_lay   =  1       ; 
    z[4].max_lay   = 38       ; 
    z[4].th_sampl  = 26.5     ; 
    z[4].th_det    =  5.0     ;
    z[4].cell_size =  10.     ; 
    z[4].mip_vis   = 875.0e-6 ; 
    z[4].e_coeff   =  35.46   ;
//---------------------------------
//                  EndCaps
//---------------------------------
    z[5].min_lay   =  1       ; 
    z[5].max_lay   = 56       ; 
    z[5].th_sampl  = 26.5     ; 
    z[5].th_det    =  5.0     ; 
    z[5].cell_size =  10.     ; 
    z[5].mip_vis   = 875.0e-6 ; 
    z[5].e_coeff   =  35.46   ;
//---------------------------------
//                  TPC
//---------------------------------
    z[6].min_lay   =  1      ; 
    z[6].max_lay   = 200     ; 
    z[6].th_sampl  = 0.5     ; 
    z[6].th_det    = 0.5     ; 
    z[6].cell_size =  6.     ; 
    z[6].mip_vis   =  1.0e-4 ; 
    z[6].e_coeff   =  1.0    ;


//         Initialize calculatable values of database
    for( int k=0 ; k < ndb ; k++ ){

//      Predicted distance to neibour for particular zone
//               1.41 is ~ sqrt(2)
      z[k].r_neibour = 2.71*((z[k].cell_size < z[k].th_sampl) ? 
			     z[k].th_sampl : z[k].cell_size);
      z[k].r_neibour = z[k].r_neibour * z[k].r_neibour;
// Volume including absorber
      z[k].cell_vol = z[k].cell_size*z[k].cell_size*z[k].th_sampl;  

//       Predicted cutoffs  for particular zone around MIP
// 0.6 -- no MIP signal less than this energy
      z[k].cut_noise = 0.5*z[k].mip_vis;   // in units of RAW data 

//       MIP peak is within 0.6 MIP and 2.0 MIP
// 2.0 -- choosen by eyes see page 267 old log book
      z[k].cut_mip  = 2.0*z[k].mip_vis;   // in units of RAW data 

// 4.7 -- choosen by eyes see page 267 in old log book
      z[k].cut_hadr = 4.7*z[k].mip_vis;   // in units of RAW data 

// MIP most probable  physical energy, i.e. including energy lost in absorber
//       Predicted MIP amplitude and density
      z[k].mip_whole = z[k].mip_vis*z[k].e_coeff;    // in [GeV]
// MIP energy density in whole cell volume
      z[k].mip_dens = z[k].mip_whole/z[k].cell_vol; // in [GeV]/[mm]^3
    }

//        ECAL and HCAL half of size in Z direction
//        to distinguish between barrel and endcap
    ecal_z_size = 2700.0;
    hcal_z_size = 2658.5;
    tpc_z_size  = 2500.0;
    tpc_inn_rad =  386.0;
    tpc_out_rad = 1626.0;
/**
    cout << "    =================================================================" 
	 << endl ;
    cout << "          Phys_Geom_Database created and initialized  "
	 <<" with "<< ndb <<" zones "<< endl;
    cout << "    =================================================================" 
	 << endl ;
 */
} 

// --------------- Destructor ----------------------------
Phys_Geom_Database::~Phys_Geom_Database()   
{
/**
    cout << "    =================================================================" 
	 << endl ;
    cout << "         Destructor called for "
	 << typeid(*this).name()
	 <<" with "<< ndb <<" zones "<< endl;
    cout << "    =================================================================" 
	 << endl ;
 */
    delete [] z;
}


   void Phys_Geom_Database::print()
{
    cout<<"==================================================================="<< endl;
    cout<<" ================ Physical Geometrical Dtatabase  ================="<< endl;
   for (int i = 0; i < ndb; i++){
    cout<<"==================================================================="<< endl;
       cout<<"                          Zone   "<< i<< endl;
       cout<<"-----------------------------------------------------------------------"<< endl;

 cout<<" Min--Max layers,    Thicknesses,           Cell size,        Volume,       R-neib  "<< endl;

       cout<<"   "<< z[i].min_lay 
	   <<'\t'<< z[i].max_lay
	   <<'\t'<< z[i].th_sampl         <<" mm, "
	   <<'\t'<< z[i].th_det           <<" mm, "
	   <<'\t'<< z[i].cell_size        <<" mm, "
	   <<'\t'<< z[i].cell_vol         <<" mm^3, "
           <<'\t'<< sqrt(z[i].r_neibour)  <<" mm"        <<endl; 
       cout<<"-----------------------------------------------------------------------"<< endl;

       cout<<"      MIP visible       Noise cut,          MIP cut,      HADR cut "<< endl
	   <<"   "<< z[i].mip_vis*1.e6  <<" [keV], "
	   <<'\t'<<'\t'<< z[i].cut_noise*1.e6 <<" [keV], "
	   <<'\t'<< z[i].cut_mip*1.e6   <<" [keV], "
	   <<'\t'<< z[i].cut_hadr*1.e6  <<" [keV]  " << endl;
       cout<<"-----------------------------------------------------------------------"<< endl;

       cout<<"     E_coeff;                   MIP whole [GeV];                 MIP density "<< endl
	   <<"   "<< z[i].e_coeff  <<" Phys./Vis. [GeV], "  
	   <<'\t'<< z[i].mip_whole <<" [GeV], "    
	   <<'\t'<< z[i].mip_dens  <<" [GeV/mm^3]" << endl <<endl;
   }
       cout<<"-----------------------------------------------------------------------"<< endl;
} // ------- end Printing


// +++++++++++ End Class Phys_Geom_Database definition  +++++++++++++++++++++++

