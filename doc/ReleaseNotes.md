# v01-17

* 2022-11-14 Thomas Madlener ([PR#35](https://github.com/iLCSoft/MarlinUtil/pull/35))
  - Remove the no longer supported gcc8 based CI workflow.

* 2022-10-20 Thomas Madlener ([PR#31](https://github.com/iLCSoft/MarlinUtil/pull/31))
  - Adding a basic Catch2 v3 based unittest setup that can be used to easily add unittests.
    - The CMake configuration allows to either build Catch2 on the fly or to discover a suitable installation of Catch2
  - Add the existing test and an example unittest to the test suit that can be run via `ctest` after building the package.
  - Disable building the tests for the coverity workflow for now to avoid polluting the output of that with Catch2 issues

* 2022-10-19 Thomas Madlener ([PR#33](https://github.com/iLCSoft/MarlinUtil/pull/33))
  - Make sure that `.ipp` files are also installed. Necessary since #30, fixes #32

* 2022-10-19 Thomas Madlener ([PR#30](https://github.com/iLCSoft/MarlinUtil/pull/30))
  - Introduce a `HelixClassT` template class and make `HelixClass` and `HellixClass_double` typedefs of this class, instead of having two separate (but practically identical) implementations that are in place currently.
    - Explicitly instantiate both versions that were in place previously to make sure things don't break downstream.
    - Mark getters as `const`

* 2022-09-19 Carl Mikael Berggren ([PR#28](https://github.com/iLCSoft/MarlinUtil/pull/28))
  The header-file has been re-organised and heavily
  commented - should work as a manual.
  
  Bug fix: a missing factor of two in propagateValErrors , the errors
  on the eigen-values (= the error on the estimated variances)
  
  New public methods:
  
  getElipsoid_r1 (_r2, _r3)
  getElipsoid_r1Error (_r2, _r3)
  getElipsoid_r_ave
  * getElipsoid_r_forw (_back)
  getElipsoid_vol
  getElipsoid_density
  getLongitudinalElipsis_eccentricity
  getTransverseElipsis_eccentricity
  * getMaxDist
  * getElipsoid_FractionInside
  
  Transformations:
  
  TransformPointToEigenSyst
  * TransformToEigenSyst
  * TransformAlongDirection (2 versions)
  
  (The ones with a * in front are only useful with Reco-input, since
   they need to have access to the CaloHits)
  getters of transformed properties:
  
  * get_x_trans (_y_, _z_)
  * get_COG_trans
  * get_COGCov_trans
  * get_th_ref
  * get_ph_ref
  * get_xyz_ref
  
  New private methods:
  
  findMaxDist
  findFirstAndLast
  findElipsoid_FractionInside
  
  Symbols defined: _one_sigma, _CL90, _CL95, _CL99

# v01-16-02

* 2022-06-28 Thomas Madlener ([PR#27](https://github.com/iLCSoft/MarlinUtil/pull/27))
  - Make the doxygen cmake configuartion work with newer versions of CMake (>= 3.17)

# v01-16-01

* 2021-09-29 Bohdan Dudar ([PR#23](https://github.com/iLCSoft/MarlinUtil/pull/23))
  - Add generic `MarlinUtil::getDetData` for getting DDRec detector extension data from dd4hep.
  - Add dedicated `getVXDData`, `getSITData`, `getFTDData`, `getTPCData`, `getSETData` in `MarlinUtil::ILD` namespace which return DDRec extensions of corresponding detector elements of the ILD detector.
  - Fix warnings about catching exceptions by value in `TrueJet_Parser`

* 2021-08-23 Andre Sailer ([PR#20](https://github.com/iLCSoft/MarlinUtil/pull/20))
  - CI: build against LCG_99python2 gcc8 and LCG_100 gcc10, clang11

# v01-16

* 2021-06-15 Thomas Madlener ([PR#18](https://github.com/iLCSoft/MarlinUtil/pull/18))
  - Move `TrueJet_Parser` utility class from MarlinReco to MarlinUtil.
    - Make it possible to use this in analysis code outside of MarlinReco

* 2021-06-09 Thomas Madlener ([PR#17](https://github.com/iLCSoft/MarlinUtil/pull/17))
  - Migrate the CI setup to use github actions instead of travis

* 2020-04-12 Frank Gaede ([PR#15](https://github.com/iLCSoft/MarlinUtil/pull/15))
  - make compatible w/ c++17 for macos/clang
        - patch provided by K.Fujii

# v01-15-01

* 2019-08-26 Andre Sailer ([PR#13](https://github.com/iLCSoft/MarlinUtil/pull/13))
  - MarlinUtilConfig: explicitly add DD4hep dependency so that dependent packages resolve the DD4hep::DDCore etc. libraries

* 2018-04-20 Erica Brondolin ([PR#11](https://github.com/iLCSoft/MarlinUtil/pull/11))
  - This PR is related to a detector issue discussed in: https://github.com/AIDASoft/DD4hep/issues/356
  Given that the return value of det.field().isValid() is always true for design, a different approach on the check of the geometry validity is applied.
  - This PR must be integrated with the following PR in AIDASoft/DD4hep: https://github.com/AIDASoft/DD4hep/pull/367

# v01-15

# v01-14

* 2017-05-03 Frank Gaede ([PR#3](https://github.com/iLCSoft/MarlinUtil/pull/3))
  - fix -Wabsolute-value in MarlinCED.cc

* 2017-05-03 Andre Sailer ([PR#2](https://github.com/iLCSoft/MarlinUtil/pull/2))
  - Fix all warnings in MarlinUtil for gcc and llvm, and cmake
      - use vectors instead of arrays to simplify memory handling in many classes, shadow warnings, unused parameters, unused member variables (LLVM), initialising members
  - Fix bug in copy constructor of SimpleHelix if the LCErrorMatrix argument was not NULL (https://github.com/iLCSoft/MarlinUtil/commit/4a65837e720e5ba001f3031622c1eeb6cc20be3e) the LCErrorMatrix member was never properly initialised
  - Separated ANN code into separate library to ignore all warnings, as this is exernal code (https://github.com/iLCSoft/MarlinUtil/commit/773a4d1dfd4ea7ad468d705a6c580e7a2890e165)

* 2017-05-05 Andre Sailer ([PR#4](https://github.com/iLCSoft/MarlinUtil/pull/4))
  - Add GeometryUtil file: implemented MarlinUtil::getBFieldInZ and MarlinUtil::getDetectorExtension
    - getBFieldInZ returns bfield at (0 0 0) in z direction in Tesla from GEAR or DD4hep automagically
    - getDetectorExtension returns DDRec detector extension for given flags
  - HelixClass[_double]::getDistanceToPoint, first argument is now const to allow passing return values from lcio functions. Fully backward compatible

* 2017-05-04 Andre Sailer ([PR#5](https://github.com/iLCSoft/MarlinUtil/pull/5))
  - Fixed warnings for gcc49
  - add Werror to CI configuration, no longer accepting PRs creating warnings

* 2017-06-20 Shaojun Lu ([PR#10](https://github.com/iLCSoft/MarlinUtil/pull/10))
  - Update DD4hep::DDRec::LayeredCalorimeterData to dd4hep::rec::LayeredCalorimeterData for backward compatible.
  - Adapt namespaces to changes in DD4hep

* 2017-06-03 Andre Sailer ([PR#8](https://github.com/iLCSoft/MarlinUtil/pull/8))
  - MarlinUtil::getAbsMomentum: bugfix: use `delete[]` instead of `delete` as getMomentum uses new[]

* 2017-05-29 Andre Sailer ([PR#7](https://github.com/iLCSoft/MarlinUtil/pull/7))
  - MarlinUtil::getAbsMomentum: fix memory leak
  - ClusterExtended::setHelix: use const reference instead of object as argument for function
  - HelixClass[_double]::getDistanceToHelix: fix out-of-bounds access of array
  - DDMarlinCED::*ParameterConversion: initialise isBarrel for CEDGeoTubeParams
  - DDMarlinCED:: draw: getchar returns int, value otherwise truncated
  - DDMarlinCED::kbhit: initialise fd_set

* 2017-05-22 Andre Sailer ([PR#6](https://github.com/iLCSoft/MarlinUtil/pull/6))
  - Changed `class MarlinUtil` to `namespace MarlinUtil`
  - MarlinUtilConfig: include the dependency on DD4hep and ROOT

# v01-13

# v01-12-01
  - patch release
  - replace LayeredCalorimeter::Layer::thickness
    ( needed for DD4hep v00-18 )

# v01-12

F. Gaede
   - protect against missing SurfaceManager object in DD4hep model
   - updated to use gsl 2.1
   - made compatible with c++11
   - removed -ansi -pedantic -Wno-long-long 

# v01-11

 M.Berggren
  - Added WeightedPoints3D class, a class that analyses a set of points
    in 3 dimensions (possibly with weights), and calculates things like
    covariance matrix, averge position, eigen-vectors/values ...

 F.Gaede
   M /MarlinUtil/trunk/source/src/ClusterShapes.cc
 - merge back changes from version copied to
   MarlinReco for HLR week at Desy:
    - added dummy error functions to ClusterShapes for MarlinPandora development (J.List)
    - cluster shape variables added for mu-pi seaparation (H.Sert)

 T.Quast (RWTH Aachen)
 
   - added DDMarlinCED.cc
     utility and detector draw routines for CEDViewer


# v01-10

    - added helper function in CHT (J.List)
      CHT::CaloType caloTypeFromString(const std::string& name) 

    - fixes in MarlinCED (F.Gaede)
      - fixed warning -Wc++11-narrowing
      - draw inner and outer edge of FTD disks


# v01-09
     - added functions for cluster shape calculation (M. Kurata) 
     - mad compatible for drawing CLIC detector (F.Gaede/M.Petric)

# v01-08-01
     - draw all layers for SIT and SET
       -> draw two layers per double layer

# v01-08
     - improved FPCCD classes (T.Mori)
     - added HelixClass_double 


# v01-07-01
    - added wrapper function (extern "C")
      void draw_helix( ... ) to MarlinCED::drawHelix(...)
       

# v01-07
    - updated code for the ANN library to version 1.1.2 
      (see: http://www.cs.umd.edu/~mount/ANN/ )

    - added helper function to library (extern "C")
      void drawDetectorFromGearFile( const char* fname ) ;

# v01-06-01
    - adopted for multi-module support in MarlinCED for the TPC detector


# v01-06
		- Akiya Miyamoto:
	   implemented flags to turn On and Off geometry of EcalBarrel, EcalEndcap, HcalBarrel, HcalRing, HcalEndcap individually.


# v01-05-03

    - fixed drawing of Hcal endcap (outer radius in gear file)
    - fixed uninitialised variables


# v01-05-02
    
    - bug fix MarlinCED: high pt track straight line was sometimes drawn in wrong direction
    - fixed warnings and made compatible w/ clang++


# v01-05-01

   bug fixes:

       /MarlinUtil/trunk/source/src/MarlinCED.cc:

            - updated for new GEAR attribute getNSensors
            - updated for new definition of FDT sensors (max #2)

            - changed vertex draw method from ced_gebox_r to ced_gebox_r_solid to ensure transparency
       

       /MarlinUtil/trunk/source/src/CalorimeterHitType.cc:

            - changed order of ring and endcap in order to get HcalEndCapRingsCollection as ring layout (reported by G.Grenier)


# v01-05

    - added drawing of the SET
    - added picking print funtions for TrackerHitPlane and TrackerHitZCylinder
    - changed drawing of SIT to 3D shape with surface 
    - modified to match change of MAX_LAYER to CED_MAX_LAYER in ced_cli.h (requires CED 1.4)
    - added plannar SIT and protected against Gear exceptions for missing parameters


# v01-04

    New features:

        - MarlinCED: 
	        - added more subdetectors: FTD, LHCal, LCal, Coil, Yoke, SIT

	        - changed detector colors to match the default colors used in Mokka
     
            - detector components have now layers and layer description

            - Cylinder: inner edges can be different from outeredges, and the inner
                cylinder can be rotated with an other angle as the outer cylinder shape

            - made compatible with CLIC (no LHCal and 11 FTD layers):
               LHCal is now optional and the FTD is drawn dynamically from gear....


         To use the new CED features use CED_GeoTube instead of CED_GeoCylinder.


        - added ILDCellIDEncoding which is needed for development of new tracking
          system, could be removed at a later stage of development

        - header files now installed into subdirectory (marlinutil) to avoid name clashes
            still backwards compatible

        - moved classes FPCCDData and FPCCDPixelHit from MarlinReco
          to MarlinUtil to avoid dependency of Overlay on MarlinReco

    Fixes:
    
        - MarlinCED:
            - old ced_hit_id is now deprecated (the old ced_hit_id used one int to store type and layer number).
	            The new one have one argument more, (layer and type are now 2 fields)
            - Added try/catch blocks to prevent gear::UnknownParameterException, while accessing undefined detector components.
            - fixed drawGEARDetector
            - fixed length of VXD detector ( had used half size only !)
            - allow to toggle visibility of VXD
            - made compatible with new gear (used for ILD_01):
            - use optionally FTDParameters
            - use size of SIT arrays (will have to be changed to ZPlanarLayout)
            - optionally use new coil parameter names



# v01-03

    - Changes (A. Sailer):
        
        - Rewrote the fast HelixClass::getDistanceToPoint

        - Remove unnecessary code from HelixClass:getDistanceToPoint,
            add second getDistanceToPoint function, that stops if distXY is
            larger than given distance, keeps cpu intensive functions from being run.


# v01-02

    - Changes (H. Hoelbe, DESY):
        - Implement remote access to CED:

            * Now you are able to connect to an user defined host and portnumber.

            * Example: To run CED on host "foo" port 7777 and Marlin on host "bar":
                      ssh foo
                      export CED_PORT="7777"
                      glced --trust bar

                      ssh bar
                      export CED_HOST="foo"
                      export CED_PORT="7777"
                      Marlin viewer.xml

    - Simplified cmake scripts
        Improved dependencies:
            - removed dependencies GEAR, LCIO and streamlog (now exported through Marlin)
            - exchanged CMakeModules dependency with new package ILCUTIL


# v01-01

       Changes (H. Hoelbe, DESY):
    -  New methods: MarlinCED::set_layer_description and
       MarlinCED::add_layer_description which stores the 
       layer descriptions and call ced_describe_layer only 
       one time per layer.
    -  Add usleep to the mainloop to avoid busy waiting
       (speed up CED)
    - CalorimeterHitType - added functions:
         CHT::Layout layoutFromString(const std::string& name) ;
         CHT::CaloID caloIDFromString(const std::string& name) ;
      (moved from MarlinReco/LDCCaloDigi/src/CHT_helper)
    - fixed: removed gsl include statements from public header files

# v01-00

     - added support for picking to MarlinCED 
  
     - made compliant w/ gcc 4.4.3 
     
     - made compliant w/ MacOSX


# v00-14-01

    * bug fix: incorrect library version numbers

# v00-14

        * New helper functions for MarlinCED< print(andDraw)MCFamily
        * Header file changes associated with new picking/client-server (as in MarlinCED.cc)
        * src/MarlinCED.cc: SM-H: Added Picking code 
        * Added a (as yet commented out) version of client server communication
        * Changed references to pi = acos(-1) to pi = M_PI
          (Stewart Martin-Haugh)        

# v00-13
    * added helper class CalorimeterHitType for encoding/decoding CalorimeterHit      types (for the ILD detector)

# v00-12
    * added drawing of the VTX detector in drawGEARDetector() [S.Aplin]
    * improved drawing of neutral particles in drawHelix (straight lines)
    * made cmake 2.6 compliant
    * added LIBRARY_DIRS
    * added 32 bit compatibility build option

# v00-11
    * bug fix MarlinCED.h (removed unneeded "#include ced.h" statement)
    * bug fix in CED dependency:
        -> MarMarlinUtilConfig.cmake.in now exports ced_cli.h
        -> needs CED >= v00-05
 
# v00-10
    * removed duplicated code from CED 
      -> made dependant on CED ( >= v00-04-01 ) 
     

# v00-09

    * src/TrackPair.cc: 
      New class used to associate tracks coming from neutral vertices
      (needed by MarlinReco/Tracking/V0Finder) 

# v00-08

    * Updated HelixClass. (A.Raspereza)
      added new  method added which identifies neutral vertices

    * cleaned up code and made gcc4.3 compliant (S.Aplin)

# v00-07
   * removed optional AIDA dependency ( SelectEvents.cc )  J.Engels

   * fixed compiler warnings  J.Samson

   * removed  deprecated HelixClass::getDistanceToHelix()

   * modfied MarlinCED::newEvent to not draw detector for modelID=-1

# v00-06
	* MarlinCED:  modified drawGEARDetector(): incl. Hcal ring/
	  removed overlap Ecal/Hcal (M.Morgunov)

        * Added Circle class: finds a circle from 3 points. Used by TPCDigi. 
          (S.Aplin)

        * improved cmake scriots (J.Engels)

# v00-05
    - improved versions of Track utility classes (A.Raspereza)

    - made compatible with CLHEP 2.0.3.x
        - removed HepPDT dependency

    - added code for writing millipede files 

    - bug fixes:
        -  Doxyfile predefined USE_CLHEP
	-  patched for 64 bit (M.Killenberg)
        - ...

# v00-04
      - cmake now default build tool (J.Engels)
         -> see README for instructions

      - TrackwiseClusters (O.Wendt) 
          improved clustering in 'trajectory coordinates': 
          path-length and distance to trajectory


      - MarlinCED / ced_cli.cc (A.Bulgheroni)
          added drawGEARTelescope
	  added geo_boxes

      - ClusterShapes (O.Wendt) Fixed a bug in FitHelix(int max_iter, int
	  status_out, int parametrisation,float* parameter, float*
	  dparameter, float& chi2,float& distmax, int direction), thanks to
	  Hengne Li @ LAL

      - additional bug fixes 
          + compliance with SL4 (gcc3.4/gcc4)
          + ...

      see Changelog for details 

# v00-03
  - suport for Full Tracking  (A.Raspereza)

  -  suport for TrackBasedPFA (O.Wendt)

  - suport for PhotonFinderKit (P. Krstonosic)

  - first version of Trajectory implementation
    SimpleLine and SimpleHelix (Th. Kraemer)

  - included the ANN (Aproximate Nearest Neighbor) library
    [ see: http://www.cs.umd.edu/~mount/ANN ]

  - added cmake support (epxerimental)
 

 for details see the ChangeLog file


