#ifndef WeightedPoints3D_h
#define WeightedPoints3D_h


#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include<vector>

/** these are the numbers to scale the three standard deviations by to get probability 
 * inside a 3dim Gaussian equal to 0.64, 0.90, 0.95 and 0.99, respectively.
 */

#define  _one_sigma 1.56 
#define  _CL90  2.12
#define  _CL95  2.39 
#define  _CL99  2.93

/**
 *    Utility class to derive properties of a set of weighted points in 3D, 
 *    such as centre of gravity, eigen-values and eigen-vectors, and their
 *    uncertainties.
 *
 *    @authors M. Berggren (DESY)
 *    @version $Id: WeightedPoints3D.h,v 1.0
 *
 */
class WeightedPoints3D {

public:

  /**
   *    Constructor
   *    @param nhits : number of hits in the cluster
   *    @param a     : amplitudes of elements ('cells') of the cluster. Stored in 
   *                   an array, with one entry for each element ('cell'). Each entry 
   *                   is depending on coordinates x,y,z (Cartesian), which are stored 
   *                   in the arrays x,y,z.
   *    @param x,y,z : array of coordinates corresponding to the array of amplitudes a.
   *
   *
   */
  WeightedPoints3D(int nhits, double* a, double* x, double* y, double* z);
  WeightedPoints3D(const std::vector<double> &cog,
		   const std::vector<double> &cov, const std::vector<double> &mayor_axis_error =  std::vector<double>() , 
		   int npnt = 0 ,double wgtsum =0.0 , double wgt2sum=0.0 , double wgt4sum=0.0 );

  /**
   *    Destructor
   */
  ~WeightedPoints3D();

  /**
   *    Defining position errors. (Currently not used) 
   */
  void setErrors(double *ex, double* ey, double *ez);



  /**
   * returns the number of elements of the set of points
   */
  int getNumberOfPoints();

  /**
   * returns the summed weight for the whole set of points.
   * Applied to calorimeter clusters this is usually, but not
   * always, the reconstructed energy of the cluster. It might
   * differ if down-stream corrections (crack-corrections etc.)
   * were applied.
   */
  double getTotalWeight();

  /**
   * returns the sum of weight^2 for the whole set of points.
   */
  double getTotalSquaredWeight();

  /**
   * returns the sum of weight^4 for the whole set of points.
   */
  double getTotalQuarticWeight();

   /**
   * returns an array, which represents a vector from the origin of the
   * coordinate system, i. e. the nominal IP, to the average (aka the centre of 
   * gravity) of the set of points, with the weighs of each point taken into account.
   */
  double* getCentreOfGravity();
 
  /**
   * returns the covariance matrix between the three coordinates of the
   * centre of gravity. Returned as a vector storing the matrix in row-major order.
   */

  double* getCentreOfGravityErrors();

  /** US spelling of getCentreOfGravity */
  inline double* getCenterOfGravity() { return getCentreOfGravity() ; }
  inline double* getCenterOfGravityErrors() { return getCentreOfGravityErrors() ; }
  
  /**
   * Array of the eigen-values of the covariance matrix, sorted in ascending size.
   * The eigen-values are the _variances_ of the marginal distribution of the weighted
   * points, projected on the corresponding eigen-vectors
   */
  double* getEigenVal();

  /**
   * Returns a vector containing the variances of the three eigen-values. NB that
   * if the points are not available (i.e. if the input is a DST-file), the errors
   * on the eigen-values are approximative.
   */

  double* getEigenValErrors();

  /**
   * array of the three main axes of the covariance matrix, starting
   * with the axis corresponding to the smallest variance.
   * The main principal axis is thus the last one. All axes are normalised to a length 
   * of 1, and they form a right-handed system (ie. e_3 = e_1 X e_2).
   * the eigen-values are the columns of the matrix, which is returned as a vector
   * storing the matrix in row-major order.
   * 
   */
  double* getEigenVecCartesian();

  /**
   * The same in polar. Each eigen-vector is represented by two angles in the order
   * (theta, phi). Also returned as a vector storing the matrix in row-major order.
   */  
  double* getEigenVecPolar();
 
  /**
   *  The covariance matrix of the eigen-vectors, one matrix/eigen-vector, so it
   *  is a (2/3) x (2/3) x 3 matrix (2 for polar, 3 for Cartesian), which is
   *  returned as a vector in depth->row major order.
   *
   * NB that if the points are not available (i.e. if the input is a DST-file), the errors
   * on the first and second eigen-vector are approximative, while the third (largest)
   * one is exact.
   */
  double* getEigenVecCartesianErrors();
  double* getEigenVecPolarErrors();

   
  /**
   * distance to the centre of gravity measured from IP
   * (absolute value of the vector to the centre of gravity)
   */
  inline double radius() { return _radius; }

  /**
   * Transform a point point-distribution to the Eigen system (i.e. a system
   * with the origin at the C.O.G, the z-axis along the longest eigen-vector.
   */
  double*  TransformPointToEigenSyst(double* xyz) ;

  // stolen from cluster-shapes:

  /* float getEmax(float* xStart, int& index_xStart, float* X0, float* Rm); */


  /* //shower max of the hits from the shower start hit */
  /* float getsmax(float* xStart, int& index_xStart, float* X0, float* Rm); */

  /* //radius where 90% of the cluster energy exists */
  /* float getxt90(float* xStart, int& index_xStart, float* X0, float* Rm); */

  /* //length where less than 20% of the cluster energy exists */
  /* float getxl20(float* xStart, int& index_xStart, float* X0, float* Rm); */

  /*  //Mean of the radius of the hits */
  /* float getRhitMean(float* xStart, int& index_xStart, float* X0, float* Rm); */

  /* //RMS of the radius of the hits */
  /* float getRhitRMS(float* xStart, int& index_xStart, float* X0, float* Rm); */

  /**
   * medium spatial semi-axis lengths of the ellipsoid derived
   * by the covariance matrix, ie.the eigen-values from largest to smallest.
   * The values can be scaled with cl3d such that, in the Gaussian case, the probability
   * the probability to be inside the ellipsoid is the generic CLs (1 sigma(~64), 90, 95 or 
   * 99% - use the symbols  _one_sigma, _CL90, _CL95 or _CL99 for that. As can be seen,
   * the default is no scaling.
   *
   * (Note that the eigen-values corresponds to the *variance*, r1-r3 to the *standard deviation*)
   */
  double getElipsoid_r1(double cl3d = 1.0 );
  double getElipsoid_r2(double cl3d = 1.0 ) ;
  double getElipsoid_r3(double cl3d = 1.0 ) ;

  /**
   *  errors on the same
   *
   * NB that if the points are not available (i.e. if the input is a DST-file), the errors
   * are approximative.
   */
  double getElipsoid_r1Error(double cl3d = 1.0) ;
  double getElipsoid_r2Error(double cl3d = 1.0) ;
  double getElipsoid_r3Error(double cl3d = 1.0) ;


 
  /**
   * volume of the ellipsoid
   */
  double getElipsoid_vol(double cl3d = 1.0) ;

  /**
   * average radius of the ellipsoid (cubic root of volume)
   */
  double getElipsoid_r_ave(double cl3d = 1.0) ;


  /**
   * density of the ellipsoid defined by: totAmpl/vol
   */
  double getElipsoid_density(double cl3d = 1.0);


  /**
   * Longitudinal eccentricity of the ellipsis defined by: 
   * sqrt(1-((r2+r3)/2)^2/r1^2) (I.e 0 for a circle, close 
   * to 1 for a very elongated ellipsis)
   */
  double getLongitudinalElipsis_eccentricity();

  /**
   * Transverse eccentricity of the ellipsis defined by: 
   * sqrt(1-r2^2/r3^2)(I.e 0 for a circle, close 
   * to 1 for a very elongated ellipsis)
   */
  double getTransverseElipsis_eccentricity();

  /*************
   * The following methods only works if the points are available,
   * that is, if the input data is a Rec-file, and the object
   * has been instantiated with the second above c'tor.
   ************
   */

/**
   * 'mean' width of the cluster perpendicular to the main 
   * principal axis, defined as: 
   * width := sqrt( Sum(wgt_i * d^2 )/Sum(wgt)  )
   * where d[i] is the distance of the i-th point to the main
   * principal axis. 
   */
  double getWidth();

  /**
   * Get the hit furthest away from the main shower axis.
  */

  double getMaxDist();

  /**
   * distances from the centre of gravity to the projection
   * on the main axis of first (getElipsoid_r_back) and last
   * (getElipsoid_r_forw) points in the set. The eigen-vector
   * of the largest eigen-value is oriented such that it points
   * outwards, the projection on the axis of first point is
   * the one closest to the IP, the last point the one furtherest 
   * away. The distance between the first and last points (projected 
   * on the axis) is therefore 
   * getElipsoid_r_forw()+getElipsoid_r_back() 
   */
  double getElipsoid_r_forw();
  double getElipsoid_r_back() ;

  /**
   * fraction of weight ( ~energy) inside ellipsoid, with axes scaled by cl3d.
   * For a 3D-Gaussian, one would expect 0.95 for cl3d = _CL95, etc. so this
   * method can be used to gauge the "non-gaussianity" of the distribution of
   * the points.
   */
  double getElipsoid_FractionInside(double cl3d = 1.0);


  /** 
   * Functions to change coordinates of the hits: To the eigen-system, along some
   * direction, with origin on the hit closest to the main axis, or with an
   * arbitrary point as origin (e.g.: set tht, pht, and xyz_start as the track-state
   * at the calorimeter, and check how well this matches the cluster)
   */

  void  TransformToEigenSyst() ;
  void  TransformAlongDirection(double tht , double pht , int on_axis = 1 ) ;
  void  TransformAlongDirection(double tht , double pht , double* xyz_start , int on_axis = 1) ;

  /**
   * getters for the results of the transformations
   */
  inline std::vector<double> get_x_trans() {return _x_trans ;};
  inline std::vector<double> get_y_trans() {return _y_trans ;};
  inline std::vector<double> get_z_trans() {return _z_trans ;};
  inline double* get_COG_trans() {return _COG_trans;};
  inline double* get_COGCov_trans() {return _COGCov_trans[3];};
  inline double get_th_ref(){return _th_ref;};
  inline double get_ph_ref(){return _ph_ref;};
  inline double* get_xyz_ref(){return _xyz_ref;};


protected:
 //private:
  int _nPoints;

  std::vector<double> _Wgt;
  std::vector<double> _x;
  std::vector<double> _y;
  std::vector<double> _z;
  std::vector<double> _ex;
  std::vector<double> _ey;
  std::vector<double> _ez;
  std::vector<double> _xl;
  std::vector<double> _xt;
  std::vector<double> _x_trans;
  std::vector<double> _y_trans;
  std::vector<double> _z_trans;
  
  int _current_transformation; // interpretation of transformed variables:
                               // 0 : none, ie. _x_trans = _x etc.
                               // 1 : eigen-system
                               // 2 : as 1, but reference point is the point closest to the origin (rather than the COG)
                               // 3 : as 2, but the reference point is projected unto the major axis.
                               // 4 : as 2, but the reference point is externally defined.
                               // 5 : as 4, but the reference point is projected unto the major axis.
                               // 6 : A general transformation, with z-axis direction given by _th_ref and _ph_ref,
                               //     reference point by _xyz_ref
  double _th_ref;
  double _ph_ref;
  double _xyz_ref[3];
  int _ind_first;


private:
  int   _ifNotPointsGiven;
  int   _ifNotCOG;
  int   _ifNotCOGErrors;
  double _sumWgt;
  double _wgt_sqr_sum ;
  double _wgt_4_sum ;
  double _radius;
  double _xgr;
  double _ygr;
  double _zgr;
  double _COG[3];
  double _COGCov[3][3];
  double _COGCovCov[6][6];
  double _COG_trans[3];
  double _COGCov_trans[3][3];
  double _RotToEigen[3][3];

  int   _ifNotWidth;
  double _Width;
  int   _ifNotMaxDist;
  double _MaxDist;
  double xyz_eigen[3];

  int   _ifNotEigenSolved;
  int   _ifNotEigenToPolarDone ;
  int   _ifNotCovErrors;
  int   _ifNotVecErrorsPolarPropagated;
  int   _ifNotVecErrorsCartesianPropagated;
  int   _ifNotValErrorsPropagated;
  int   _ifNot_dVec_dCov;
  int   _ifNot_fourthmom;
  int   _last_evec_to_error_propagate;
  double _EigenVal[3];
  // double _EigenValCov[3][3];  // Is this really not used: Check F08 code !
  double _EigenVec[3][3];
  double _EigenVecAngle[2][3];
  // double _EigenVecCov[9][9]; // Is this really not used: Check F08 code !
  double _theta_phi_cov[2][2][3];
  double _xyz_cov[3][3][3];
  double _var_lam[3];
  double _nfact1;
  double _nfact2;
  double _nfact3;
  double _nfact4;
  double _fourth_mom[3][3][3][3];
  double dang[6][2][3]; 
  double _r1           ;  // Cluster spatial axis length -- the largest
  double _r2           ;  // Cluster spatial axis length -- less
  double _r3           ;  // Cluster spatial axis length -- less
  double _vol          ;  // Cluster ellipsoid volume
  double _r_ave        ;  // Cluster average radius  (cubic root)
  double _density      ;  // Cluster density
  double _logitudinaleccentricity ;  // Cluster Eccentricity longitudinally
  double _transverseeccentricity ;  // Cluster Eccentricity transversally
  int _ifNotfirst_and_last_found ;
  double _r1_forw      ;  
  double _r1_back      ;
  double _FractionInside ; 


  void  findCOG();
  void  findCOGErrors();
  double findDistance(int i);
  void  solveEigenValEq();
  void  findCovErrors();
  void  findWidth();
  void  findMaxDist();
  void  propagateVecErrorsPolar();
  void  propagateVecErrorsCartesian();
  void  propagateValErrors();
  void  find_dVec_dCov();
  void  findFirstAndLast() ;
  void  findElipsoid_FractionInside(double cl3d);
};


#endif
