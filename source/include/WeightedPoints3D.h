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


/**
 *    Utility class to derive properties of a set od weighted points in 3D, 
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
   * coordiante system, i.\ e.\ IP, to average (aka the centre of gravity) of the 
   * set of points, with the weighs of each point taken into account.
   */
  double* getCentreOfGravity();
 
  /**
   * returns the covarinace matrix between the three coordinates of the
   * centre of gravity. Returned as a vector storing the matrix in row-major order.
   */

  double* getCentreOfGravityErrors();

  /** US spelling of getCentreOfGravity */
  inline double* getCenterOfGravity() { return getCentreOfGravity() ; }
  inline double* getCenterOfGravityErrors() { return getCentreOfGravityErrors() ; }
  
  /**
   * array of the eigen-values of the covariance matrix, sorted in ascending size.
   */
  double* getEigenVal();

  /**
   * Retuns a vector containing the variances of the three eigen-values.
   */

  double* getEigenValErrors();

  /**
   * array of the three main axes of the covariance matrix, starting
   * with the axis corresponding to the smallest covariance.
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
   *  is a (2/3) x (2/3) x 3 matrix (2 for polar, 3 for cartesian), which is
   *  returned as a vector in depth->row major order.
   */
  double* getEigenVecCartesianErrors();
  double* getEigenVecPolarErrors();

  /**
   * 'mean' width of the cluster perpendicular to the main 
   * principal axis, defined as: 
   * width := sqrt( 1/Sum(wgt) * Sum(wgt_i * d^2 ))
   * where d[i] is the distance of the i-th point to the main
   * principal axis.
   */
  double getWidth();

 
/**
   * distance to the centre of gravity measured from IP
   * (absolut value of the vector to the centre of gravity)
   */
  inline double radius() { return _radius; }


private:

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
  std::vector<double> _t;
  std::vector<double> _s;

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


  int   _ifNotWidth;
  double _Width;


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
  double _EigenVec[3][3];
  double _EigenVecAngle[2][3];
  double _theta_phi_cov[2][2][3];
  double _xyz_cov[3][3][3];
  double _var_lam[3];
  double _nfact1;
  double _nfact2;
  double _nfact3;
  double _nfact4;
  double _fourth_mom[3][3][3][3];
  double dang[6][2][3]; 



  void  findCOG();
  void  findCOGErrors();
  double findDistance(int i);
  void  solveEigenValEq();
  void  findCovErrors();
  void  findWidth();
  void  propagateVecErrorsPolar();
  void  propagateVecErrorsCartesian();
  void  propagateValErrors();
  void  find_dVec_dCov();

};


#endif
