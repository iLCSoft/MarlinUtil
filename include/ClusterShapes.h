#ifndef ClusterShapes_h
#define ClusterShapes_h


#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_rng.h>




/**
 *    Utility class to derive properties of clusters (that means sets of points),
 *    such as centre of gravity, axes of inertia and so on.
 *
 *    @authors V. Morgunov (ITEP/DESY), A. Raspereza (DESY), O. Wendt (DESY)
 *    @version $ld: $
 *
 */
class ClusterShapes {

 public:

  /**
   *    Constructor:
   *    nhits : number of hits in the cluster
   *    a     : amplitudes of elements ('cells') of the cluster. Stored in an array,
   *            with one entry for each element ('cell'). Each entry is depending on
   *            coordinates x,y,z (Cartesian), which are stored in the arrays x,y,z.
   *    x,y,z : array of coordinates corresponding to the array of amplitudes a.
   *
   *
   */
  ClusterShapes(int nhits, float* a, float* x, float* y, float* z);

  ~ClusterShapes();


  /**
   * returns the number of elements of the cluster
   */
  int getNumberOfHits();

  /**
   * returns the accumulated amplitude for the whole cluster
   */
  float getTotalAmplitude();

  /**
   * returns an 'vector' from the origin to the centre of gravity
   * (weighted with the amplitudes per element) of the cluster
   */
  float* getCentreOfGravity();

  /**
   * array of the inertias of mass ('amplitudes') corresponding
   * to the three main axes of inertia. The array is sorted in
   * ascending order.
   */
  float* getEigenValInertia();

  /**
   * array of the three main axes of inertia (9 entries) starting
   * with the axis corresponding to the smallest inertia of mass 
   * (main principal axis). All axes are normalised to a length 
   * of 1.
   */
  float* getEigenVecInertia();

  /**
   * 'mean' width of the cluster perpendicular to the main 
   * principal axis, defined as: 
   * width := sqrt( 1/TotalAmplitude * Sum(a[i]*d[i]*d[i]) ),
   * where d[i] is the distance of the i-th point to the main
   * principal axis.
   */
  float getWidth();

  /**
   * performs a least square fit on the shape of an electro-
   * magnetic-shower, which is defined as:
   * A[i] = a * xl[i]^b * exp(-c*xl[i]) * exp(-d*xt[i]),
   * where A[i] is the array of amplitudes, xl[i] is the 
   * coordinate of the actual point along the main principal 
   * axis and xt[i] the coordinate perpendicular to it. a,b,c,d  
   * are the parameters to be fitted.    
   * The method returns the chi2 and the parameters a,b,c,d of
   * the fit as well as xStart, which is an 3-dim array to the 
   * point closest to IP.
   * The return value of the method itself is not used at the
   * moment (always returns 0).
   */
  int Fit3DProfile(float& chi2, float& a, float& b, float& c, float& d, float * xStart);

  /**
   * returns the chi2 of the fit in the method Fit3DProfile
   * for a given set of parameters a,b,c,d
   */
  float getChi2Fit3DProfile(float a, float b, float c, float d);

  /**
   * performs a least square fit on a helix path in space, which
   * which is defined as (Cartesian coordiantes):
   * x[i] = x0 + R*cos(b*z[i] + phi0)
   * y[i] = y0 + R*sin(b*z[i] + phi0)
   * z[i] = z[i],
   * where x0,y0,R,b and phi0 are the parameters to be fitted and
   * x[i],y[i],z[i] are the (Cartesian) coordiantes of the space
   * points.
   * The following output/input parameters are returned/needed:
   *
   * OUTPUTS:
   * method itself : returns 1 if an error occured and 0 if not
   * parameter  : array of parameters to be fitted (defined
   *              as : parameter[5] = {x0,y0,R,b,phi0}
   * dparameter : error on the parameters, that means: 
   *              dparameter[i] = sqrt( CovarMatrix[i][i] )
   * chi2       : chi2 of the fit
   * distmax    : maximal distance between the points x[i],y[i]
   *              z[i] and the fitted function
   *
   * INPUTS:
   * parametrisation : 1 for first and 2 for second parametrisation
   * max_iter   : maximal number of iterations, which should be 
   *              performed in the fit
   * status_out : if set to 1, only the initial parameters of
   *              the fit are calculated and are stored in
   *              parameter. The entries of dparameter are
   *              set to 0.0
   */

  // CHANGE DOCUMENTATION!!!!
 
  int FitHelix(int max_iter, int status_out, int parametrisation,
	       float* parameter, float* dparameter, float& chi2, float& distmax);

  /**
   * distance to the centre of gravity measured from IP
   * (absolut value of the vector to the centre of gravity)
   */
  inline float radius() { return _radius; }

  /**
   * largest spatial axis length of the ellipsoid derived
   * by the inertia tensor (by their eigenvalues and eigen-
   * vectors)
   */
  inline float getElipsoid_r1() { return _r1; }

  /**
   * medium spatial axis length of the ellipsoid derived
   * by the inertia tensor (by their eigenvalues and eigen-
   * vectors)
   */
  inline float getElipsoid_r2() { return _r2; }

  /**
   * smallest spatial axis length of the ellipsoid derived
   * by the inertia tensor (by their eigenvalues and eigen-   
   * vectors)
   */
  inline float getElipsoid_r3() { return _r3; }

  /**
   * volume of the ellipsoid
   */
  inline float getElipsoid_vol() { return _vol; }

  /**
   * average radius of the ellipsoid (qubic root of volume)
   */
  inline float getElipsoid_r_ave() { return _r_ave; }

  /**
   * density of the ellipsoid defined by: totAmpl/vol
   */
  inline float getElipsoid_density() { return _density; }

  /**
   * eccentricity of the ellipsoid defined by: 
   * Width/r1
   */
  inline float getElipsoid_eccentricity() { return _eccentricity; }

  /**
   * distance from centre of gravity to the point most far 
   * away from IP projected on the main principal axis
   */
  inline float getElipsoid_r_forw() { return _r1_forw; }

  /**
   * distance from centre of gravity to the point nearest 
   * to IP projected on the main principal axis    
   */
  inline float getElipsoid_r_back()       { return _r1_back; }





 private:

  int _nHits;

  float* _aHit;
  float* _xHit;
  float* _yHit;
  float* _zHit;

  int   _ifNotGravity;
  float _totAmpl;
  float _radius;
  float _xgr;
  float _ygr;
  float _zgr;
  float _analogGravity[3];

  int   _ifNotWidth;
  float _analogWidth;

  int   _ifNotInertia;
  float _ValAnalogInertia[3];
  float _VecAnalogInertia[9];

  int   _ifNotElipsoid;
  float _r1           ;  // Cluster spatial axis length -- the largest
  float _r2           ;  // Cluster spatial axis length -- less
  float _r3           ;  // Cluster spatial axis length -- less
  float _vol          ;  // Cluster ellipsoid volume
  float _r_ave        ;  // Cluster average radius  (qubic root)
  float _density      ;  // Cluster density
  float _eccentricity ;  // Cluster Eccentricity
  float _r1_forw      ;
  float _r1_back      ;

  void  findElipsoid();
  void  findGravity();
  void  findInertia();
  void  findWidth();
  float findDistance(int i);
  float vecProduct(float * x1, float * x2);
  float vecProject(float * x, float * axis);
  float DistanceHelix(float x, float y, float z, float X0, float Y0, float R0, float bz, float phi0);


  // private methods for non-linear, multidim. fitting (helix)
  // static int functParametrisation1(const gsl_vector* par, void* data, gsl_vector* f);
  // static int dfunctParametrisation1(const gsl_vector* par, void* d, gsl_matrix* J);
  // static int fdfParametrisation1(const gsl_vector* par, void* d, gsl_vector* f, gsl_matrix* J);


};

#endif
