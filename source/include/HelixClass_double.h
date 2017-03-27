#ifndef HelixClass_double_H
#define HelixClass_double_H 1
#include <vector>
#include <cmath>

/**
 *    Utility class to manipulate with different parameterisations <br>
 *    of helix. Helix can be initialized in a three different <br>
 *    ways using the following public methods : <br>
 *    1) Initialize_VP(float * pos, float * mom, float q, float B) : <br>
 *       initialization of helix is done using <br>
 *       - position of the reference point : pos[3]; <br>
 *       - momentum vector at the reference point : mom[3];<br>
 *       - particle charge : q;<br>
 *       - magnetic field (in Tesla) : B;<br>
 *    2) Initialize_BZ(float xCentre, float yCentre, float radius, <br>
 *				   float bZ, float phi0, float B, float signPz,<br>
 *				   float zBegin):<br>
 *       initialization of helix is done according to the following<br>
 *       parameterization<br>
 *       x = xCentre + radius*cos(bZ*z + phi0)<br>
 *       y = yCentre + radius*sin(bZ*z + phi0)<br>
 *       where (x,y,z) - position of point on the helix<br>
 *       - (xCentre,yCentre) is the centre of circumference in R-Phi plane<br>
 *       - radius is the radius of circumference<br>
 *       - bZ is the helix slope parameter<br>
 *       - phi0 is the initial phase of circumference<br>
 *       - B is the magnetic field (in Tesla)<br>
 *       - signPz is the sign of the z component of momentum vector<br>
 *       - zBegin is the z coordinate of the reference point;<br>
 *    3) void Initialize_Canonical(float phi0, float d0, float z0, float omega,<br> 
 *			      float tanlambda, float B) :<br>
 *    canonical (LEP-wise) parameterisation with the following parameters<br>
 *       - phi0 - Phi angle of momentum vector at the point of<br>
 *         closest approach to IP in R-Phi plane;
 *       - d0 - signed distance of closest approach to IP in R-Phi plane;<br>
 *       - z0 - z coordinate of the point of closest approach in R-Phi plane;<br>
 *       - omega - signed curvature;<br>
 *       - tanlambda - tangent of dip angle;<br>
 *       - B - magnetic field (in Tesla);<br>
 *    A set of public methods (getters) provide access to <br>
 *    various parameters associated with helix. Helix Class contains<br>
 *    also several utility methods, allowing for calculation of helix<br>
 *    intersection points with planes parallel and perpendicular to <br>
 *    z (beam) axis and determination of the distance of closest approach<br>
 *    from arbitrary 3D point to the helix. <br>
 *    @author A. Raspereza (DESY)<br>
 *    @version $Id: HelixClass_double.h,v 1.16 2008-06-05 13:47:18 rasp Exp $<br>
 *
 */

#include "LineClass.h"
class HelixClass_double;

class HelixClass_double {
 public:

/**
 *  Constructor. Initializations of constants which are used
 *  to calculate various parameters associated with helix.
 */ 
    HelixClass_double();
/**
 *  Destructor. 
 */
    ~HelixClass_double();
/**
 *   Initialization of helix using <br>
 *     - position of the reference point : pos[3]; <br>
 *     - momentum vector at the reference point : mom[3];<br>
 *     - particle charge : q;<br>
 *     - magnetic field (in Tesla) : B<br>
 */  
    void Initialize_VP(double * pos, double * mom, double q, double B);

/**
 *   Initialization of helix according to the following<br>
 *   parameterization<br>
 *   x = xCentre + radius*cos(bZ*z + phi0)<br>
 *   y = yCentre + radius*sin(bZ*z + phi0)<br>
 *     where (x,y,z) - position of point on the helix<br>
 *     - (xCentre,yCentre) is the centre of circumference in R-Phi plane<br>
 *     - radius is the radius of circumference<br>
 *     - bZ is the helix slope parameter<br>
 *     - phi0 is the initial phase of circumference<br>
 *     - B is the magnetic field (in Tesla)<br>
 *     - signPz is the sign of the z component of momentum vector<br>
 *     - zBegin is the z coordinate of the reference point<br>
 */  
    void Initialize_BZ(double xCentre, double yCentre, double radius, 
				   double bZ, double phi0, double B, double signPz,
				   double zBegin);
/**
 *  Canonical (LEP-wise) parameterisation with the following parameters<br>
 *     - phi0 - Phi angle of momentum vector at the point of<br>
 *       closest approach to IP in R-Phi plane;
 *     - d0 - signed distance of closest approach in R-Phi plane;<br>
 *     - z0 - z coordinate of the point of closest approach to IP 
 *       in R-Phi plane;<br>
 *     - omega - signed curvature;<br>
 *     - tanlambda - tangent of dip angle;<br>
 *     - B - magnetic field (in Tesla)<br>
 */  
    void Initialize_Canonical(double phi0, double d0, double z0, double omega, 
			      double tanlambda, double B);
    /**
     *  Returns momentum of particle at the point of closest approach <br>
     *  to IP <br>
     */
    const double * getMomentum();

    /**
     *  Returns reference point of track <br>     
     */
    const double * getReferencePoint();

    /**
     *  Returns Phi angle of the momentum vector <br>
     *  at the point of closest approach to IP <br>     
     */
    double getPhi0();

    /**
     *  Returns signed distance of closest approach <br>
     *  to IP in the R-Phi plane <br>     
     */
    double getD0();

    /**
     *  Returns z coordinate of the point of closest 
     *  approach to IP in the R-Phi plane <br>     
     */
    double getZ0();

    /**
     *  Returns signed curvature of the track <br>     
     */
    double getOmega();

    /**
     *  Returns tangent of dip angle of the track <br>     
     */
    double getTanLambda();

    /**
     *  Returns transverse momentum of the track <br>     
     */
    double getPXY();


    /**
     *  Returns x coordinate of circumference
     */
    double getXC();

    /**
     *  Returns y coordinate of circumference
     */
    double getYC();


     /**
     *  Returns radius of circumference
     */
    double getRadius();   


    /**
     *  Returns helix intersection point with the plane <br>
     *  parallel to z axis. Plane is defined by two coordinates <br>
     *  of the point belonging to the plane (x0,y0) and normal <br>
     *  vector (ax,ay).  ref[3] is the reference point of the helix. <br>
     *  point[3] - returned vector holding the coordinates of <br>
     *  intersection point <br>     
     */
    double getPointInXY(double x0, double y0, double ax, double ay, 
			      double * ref , double * point);

    /**
     *  Returns helix intersection point with the plane <br>
     *  perpendicular to z axis. Plane is defined by z coordinate <br>
     *  of the plane. ref[3] is the reference point of the helix. <br>
     *  point[3] - returned vector holding the coordinates of <br>
     *  intersection point <br>     
     */
    double getPointInZ(double zLine, double * ref, double * point);

    /**
     * Return distance of the closest approach of the helix to <br>
     * arbitrary 3D point in space. xPoint[3] - coordinates of <br>
     * space point. Distance[3] - vector of distances of helix to <br> 
     * a given point in various projections : <br>
     * Distance[0] - distance in R-Phi plane <br>
     * Distance[1] - distance along Z axis <br>
     * Distance[2] - 3D distance <br> 
     */
    double getDistanceToPoint(double * xPoint, double * Distance);

    /**
     * Return distance of the closest approach of the helix to <br>
     * arbitrary 3D point in space. xPoint[3] - coordinates of <br>
     * space point. distCut - limit on the distance between helix <br> 
     * and the point to reduce calculation time <br>
     * If R-Phi is found to be greater than distCut, rPhi distance is returned <br>
     * If the R-Phi distance is not too big, than the exact 3D distance is returned <br>
     * This function can be used, if the exact distance is not always needed <br>
     */
    double getDistanceToPoint(const double* xPoint, double distCut);
    double getDistanceToPoint(const std::vector<double>& xPoint, double distCut);

    /**
     * This method calculates coordinates of both intersection <br>
     * of the helix with a cylinder. <br>
     * Rotational symmetry with respect to z axis is assumed,  <br>
     * meaning that cylinder axis corresponds to the z axis <br>
     * of reference frame. <br>
     * Inputs : <br> 
     * Radius - radius of cylinder. <br>
     * ref[3] - point of closest approach to the origin of the helix. <br>
     * Output : <br>
     * point[6] - coordinates of intersection point. <br>
     * Method returns also generic time, defined as the <br>
     * ratio of helix length from reference point to the intersection <br>
     * point to the particle momentum <br>
     */
    double getPointOnCircle(double Radius, double * ref, double * point);

    /** Returns distance between two helixes <br>
     * Output : <br> 
     * pos[3] - position of the point of closest approach <br>
     * mom[3] - momentum of V0 <br>
     */
    double getDistanceToHelix(HelixClass_double * helix, double * pos, double * mom);

    /**
     * Set Edges of helix 
     */
    void setHelixEdges(double * xStart, double * xEnd);

    /**
     * Returns starting point of helix
     */
    double * getStartingPoint() {return _xStart;}

    /**
     * Returns endpoint of helix
     */
    double * getEndPoint() {return _xEnd;}
    
    /**
     * Returns BZ for the second parameterization
     */
    double getBz();

    /**
     * Returns Phi for the second parameterization
     */
    double getPhiZ();

    /**
     * Returns extrapolated momentum
     */
    void getExtrapolatedMomentum(double * pos, double * momentum);

    /**
     * Returns charge 
     */
    double getCharge();

 private:    
    double _momentum[3]; // momentum @ ref point 
    double _referencePoint[3]; // coordinates @ ref point
    double _phi0=0; // phi0 in canonical parameterization
    double _d0=0;   // d0 in canonical parameterisation
    double _z0=0;   // z0 in canonical parameterisation
    double _omega=0; // signed curvuture in canonical parameterisation
    double _tanLambda=0; // TanLambda
    double _pxy=0; // Transverse momentum
    double _charge=0; // Particle Charge
    double _bField=0; // Magnetic field (assumed to point to Z>0)
    double _radius=0; // radius of circle in XY plane
    double _xCentre=0; // X of circle centre
    double _yCentre=0; // Y of circle centre
    double _phiRefPoint=0; // Phi w.r.t. (X0,Y0) of circle @ ref point
    double _phiAtPCA=0; // Phi w.r.t. (X0,Y0) of circle @ PCA
    double _xAtPCA=0; // X @ PCA
    double _yAtPCA=0; // Y @ PCA
    double _pxAtPCA=0; // PX @ PCA
    double _pyAtPCA=0; // PY @ PCA
    double _phiMomRefPoint=0; // Phi of Momentum vector @ ref point
    double _const_2pi=2.0*M_PI; // 2*PI
    double _const_pi2=0.5*M_PI; // PI/2
    double _FCT=2.99792458E-4; // 2.99792458E-4
    double _xStart[3]; // Starting point of track segment
    double _xEnd[3]; // Ending point of track segment

    double _bZ=0;
    double _phiZ=0;

};


#endif
