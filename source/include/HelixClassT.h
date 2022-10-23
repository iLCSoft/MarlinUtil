#ifndef HELIXCLASST_H
#define HELIXCLASST_H 1

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
 *    @version $Id: HelixClassT.h,v 1.16 2008-06-05 13:47:18 rasp Exp $<br>
 *
 */

#include "LineClass.h"

template<typename FloatT>
class HelixClassT {
 public:
    /**
     * The floatint point type used internally. Useful for generic programming
     */
    using float_type = FloatT;
/**
 *  Constructor. Initializations of constants which are used
 *  to calculate various parameters associated with helix.
 */
    HelixClassT() = default;
/**
 *  Destructor.
 */
    ~HelixClassT() = default;
/**
 *   Initialization of helix using <br>
 *     - position of the reference point : pos[3]; <br>
 *     - momentum vector at the reference point : mom[3];<br>
 *     - particle charge : q;<br>
 *     - magnetic field (in Tesla) : B<br>
 */
    void Initialize_VP(FloatT * pos, FloatT * mom, FloatT q, FloatT B);

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
    void Initialize_BZ(FloatT xCentre, FloatT yCentre, FloatT radius,
				   FloatT bZ, FloatT phi0, FloatT B, FloatT signPz,
				   FloatT zBegin);
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
    void Initialize_Canonical(FloatT phi0, FloatT d0, FloatT z0, FloatT omega,
			      FloatT tanlambda, FloatT B);

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
 *	   - pos - reference point at which the helix parameteres are given
 *	     (by default PCA is used as the reference point);<br>
 */
  	void Initialize_Canonical(FloatT phi0, FloatT d0, FloatT z0,
  				  FloatT omega, FloatT tanLambda, FloatT B, FloatT * pos);

    /**
     *  Returns momentum of particle at the point of closest approach <br>
     *  to IP <br>
     */
    const FloatT* getMomentum() const { return _momentum; }

    /**
     *  Returns reference point of track <br>
     */
    const FloatT* getReferencePoint() const { return _referencePoint; }

    /**
     *  Returns Phi angle of the momentum vector <br>
     *  at the point of closest approach to IP <br>
     */
    FloatT getPhi0() const;

    /**
     *  Returns signed distance of closest approach <br>
     *  to IP in the R-Phi plane <br>
     */
    FloatT getD0() const { return _d0; }

    /**
     *  Returns z coordinate of the point of closest
     *  approach to IP in the R-Phi plane <br>
     */
    FloatT getZ0() const { return _z0; }

    /**
     *  Returns signed curvature of the track <br>
     */
    FloatT getOmega() const { return _omega; }

    /**
     *  Returns tangent of dip angle of the track <br>
     */
    FloatT getTanLambda() const { return _tanLambda; }

    /**
     *  Returns transverse momentum of the track <br>
     */
    FloatT getPXY() const { return _pxy; }


    /**
     *  Returns x coordinate of circumference
     */
    FloatT getXC() const { return _xCentre; }

    /**
     *  Returns y coordinate of circumference
     */
    FloatT getYC() const { return _yCentre; }

    /**
     *  Returns x coordinate of point of closest approach
     */
    FloatT getXPCA() const { return _xAtPCA; }

    /**
     *  Returns y coordinate of point of closest approach
     */
    FloatT getYPCA() const { return _yAtPCA; }


     /**
     *  Returns radius of circumference
     */
    FloatT getRadius() const { return _radius; }


    /**
     *  Returns helix intersection point with the plane <br>
     *  parallel to z axis. Plane is defined by two coordinates <br>
     *  of the point belonging to the plane (x0,y0) and normal <br>
     *  vector (ax,ay).  ref[3] is the reference point of the helix. <br>
     *  point[3] - returned vector holding the coordinates of <br>
     *  intersection point <br>
     */
    FloatT getPointInXY(FloatT x0, FloatT y0, FloatT ax, FloatT ay,
			      FloatT * ref , FloatT * point) const;

    /**
     *  Returns helix intersection point with the plane <br>
     *  perpendicular to z axis. Plane is defined by z coordinate <br>
     *  of the plane. ref[3] is the reference point of the helix. <br>
     *  point[3] - returned vector holding the coordinates of <br>
     *  intersection point <br>
     */
    FloatT getPointInZ(FloatT zLine, FloatT * ref, FloatT * point) const;

    /**
     * Return distance of the closest approach of the helix to <br>
     * arbitrary 3D point in space. xPoint[3] - coordinates of <br>
     * space point. Distance[3] - vector of distances of helix to <br>
     * a given point in various projections : <br>
     * Distance[0] - distance in R-Phi plane <br>
     * Distance[1] - distance along Z axis <br>
     * Distance[2] - 3D distance <br>
     */
    FloatT getDistanceToPoint(FloatT const* xPoint, FloatT * Distance) const;

    /**
     * Return distance of the closest approach of the helix to <br>
     * arbitrary 3D point in space. xPoint[3] - coordinates of <br>
     * space point. distCut - limit on the distance between helix <br>
     * and the point to reduce calculation time <br>
     * If R-Phi is found to be greater than distCut, rPhi distance is returned <br>
     * If the R-Phi distance is not too big, than the exact 3D distance is returned <br>
     * This function can be used, if the exact distance is not always needed <br>
     */
    FloatT getDistanceToPoint(const FloatT* xPoint, FloatT distCut) const;
    FloatT getDistanceToPoint(const std::vector<FloatT>& xPoint, FloatT distCut) const;

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
    FloatT getPointOnCircle(FloatT Radius, FloatT * ref, FloatT * point) const;

    /** Returns distance between two helixes <br>
     * Output : <br>
     * pos[3] - position of the point of closest approach <br>
     * mom[3] - momentum of V0 <br>
     */
    FloatT getDistanceToHelix(HelixClassT * helix, FloatT * pos, FloatT * mom) const;

    /**
     * Set Edges of helix
     */
    void setHelixEdges(FloatT * xStart, FloatT * xEnd);

    /**
     * Returns starting point of helix
     */
    const FloatT* getStartingPoint() const {return _xStart;}

    /**
     * Returns endpoint of helix
     */
    const FloatT* getEndPoint() const {return _xEnd;}

    /**
     * Returns BZ for the second parameterization
     */
    FloatT getBz() const { return _bZ; }

    /**
     * Returns Phi for the second parameterization
     */
    FloatT getPhiZ() const { return _phiZ; }

    /**
     * Returns extrapolated momentum
     */
    void getExtrapolatedMomentum(FloatT * pos, FloatT * momentum) const;

    /**
     * Returns charge
     */
    FloatT getCharge() const { return _charge; }

 private:
    FloatT _momentum[3]; // momentum @ ref point
    FloatT _referencePoint[3]; // coordinates @ ref point
    FloatT _phi0=0.0; // phi0 in canonical parameterization
    FloatT _d0=0.0;   // d0 in canonical parameterisation
    FloatT _z0=0.0;   // z0 in canonical parameterisation
    FloatT _omega=0.0; // signed curvuture in canonical parameterisation
    FloatT _tanLambda=0.0; // TanLambda
    FloatT _pxy=0.0; // Transverse momentum
    FloatT _charge=0.0; // Particle Charge
    FloatT _bField=0.0; // Magnetic field (assumed to point to Z>0)
    FloatT _radius=0.0; // radius of circle in XY plane
    FloatT _xCentre=0.0; // X of circle centre
    FloatT _yCentre=0.0; // Y of circle centre
    FloatT _phiRefPoint=0.0; // Phi w.r.t. (X0,Y0) of circle @ ref point
    FloatT _phiAtPCA=0.0; // Phi w.r.t. (X0,Y0) of circle @ PCA
    FloatT _xAtPCA=0.0; // X @ PCA
    FloatT _yAtPCA=0.0; // Y @ PCA
    FloatT _pxAtPCA=0.0; // PX @ PCA
    FloatT _pyAtPCA=0.0; // PY @ PCA
    FloatT _phiMomRefPoint=0.0; // Phi of Momentum vector @ ref point
    constexpr static double _const_2pi=2.0*M_PI; // 2*PI
    constexpr static double _const_pi2=0.5*M_PI; // PI/2
    double _FCT=2.99792458E-4; // 2.99792458E-4
    FloatT _xStart[3]; // Starting point of track segment
    FloatT _xEnd[3]; // Ending point of track segment

    FloatT _bZ=0.0;
    FloatT _phiZ=0.0;

};

#include "HelixClassT.ipp"

#endif
