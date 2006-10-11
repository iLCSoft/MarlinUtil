#ifndef LCPlane3D_H
#define LCPlane3D_H 1

// #include "CLHEP/Vector/ThreeVector.h"
#include <LCGeometryTypes.h>

/** Definition of a LCPlane3D describing a geometrical plane in 3D space.
 *  @author T.Kraemer, DESY
 *  @version $Id: LCPlane3D.h,v 1.1 2006-10-11 16:03:24 tkraemer Exp $
 */

class LCPlane3D {

public:

  /**
   * Constructor from four numbers - creates plane a*x+b*y+c*z+d=0. 
   * @param a
   * @param b
   * @param c
   * @param d
   */
  LCPlane3D(double a = 0, double b = 0, double c = 1, double d = 0) ;

  /**
   * Constructor from normal and point. 
   * @param normal vector pointing in the direction of the normal. 
   *               This vector does not have to be normalised.
   * @param point Point on the plane.
   */
  LCPlane3D(LCVector3D normal, LCVector3D point) ;

  /**
   * Constructor from three different points. 
   * @param point1 Point on the plane.
   * @param point2 Point on the plane.
   * @param point3 Point on the plane.
   */
  LCPlane3D(LCVector3D point1, LCVector3D point2, LCVector3D point3) ;

  /** Copy constructor.
   * @param plane plane is an other LCPlane3D.
   */
  LCPlane3D(const LCPlane3D & plane) ;

  /**
   * Destructor. */
  ~LCPlane3D() {}

  /**
   * Assignment. */

  LCPlane3D & operator=(const LCPlane3D & rhs) ;

  /**
   * Returns the a-coefficient in the plane equation: a*x+b*y+c*z+d=0. */
  double a() const ; 

  /**
   * Returns the b-coefficient in the plane equation: a*x+b*y+c*z+d=0. */
  double b() const ; 

  /**
   * Returns the c-coefficient in the plane equation: a*x+b*y+c*z+d=0. */
  double c() const ; 

  /**
   * Returns the free member of the plane equation: a*x+b*y+c*z+d=0. */
  double d() const ;

  /**
   * Returns normal. */
  LCVector3D normal() const ;

  /**
   * Normalization. */
  LCPlane3D & normalize() ;

  /**
   * Distance of a point to the plane. 
   * @param point point is a point in space
   */
  double distance(const LCVector3D & point) const ;  

  /**
   * Projection of a point on to the plane. 
   * @param point point is a point in space.
   */
  LCVector3D projectPoint(const LCVector3D & point) const ;

  /**
   * Projection of the origin onto the plane. */
  LCVector3D projectPoint() const ;

  /**
   * Test for equality. */
  bool operator==(const LCPlane3D & plane) const ; 

  /**
   * Test for inequality. */
  bool operator!=(const LCPlane3D & plane) const ; 

protected:
  double _a, _b, _c, _d;
};

#endif /* ifndef LCPlane3D_H */
