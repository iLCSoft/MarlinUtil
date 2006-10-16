#ifndef LCLine3D_H
#define LCLine3D_H 1

// #include "CLHEP/Vector/ThreeVector.h"
#include <LCGeometryTypes.h>

/** Definition of a LCLine3D describing a geometrical line in 3D space.
 *  @author T.Kraemer, DESY
 *  @version $Id: LCLine3D.h,v 1.2 2006-10-16 15:38:05 tkraemer Exp $
 */

class LCLine3D {

public:

  /**
   * Constructor from a point and a direction.
   * @param point Point is a point of the line
   * @param direction Direction is the directional vector of the line.
   */
  LCLine3D(LCVector3D point, LCVector3D direction) ;

  /** Copy constructor.
   * @param line line is an other LCLine3D.
   */
  LCLine3D(const LCLine3D & line) ;

  /**
   * Destructor. */
  ~LCLine3D() {}

  /**
   * Assignment. */
  LCLine3D & operator=(const LCLine3D & rhs) ;

  /**
   * Position is the point of the line after a distance s. 
   * Is is given with respect to the point of closes approach to the origen of
   * the coordinate system. 
   * @param s s is the path length along the line */
  LCVector3D position(double s = 0) const ;

  /** Direction of the line 
   */
  LCVector3D direction() const ;

  /**
   * Distance of a point to the line. 
   * @param point point is a point in space
   */
  double distance(const LCVector3D & point) const ;  

  /**
   * Projection of a point on to the line. 
   * @param point point is a point in space.
   */
  double projectPoint(const LCVector3D & point) const ;

  /**
   * Test for equality. */
  bool operator==(const LCLine3D & rhs) const ; 

  /**
   * Test for inequality. */
  bool operator!=(const LCLine3D & rhs) const ; 

protected:

  LCVector3D _point;
  LCVector3D _direction;

};

#endif /* ifndef LCLine3D_H */
