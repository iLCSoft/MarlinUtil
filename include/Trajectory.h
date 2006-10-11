#ifndef Trajectory_H
#define Trajectory_H 1

#include "LCGeometryTypes.h"

/** Abstract trajectory interface describing a geometrical path in 3D space.
 *  @author F.Gaede, DESY
 *  @version $Id: Trajectory.h,v 1.2 2006-10-11 09:01:31 tkraemer Exp $
 */

class Trajectory {
  
public:
  
  /** Position at path length s - s==0 corresponds to P.C.A to the origin.
   *  @param s      path length
   *  @param errors 3x3 matrix, return argument - not computed if NULL
   */
  virtual LCVector3D getPosition(double s, LCErrorMatrix* errors=0) const = 0;
  
  /** Direction at path length s, i.e. (dx/ds,dy/ds,dz/ds) 
   *  @param s      path length
   *  @param errors 3x3 matrix, return argument - not computed if NULL
   */
  virtual LCVector3D getDirection(double s,  LCErrorMatrix* errors=0) const = 0;
  
  /** Full covariance Matrix of x,y,z,px,py,pz   
   *  @param s      path length
   */
  virtual LCErrorMatrix getCovarianceMatrix( double s) const = 0;
  
  
  /** Pathlength at point on trajectory closest to given position.  
   *  In order to get the distance use for example:  <br>  
   *     LCVector3D pt = t.getPosition( t.getPathAtClosestPoint( p ) ) ; <br>
   *     double d = LCVector3D( pt - p ).mag()  ; <br> 
   */
  
  virtual double getPathAt(const LCVector3D position ) const = 0;
  
  
  /*----------------------------------------------------------------------*/
  
  /** Pathlength at closest intersection point with plane - undefined 
   *  if pointExists==false. 
   */
  virtual double getIntersectionWithPlane( LCPlane3D p, bool& pointExists) const = 0 ;
  
  
  /** Pathlength at closest intersection point with cylinder - undefined 
   *  if pointExists==false. 
   * @param center center point of cylinder
   * @param axis   axis, mag( axis) is cylinder half length
   * @param radius 
   */

  virtual  double getIntersectionWithCylinder(LCVector3D center, 
					      LCVector3D axis, 
					      double radius, 
					      bool & pointExists) const = 0;
  
  
}; // class 



/** Physical trajectory describing a (charged) particle's  path in a B 
 *  field and material. 
 *  @author F.Gaede, DESY
 *  @version $Id: Trajectory.h,v 1.2 2006-10-11 09:01:31 tkraemer Exp $
 */

class PhysicalTrajectory : public Trajectory{
  
  /** Particle's momentum at path length s. Implementations will have to  have knowledge
   *  about the particle type, B-field and material.
   */
  virtual LCLorentzVector get4Momentum( double s ) const = 0;
};


#endif /* ifndef Trajectory_H */
