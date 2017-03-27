#ifndef SimpleLine_H
#define SimpleLine_H 1

#include "Trajectory.h"
#include <LCLine3D.h>

/** Simple line trajectory.
 *  @author F.Gaede, T.Kraemer, DESY
 *  @version $Id: SimpleLine.h,v 1.4 2006-10-25 09:21:37 tkraemer Exp $
 */

class SimpleLine : public Trajectory {

protected:
  SimpleLine() {} 

  LCLine3D _line{};

public:

  virtual ~SimpleLine() {} 
 
  /** Construct Line from reference point and direction.
   */
  SimpleLine( LCVector3D ref , LCVector3D direction ) ;
  
  /** Position at path length s - s==0 corresponds to P.C.A to the origin.
   *  @param s      path length
   *  @param errors return argument - not computed if NULL
   */
  virtual LCVector3D getPosition(double s, LCErrorMatrix* errors=0) const ;
  
  /** Direction at path length s, i.e. (dx/ds,dy/ds,dz/ds) 
   *  @param s      path length
   *  @param errors return argument - not computed if NULL
   */
  virtual LCVector3D getDirection(double s,  LCErrorMatrix* errors=0) const ;
  
  /** Full covariance Matrix of x,y,z,px,py,pz   
   *  @param s      path length
   */
  virtual LCErrorMatrix getCovarianceMatrix( double s) const ;
  
  
  /** Pathlength at point on trajectory closest to given position.  
   *  In order to get the distance use for example:  <br>  
   *     LCVector3D pt = t.getPosition( t.getPathAtClosestPoint( p ) ) ; <br>
   *     double d = LCVector3D( pt - p ).mag()  ; <br> 
   */
  virtual double getPathAt(const LCVector3D position ) const ;
  
  
  /*----------------------------------------------------------------------*/
  
  /** Pathlength at closest intersection point with plane - undefined 
   *  if pointExists==false. 
   */
  virtual double getIntersectionWithPlane( LCPlane3D p, bool& pointExists) const  ;
  
  
  /** Pathlength at closest intersection point with cylinder - undefined 
   *  if pointExists==false. 
   * @param cylinder cylinder object to intersect with
   */
  virtual  double getIntersectionWithCylinder(const LCCylinder & cylinder,
                                              bool & pointExists) const ;

}; // class 



// /** Physical trajectory describing a (charged) particle's  path in a B 
//  *  field and material. 
//  *  @author F.Gaede, DESY
//  *  @version $Id: SimpleLine.h,v 1.4 2006-10-25 09:21:37 tkraemer Exp $
//  */

// class PhysicalSimpleLine : public SimpleLine{
  
//   /** Particle's momentum at path length s. Implementations will have to  have knowledge
//    *  about the particle type, B-field and material.
//    */
//   virtual LCLorentzVector get4Momentum( double s ) const ;
// }


#endif /* ifndef SimpleLine_H */
