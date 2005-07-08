#ifndef HELIXAR_H
#define HELIXAR_H 1

class HelixClass {
 public:
    HelixClass();
    ~HelixClass();
    void Initialize_VP(float * pos, float * mom, float q, float B);
    void Initialize_Canonical(float phi0, float d0, float z0, float omega, float tanlambda, float B);
    const float * getMomentum();
    const float * getReferencePoint();
    float getPhi0();
    float getD0();
    float getZ0();
    float getOmega();
    float getTanLambda();
    float getPXY();
    float getPointInXY(float x0, float y0, float ax, float ay, 
			      float * ref , float * point);
    float getPointInZ(float zLine, float * ref, float * point);



 private:    
    float _momentum[3]; // momentum @ ref point 
    float _referencePoint[3]; // coordinates @ ref point
    float _phi0; // phi0 in canonical parameterization 
    float _d0;   // d0 in canonical parameterisation
    float _z0;   // z0 in canonical parameterisation
    float _omega; // signed curvuture in canonical parameterisation
    float _tanLambda; // TanLambda 
    float _pxy; // Transverse momentum
    float _charge; // Particle Charge
    float _bField; // Magnetic field
    float _radius; // radius of circle in XY plane
    float _xCentre; // X of circle centre
    float _yCentre; // Y of circle centre
    float _phiRefPoint; // Phi w.r.t. (X,Y) of circle centre @ ref point
    float _phiAtPCA; // Phi @ PCA (Point of closest approach in XY plane)
    float _xAtPCA; // X @ PCA
    float _yAtPCA; // Y @ PCA
    float _pxAtPCA; // PX @ PCA
    float _pyAtPCA; // PY @ PCA
    float _phiMomRefPoint; // Phi of Momentum vector @ ref point
    float _const_pi; // PI
    float _const_2pi; // 2*PI
    float _const_pi2; // PI/2    
    float _FCT; // 2.99792458E-4

};


#endif
