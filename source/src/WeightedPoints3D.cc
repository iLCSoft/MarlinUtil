

/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */

#include "WeightedPoints3D.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>







// ##########################################
// #####                                #####
// #####   Constructor and Destructor   #####
// #####                                #####
// ##########################################

//=============================================================================

// alternative c'tor if cov known - no points available ....
WeightedPoints3D::WeightedPoints3D(const std::vector<double> &cog,
				   const std::vector<double> &cov, const std::vector<double> &mayor_axis_error, int npnt , 
				   double wgtsum  , double wgt2sum , double wgt4sum ):

  _nPoints(npnt),
  _Wgt(),
  _x  (),
  _y  (),
  _z  (),
  _ex (),
  _ey (),
  _ez (),
  _xl (),
  _xt (),
  _x_trans  (),
  _y_trans  (),
  _z_trans  (),
  _current_transformation (0),
  _th_ref (0.0),
  _ph_ref (0.0),
  _ind_first (0),


  _ifNotPointsGiven(1),

  _ifNotCOG(0),
  _ifNotCOGErrors(0),

  _sumWgt(0.0),
  _wgt_sqr_sum (0.0),
  _wgt_4_sum (0.0),
  _radius(0.0),
  _xgr(0.0),
  _ygr(0.0),
  _zgr(0.0),

  _ifNotWidth(1),
  _Width(0.0),
  _ifNotMaxDist(1),
  _MaxDist(0.0),
  _ifNotEigenSolved (1),
  _ifNotEigenToPolarDone (1),
  _ifNotCovErrors(1),
  _ifNotVecErrorsPolarPropagated( 1),
  _ifNotVecErrorsCartesianPropagated( 1),
  _ifNotValErrorsPropagated( 1),
  _ifNot_dVec_dCov( 1),
  _ifNot_fourthmom ( 1),
  _last_evec_to_error_propagate ( 2),

  _nfact1(0.0),
  _nfact2(0.0),
  _nfact3(0.0),
  _nfact4(0.0),
  _r1(0.0),
  _r2(0.0),
  _r3(0.0),
  _vol(0.0),
  _r_ave(0.0),
  _density(0.0),
  _logitudinaleccentricity(0.0),
  _transverseeccentricity(0.0),
  _ifNotfirst_and_last_found (1),
  _r1_forw(0.0),
  _r1_back(0.0),
  _FractionInside(0.0)
{

  //  _ifNotEigensystem ( 1)

  for ( int i=0 ; i < 3 ; i++) {
    _COG[i] = cog[i];
  }
  if ( cov.size() == 6 ) {
    int kkk=0;
    for ( int i=0 ; i < 3 ; i++) {
      for ( int j=0 ; j <=i ; j++) {
        _COGCov[i][j] = cov[kkk]; 
        _COGCov[j][i] = cov[kkk]; kkk++;
      }
    }
  } else if (  cov.size() == 9 ) {
    for ( int i=0 ; i < 3 ; i++) {
      for ( int j=0 ; j < 3 ; j++) {
        _COGCov[i][j] = cov[j+i*3];
      }
    }
  } 
  if ( wgtsum > 0.0 && wgt2sum > 0.0 ) {
    if ( cov.size() == 6 ) {
      for ( int i=0 ; i < 3 ; i++) {
        for ( int j=0 ; j < 3 ; j++) {
          _COGCov[i][j] =  _COGCov[i][j]*wgtsum*wgtsum/wgt2sum;
        }
      }
    }
    if ( _nPoints > 1 ) {
      _nfact3 = pow(1.0*_nPoints,4.)*(wgt2sum*wgt2sum-wgt4sum)/(pow(wgtsum,4.0)*_nPoints*(_nPoints-1));
    }
  }
  _sumWgt = wgtsum ; _wgt_sqr_sum =  wgt2sum; _wgt_4_sum =  wgt4sum;
  
  for ( int k=0 ; k < 3; k++) {
    for ( int i=0 ; i < 2 ; i++) {
      for ( int j=0 ; j < 2 ; j++) {
        _theta_phi_cov[i][j][k] = 0.0 ;
      }
    }
  }
  if ( mayor_axis_error.size() == 3 ) {
    int kkk=0;
    for ( int i=0 ; i < 2 ; i++) {
      for ( int j=0 ; j <=i ; j++) {
        _theta_phi_cov[i][j][2] =  mayor_axis_error[kkk]; 
        _theta_phi_cov[j][i][2] =  mayor_axis_error[kkk]; kkk++;
      }
    }  
  } else if ( mayor_axis_error.size() == 4 ) {
    for ( int i=0 ; i < 2 ; i++) {
      for ( int j=0 ; j < 2 ; j++) {
	_theta_phi_cov[i][j][2] =  mayor_axis_error[j+i*2];
      }
    }
  } else if ( mayor_axis_error.size() == 0 ) {
    _last_evec_to_error_propagate = 3;
  }  
}
//
WeightedPoints3D::WeightedPoints3D(int nhits, double* a, double* x, double* y, double* z):
  _nPoints(nhits),
  _Wgt(nhits, 0.0),
  _x  (nhits, 0.0),
  _y  (nhits, 0.0),
  _z  (nhits, 0.0),
  _ex (nhits, 1.0),
  _ey (nhits, 1.0),
  _ez (nhits, 1.0),
  _xl (nhits, 0.0),
  _xt (nhits, 0.0),
  _x_trans  (nhits, 0.0),
  _y_trans  (nhits, 0.0),
  _z_trans  (nhits, 0.0),
  _current_transformation (0),
  _th_ref (0.0),
  _ph_ref (0.0),
  _ind_first (0),

  _ifNotPointsGiven(0),

  _ifNotCOG     ( 1),
  _ifNotCOGErrors ( 1),
  _sumWgt(0.0),
  _wgt_sqr_sum(0.0),
  _wgt_4_sum(0.0),
  _radius(0.0),
  _xgr(0.0),
  _ygr(0.0),
  _zgr(0.0),

  _ifNotWidth(1),
  _Width(0.0),
  _ifNotMaxDist(1),
  _MaxDist(0.0),
  _ifNotEigenSolved (1),
  _ifNotEigenToPolarDone (1),
  _ifNotCovErrors( 1),
  _ifNotVecErrorsPolarPropagated( 1),
  _ifNotVecErrorsCartesianPropagated( 1),
  _ifNotValErrorsPropagated( 1),
  _ifNot_dVec_dCov( 1),
  _ifNot_fourthmom ( 1),
  //  _ifNotEigensystem ( 1),
  _last_evec_to_error_propagate ( 3),
  _nfact1(0.),
  _nfact2 (0.),
  _nfact3(0.),
  _nfact4(0.),
  _r1(0.0),
  _r2(0.0),
  _r3(0.0),
  _vol(0.0),
  _r_ave(0.0),
  _density(0.0),
  _logitudinaleccentricity(0.0),
  _transverseeccentricity(0.0),
  _ifNotfirst_and_last_found (1),
  _r1_forw(0.0),
  _r1_back(0.0),
  _FractionInside(0.0)
  
{

 
  //if ( _nPoints < 3 ) { throw int (); }


  for ( int i=0 ; i<3 ; i++ ) { _xyz_ref[i] = 0.0 ; }

  for (int i(0); i < nhits; ++i) {
    _Wgt[i] = a[i];
    _x[i] = x[i];
    _y[i] = y[i];
    _z[i] = z[i];
    _x_trans[i] = x[i];
    _y_trans[i] = y[i];
    _z_trans[i] = z[i];
  }

}


//=============================================================================

WeightedPoints3D::~WeightedPoints3D() {

}

//=============================================================================

void WeightedPoints3D::setErrors(double *ex, double *ey, double *ez) {

  if ( _ifNotPointsGiven == 1 ) return;
  for (int i=0; i<_nPoints; ++i) {
    _ex[i] = ex[i];
    _ey[i] = ey[i];
    _ez[i] = ez[i];
  }	

}




// ##########################################
// #####                                #####
// #####        public methods          #####
// #####                                #####
// ##########################################

//=============================================================================

int WeightedPoints3D::getNumberOfPoints() {
  return _nPoints;
}

//=============================================================================

double WeightedPoints3D::getTotalWeight() {
  if (_ifNotCOG == 1) findCOG();
  return _sumWgt;
}

//=============================================================================

double WeightedPoints3D::getTotalSquaredWeight() {
  if (_ifNotCOG == 1) findCOG();
  return _wgt_sqr_sum ; 
}
//=============================================================================

double WeightedPoints3D::getTotalQuarticWeight() {
  if (_ifNotCOG == 1) findCOG();
  return _wgt_4_sum ; 
}

//=============================================================================

double* WeightedPoints3D::getCentreOfGravity() {
  if (_ifNotCOG == 1) findCOG() ;
  return &_COG[0] ;
}

double* WeightedPoints3D::getCentreOfGravityErrors() {
  if (_ifNotCOGErrors == 1) findCOGErrors() ;
  return &_COGCov[0][0] ;
}


//=============================================================================

double* WeightedPoints3D::getEigenVal() {
  if (_ifNotEigenSolved == 1) solveEigenValEq();
  return &_EigenVal[0] ;
}
double* WeightedPoints3D::getEigenValErrors() {
  if (_ifNotValErrorsPropagated == 1) propagateValErrors();
  return &_var_lam[0];
}

//=============================================================================

double* WeightedPoints3D::getEigenVecCartesian() {
  if (_ifNotEigenSolved == 1) solveEigenValEq();
  return &_EigenVec[0][0] ;
}
double* WeightedPoints3D::getEigenVecPolar() {
  if (_ifNotEigenSolved == 1) solveEigenValEq();
  if (_ifNotEigenToPolarDone == 1) {
    for ( int k=0 ; k < 3 ; k++ ) {
      double ev_T_main = sqrt(_EigenVec[0][k]*_EigenVec[0][k] + _EigenVec[1][k]*_EigenVec[1][k]);
      _EigenVecAngle[0][k]=atan2( ev_T_main,_EigenVec[2][k] );
      _EigenVecAngle[1][k]=atan2(_EigenVec[1][k],_EigenVec[0][k]);
    }
    _ifNotEigenToPolarDone = 0;
  }
  return &_EigenVecAngle[0][0] ;
}
double* WeightedPoints3D::getEigenVecCartesianErrors() {
  if (_ifNotVecErrorsCartesianPropagated == 1) propagateVecErrorsCartesian();
  return &_xyz_cov[0][0][0];
}
double* WeightedPoints3D::getEigenVecPolarErrors() {
  if (_ifNotVecErrorsPolarPropagated == 1) propagateVecErrorsPolar();
  return &_theta_phi_cov[0][0][0];
}
double WeightedPoints3D::getElipsoid_r1(double cl3d ) {
  if (_ifNotEigenSolved == 1) solveEigenValEq();
  return sqrt(_EigenVal[2])*cl3d;
}
double WeightedPoints3D::getElipsoid_r2(double cl3d) {
  if (_ifNotEigenSolved == 1) solveEigenValEq();
  return sqrt(_EigenVal[1])*cl3d;
}
double WeightedPoints3D::getElipsoid_r3(double cl3d) {
  if (_ifNotEigenSolved == 1) solveEigenValEq();
  return sqrt(_EigenVal[0])*cl3d;
}
double WeightedPoints3D::getElipsoid_r1Error(double cl3d) {
  if (_ifNotEigenSolved == 1) solveEigenValEq();
  if (_ifNotValErrorsPropagated == 1) propagateValErrors();
  return sqrt(_var_lam[2])/(2.0*sqrt(_EigenVal[2]))*cl3d;
}
double WeightedPoints3D::getElipsoid_r2Error(double cl3d) {
  if (_ifNotEigenSolved == 1) solveEigenValEq();
  if (_ifNotValErrorsPropagated == 1) propagateValErrors();
  return sqrt(_var_lam[1])/(2.0*sqrt(_EigenVal[1]))*cl3d;
}
double WeightedPoints3D::getElipsoid_r3Error(double cl3d) {
  if (_ifNotEigenSolved == 1) solveEigenValEq();
  if (_ifNotValErrorsPropagated == 1) propagateValErrors();
  return sqrt(_var_lam[0])/(2.0*sqrt(_EigenVal[0]))*cl3d;
}
double WeightedPoints3D::getElipsoid_vol(double cl3d) {
  if (_ifNotEigenSolved == 1) solveEigenValEq();
  _vol = 4.*M_PI* 
         sqrt(_EigenVal[2] *  _EigenVal[1] *  _EigenVal[0])* 
         pow(cl3d/2.0,3.) /3.;  
  return  _vol ;
}
double WeightedPoints3D::getElipsoid_r_ave(double cl3d) {
  _r_ave = pow(getElipsoid_vol(cl3d),1./3.);
  return   _r_ave ;
}
double WeightedPoints3D::getElipsoid_density(double cl3d) {
  if (_ifNotEigenSolved == 1) solveEigenValEq();
  _density = _sumWgt/getElipsoid_vol(cl3d);
  return  _density ;
}
double WeightedPoints3D::getLongitudinalElipsis_eccentricity(){
  if (_ifNotEigenSolved == 1) solveEigenValEq();
  _logitudinaleccentricity = sqrt(1.0-pow(((sqrt(_EigenVal[0])+sqrt(_EigenVal[1]))/2.),2)  /_EigenVal[2]);
  return   _logitudinaleccentricity ;
}
double WeightedPoints3D::getTransverseElipsis_eccentricity(){
  if (_ifNotEigenSolved == 1) solveEigenValEq();
  _transverseeccentricity = sqrt(1.0-_EigenVal[0]/_EigenVal[1]);
  return   _transverseeccentricity ;
}
double WeightedPoints3D::getElipsoid_r_forw() { 
  if ( _ifNotfirst_and_last_found ) { findFirstAndLast() ;}
  return _r1_forw; 
}
double WeightedPoints3D::getElipsoid_r_back() { 
  if ( _ifNotfirst_and_last_found ) { findFirstAndLast() ;}
  return _r1_back; 
}

//=============================================================================

double WeightedPoints3D::getWidth() {
  if (_ifNotWidth == 1) findWidth();
  return _Width;
}

double WeightedPoints3D::getMaxDist() {
  if (_ifNotMaxDist == 1) findMaxDist();
  return _MaxDist;
}

double WeightedPoints3D::getElipsoid_FractionInside(double cl3d) {
  if ( _ifNotPointsGiven == 0 ) { findElipsoid_FractionInside(cl3d); }
  return _FractionInside ;
}

// ##########################################
// #####                                #####
// #####        private methods         #####
// #####                                #####
// ##########################################

//=============================================================================


//=============================================================================

void WeightedPoints3D::findCOG() {

  //  _sumWgt = sum(_Wgt)
  //  forall (j=1:3) _COG(j)   = sum(xyz(j,:)*_Wgt)/_sumWgt 

  if ( _ifNotPointsGiven == 1 ) return;
  _wgt_sqr_sum = 0.0;
  _wgt_sqr_sum = 0.0;
  _wgt_4_sum = 0.0;
  _sumWgt = 0. ;
  for (int i(0); i < 3; ++i) {
    _COG[i] = 0.0 ;
  }
  for (int i(0); i < _nPoints; ++i) {
    _sumWgt+=_Wgt[i] ;
    _wgt_sqr_sum+=_Wgt[i]*_Wgt[i] ;
    _wgt_4_sum+=_Wgt[i]*_Wgt[i]*_Wgt[i]*_Wgt[i] ;
    _COG[0]+=_Wgt[i]*_x[i] ;
    _COG[1]+=_Wgt[i]*_y[i] ;
    _COG[2]+=_Wgt[i]*_z[i] ;
  }
  for (int i(0); i < 3; ++i) {
    _COG[i]/=_sumWgt ;
  }
  _xgr = _COG[0];
  _ygr = _COG[1];
  _zgr = _COG[2];
  _ifNotCOG = 0;
} // findCOG

//=============================================================================

void WeightedPoints3D::findCOGErrors() {
  
  //  nWgt = Wgt)*_nPoints/sum(Wgt)
  //  _nfact1=(_nPoints    - sum(nWgt**2) / _nPoints )/(_nPoints-1.) ;
  //  forall (i=1:3, j=1:3) _COGCov(i,j) = sum((xyz(i,:)-_COG(i))*(xyz(j,:)-_COG(j))*nWgt)/(_nfact1*(_nPoints-1.))
  
 for ( int j = 0 ; j < 3 ; j++ ) {
   for ( int k = 0 ; k < 3 ; k++ ) {
     _COGCov[j][k] =  0.0;
   }
  } 
  if ( _ifNotPointsGiven == 1 ) return;


  double s[3][3][3][3];

  if (_ifNotCOG == 1) findCOG();

  for ( int j = 0 ; j < 3 ; j++ ) {
    for ( int k = 0 ; k < 3 ; k++ ) {
      _COGCov[j][k] = 0;
      for ( int l = 0 ; l < 3 ; l++ ) {
        for ( int m = 0 ; m < 3 ; m++ ) {
          s[j][k][l][m] = 0;
          _fourth_mom[j][k][l][m] = 0;
        }
      }
    }
  } 

  if ( _nPoints < 2 ) return;

  double nwgt_sqr_sum = 0.0;
  double wgt_cube_sum = 0.0;
  double wgt_4_sum = 0.0;
  double pos[3];
  std::vector<double> nWgt;
  for (int i(0); i < _nPoints; ++i) {
    nWgt.push_back(_Wgt[i]* _nPoints/_sumWgt);
  }
  for (int i(0); i < _nPoints; ++i) {
    nwgt_sqr_sum+=nWgt[i]*nWgt[i];
    wgt_cube_sum+=nWgt[i]*nWgt[i]*nWgt[i];
    wgt_4_sum+=nWgt[i]*nWgt[i]*nWgt[i]*nWgt[i];
    // pos = { _x[i] - _COG[0], _y[i] - _COG[1],  _z[i] - _COG[2] };
    pos[0] = _x[i] - _COG[0];
    pos[1] = _y[i] - _COG[1];
    pos[2] = _z[i] - _COG[2];
    for ( int j = 0 ; j < 3 ; j++ ) {
      for ( int k = 0 ; k < 3 ; k++ ) {
        _COGCov[j][k] += double(nWgt[i]*pos[j]*pos[k]);
        for ( int l = 0 ; l < 3 ; l++ ) {
          for ( int m = 0 ; m < 3 ; m++ ) {
            s[j][k][l][m] += nWgt[i]*pos[j]*pos[k]*pos[l]*pos[m];
          }
        }
      }
    } 
  }
  double rPoints=1.0*_nPoints;
  _nfact1=(rPoints   - nwgt_sqr_sum/ rPoints )/(rPoints-1.) ;
  
   for ( int j = 0 ; j < 3 ; j++ ) {
     for ( int k = 0 ; k < 3 ; k++ ) {
       _COGCov[j][k] =  _COGCov[j][k]/(_nfact1*(rPoints-1.)) ;
     }
   } 

  //  fact5=(2*( sum(nWgt)**2)*(sum(nWgt**2))-2*( sum(nWgt)*(sum(nWgt**3)) - 3*(sum(nWgt**2)**2)+3*sum(nWgt**4))/(2.0* sum(nWgt)**3)
  //  fact6=(( sum(nWgt)**4-4*( sum(nWgt)**2)*(sum(nWgt**2))+6*( sum(nWgt)*(sum(nWgt**3))-3*sum(nWgt**4))/( sum(nWgt)**3))/(_nPoints-4)
  //  forall (i=1:3, j=1:3, m=1:3, l=1:3) &
  //    fourth_mom(i,j,m,l) =  (sum((xyz(i,:)-_COG(i))*(xyz(j,:)-_COG(j))*(xyz(m,:)-_COG(m))*(xyz(l,:)-_COG(l))*nWgt) - &
  //                           (2.0*fact5)* ( _COGCov(i,j)* _COGCov(m,l) +  &
  //                           _COGCov(i,m)* _COGCov(j,l) +  _COGCov(i,l)* _COGCov(j,m)))/(fact6*(nPoints-4.))
  //
  //    _nfact2= sum(nWgt)/_nPoints 
  //    _nfact3=(sum(nWgt**2)**2-sum(nWgt**4))/(_nPoints*(_nPoints-1))
  //    _nfact4=(2*sum(nWgt**4) + (sum(nWgt)**2)*sum(nWgt**2) - sum(nWgt**2)**2 - 2* ( sum(nWgt)*sum(nWgt**3)) / &
  //                  ((_nPoints*(_nPoints-1)*(_nPoints-2)))

  double fact5 = 0. ;
  double fact6 = 0.;
  _nfact2= rPoints/rPoints  ;
  _nfact3=(nwgt_sqr_sum*nwgt_sqr_sum - wgt_4_sum )/(rPoints*(rPoints-1));

  if ( _nPoints <= 4 ) return ;
  fact5=( 2*rPoints*rPoints*nwgt_sqr_sum - 2*rPoints*wgt_cube_sum - 3*nwgt_sqr_sum*nwgt_sqr_sum + 3*wgt_4_sum )/(2.0*rPoints*rPoints*rPoints);
  fact6=(( rPoints*rPoints*rPoints*rPoints - 4*rPoints*rPoints*nwgt_sqr_sum + 6*rPoints*wgt_cube_sum - 3*wgt_4_sum)/( rPoints*rPoints*rPoints))/(rPoints-4);
  for ( int j = 0 ; j < 3 ; j++ ) {
    for ( int k = 0 ; k < 3 ; k++ ) {
      for ( int l = 0 ; l < 3 ; l++ ) {
        for ( int m = 0 ; m < 3 ; m++ ) {
          _fourth_mom[j][k][l][m] = (s[j][k][l][m] - 2*fact5*( _COGCov[j][k]* _COGCov[m][l] +  
                       _COGCov[j][m]* _COGCov[k][l] +  
                       _COGCov[j][l]* _COGCov[k][m]))/(fact6*(rPoints-4.));
        }
      }
    }
  } 
  _nfact4=(2*wgt_4_sum  + rPoints*rPoints*nwgt_sqr_sum  - nwgt_sqr_sum*nwgt_sqr_sum - 2* rPoints*wgt_cube_sum) / 
    (rPoints*(rPoints-1)*(rPoints-2));

  _ifNotCOGErrors = 0;
  _ifNot_fourthmom = 0;
} // findCOGErrors
//=============================================================================

void WeightedPoints3D::solveEigenValEq() {

  //
  //  CALL EISRS1 (3,3,_COGCov(1,1),_ValCov,EigenVec,error,wsp)
  //  CALL Cross(EigenVec(1,1),EigenVec(1,2),ocrt)
  //  EigenVec(:,3) =  dot_product( ocrt,  EigenVec(:,3)) *  EigenVec(:,3)
  //  IF ( dot_product(_COG,EigenVec(:,3)) < 0 ) EigenVec(:,2:3) = -EigenVec(:,2:3) 
  //  _radius=sqrt(sum(cog**2))
  //  forall (j=1:3) EigenVecAngle(1,j)=atan2(sqrt(sum(EigenVec(1:2,j)**2)),EigenVec(3,j))
  //  EigenVecAngle(2,:) =    phi_v  = atan2(EigenVec(2,:),EigenVeC(1,:))
  //



  if (_ifNotCOG == 1) findCOG();
  if (_ifNotCOGErrors == 1) findCOGErrors();

  for (int i(0); i < 3; i++) {
    _EigenVal[i] = 0.0;
    for (int j(0); j < 3; j++) {
      _EigenVec[i][j] = 0.0;
    }
  }

  if ( _nPoints < 2 ) return ;

  double cov[3][3];
  for ( int j = 0 ; j < 3 ; j++ ) {
    for ( int k = 0 ; k < 3 ; k++ ) {
      cov[j][k]= _COGCov[j][k] ;
    }
  }
  gsl_matrix_view aMatrix = gsl_matrix_view_array((double*)cov,3,3);
  gsl_vector* aVector = gsl_vector_alloc(3);
  gsl_matrix* aEigenVec = gsl_matrix_alloc(3,3);
  gsl_eigen_symmv_workspace* wa = gsl_eigen_symmv_alloc(3);

  gsl_eigen_symmv(&aMatrix.matrix,aVector,aEigenVec,wa);
  gsl_eigen_symmv_free(wa);
  gsl_eigen_symmv_sort(aVector,aEigenVec,GSL_EIGEN_SORT_ABS_ASC);

  for (int i(0); i < 3; i++) {
    _EigenVal[i] = gsl_vector_get(aVector,i);
    for (int j(0); j < 3; j++) {
      _EigenVec[i][j] = gsl_matrix_get(aEigenVec,i,j);
    }
  }
  double cross_0_1[3];
  cross_0_1[0]= _EigenVec[1][0]*_EigenVec[2][1] - _EigenVec[2][0]*_EigenVec[1][1] ; 
  cross_0_1[1]=-_EigenVec[0][0]*_EigenVec[2][1] + _EigenVec[2][0]*_EigenVec[0][1] ; 
  cross_0_1[2]= _EigenVec[0][0]*_EigenVec[1][1] - _EigenVec[1][0]*_EigenVec[0][1] ; 
  double major_dot_cross =  cross_0_1[0]*_EigenVec[0][2]+cross_0_1[1]*_EigenVec[1][2]+cross_0_1[2]*_EigenVec[2][2];
  for ( int i=0 ; i<3 ; i++ ) { _EigenVec[i][2]=major_dot_cross*_EigenVec[i][2];}

  _radius = 0.;

  for (int i(0); i < 3; ++i) { 
    _radius += _COG[i]*_COG[i];
  }
  _radius = sqrt(_radius);

  double dot_major_cog=_COG[0]*_EigenVec[0][2]+_COG[1]*_EigenVec[1][2]+_COG[2]*_EigenVec[2][2];
  if (dot_major_cog < 0. ) {
    for (int i=1; i < 3; ++i){ // from 3, since we want to keep the system right-handed !
      for (int j=0; j < 3; ++j){ // from 3, since we want to keep the system right-handed !
	_EigenVec[j][i] = - _EigenVec[j][i];
      }
    }
  }

  // rotation = transpose of eigen-vectors:

  for ( int i=0 ; i < 3 ; i++ ) {
    for ( int j=0 ; j < 3 ; j++ ) {
      _RotToEigen[i][j] = _EigenVec[j][i];
    }
  }


  _ifNotEigenSolved = 0;

  gsl_vector_free(aVector);
  gsl_matrix_free(aEigenVec);


} // solveEigenValEq

//=============================================================================

void WeightedPoints3D::findCovErrors() {

  // ind = reshape ( [ 1 ,2 ,3, 1, 1, 2 ,&
  //                      1 ,2 ,3, 2, 3, 3 ] , [ 6 , 2 ] )
  //  forall (i=1:6, j=1:6) &
  //    _COGCovCov(i,j)=( (nfact4+(nfact3-nfact4)/(ntot-1))*(forth_moms(ind(i,1),ind(i,2),ind(j,1),ind(j,2)) - & 
  //                                     _COGCov(ind(i,1),ind(i,2))*_COGCov(ind(j,1),ind(j,2))) +  &
  //                                     nfact3*(_COGCov(ind(i,1),ind(j,2))*_COGCov(ind(j,1),ind(i,2)) + & 
  //                                     (_COGCov(ind(i,1),ind(j,1))*_COGCov(ind(i,2),ind(j,2))))/(nPoints-1)) / &
  //                                              (nPoints*((nfact1*nfact2)**2))
    
    for ( int i=0 ; i < 6 ; i++ ) {
      for ( int j=0 ; j < 6 ; j++ ) {
        _COGCovCov[i][j]=0.0;
      }
    }
    if ( _nPoints <= 1 || _nfact3 <= 0. ) return;
    if (_ifNotCOGErrors == 1) findCOGErrors() ;

    int ind[6][2];
    ind[0][0] = 0 ;
    ind[1][0] = 1 ;
    ind[2][0] = 2 ;
    ind[3][0] = 0 ;
    ind[4][0] = 0 ;
    ind[5][0] = 1 ;
    ind[0][1] = 0 ;
    ind[1][1] = 1 ;
    ind[2][1] = 2 ;
    ind[3][1] = 1 ;
    ind[4][1] = 2 ;
    ind[5][1] = 2 ;


    for ( int i=0 ; i < 6 ; i++ ) {
      for ( int j=0 ; j < 6 ; j++ ) {
  
        double vv=_nfact3*(_COGCov[ind[i][0]][ind[j][1]]*_COGCov[ind[j][0]][ind[i][1]] +  
			   _COGCov[ind[i][0]][ind[j][0]]*_COGCov[ind[i][1]][ind[j][1]])/(_nPoints-1);
	if ( _ifNot_fourthmom == 0 && _nPoints > 4) {
          _COGCovCov[i][j]=( (_nfact4+(_nfact3-_nfact4)/(_nPoints-1))*(_fourth_mom[ind[i][0]][ind[i][1]][ind[j][0]][ind[j][1]] -  
				  _COGCov[ind[i][0]][ind[i][1]]*_COGCov[ind[j][0]][ind[j][1]])  
		     +     vv ) / ( _nPoints*((_nfact1*_nfact2)*(_nfact1*_nfact2)));
	  if ( _COGCovCov[i][j] != _COGCovCov[i][j] ) { 
             std::cout << " _COGCovCov " << i << " " << j << " is NaN " << std::endl ; 
	     std::cout << _nPoints << " " << _nfact1 << " " << _nfact2 << " " << _nfact1 << " " << _nfact3 << " " << _nfact4 << std::endl;
	     std::cout << vv << " " << _fourth_mom[ind[i][0]][ind[i][1]][ind[j][0]][ind[j][1]]  << " " <<  _COGCov[ind[i][0]][ind[i][1]] << " " 
                       << _COGCov[ind[j][0]][ind[j][1]] << std::endl;
          }
        } else {
          // approximation, exact if the distribution is a 3D gaussian
          _COGCovCov[i][j]= vv ; 
        }
      }
    }
 _ifNotCovErrors = 0;
    
  }  //findCovErrors


void WeightedPoints3D::find_dVec_dCov(){

  //    theta=_EigenVecAngle(1,3) ; phi=_EigenVecAngle(2,3)
  //    call cross([-sin(phi)*sin(theta) ,  cos(phi)*sin(theta) , 0.]  , ev(:,1)/sin(theta), rXu)
  //    sin_psi=dot_product(rXu,ev(:,3))
  //    cos_psi= dot_product([-sin(phi)*sin(theta) ,  cos(phi)*sin(theta) , 0.]  , ev(:,1))/sin(theta)
  //    u=EigenVec(:,3)
  //    dcu_mat=reshape ( [ u(1) , 0.   , 0.    , &
  //                        0.   , u(2) , 0.    , &
  //                        0.   , 0.   , u(3)  , &
  //                        u(2) , u(1) , 0.    , &
  //                        u(3) , 0.   , u(1)  , &
  //                        0.   , u(3) , u(2)  ], [ 3 , 6 ])
  //    forall (j=1:2,i = 1:6) dpv(j,i)=dot_product(ev(:,j),dcu_mat(:,i))
  //    forall (i = 1:6) &
  //       dang(i,:)=matmul(transpose(reshape (&
  //                 [ dpv(1,i)/(sin(theta)*(lam(3)-lam(1))) , -dpv(2,i)/(sin(theta)*(lam(3)-lam(2))) , &
  //                   dpv(2,i)/((lam(3)-lam(2)))            , dpv(1,i)/((lam(3)-lam(1)))], [2,2]) ), &
  //                   [ cos_psi , sin_psi ] )


  int ind[3];
  for ( int k=0 ; k < 3 ; k++ ) {
    for ( int i=0 ; i<6 ; i++ ) {   
      dang[i][0][k]=0.0;
      dang[i][1][k]=0.0;
    }
  } 
  if ( _nPoints <= 2 ) return ;
  for ( int k=0 ; k < 3 ; k++ ) {
    if ( k == 2 ) { ind[0]=0 ; ind[1]=1 ; ind[2]=2 ; }
    if ( k == 1 ) { ind[0]=2 ; ind[1]=0 ; ind[2]=1 ; }
    if ( k == 0 ) { ind[0]=1 ; ind[1]=2 ; ind[2]=0 ; }
    if ( _EigenVal[ind[2]] <= 0.0 || ( _EigenVal[ind[0]] <= 0.0 &&  _EigenVal[ind[1]]  <= 0.0 ) ) continue ; 
    double theta=_EigenVecAngle[0][ind[2]];
    double phi=_EigenVecAngle[1][ind[2]];
    double cross_0_1[3];
    double r[3] ;
    r[0]=-sin(phi)*sin(theta)/sin(theta);
    r[1]=cos(phi)*sin(theta)/sin(theta);
    r[2]=0.;
    cross_0_1[0]= r[1]*_EigenVec[2][ind[0]]  - r[2]*_EigenVec[1][ind[0]] ;
    cross_0_1[1]=-r[0]*_EigenVec[2][ind[0]]  + r[2]*_EigenVec[0][ind[0]] ;
    cross_0_1[2]= r[0]*_EigenVec[1][ind[0]]  - r[1]*_EigenVec[0][ind[0]] ;
    double sin_psi = cross_0_1[0]*_EigenVec[0][ind[2]]+cross_0_1[1]*_EigenVec[1][ind[2]]+cross_0_1[2]*_EigenVec[2][ind[2]];
    double cos_psi = r[0]*_EigenVec[0][ind[0]] + r[1]*_EigenVec[1][ind[0]] ;
      //sqrt( 1-sin_psi*sin_psi );
    double dcu_mat[3][6];
    for ( int i = 0 ; i <3 ; i++ ) {
      for ( int j = 0 ; j <6 ; j++ ) {
        dcu_mat[i][j] = 0.0 ;
      }
    }
    dcu_mat[0][0]=_EigenVec[0][ind[2]]; 
    dcu_mat[1][1]=_EigenVec[1][ind[2]]; 
    dcu_mat[2][2]=_EigenVec[2][ind[2]]; 
    dcu_mat[0][3]=_EigenVec[1][ind[2]]; 
    dcu_mat[1][3]=_EigenVec[0][ind[2]]; 
    dcu_mat[0][4]=_EigenVec[2][ind[2]]; 
    dcu_mat[2][4]=_EigenVec[0][ind[2]]; 
    dcu_mat[1][5]=_EigenVec[2][ind[2]]; 
    dcu_mat[2][5]=_EigenVec[1][ind[2]]; 

    double dpv[2][6] ;
    for ( int j = 0 ; j <2 ; j++ ) {
      for ( int i = 0 ; i <6 ; i++ ) {
        dpv[j][i]=dcu_mat[0][i]*_EigenVec[0][ind[j]]+dcu_mat[1][i]*_EigenVec[1][ind[j]]+dcu_mat[2][i]*_EigenVec[2][ind[j]];
      }
    }

    double rv[2];
    double rmat[2][2];
    rv[0]= cos_psi ;
    rv[1]= sin_psi;
    // sqrt( 1-cos_psi*cos_psi );

    for ( int i=0 ; i<6 ; i++ ) {   
      rmat[0][0]=dpv[0][i]/(sin(theta)*(_EigenVal[ind[2]]-_EigenVal[ind[0]]));
      rmat[0][1]=-dpv[1][i]/(sin(theta)*(_EigenVal[ind[2]]-_EigenVal[ind[1]]));
      rmat[1][0]=dpv[1][i]/(_EigenVal[ind[2]]-_EigenVal[ind[1]]);
      rmat[1][1]=dpv[0][i]/(_EigenVal[ind[2]]-_EigenVal[ind[0]]);
      dang[i][0][k]=rmat[0][0]*rv[0] + rmat[0][1]*rv[1] ;
      dang[i][1][k]=rmat[1][0]*rv[0] + rmat[1][1]*rv[1] ;
      if ( dang[i][0][k] != dang[i][0][k] || dang[i][1][k] !=  dang[i][1][k] ) { 
        std::cout << " dang terms 0 " << i << " " << k << std::endl ; 
        std::cout << " dang terms 1 " << _nPoints << " " << rmat[0][0] << " " << rmat[0][1] << " " << rmat[1][0] << " " <<  rmat[1][1] << std::endl ; 
        std::cout << " dang terms 2 " << rv[0]      << " " << rv[1]      << " " << sin(theta) << std::endl ; 
        std::cout << " dang terms 3 " <<  dpv[0][i] << " " <<  dpv[1][i] << std::endl ; 
        std::cout << " dang terms 4 " <<  _EigenVal[ind[0]] << " "  <<  _EigenVal[ind[1]] << " " <<  _EigenVal[ind[2]] << " " << std::endl ; 
      }
    }
  }
   _ifNot_dVec_dCov = 0;
} //find_dVec_dCov 

void WeightedPoints3D::propagateVecErrorsPolar(){

  //    theta_phi_cov = matmul(transpose(dang(:,1:2)),matmul(_COGCovCov,dang(:,1:2)))

  for ( int l=0 ; l < _last_evec_to_error_propagate ; l++ ) {
    for ( int i=0 ; i < 2 ; i++ ) {
      for ( int j=0 ; j < 2 ; j++ ) { 
        _theta_phi_cov[i][j][l] = 0.0;
      }
    }
  }
  if (_ifNotEigenSolved == 1) solveEigenValEq();
  if (_ifNot_dVec_dCov == 1) find_dVec_dCov();
  if (_ifNotCovErrors == 1) findCovErrors();

  double vv_dang[6][2];

  for ( int l=0 ; l < _last_evec_to_error_propagate ; l++ ) {
    for ( int i=0 ; i < 6 ; i++ ) {
      for ( int j=0 ; j < 2 ; j++ ) {
        vv_dang[i][j]=0.0;
        for ( int k=0 ; k < 6 ; k++ ) {
          vv_dang[i][j]+=_COGCovCov[i][k]*dang[k][j][l];
	  if ( dang[k][j][l] != dang[k][j][l] ) { std::cout << " dang " << k << " " << j << " " << l << " is NaN " <<  dang[k][j][l] << std::endl ; }
	  if ( j == 0 ) {if ( _COGCovCov[i][k] != _COGCovCov[i][k] ) { std::cout << " _COGCovCov " << i << " " << k << " is NaN " << std::endl ; }} 
        }
	if ( vv_dang[i][j] != vv_dang[i][j] ) { std::cout << " _vv_dang " << i << " " << j << " is NaN " << std::endl ; }
      }
    }
    for ( int i=0 ; i < 2 ; i++ ) {
      for ( int j=0 ; j < 2 ; j++ ) {
        _theta_phi_cov[i][j][l]=0.0;
        for ( int k=0 ; k < 6 ; k++ ) {
          _theta_phi_cov[i][j][l]+=dang[k][i][l]*vv_dang[k][j];
        }
	if ( _theta_phi_cov[i][j][l] != _theta_phi_cov[i][j][l]  ) { std::cout << " _theta_phi_cov " << i << " " << j << " " << l << " is NaN " << std::endl ; }
      }
    }
  }
 _ifNotVecErrorsPolarPropagated = 0;
    
}  //propagateVecErrorsPolar
void WeightedPoints3D::propagateVecErrorsCartesian(){

  //    theta_phi_cov = matmul(transpose(dang(:,1:2)),matmul(_COGCovCov,dang(:,1:2)))

  for ( int l=0 ; l < _last_evec_to_error_propagate ; l++ ) {
    for ( int i=0 ; i < 2 ; i++ ) {
      for ( int j=0 ; j < 2 ; j++ ) { 
        _xyz_cov[i][j][l] = 0.0;
      }
    }
  }
  if (_ifNotVecErrorsPolarPropagated == 1) propagateVecErrorsPolar();
  if (_ifNotEigenSolved == 1) solveEigenValEq();
  if (_ifNotCovErrors == 1) findCovErrors();
  for ( int l=0 ; l < 3 ; l++ ) {

 
     // Dmat=RESHAPE([ -sin(phi)*sin(theta) , -cos(phi)*sin(theta) , 0. ,&
     //              cos(phi)*cos(theta) ,  sin(phi)*cos(theta) , -sin(theta) ] , [3,2])
     //  Cxyz = MatMul( Dmat , MatMul (theta_phi_cov, Transpose(Dmat))) 

    double theta=_EigenVecAngle[0][l];
    double phi=_EigenVecAngle[1][l];
    double Dmat[3][2] ;
    Dmat[0][0] =  -sin(phi)*sin(theta) ; 
    Dmat[1][0] =  -cos(phi)*sin(theta) ;
    Dmat[2][0] =  0.0;
    Dmat[0][1] =  cos(phi)*cos(theta)  ;
    Dmat[1][1] =  sin(phi)*cos(theta) ;
    Dmat[2][1] =  -sin(theta) ;
    double thph_dmat [2][3];
    for ( int i =0 ; i < 3 ; i++ ) {
      for ( int j=0 ; j < 2 ; j++ ) {
        thph_dmat [i][j] = 0. ;
        for (int k=0 ; k < 2 ; k++ ) {
          thph_dmat [i][j] += _theta_phi_cov[i][k][l]*Dmat[j][k];
        }
      }
    }
    for ( int i =0 ; i < 3 ; i++ ) {
      for ( int j=0 ; j < 3 ; j++ ) {
        _xyz_cov[i][j][l] = 0. ;
        for (int k=0 ; k < 2 ; k++ ) {
	  _xyz_cov[i][j][l] += Dmat[i][k]*thph_dmat[k][j] ;
        }
      }
    }
  }
 _ifNotVecErrorsCartesianPropagated = 0;
    
}  //propagateVecErrorsCartesian

void WeightedPoints3D::propagateValErrors(){

    //    d_lam=Reshape([u(1)**2,u(2)**2,u(3)**2,2*u(1)*u(2),2*u(1)*u(3),2*u(2)*u(3)],[6,1])
    //    lam_var_mat = matmul(transpose(d_lam),matmul(vv,d_lam)) ; lam_var=lam_var_mat(1,1)

  for ( int j=0 ; j < 3 ; j++ ) {
    _var_lam[j] = 0.0;
  }
  if (_ifNotEigenSolved == 1) solveEigenValEq();
  if (_ifNotCovErrors == 1) findCovErrors();

  double vv_dlam[6];
  double d_lam[6];
  for ( int j=0 ; j < 3 ; j++ ) {
    for ( int i=0 ; i <3 ; i++ ) {
      d_lam[i]= _EigenVec[i][j]*_EigenVec[i][j];
    }
    d_lam[3]= 2.0*_EigenVec[0][j]*_EigenVec[1][j];
    d_lam[4]= 2.0*_EigenVec[0][j]*_EigenVec[2][j];
    d_lam[5]= 2.0*_EigenVec[1][j]*_EigenVec[2][j];
    for ( int i=0 ; i < 6 ; i++ ) {
      vv_dlam[i]=0.0;
      for ( int k=0 ; k < 6 ; k++ ) {
        vv_dlam[i]+=_COGCovCov[i][k]*d_lam[k];
      }
    }
    _var_lam[j]=0.0;
    for ( int k=0 ; k < 6 ; k++ ) {
      _var_lam[j]+=d_lam[k]*vv_dlam[k];
    }
  }
  _ifNotValErrorsPropagated = 0;
    
}  //propagateValErrors

void WeightedPoints3D::findWidth() {

  // forall (i=1,_nPoints) _Width=findDistance(i)**2 +  _Width
  // _Width=sqrt(_Width/ _sumWgt)

  _Width  = 0.;
  if ( _ifNotPointsGiven == 1 ) return;

  double dist = 0.0;
  if (_ifNotEigenSolved == 1)  solveEigenValEq() ;
  _Width  = 0.0 ;
  for (int i(0); i < _nPoints; ++i) {
    dist = findDistance(i) ;
    _Width+=_Wgt[i]*dist*dist ;
  }
  _Width  = sqrt(_Width / _sumWgt) ;
  _ifNotWidth = 0 ;
} // findWidth

void WeightedPoints3D::findMaxDist() {

  // forall (i=1,_nPoints) _Width=findDistance(i)**2 +  _Width
  // _Width=sqrt(_Width/ _sumWgt)

  _MaxDist= 0.;
  if ( _ifNotPointsGiven == 1 ) return;

  double dist = 0.0;
  if (_ifNotEigenSolved == 1)  solveEigenValEq() ;
  _MaxDist = 0.0 ;
  for (int i(0); i < _nPoints; ++i) {
    dist = findDistance(i) ;
    if ( dist > _MaxDist ) {
      _MaxDist = dist;
    }
  }  
  _ifNotMaxDist = 0 ;
} // findMaxDist

//=============================================================================

double WeightedPoints3D::findDistance(int i) {

  // f = sqrt(sum((xyz(:,i)-_COG)**2)-dot_product((xyz(:,i)-_COG),  _EigenVector(:,3))**2)

  double cx = 0.0;
  double cy = 0.0;
  double cz = 0.0;
  double dx = 0.0;
  double dy = 0.0;
  double dz = 0.0;
  cx = _EigenVec[0][2] ;
  cy = _EigenVec[1][2] ;
  cz = _EigenVec[2][2] ;
  dx = _COG[0] - _x[i] ;
  dy = _COG[1] - _y[i] ;
  dz = _COG[2] - _z[i] ;
  double tx = cy*dz - cz*dy ;
  double ty = cz*dx - cx*dz ;
  double tz = cx*dy - cy*dx ;
  double tt = sqrt(tx*tx+ty*ty+tz*tz) ;
  double ti = sqrt(cx*cx+cy*cy+cz*cz) ;
  double f = tt / ti ;
  return f ;
} // findDistance

double*  WeightedPoints3D::TransformPointToEigenSyst(double* xyz){
  if (_ifNotEigenSolved == 1 )  solveEigenValEq();
  //  double xyz_eigen[3] ;
  for ( int j = 0 ; j <3 ; j++ ) { 
    xyz_eigen [j] = 0.;
    for ( int k = 0 ; k <3 ; k++ ) { 
      xyz_eigen[j]+=_RotToEigen[j][k]*(xyz[k] - _COG[k]);
    }
  
  } 
  return &xyz_eigen[0];
}

void  WeightedPoints3D::TransformToEigenSyst(){
  /* Transform all points - in _(x,y,z)_trans -to the eigen-system, 
   * the origin - in _COG_trans - to (0,0,0), the covariance -  in _COGCov_trans - to 
   * diagonal (with the eigen values along the diagonal.
   *  _xyz_ref, _th_ref, and _ph_ref contains the position and direction
   *  of this sytem wrt the lab-frame,
   */
  if (_ifNotEigenSolved == 1 )  solveEigenValEq();
  if ( _current_transformation == 1 ) return;
  // COG = null-vector:
  for ( int i=0 ; i < 3 ; i++ ) {
    _COG_trans[i] = 0.0;
  }
  // covariance matrix = diag(eigen-values):
  for ( int kkk=0 ; kkk < 3 ; kkk++ ) {
    for ( int lll=0 ; lll<3 ; lll++) {
      _COGCov_trans[lll][kkk] = 0. ;
      if ( kkk == lll ) { _COGCov_trans[lll][kkk]  = _EigenVal[lll] ;}
    }
  }
  // transform all points:
  //double r_mod_trans[3] ;
  // double r_mod[3] ;
  double r[3] ;
  for ( int i = 0 ; i < _nPoints ; i++ ) {
    r[0] = _x[i] ; r[1] = _y[i] ; r[2] = _z[i] ;
    double* r_mod_trans=TransformPointToEigenSyst(r);
    // for ( int j = 0 ; j <3 ; j++ ) { r_mod[j]=r[j] - _COG[j] ; }
    // for ( int j = 0 ; j <3 ; j++ ) { 
    //   r_mod_trans[j] = 0.;
    //   for ( int k = 0 ; k <3 ; k++ ) { 
    //     r_mod_trans[j]+=_RotToEigen[j][k]* r_mod[k];
    //   }
  
    //    } 
    _x_trans[i] = r_mod_trans[0] ;
    _y_trans[i] = r_mod_trans[1] ;
    _z_trans[i] = r_mod_trans[2] ;
  }
  _ind_first = 0;
  _current_transformation = 1;
  _th_ref= _EigenVecAngle[0][2];
  _ph_ref= _EigenVecAngle[1][2];
  for ( int j = 0 ; j <3 ; j++ ) { _xyz_ref[j] = _COG[j] ; }
}
void  WeightedPoints3D::TransformAlongDirection(double tht , double pht , int on_axis ){
  /* transform all points to a system along the direction given by tht and pht. The origin
   * will be the point that is closest to the third eigen-vector, either (on_axis == 1)
   * projected onto the third eigen-vector, or  (on_axis == 0) the point itself.
   */
  if ( _current_transformation != 1 ) TransformToEigenSyst();
  double r2min=100000.;
  int imin = 0 ;
  for ( int i = 0 ; i < _nPoints ; i++ ) {
    if ( _x_trans[0]*_x_trans[0] + _x_trans[1]*_x_trans[1] + _x_trans[2]*_x_trans[2] < r2min ) {
      r2min =  _x_trans[0]*_x_trans[0] + _x_trans[1]*_x_trans[1] + _x_trans[2]*_x_trans[2] ;
      imin = 1 ;
    }
  }
  double xyz[3] ;
  // xyz[0] = _x_trans[imin] + _COG[0] ;
  // xyz[1] = _y_trans[imin] + _COG[1]  ;
  // xyz[2] = _z_trans[imin] + _COG[2]  ;
  xyz[0] = _x[imin] ;
  xyz[1] = _y[imin] ;
  xyz[2] = _z[imin] ;
  TransformAlongDirection(tht , pht , xyz , on_axis) ;
  if ( on_axis != 1 ) _current_transformation = 2 ;
  if ( on_axis == 1 ) _current_transformation = 3 ;
  _th_ref= _EigenVecAngle[0][2];
  _ph_ref= _EigenVecAngle[1][2];
  // for ( int j = 0 ; j <3 ; j++ ) { _xyz_ref[j] = xyz[j] ; }
  _xyz_ref[0] =_x_trans[imin] ;  
  _xyz_ref[1] =_y_trans[imin] ;  
  _xyz_ref[2] =_z_trans[imin] ; 
  _ind_first = imin ;  
}
void  WeightedPoints3D::TransformAlongDirection(double tht , double pht , double* xyz_start, int on_axis  ){
  /* As the previous, except that the reference point is also input, rather than the
   * point closest to the axis
   */
                // transform cluster to syst along track direction, then check transversal pos wrt track
  
                // Fortran:
                //
                // tht = atan(1./stat%tanlambda) 
                // pht = stat%phi
                // erot1 = reshape([cos(tht)*cos(pht) , -sin(pht), sin(tht)*cos(pht)  , &
                //                  cos(tht)*sin(pht) , cos(pht) , sin(tht)*sin(pht)  , &
                //                  -sin(tht)         ,      0.  , cos(tht) ] , [3,3])
                // cog_rot_trans =  matmul(erot,clu%getPosition-xyz_calo)
                // seen_covmat_esys= matmul(erot,matmul(clu%seen_covmat,transpose(erot)) )
  
                // C++ :
  if (_ifNotCOG == 1) findCOG();
  if (_ifNotCOGErrors == 1) findCOGErrors() ;
  double xyz[3];
  if ( on_axis == 1 ) {
    double xp[3];
    for ( int i = 0 ; i <3 ; i++ ) { xp[i] = xyz_start[i] - _EigenVec[i][2];}
    double d=xp[0]*_EigenVec[0][2]+xp[1]*_EigenVec[1][2]+xp[2]*_EigenVec[2][2];
    for ( int i = 0 ; i <3 ; i++ ) { xyz[i]=_COG[i]+d*_EigenVec[i][2] ; }
  } else {
    for ( int i = 0 ; i <3 ; i++ ) { xyz[i]=xyz_start[i];}
  }

  if ( abs( tht - _EigenVecAngle[0][2] ) < 0.0001 &&  abs( pht - _EigenVecAngle[1][2] ) < 0.0001 ) {

    TransformToEigenSyst();

    double* xyz_t=TransformPointToEigenSyst(xyz);
    for ( int kkk=0 ; kkk<3 ; kkk++) {
      _COG_trans[kkk] =  _COG_trans[kkk]-xyz_t[kkk];
    }
    for ( int i = 0 ; i < _nPoints ; i++ ) {
      _x_trans[i] = _x_trans[i]  - xyz_t[0] ; 
      _y_trans[i] = _y_trans[i]  - xyz_t[1] ; 
      _z_trans[i] = _z_trans[i]  - xyz_t[2] ; 
    }
    _ind_first = 0;
    if ( on_axis != 1 ) _current_transformation = 4 ;
    if ( on_axis == 1 ) _current_transformation = 5 ;
    _th_ref= _EigenVecAngle[0][2];
    _ph_ref= _EigenVecAngle[1][2];
    for ( int j = 0 ; j <3 ; j++ ) { _xyz_ref[j] = xyz_start[j] ; }

  } else {
 
    if ( tht < 0.0 ) { tht = 3.14159 + tht ;}
    float erot[3][3] ;
    erot[0][0] = cos(tht)*cos(pht) ;
    erot[1][0] = -sin(pht) ;
    erot[2][0] =  sin(tht)*cos(pht) ;
    erot[0][1] = cos(tht)*sin(pht) ;
    erot[1][1] =  cos(pht) ;
    erot[2][1] = sin(tht)*sin(pht) ;
    erot[0][2] =  -sin(tht)  ;
    erot[1][2] = 0. ;
    erot[2][2] = cos(tht) ;
    //              float seen_covmat_esys[3][3] ;
    float temp[3][3];
    for ( int kkk=0 ; kkk < 3 ; kkk++ ) {
      for ( int lll=0 ; lll<3 ; lll++) {
        temp[lll][kkk] = 0. ;
        for ( int nnn=0 ; nnn<3 ; nnn++) {
          temp[lll][kkk] +=   _COGCov[lll][nnn]*erot[kkk][nnn];
        }
      }
    }
    for ( int kkk=0 ; kkk < 3 ; kkk++ ) {
      for ( int lll=0 ; lll<3 ; lll++) {
        _COGCov_trans[lll][kkk] = 0. ;
        for ( int nnn=0 ; nnn<3 ; nnn++) {
          _COGCov_trans[lll][kkk] +=  erot[lll][nnn]*temp[nnn][kkk];
        }
      }
    }
  
   //  float xyz_rot[3]; 
   // for ( int kkk=0 ; kkk<3 ; kkk++) {
   //    xyz_rot[kkk] = 0.0 ;
   //    for ( int lll=0 ; lll<3 ; lll++) {
   //      xyz_rot[kkk] += (xyz[lll]-_COG[lll])*erot[kkk][lll];
   //    }
   //  }
    
    //              float cog_rot_trans[3] ;
  
    float cog_mod[3] ;
    for ( int kkk=0 ; kkk<3 ; kkk++) {
      cog_mod[kkk] =  _COG[kkk]-xyz[kkk];
    }
    for ( int kkk=0 ; kkk<3 ; kkk++) {
      _COG_trans[kkk] = 0.0 ;
      for ( int lll=0 ; lll<3 ; lll++) {
        _COG_trans[kkk] += cog_mod[lll]*erot[kkk][lll];
      }
    }
  
    double r_mod_trans[3] ;
    double r_mod[3] ;
    double r[3] ;
    for ( int i = 0 ; i < _nPoints ; i++ ) {
      r[0] = _x[i] ; r[1] = _y[i] ; r[2] = _z[i] ;
      for ( int j = 0 ; j <3 ; j++ ) { r_mod[j]=r[j] - xyz[j] ; }
      for ( int j = 0 ; j <3 ; j++ ) { 
        r_mod_trans[j] = 0.;
        for ( int k = 0 ; k <3 ; k++ ) { 
          r_mod_trans[j]+=erot[j][k]* r_mod[k];
        }
    
      } 
      _x_trans[i] = r_mod_trans[0] ;
      _y_trans[i] = r_mod_trans[1] ;
      _z_trans[i] = r_mod_trans[2] ;
    }
    _ind_first = 0;
    _current_transformation = 6 ;
    _th_ref= tht;
    _ph_ref= pht;
    for ( int j = 0 ; j <3 ; j++ ) { _xyz_ref[j] = xyz_start[j] ; }
  }
}
void WeightedPoints3D::findFirstAndLast() {
  if ( _ifNotPointsGiven == 0 ) {
    if (_ifNotEigenSolved == 1) solveEigenValEq();
    double r_hit_max, d_begn, d_last, r_max, proj;
    double d[3] ;
    r_hit_max = -100000.;
    d_begn    =  100000.;
    d_last    = -100000.;
    for (int i(0); i < _nPoints; ++i) {
      d[0] = _x[i] - _COG[0];
      d[1] = _y[i] - _COG[1];
      d[2] = _z[i] - _COG[2];
      r_max = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);;
      if(r_max > r_hit_max) r_hit_max = r_max;
      proj = 0.;
      for ( int j=0 ; j<3 ; j++ ) { proj += d[j]*_EigenVec[j][2] ; }
      if(proj < d_begn)
        d_begn = proj;
      if(proj > d_last)
        d_last = proj;
    }
    _r1_forw = abs(d_last);
    _r1_back = abs(d_begn);
  } else {
    _r1_forw = 0. ;
    _r1_back = 0. ;
  }
  _ifNotfirst_and_last_found = 0;
}
void WeightedPoints3D::findElipsoid_FractionInside(double cl3d) {
  double wgt_inside = 0.0;
  if ( _ifNotPointsGiven == 0 ) {
    TransformToEigenSyst();
    for ( int i = 0 ; i < _nPoints ; i++ ) {
      if ( pow(_x_trans[i]/cl3d,2)/_EigenVal[0] +  
           pow(_y_trans[i]/cl3d,2)/_EigenVal[1] +   
           pow(_z_trans[i]/cl3d,2)/_EigenVal[2] < 1.0 ) {
        wgt_inside+=_Wgt[i];
      } 
    }
  }
  _FractionInside = wgt_inside/_sumWgt;
}

