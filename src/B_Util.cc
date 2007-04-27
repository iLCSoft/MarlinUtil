#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <algorithm>

#include <B_Util.h>

using namespace std;

//============================================================================
Point3D Projection(Point3D &p,Point3D &p0,Point3D &p1){
//============================================================================
  Vector3D v(p0,p1);
  Vector3D w(p0,p);
  Point3D pr=v.linear(w.dot(v)/v.dot(v),p0);
  return pr;
}
//============================================================================
bool Point3D::IsBetween(Point3D &p0,Point3D &p1){
//============================================================================
  if((x>max(p0.x,p1.x)) || (x<min(p0.x,p1.x))||
     (y>max(p0.y,p1.y)) || (y<min(p0.y,p1.y))||
     (z>max(p0.z,p1.z)) || (z<min(p0.z,p1.z)))
    return false;
  return true;
}
//============================================================================
double angle_dist(double t1,double p1,double t2,double p2){
//============================================================================
  double angle=sin(t1)*cos(p1)*sin(t2)*cos(p2)+
               sin(t1)*sin(p1)*sin(t2)*sin(p2)+
               cos(t1)*cos(t2);
  if(angle>1.0)
    return 0.0;
  return fabs(acos(angle));
}
//============================================================================
void Order_Max(int n,double *a,int m, double *b, int *ia){
//============================================================================
/*      Utility routine to find M smallest values                    
 from N values and thier addresses in the initial array
 Can be used if M << N in opposite case better use SORTZV (CERN-lib)                
    Author    V.L. Morgunov, V.A. Ivanova  created   04-Sep-2001
    rewitten to C++ at 30-Nov-2006
 INPUT  n -- initial array size                               
        a -- Array with tested numbers                        
        m -- Number of needed min numbers                     
 OUTPUT
        b -- Sorted array                                     
        ia -- addresses of sorted numbers in initial array    
*/
  for(int l = 0; l < m; l++)
    b[l] = 1.e10;

  int key = 0;
  int i,j,k;
  double aa;
  for(i = 0; i < n; i++){
    aa = a[i];
    for(j = 0; j < m; j++){
      if(aa > b[j]){
	for(k=key; k > j; k--){
	   b[k] =  b[k-1];
	  ia[k] = ia[k-1];
	}
	key = key + 1;
	if(key >= m-1) 
	  key = m-1;
	b[j] = aa;
	ia[j] = i;
	break;
      }
    }
  }
} 
