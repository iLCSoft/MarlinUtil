/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifndef NNClusters_h
#define NNClusters_h 1

/** File with various classes for a generic nearest neighbour type clustering.
 *
 *  @author F.Gaede (DESY)
 *  @version $Id: NNClusters.h,v 1.5 2007-06-05 15:35:49 engels Exp $
 */

#include <list>
#include <vector>

#include "lcio.h"
#include "EVENT/LCObject.h"
#include "EVENT/LCCollection.h"
#include "IMPL/ClusterImpl.h"

#include "ClusterShapes.h"
#include "CLHEP/Vector/ThreeVector.h"

// fix for transition from CLHEP 1.8 to 1.9
namespace CLHEP{}
using namespace CLHEP ;

// forward declarations:

template <class U>
class GenericCluster ;


/** Simple nearest neighbour (NN) clustering algorithm. Users have to provide an input iterator of
 *  GenericHit objects and an output iterator for the clusters found. The predicate has to have 
 *  a method with the following signature: bool mergeHits( GenericHit<T>*, GenericHit<T>*) where
 *  T is the original (LCIO) type of the hit objects. All pairs of hits for which this method returns
 *  'true' will be merged into one output cluster  - all other pairs of hits will be in distinct 
 *  clusters.
 * 
 *  $see GenericCluster
 *  @see GenericHit
 *  @see NNDistance
 * 
 *  @author F.Gaede (DESY)
 *  @version $Id: NNClusters.h,v 1.5 2007-06-05 15:35:49 engels Exp $
 */

template <class In, class Out, class Pred > 
void cluster( In first, In last, Out result, Pred* pred ) {

  typedef typename In::value_type GenericHitPtr ;
  typedef typename Pred::hit_type HitType ;

  typedef std::vector< GenericCluster<HitType >* >  ClusterList ;

  ClusterList tmp ; 
  tmp.reserve( 1024 ) ;
  
//   int i(0),j(0) ;

  while( first != last ) {

//     j=i+1 ;

    for( In other = first+1 ; other != last  ; other ++ ) {
      
//       std::cout << "  in nested loop " << i << " , " << j << std::endl ;

      if( pred->mergeHits( (*first) , (*other) ) ) {
	
 	if( (*first)->second == 0 && (*other)->second == 0 ) {  // no cluster exists
	  
 	  GenericCluster<HitType >* cl = new GenericCluster<HitType >( (*first) ) ;

	  cl->addHit( (*other) ) ;
	  
	  tmp.push_back( cl ) ;

	}
 	else if( (*first)->second != 0 && (*other)->second != 0 ) { // two clusters
	  
// 	  if(  (*first)->second == (*other)->second ) 
// 	    std::cout << " Merging identical clusters !? " << std::endl ;

	  if(  (*first)->second != (*other)->second )  // this is a bug fix for old gcc 3.2 compiler !
	    (*first)->second->mergeClusters( (*other)->second ) ;
	  
	} else {  // one cluster exists
	  
 	  if( (*first)->second != 0 ) {
	    
	    (*first)->second->addHit( (*other)  ) ;
	    
	  } else {                           
	    
	    (*other)->second->addHit( (*first)  ) ;
	  }
	}
	
      } // dCut 
//       ++j ;
    }
//     ++i ;
    ++first ;
  }

  // remove empty clusters 
  //   std::remove_copy_if( tmp.begin() , tmp.end() , result , &empty_list< GenericCluster<HitType > > ) ;

  for( typename ClusterList::iterator i = tmp.begin(); i !=  tmp.end() ; i++ ){

    if( (*i)->size() > 0 ) {

      result++ = *i ;
    } 
    else {  delete *i ; }
  }
}

/** Templated class for generic hit type objects that are to be clustered with
 *  an NN-like clustering algorithm. Holds a pointer to a generalized cluster
 *  object that is templated with the same type. 
 *
 *  @see GenericCluster
 *  @author F.Gaede (DESY)
 *  @version $Id: NNClusters.h,v 1.5 2007-06-05 15:35:49 engels Exp $
 */
template <class T>
class GenericHit : public  std::pair< T*, GenericCluster<T>* >{

  typedef T value_type ;
  typedef std::pair< T*, GenericCluster<T>* > Pair ;

public:

  /** Default c'tor takes a pointer to the original hit type object. The optioal index can be used to 
   *  code nearest neighbour bins, e.g. in z-coordinate to speed up the clustering process.
   */
  GenericHit(T* hit, int index0 = 0 ) : Index0( index0 )   {  
    Pair::first = hit ; 
    Pair::second = 0 ;
  }

  /** C'tor that also takes a pointer to the cluster this hit belongs to - in case seed hits/clusters are used.
   */
  GenericHit(T* hit ,  GenericCluster<T>* cl , int index0 = 0) : Index0( index0 ) {
    Pair::first =  hit ;
    Pair::second = cl ;
  }
  
  /** Index that can be used to code nearest neighbour bins, e.g. in z-coordinate
   *  to speed up the clustering process.
   */
  int Index0 ;
  
protected:
  /** Don't allow default c'tor w/o hit */
  GenericHit() ;
} ;


/** Templated class for generic clusters  of GenericHits that are clustered with
 *  an NN-like clustering algorithm. Effectively this is just a list of hits.
 * 
 *  @see GenericHit
 *  @author F.Gaede (DESY)
 *  @version $Id: NNClusters.h,v 1.5 2007-06-05 15:35:49 engels Exp $
 */
template <class T >
class GenericCluster : public std::list< GenericHit<T> * > {

public :

  /** C'tor that takes the first hit */
  GenericCluster( GenericHit<T>* hit) {
    addHit( hit ) ;
  }
  
  /** Add a hit to this cluster - updates the hit's pointer to cluster */
  void addHit( GenericHit<T>* hit ) {
    
    hit->second = this ;
    this->push_back( hit ) ;

  }

  /** Merges all hits from the other cluster cl into this cluster */
  void mergeClusters( GenericCluster<T>* cl ) {
    
    for( typename GenericCluster<T>::iterator it = cl->begin() ; it !=  cl->end() ; it++ ){
      (*it)->second = this  ;
    }
    this->merge( *cl ) ;
  }
} ;


/** Helper vector of GenericHit<T> taking care of memory management, i.e. deletes all
 *  GenericHit<T> objects when it goes out of scope.
 */
template <class T> 
class GenericHitVec : public std::vector< GenericHit<T>* > {
  typedef std::vector< GenericHit<T>* > Vector ;
public:
  ~GenericHitVec() {
    for( typename GenericHitVec::iterator i = Vector::begin() ; i != Vector::end() ; i++) delete *i ;
  }
};


/** Helper method that copies all hit pointers from an LCIO collection that fullfill the predicate to
 *  a GenericHitVec. The predicate can either be a bool funtion or functor that takes a T*, e.g.
 *  @see EnergyCut
 */
template <class T, class Pred> 
void addToGenericHitVec(GenericHitVec<T>& v, LCCollection* col, Pred pred ){

  for( int i=0 ; i < col->getNumberOfElements() ; i++ ){

    T* hit = dynamic_cast<T*>( col->getElementAt( i) ) ;

    if( pred( hit ) ){

      v.push_back( new GenericHit<T>( hit ) ) ;
    }
  }
} 

/** Same as addToGenericHitVec(GenericHitVec<T>& v, LCCollection* col, Pred pred ) except that an additional 
 *  order function/functor can be given that defines the index of the hit, e.g.
 *  @see ZIndex.
 */
template <class T, class Pred, class Order> 
void addToGenericHitVec(GenericHitVec<T>& v, LCCollection* col, Pred pred , Order order ){

  for( int i=0 ; i < col->getNumberOfElements() ; i++ ){

    T* hit = dynamic_cast<T*>( col->getElementAt( i) ) ;

    if( pred( hit ) ){

      v.push_back( new GenericHit<T>( hit , order(hit) ) ) ;
    }
  }
} 

/** Helper vector of GenericCluster<T> taking care of memory management.
 */
template <class T> 
class GenericClusterVec : public std::list< GenericCluster<T>* > {
  typedef std::list< GenericCluster<T>* > List ;
public:
  ~GenericClusterVec() {
    for( typename GenericClusterVec::iterator i = List::begin() ; i != List::end() ; i++) delete *i ;
  }
};


/** Simple predicate class for NN clustering. Requires 
 *  PosType* HitClass::getPosition(), e.g for CalorimeterHits use: <br>
 *  NNDistance<CalorimeterHit,float> dist( myDistCut ) ; 
 */
template <class HitClass, typename PosType > 
class NNDistance{
public:

  /** Required typedef for cluster algorithm 
   */
  typedef HitClass hit_type ;

  /** C'tor takes merge distance */
  NNDistance(float dCut) : _dCutSquared( dCut*dCut ) , _dCut(dCut)  {} 


  /** Merge condition: true if distance  is less than dCut given in the C'tor.*/ 
  inline bool mergeHits( GenericHit<HitClass>* h0, GenericHit<HitClass>* h1){
    
    if( std::abs( h0->Index0 - h1->Index0 ) > 1 ) return false ;

    const PosType* pos0 =  h0->first->getPosition() ;
    const PosType* pos1 =  h1->first->getPosition() ;
    
    return 
      ( pos0[0] - pos1[0] ) * ( pos0[0] - pos1[0] ) +
      ( pos0[1] - pos1[1] ) * ( pos0[1] - pos1[1] ) +
      ( pos0[2] - pos1[2] ) * ( pos0[2] - pos1[2] )   
      < _dCutSquared ;
  }
  
protected:
  NNDistance() ;
  float _dCutSquared ;
  float _dCut ;
} ;


/** Simple predicate class for applying an energy cut to the objects of type T.
 *  Requires float/double T::getEnergy().
 */
template <class T>
class EnergyCut{
public:
  EnergyCut( double eCut ) : _eCut( eCut ) {}  

  inline bool operator() (T* hit) {  return hit->getEnergy() > _eCut ; }

protected:

  EnergyCut() {} ;
  double _eCut ;
} ;


/** Simple predicate class for computing an index from N bins of the z-coordinate of LCObjects
 *  that have a float/double* getPostion() method.
 */

template <class T, int N>
class ZIndex{
public:
  /** C'tor takes zmin and zmax - NB index can be negative and larger than N */
  ZIndex( float zmin , float zmax ) : _zmin( zmin ), _zmax( zmax ) {}  

  inline int operator() (T* hit) {  

    return (int) std::floor( ( hit->getPosition()[2] - _zmin ) / ( _zmax - _zmin ) * N ) ; 
  }

protected:

  ZIndex() {} ;
  float _zmin ;
  float _zmax ;
} ;


/** Helper class that creates an lcio::Cluster from a generic cluster with hit types that have a 
 *  getPosition() and a getEnergy() method.
 */
template <class T>
struct LCIOCluster{

  inline lcio::Cluster* operator() (GenericCluster<T>* c) {  
    
    ClusterImpl* clu = new ClusterImpl ;
    
    unsigned n = c->size() ;
    unsigned i=0 ;

    float a[n], x[n], y[n], z[n] ;

    for( typename GenericCluster<T>::iterator hi = c->begin(); hi != c->end() ; hi++) {
      
      T* hit = (*hi)->first ;  

      a[i] = hit->getEnergy() ;
      x[i] = hit->getPosition()[0] ;
      y[i] = hit->getPosition()[1] ;
      z[i] = hit->getPosition()[2] ;

      clu->addHit( hit , a[i] ) ;
      
      ++i ;
    }
    
    ClusterShapes cs( n,a,x,y,z) ;

    clu->setEnergy( cs.getTotalAmplitude()  ) ;
    clu->setPosition( cs.getCenterOfGravity() ) ;

    // direction of cluster's PCA
    float* d = cs.getEigenVecInertia() ;

    Hep3Vector v( d[0], d[1], d[2] ) ;
    
    clu->setITheta( v.theta() )  ;
    clu->setIPhi(   v.phi() ) ;
  
    std::vector<float> param(5) ;

//     float* ev = cs.getEigenValInertia() ;
//     param[0] = ev[0] ;
//     param[1] = ev[1] ;
//     param[2] = ev[2] ;

    param[0] = cs.getElipsoid_r1() ;
    param[1] = cs.getElipsoid_r2() ;
    param[2] = cs.getElipsoid_r3() ;
    param[3] = cs.getElipsoid_vol() ;
    param[4] = cs.getWidth() ;

    clu->setShape( param ) ;

    return clu ;

  }

} ;



#endif


