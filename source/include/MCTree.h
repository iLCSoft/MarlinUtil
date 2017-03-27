#ifndef MCTree_h
#define MCTree_h

#include <stdio.h>
#include <iostream>

#include <string>
#include <vector>
#include <map>
#include <set>
#include <stack>
#include <sstream>
#include <queue>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/MCParticleImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <iomanip>
#include <lcio.h>

using namespace lcio;

/**
 *    Utility class for MC particles.
 *
 *    @authors P. Krstonosic (DESY), O. Wendt (DESY)
 *    @version $Id: MCTree.h,v 1.2 2005-08-07 16:16:08 gaede Exp $
 *
 */

class MCTree {

 public:

  /**
   * An object of this class is constructed with a pointer to a LCCollection of MC particles
   */
  MCTree(LCCollection* col);

  ~MCTree();


  /**
   * Prints the whole MC Tree on the standard output. The following print options are possible:
   * opt = -1  : print MC Tree only
   * opt =  0  : print status bit field additionally
   * opt =  1  : print four vector of momentum and mass additionally
   * opt =  2  : print energy and mass additionally
   */
   void print(int opt);

 private:
  
  LCCollection* _col;
  MCTree(const MCTree&) = default;
  MCTree& operator=(const MCTree&) = default;

  int  printShortMCInfo(const MCParticle* part, const unsigned int pdg_size, 
			const int options );
  int length_of_int( const int&  to_convert);
  std::string pdg_to_string( const int&  to_convert , const unsigned int&  max_size);
  std::string adjust_position( const unsigned int&  pdg_size, const unsigned int& index);
  std::string adjust_position1(const  unsigned int& n_blanks);

};

#endif
