#ifndef MarlinUtil_h
#define MarlinUtil_h 1

#include <iostream>
#include <cmath>

#include <lcio.h>
//#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>


using namespace lcio;


class MarlinUtil {


 public:




 
  /** Function to get the accumulated sum of the energy per event and the number of particles within different categories at IP. The return values are given in the array accumulatedEnergies of size 20 with the following content. Only MC particles with generator status 1 are considered:
   *
   *  accumulatedEnergies[0]  : energy of MC particles which is possible to measure
   *  accumulatedEnergies[1]  : energy lost in tube
   *  accumulatedEnergies[2]  : energy of neutrinos
   *  accumulatedEnergies[3]  : energy of electrons 
   *  accumulatedEnergies[4]  : energy of muons
   *  accumulatedEnergies[5]  : energy of photons 
   *  accumulatedEnergies[6]  : energy of pi0s
   *  accumulatedEnergies[7]  : energy of long lived hadrons
   *  accumulatedEnergies[8]  : energy of short lived hadrons
   *  accumulatedEnergies[9]  : energy of charged hadrons
   *  accumulatedEnergies[10] : number of MC particles which are possible to measure
   *  accumulatedEnergies[11] : number of MC particles lost in tube
   *  accumulatedEnergies[12] : number of neutrinos
   *  accumulatedEnergies[13] : number of electrons 
   *  accumulatedEnergies[14] : number of muons
   *  accumulatedEnergies[15] : number of photons
   *  accumulatedEnergies[16] : number of pi0s
   *  accumulatedEnergies[17] : number of long lived hadrons
   *  accumulatedEnergies[18] : number of short lived hadrons
   *  accumulatedEnergies[19] : number of charged hadrons
   */   
  static void getMC_Balance(LCEvent* evt, double* accumulatedEnergies);
  

};

#endif
