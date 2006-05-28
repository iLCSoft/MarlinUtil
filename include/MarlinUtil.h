#ifndef MarlinUtil_h
#define MarlinUtil_h 1

#include <iostream>
#include <fstream>

#include <cmath>
#include <string>
#include <vector>

#include <CLHEP/HepPDT/TableBuilder.hh>
#include <CLHEP/HepPDT/ParticleDataTable.hh>
#include <CLHEP/HepPDT/TempParticleData.hh>


#include <lcio.h>
//#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>

#include "HelixClass.h"

using namespace lcio;


class MarlinUtil {


 public:
  static void printTrack(Track* track, double BField);
  static void printCluster(Cluster* cluster);
  static void printRecoParticle(ReconstructedParticle* recoParticle, double BField);
  static void printMCParticle(MCParticle* MCP, bool printDaughters = false);

  /** Function to get the accumulated sum of the energy per event and the number of particles within different categories at IP. The return values are given in the array accumulatedEnergies of size 21 with the following content. Only MC particles with generator status 1 are considered:
   *
   *  accumulatedEnergies[0]  : energy of MC particles which is possible to measure in calorimeter (This means muons are presented as 1.6 GeV each)
   *  accumulatedEnergies[1]  : energy lost in tube
   *  accumulatedEnergies[2]  : energy of neutrinos
   *  accumulatedEnergies[3]  : energy of electrons 
   *  accumulatedEnergies[4]  : energy of muons
   *  accumulatedEnergies[5]  : energy of photons 
   *  accumulatedEnergies[6]  : energy of pi0s
   *  accumulatedEnergies[7]  : energy of long lived neutral hadrons
   *  accumulatedEnergies[8]  : energy of short lived neutral hadrons
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
   *  accumulatedEnergies[20] : energy of MC particles which is possible to measure (real sum [see 0])
   */   
  static void getMC_Balance(LCEvent* evt, double* accumulatedEnergies);
  

  static std::string getMCName(int PDGCode);
  static int getPDGCode(std::string name);
  static bool DecayChainInTree(std::vector<int> DecayChannel, LCEvent* evt);

};

#endif
