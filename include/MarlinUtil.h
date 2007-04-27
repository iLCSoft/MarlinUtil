#ifndef MarlinUtil_h
#define MarlinUtil_h 1

#include <iostream>
#include <fstream>

#include <cmath>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>

#include <CLHEP/HepPDT/TableBuilder.hh>
#include <CLHEP/HepPDT/ParticleDataTable.hh>
#include <CLHEP/HepPDT/TempParticleData.hh>


#include <lcio.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/ReconstructedParticle.h>

#include "HelixClass.h"

#include "csvparser.h"


using namespace lcio;


class MarlinUtil {


 public:

  static void printMCParticle(MCParticle* MCP, bool printDaughters = false);
  static std::string getMCName(int PDGCode);
  static int getPDGCode(std::string name);
  static MCParticleVec getAllMCParents(MCParticle* mcPart );
  static MCParticleVec getAllMCDaughters(MCParticle* mcPart);
  static bool isDaughterOf( MCParticle* daughter, MCParticle* parent );
  static bool DecayChainInTree(std::vector<int> DecayChannel, LCEvent* evt);

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

  static void printTrack(Track* track, double bField=4.0);
  static const double* getMomentum(Track* track, double bField=4.0);
  static const double getAbsMomentum(Track* track, double bField=4.0);
  static void printCluster(Cluster* cluster);
  static void printRecoParticle(ReconstructedParticle* recoParticle, double bField=4.0);
  static int countAllSimTrackerHits(LCEvent* evt,MCParticle* MCP);
  static int countAllSimCalorimeterHits(LCEvent* evt,MCParticle* MCP,double& accumulatedSimCaloEnergy);
  static double getEnergyDepositedInFullCalorimeter(LCEvent* evt);


};




// FIXME: the whole MC particle related stuff should be placed somewhere else (see MarlinUtil as well)


class MCParticleHelper {

 public: 

  MCParticleHelper();
  std::string getMCCharge(int PDGCode);


 private: 

  std::vector<int> _pdgCodesMCParticles;
  std::vector<double> _massMCParticles;
  std::vector<std::string> _nameMCParticles;
  std::vector<std::string> _chargeMCParticles;


};

#endif
