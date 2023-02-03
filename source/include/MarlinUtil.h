#ifndef MarlinUtil_h
#define MarlinUtil_h 1

#include <lcio.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/ReconstructedParticle.h>

#include <string>
#include <vector>
#include <fstream>


namespace MarlinUtil {


  void printMCParticle(lcio::MCParticle* MCP, bool printDaughters = false);
  std::string getMCName(int PDGCode);
  int getPDGCode(std::string name);
  lcio::MCParticleVec getAllMCParents(lcio::MCParticle* mcPart );
  lcio::MCParticleVec getAllMCDaughters(lcio::MCParticle* mcPart);
  bool isDaughterOf( lcio::MCParticle* daughter, lcio::MCParticle* parent );
  bool DecayChainInTree(std::vector<int> DecayChannel, lcio::LCEvent* evt);

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
  void getMC_Balance(lcio::LCEvent* evt, double* accumulatedEnergies);

  void printTrack(lcio::Track* track, double bField=4.0);
  const double* getMomentum(lcio::Track* track, double bField=4.0);
  double getAbsMomentum(lcio::Track* track, double bField=4.0);
  void printCluster(lcio::Cluster* cluster);
  void printRecoParticle(lcio::ReconstructedParticle* recoParticle, double bField=4.0);
  int countAllSimTrackerHits(lcio::LCEvent* evt,lcio::MCParticle* MCP);
  int countAllSimCalorimeterHits(lcio::LCEvent* evt,lcio::MCParticle* MCP,double& accumulatedSimCaloEnergy);
  double getEnergyDepositedInFullCalorimeter(lcio::LCEvent* evt);


  /** Return track weight contribution encoded as trackwgt = (int(wgt)%10000)/1000. by
   * the MarlinReco/Analysis/RecoMCTruthLink/include/RecoMCTruthLinker.h for the
   * ReconstructedParticle-MCParticle (or vise-versa) type of relations.
   * 
   */
  inline float getTrackWeight(float encodedWeight){
      return float( int(encodedWeight) % 10000 ) / 1000.f;
  }


  /** Return cluster weight contribution encoded as clusterwgt = (int(wgt)/10000)/1000. by
   * the MarlinReco/Analysis/RecoMCTruthLink/include/RecoMCTruthLinker.h for the
   * ReconstructedParticle-MCParticle (or vise-versa) type of relations.
   */
  inline float getClusterWeight(float encodedWeight){
      return float( int(encodedWeight) / 10000 ) / 1000.f;
  }


  /** Comparator function to compare weights with track weight encoding trackwgt =
   * (int(wgt)%10000)/1000. set by the MarlinReco/Analysis/RecoMCTruthLink/include/RecoMCTruthLinker.h 
   * for the ReconstructedParticle-MCParticle (or vise-versa) type of relations.
   */
  inline bool compareTrackWeights(float a, float b) {
    return getTrackWeight(a) < getTrackWeight(b);
  }

  /** Comparator function to compare weights with cluster weight encoding clusterwgt = (int(wgt)/10000)/1000.
   * set by the MarlinReco/Analysis/RecoMCTruthLink/include/RecoMCTruthLinker.h 
   * for the ReconstructedParticle-MCParticle (or vise-versa) type of relations.
   */
  inline bool compareClusterWeights(float a, float b) {
    return getClusterWeight(a) < getClusterWeight(b);
  }

}




// FIXME: the whole MC particle related stuff should be placed somewhere else (see MarlinUtil as well)


class MCParticleHelper {

 public: 

  MCParticleHelper();
  std::string getMCCharge(int PDGCode);


 private: 

  std::vector<int> _pdgCodesMCParticles{};
  std::vector<double> _massMCParticles{};
  std::vector<std::string> _nameMCParticles{};
  std::vector<std::string> _chargeMCParticles{};


};

#endif
