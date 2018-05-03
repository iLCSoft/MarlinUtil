#ifndef SelectEventsByMCDecayProcessor_H
#define SelectEventsByMCDecayProcessor_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <map>

#include "IMPL/MCParticleImpl.h"
//#include "TH1.h"
//
//#ifdef MARLIN_USE_AIDA //AIDA
//#include <marlin/AIDAProcessor.h>
//#endif

using namespace lcio;
using namespace marlin;


/** SelectEventsByMCDecayProcessor Processor <br>
 *  Looks into MCParticles, selects events by decay mode(s) and, if the specified decay is found, sets Marlin::Processor::setReturnValue flag to 'true'.
 *
 *  Specify the pdg numbers of target particles in the processor parameters.
 *  All pdg numbers are absolute values - this processors doesn't differentiate between positive and negative MC pdg numbers and will not find negative input numbers.
 *  All input numbers equal to 0 will be treated as 'jokers' and will return a 'found'.
 *  Set at least one parent pdg number. If you set both parents to 0 the processor will not run and just setReturnValue(true).
 *
 *  You can select up to 2 decay chains including a parent particle, a vector of children and a vector of grandchildren.
 *  You can select if both decays must exist in one event to be tagged or only one ('both_decays_required' respectively 'and/or').
 *  If you select only one decay and set the other to 0 (or it defaults to 0), the second (0) one will always be found.
 *  In this case the 'and/or'-selection must be put to 'and' to give a meaningful result. This is the default case.
 *
 *  One decay chain is given by one mandatory parent particle (int) and optionally children and grandchildren (both vector<int>).
 *  The processor scans through the MCParticle collection and every time it finds the parent pdg number it checks if the children pdg numbers are among the daughters of the found MCParticle.
 *  If this is the case, it checks if the daughters of every found child match the corresponding grandchildren of that child.
 *
 *  If you have several children, the grandchildren need to be separated for every child.
 *  This is done by putting a '-1' in between groups of grandchildren that belong to different children.
 *  This correspondence is order-sensitive! The grandchildren of the first group are only compared to the daughters of the first child etc.
 *  This can lead to missed decays in cases of identical children that decay differently.
 *  While in reality you don't care about the ordering of those different decays, the processor does. In that case, simply add both decays in the 2 chains and select 'or'.
 *
 *  If there are fewer groups of grandchildren specified than the number of children, the missing groups are filled with 0.
 *  In the current version there are no placeholders for groups of pdg numbers.
 *
 *  Example: Look for K0s-decays to charged pions.
 *    Parent1 310
 *    Children1 211 211
 *
 *  Example: Look for H -> WW -> cx ev_e (one W to charm and s or d, the other one to electron + electron neutrino)
 *    Parent1 25
 *    Children1 24 24
 *    Grandchildren1 4 -1 11
 *    Parent2 25
 *    Children2 24 24
 *    Grandchildren2 11 -1 4
 *    both_decays_required false
 *    You only need to specify one of the decay products of each W in this case, the other follow automatically from the possible W decays.
 *    But you need to specify separately if the first or second found W decays into the charm quark.
 *    If you want to include muon and tauon decays you need to run the processor again with the according parameters.
 *
 *  @param _MCParColName - (string) name of the input MCParticle collection.
 *  default: MCParticle
 *  @param _parent1_PDG - (int) abs(PDG number) of the first parent particle, ignored if 0
 *  default: 0
 *  @param _children1_PDGVec - (vector<int>) abs(PDG numbers) of the children of the first parent particle, ignored if 0
 *  default: {0}
 *  @param _grandchildren1_PDGVec - (vector<int>) abs(PDG numbers) of the grandchildren of the first parent particle, ignored if 0
 *  use -1 as spacer between grandchildren of different children
 *  default: {0}
 *  @param _parent2_PDG - (int) abs(PDG number) of the second parent particle, ignored if 0
 *  default: 0
 *  @param _children2_PDGVec - (vector<int>) abs(PDG numbers) of the children of the second parent particle, ignored if 0
 *  default: {0}
 *  @param _grandchildren2_PDGVec - (vector<int>) abs(PDG numbers) of the grandchildren of the second parent particle, ignored if 0
 *  use -1 as spacer between grandchildren of different children
 *  default: {0}
 *  @param - _both_decays_required - (bool) flag indicating if both decays/decay chains must exist (true) or just one of them (false)
 *  default: true.
 *
 *  @author U. Einhaus, DESY
 *  @version $1$
 */
class SelectEventsByMCDecayProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new SelectEventsByMCDecayProcessor ; }
  
  
  SelectEventsByMCDecayProcessor();

  virtual ~SelectEventsByMCDecayProcessor() = default;

  SelectEventsByMCDecayProcessor(const SelectEventsByMCDecayProcessor&) = delete;
  SelectEventsByMCDecayProcessor& operator=(const SelectEventsByMCDecayProcessor&) = delete;

  
  virtual void init();
  
  virtual void processRunHeader( LCRunHeader* run );
  
  virtual void processEvent( LCEvent * evt );
  
  virtual void check( LCEvent * evt );
  
  virtual void end();

  virtual bool SearchMCParVecForPDG( MCParticleVec MCParVec, std::vector<int> PDGVec );

  virtual bool CheckEmpty( std::vector<int> PDGVec );

  virtual std::vector <std::vector<int> > DecodeGrandchildren( std::vector<int> GrandchildrenVec );


 protected:

  std::string _description = "";
  int _nRun{};
  int _nEvt{};

  std::string _MCParColName{};

  int _parent1_PDG = 0;
  int _parent2_PDG = 0;
  std::vector<int> _children1_PDGVec = {};
  std::vector<int> _children2_PDGVec = {};
  std::vector<int> _grandchildren1_PDGVec = {};
  std::vector<int> _grandchildren2_PDGVec = {};
  bool _both_decays_required = false;

  bool _skip_processor = false;
  int _nSel{};

} ;

#endif



