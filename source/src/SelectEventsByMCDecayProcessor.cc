#include "SelectEventsByMCDecayProcessor.h"
#include "marlin/Global.h"
#include "marlin/Exceptions.h"
#include "EVENT/LCCollection.h"
#include "IMPL/LCCollectionVec.h"

//#include "TCanvas.h"
//#include "TImage.h"
//#include "TStyle.h"

#include <math.h>
#include <iostream>

using namespace lcio;
using namespace marlin;


SelectEventsByMCDecayProcessor aSelectEventsByMCDecayProcessor;


SelectEventsByMCDecayProcessor::SelectEventsByMCDecayProcessor() : Processor("SelectEventsByMCDecayProcessor") {
  
  _description = "SelectEventsByMCDecayProcessor: Looks into MCParticles, selects events by decay mode(s) and sets Marlin::Processor::setReturnValue flag." ;


  registerInputCollection(LCIO::MCPARTICLE,
  			   "MCParticleCollection",
  			   "MC Particles",
  			   _MCParColName,
  			   std::string("MCParticle"));


  registerProcessorParameter("Parent1",
			     "Abs(PDG number) of first parent particle, ignored if 0, default: 0.",
			     _parent1_PDG,
			     int(0));

  registerProcessorParameter("Children1",
			     "Abs(PDG numbers) of children of first parent particle, ignored if 0, default: {0}.",
			     _children1_PDGVec,
			     std::vector<int>({0}));

  registerProcessorParameter("Grandchildren1",
			     "Abs(PDG numbers) of grandchildren of first parent particle, ignored if 0, use -1 as spacer between grandchildren of different children, default: {}.",
				 _grandchildren1_PDGVec,
			     std::vector<int>({}));

  registerProcessorParameter("Parent2",
			     "Abs(PDG number) of second parent particle, ignored if 0, default: 0.",
				 _parent2_PDG,
			     int(0));

  registerProcessorParameter("Children2",
			     "Abs(PDG numbers) of children of first parent particle, ignored if 0, default: {0}.",
				 _children2_PDGVec,
			     std::vector<int>({0}));

  registerProcessorParameter("Grandchildren2",
			     "Abs(PDG numbers) of grandchildren of first parent particle, ignored if 0, use -1 as spacer between grandchildren of different children, default: {}.",
				 _grandchildren2_PDGVec,
			     std::vector<int>({}));

    registerProcessorParameter("both_decays_required",
                 "Must both decays/decay chains exist (true) or just one of them (false), default: true.",
                 _both_decays_required,
                 bool(true));
}

void SelectEventsByMCDecayProcessor::init() {

  // usually a good idea to
  printParameters();

  _nRun = -1;
  _nEvt = 0;
  _nSel = 0;

  for (unsigned int i=0; i<_children1_PDGVec.size(); i++)
    if (_children1_PDGVec[i]<0) throw ParseException("Children1 have negative entries - only positive allowed!");
  for (unsigned int i=0; i<_children2_PDGVec.size(); i++)
    if (_children2_PDGVec[i]<0) throw ParseException("Children2 have negative entries - only positive allowed!");
  for (unsigned int i=0; i<_grandchildren1_PDGVec.size(); i++)
    if (_grandchildren1_PDGVec[i]<0 && _grandchildren1_PDGVec[i]!=-1) throw ParseException("Grandchildren1 have negative entries - only positive allowed (except for -1)!");
  for (unsigned int i=0; i<_grandchildren2_PDGVec.size(); i++)
    if (_grandchildren2_PDGVec[i]<0 && _grandchildren2_PDGVec[i]!=-1) throw ParseException("Grandchildren2 have negative entries - only positive allowed (except for -1)!");

  if (_parent1_PDG==0 && _parent2_PDG==0){
    streamlog_out(WARNING) << "No parent specified - SelectEventsByMCDecayProcessor skipped!" << std::endl;
    _skip_processor = true;
  }
//  if (_children1_PDGVec.size()==0 && _children2_PDGVec.size()==0){
//    streamlog_out(WARNING) << "No children specified - processor skipped!" << std::endl;
//    _skip_processor = true;
//  }
}

void SelectEventsByMCDecayProcessor::processRunHeader( LCRunHeader* ) {
  _nRun++ ;
} 

void SelectEventsByMCDecayProcessor::processEvent( LCEvent * evt ) {

  //streamlog_out(MESSAGE) << _nEvt << std::endl;
  if (_skip_processor) {setReturnValue(true); return;}

  LCCollection *col_mcpar{};

  try  { col_mcpar  = evt->getCollection( _MCParColName ); }
  catch(DataNotAvailableException &e)
  {
    streamlog_out(MESSAGE) << "Input collection not found - skipping event " << _nEvt << std::endl;
    setReturnValue(false);
    return;
  }
  int n_mcpar  = col_mcpar ->getNumberOfElements();

  bool parent1_found = (_parent1_PDG==0);
  bool parent2_found = (_parent2_PDG==0);
  bool parent1_decay_found = CheckEmpty(_children1_PDGVec);
  bool parent2_decay_found = CheckEmpty(_children2_PDGVec);
  bool parent1_subdecay_found = CheckEmpty(_grandchildren1_PDGVec);
  bool parent2_subdecay_found = CheckEmpty(_grandchildren2_PDGVec);

  bool decay1_total = (parent1_found && parent1_decay_found && parent1_subdecay_found);
  bool decay2_total = (parent2_found && parent2_decay_found && parent2_subdecay_found);
  bool all_total = _both_decays_required ? decay1_total && decay2_total : decay1_total || decay2_total;

  std::vector <std::vector <int> > grandchildren1_PDGVec2 = DecodeGrandchildren(_grandchildren1_PDGVec);
  std::vector <std::vector <int> > grandchildren2_PDGVec2 = DecodeGrandchildren(_grandchildren2_PDGVec);

  for (int i=0; i<n_mcpar; ++i)
  {
    streamlog_out(DEBUG) << "Starting search, i = " << i << std::endl;
    MCParticle* mcpar = dynamic_cast<MCParticle*>(col_mcpar->getElementAt(i));
    int mcpar_PDG = mcpar->getPDG();
    std::vector<int> recorded_PDG1;
    std::vector<int> recorded_PDG2;

    if ( mcpar_PDG==_parent1_PDG && !decay1_total)
    {
      streamlog_out(DEBUG) << "Parent 1 found, looking for decay chain" << std::endl;
      parent1_found = true;
      recorded_PDG1.push_back(mcpar_PDG);
      MCParticleVec children1 = mcpar->getDaughters();

      parent1_decay_found = SearchMCParVecForPDG(children1, _children1_PDGVec);
      streamlog_out(DEBUG) << "SearchMCParVecForPDG 1 left" << std::endl;
      if (parent1_decay_found) for (unsigned int j=0; j<children1.size(); j++) recorded_PDG1.push_back(children1[j]->getPDG());

      if (!parent1_subdecay_found)
      {
        // check grandchildren
        parent1_subdecay_found = true;
        for (unsigned int j=0; j<std::min(children1.size(),grandchildren1_PDGVec2.size()); j++)
        {
          MCParticleVec grandchildren1 = children1[j]->getDaughters();
          bool this_subdecay_found = SearchMCParVecForPDG(grandchildren1,grandchildren1_PDGVec2[j]);
          streamlog_out(DEBUG) << "SearchMCParVecForPDG 1-" << j+1 << " left" << std::endl;
          if (!this_subdecay_found) parent1_subdecay_found = false;
          else for (unsigned int k=0; k<grandchildren1.size(); k++) recorded_PDG1.push_back(grandchildren1[k]->getPDG());
        }
      }
    }
    streamlog_out(DEBUG) << "check1" << std::endl;
    decay1_total = (parent1_found && parent1_decay_found && parent1_subdecay_found);

    // And all again for decay (chain) 2
    if ( mcpar_PDG==_parent2_PDG && !decay2_total)
    {
      streamlog_out(DEBUG) << "Parent 2 found, looking for decay chain" << std::endl;
      parent2_found = true;
      recorded_PDG2.push_back(mcpar_PDG);
      MCParticleVec children2 = mcpar->getDaughters();

      parent2_decay_found = SearchMCParVecForPDG(children2, _children2_PDGVec);
      streamlog_out(DEBUG) << "SearchMCParVecForPDG 2 left" << std::endl;
      if (parent2_decay_found) for (unsigned int j=0; j<children2.size(); j++) recorded_PDG2.push_back(children2[j]->getPDG());

      if (!parent2_subdecay_found)
      {
        parent2_subdecay_found = true;
        for (unsigned int j=0; j<std::min(children2.size(),grandchildren2_PDGVec2.size()); j++)
        {
          MCParticleVec grandchildren2 = children2[j]->getDaughters();
          bool this_subdecay_found = SearchMCParVecForPDG(grandchildren2,grandchildren2_PDGVec2[j]);
          streamlog_out(DEBUG) << "SearchMCParVecForPDG 2-" << j+1 << " left" << std::endl;
          if (!this_subdecay_found) parent2_subdecay_found = false;
          else for (unsigned int k=0; k<grandchildren2.size(); k++) recorded_PDG2.push_back(grandchildren2[k]->getPDG());
        }
      }
    }
    decay2_total = (parent2_found && parent2_decay_found && parent2_subdecay_found);

    all_total = _both_decays_required ? decay1_total && decay2_total : decay1_total || decay2_total;
    streamlog_out(DEBUG) << "Decay chains done, found: " << all_total << std::endl;

    if (all_total)
    {
      streamlog_out(DEBUG4) << "---------------" << std::endl;
      for (unsigned int j=0; j<recorded_PDG1.size(); j++) streamlog_out(DEBUG4) << recorded_PDG1[j] << " ";
      streamlog_out(DEBUG4) << std::endl;
      for (unsigned int j=0; j<recorded_PDG2.size(); j++) streamlog_out(DEBUG4) << recorded_PDG2[j] << " ";
      streamlog_out(DEBUG4) << std::endl;
      break;
    }
  }

  setReturnValue(all_total); // from Marlin::Processor  - sets a value to be used in the steering file
  if (all_total) _nSel++;

  _nEvt++;
}


void SelectEventsByMCDecayProcessor::check( LCEvent * ) { }
  
void SelectEventsByMCDecayProcessor::end()
{
  streamlog_out(MESSAGE) << "Selected " << _nSel << " of " << _nEvt << " Events, that's " << 100.*_nSel/_nEvt << "%." << std::endl;
}

bool SelectEventsByMCDecayProcessor::SearchMCParVecForPDG( MCParticleVec MCParVec, std::vector<int> PDGVec)
{
  streamlog_out(DEBUG) << "SearchMCParVecForPDG startet \n target: ";
  for (unsigned int j=0; j<PDGVec.size(); j++) streamlog_out(DEBUG) << PDGVec[j] << " ";
  streamlog_out(DEBUG) << "\n got:    ";
  for (unsigned int j=0; j<MCParVec.size(); j++) streamlog_out(DEBUG) << MCParVec[j]->getPDG() << " ";
  streamlog_out(DEBUG) << std::endl;

  if (MCParVec.size()<PDGVec.size()) {streamlog_out(DEBUG) << "SearchMCParVecForPDG ended: false" << std::endl; return false;}

  std::vector<int> found_index = {};
  for (unsigned int j=0; j<PDGVec.size(); j++)
  {
    if (PDGVec[j]==0) continue;
    bool child_found = false;

    for (unsigned int k=0; k<MCParVec.size(); k++)
    {
      bool found_this = false;
      for (unsigned int f=0; f<found_index.size(); f++) if (f==k) found_this = true;
      if (found_this) continue;

      if (PDGVec[j]==abs(MCParVec[k]->getPDG()))
      {
        child_found = true;
        found_index.push_back(k);
        break;
      }
    }
    if (!child_found) {streamlog_out(DEBUG) << "SearchMCParVecForPDG ended: false" << std::endl; return false;}
  }
  streamlog_out(DEBUG) << "SearchMCParVecForPDG ended: true" << std::endl;
  return true;
}

bool SelectEventsByMCDecayProcessor::CheckEmpty( std::vector<int> PDGVec )
{
  bool ret = true;
  for (unsigned j=0; j<PDGVec.size(); j++) if (PDGVec[j]!=0 && PDGVec[j]!=-1) ret = false;
  return ret;
}

std::vector <std::vector<int> > SelectEventsByMCDecayProcessor::DecodeGrandchildren( std::vector<int> Grandchildren_PDGVec )
{
  std::vector <std::vector <int> > Grandchildren_PDGVec2;
  std::vector <int> subset;
  for (unsigned int j=0; j<Grandchildren_PDGVec.size(); j++)
  {
    if (Grandchildren_PDGVec[j]==0) continue;
    if (Grandchildren_PDGVec[j]==-1)
    {
      Grandchildren_PDGVec2.push_back(subset);
      subset.clear();
      continue;
    }
    subset.push_back(Grandchildren_PDGVec[j]);
  }
  Grandchildren_PDGVec2.push_back(subset);
  return Grandchildren_PDGVec2;
}
