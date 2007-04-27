#include "CutOnGEANT4Bug.h"

using namespace lcio ;
using namespace marlin ;

CutOnGEANT4Bug aCutOnGEANT4Bug ;



CutOnGEANT4Bug::CutOnGEANT4Bug() : Processor("CutOnGEANT4Bug") {

  _description = "The output condition of this processor is true as long as no track has more than a factor of k more energy deposited in the calorimeter as its energy given by momentum and mass. This should cut out events where GEANT4 produces additional energy depositions. If at least one such a track is found the return value is false. Only tracks with an energy larger than eMin are taken into account.";


  registerProcessorParameter( "eMin",
			      "minimal energy of tracks taken into account (in GeV)",
			      _eMin,
			      (double)5.0);

  registerProcessorParameter( "k",
			      "if the track has more than k times its MC energy deposited the return value is set to false",
			      _k,
			      (double)1.75);

  registerProcessorParameter( "colNameTracks" ,
			      "name of the Track collection" ,
			      _colNameTracks ,
			      std::string("Tracks") ) ;

  registerProcessorParameter( "colNameRelationTrackToMCP" , 
			      "name of the LC Relation collection between Tracks and MC particles"  ,
			      _colNameRelationTrackToMCP,
			      std::string("TrueTrackToMCP") );

  registerProcessorParameter( "colNameRelationCaloHitToSimCaloHit" , 
			      "name of the LC Relation collection between Calorimeterhits and SimCalorimeterhits"  ,
			      _colNameRelationCaloHitToSimCaloHit,
			      std::string("RelationCaloHit") );

  std::vector<float> calibrECAL;
  calibrECAL.push_back(33.0235);
  calibrECAL.push_back(93.5682);
  registerProcessorParameter("calibrCoeffECAL" , 
			     "Calibration coefficients for ECAL" ,
			     _calibrCoeffECAL,
			     calibrECAL);
  
  std::vector<float> calibrHCAL;
  calibrHCAL.push_back(21.19626);
  registerProcessorParameter("calibrCoeffHCAL" , 
			     "Calibration coefficients for HCAL" ,
			     _calibrCoeffHCAL,
			     calibrHCAL);
}


void CutOnGEANT4Bug::init() {
  
  // usually a good idea to 
  // printParameters();
  
  _nRun = 0 ;
  _nEvt = 0 ;


}


void CutOnGEANT4Bug::processRunHeader( LCRunHeader* run) {

  ++_nRun;

}


void CutOnGEANT4Bug::processEvent( LCEvent* evt ) {


  bool invalidTrackFound = false;

  try {

    std::vector< std::string >::const_iterator iter;
    const std::vector< std::string >* ColNames = evt->getCollectionNames();

    for( iter = ColNames->begin() ; iter != ColNames->end() ; iter++) {
      
      LCCollection* col = evt->getCollection( *iter ) ;
      
      if ( (col->getTypeName() == LCIO::TRACK) && (*iter == _colNameTracks) ) {

	int NTracks = col->getNumberOfElements();
	
	for(int j=0; j<NTracks; ++j){
	  
	  Track* track = dynamic_cast<Track*>(col->getElementAt(j));
	  
	  try {
	    
	    LCCollection* LCRcolTracks = evt->getCollection(_colNameRelationTrackToMCP);
	    
	    LCRelationNavigator* navTracks = new LCRelationNavigator(LCRcolTracks);
	    const LCObjectVec& relMCParticlesToTrack = navTracks->getRelatedToObjects(track); 
	    
	    if ( relMCParticlesToTrack.size() > 1 ) std::cout << "Warning: More than one MCParticle related to track." << std::endl;
    
	    MCParticle* mcpOfTrack = 0;

	    // container for the CalorimeterHits which are related to the track with hit and sub-hit energy accuracy, both run in parallel, i.e. they have the same size
	    std::vector< std::pair<CalorimeterHit*,float> > collectedCalorimeterHitsWithEnergies;
	    std::vector< std::pair<CalorimeterHit*,float> > collectedSubCalorimeterHitsWithEnergies;
	    
	    unsigned int index = 0;
	    bool alreadyCollected = false;	
	    
	    float ESumCalorimeterHits = 0.0; // accumulated energy of calorimeter hits where the track contributes
	    float ESumSubCalorimeterHits = 0.0; // accumulated energy of calorimeter hit energies, but only the part which originates from the track,i.e. sub-hit accuracy
	  
	  
	    for(unsigned int i = 0; i < relMCParticlesToTrack.size(); ++i) {
	      
	      mcpOfTrack = dynamic_cast<MCParticle*>(relMCParticlesToTrack.at(i)); 
	      
	      MCParticleVec allMCPsOfTrack = MarlinUtil::getAllMCDaughters(mcpOfTrack);
	      
	      for(unsigned int j = 0; j < allMCPsOfTrack.size(); ++j) {
		
		MCParticle* mcp = allMCPsOfTrack.at(j);
		
		std::vector< std::string >::const_iterator iter;
		const std::vector< std::string >* ColNames = evt->getCollectionNames();
		
		for( iter = ColNames->begin() ; iter != ColNames->end() ; iter++) {
		  
		  LCCollection* col = evt->getCollection( *iter ) ;
		  
		  if ( (col->getTypeName()) == LCIO::SIMCALORIMETERHIT ) {
		    
		    int nHits = col->getNumberOfElements();	 
		    
		    for(int k = 0; k < nHits; ++k) {
		      
		      SimCalorimeterHit* simCaloHit = dynamic_cast<SimCalorimeterHit*>(col->getElementAt(k));
		      
		      // debug
		      // std::cout << "simCaloHit->getNMCContributions() = " << simCaloHit->getNMCContributions() << std::endl;
		      
		      for ( int l = 0; l < simCaloHit->getNMCContributions(); ++l ) {
		      
			MCParticle* mcpOfCalo = simCaloHit->getParticleCont(l);
			
			if ( mcpOfCalo == mcp ) {
			  
			  // float ESimHit = simCaloHit->getEnergy(); // only for debugging
			  float ESimHitContribution = simCaloHit->getEnergyCont(l);
			  
			  LCCollection* LCRcolCalorimeter = evt->getCollection(_colNameRelationCaloHitToSimCaloHit);
			
			  LCRelationNavigator* navCalorimeter = new LCRelationNavigator(LCRcolCalorimeter);
			  const LCObjectVec& relCaloHitsToSimCaloHit = navCalorimeter->getRelatedFromObjects(simCaloHit);
			  
			  // there should only be one CalorimeterHit related to one SimCalorimeterHit since the CalorimeterHit consists (is related to) of several SimCalorimeterHit, 
			  // but the SimCalorimeterHit is only related to one CalorimeterHit (by the ganging in the Calorimeter digitize processor)		
			  if (relCaloHitsToSimCaloHit.size() > 1 ) {
			    
			    std::cout << "Warning: More than one (" << relCaloHitsToSimCaloHit.size() << ") CalorimeterHit related to one SimCalorimeterHit. " << std::endl;
			    
			  }
			  
			  
			  // debug
			  // std::cout << "relCaloHitsToSimCaloHit.size() = " << relCaloHitsToSimCaloHit.size() << std::endl;
			  
			  
			  for ( unsigned int m = 0; m < relCaloHitsToSimCaloHit.size(); ++m ) {
			    
			    CalorimeterHit* caloHit = dynamic_cast<CalorimeterHit*>(relCaloHitsToSimCaloHit.at(m));
			    
			    // calibration of the calorimeter
			    std::vector<float> calibration;
			    for ( unsigned int n = 0; n < _calibrCoeffECAL.size(); ++n ) calibration.push_back(_calibrCoeffECAL.at(n));
			    for ( unsigned int n = 0; n < _calibrCoeffHCAL.size(); ++n ) calibration.push_back(_calibrCoeffHCAL.at(n));
			    
			    int type   = caloHit->getType();
			    float EHit = caloHit->getEnergy();
			    // float EHitCalculated   = ESimHit*calibration.at(type); // only for debugging
			    float EHitContribution = ESimHitContribution*calibration.at(type);
			    
			    
			    // search if caloHit has already been assigned to track	  
			    for (unsigned int n = 0; n < collectedCalorimeterHitsWithEnergies.size(); ++n) {
			      
			      if (collectedCalorimeterHitsWithEnergies.at(n).first == caloHit) {
				
				index = n;
				alreadyCollected = true;
				break;
				
			      }
			      else {
				
				index = 0;
				alreadyCollected = false;
				
			      }
			      
			    
			    }
			    
			    
			    if ( !alreadyCollected ) { 		    
			      
			      std::pair<CalorimeterHit*,float> calorimeterHitWithEnergy(caloHit,EHit);
			      collectedCalorimeterHitsWithEnergies.push_back(calorimeterHitWithEnergy);
			      
			      std::pair<CalorimeterHit*,float> calorimeterSubHitWithEnergy(caloHit,EHitContribution);
			      collectedSubCalorimeterHitsWithEnergies.push_back(calorimeterSubHitWithEnergy);
			      
			    }
			    // calorimeter hit has already been assigned
			    else {
			      
			      collectedSubCalorimeterHitsWithEnergies.at(index).second += EHitContribution;
			      
			    }
		  
			    
			    // debug
			    /*
			      std::cout << "Related hit to track " << track << "  " << "MCP of track = " << MarlinUtil::getMCName(mcp->getPDG()) 
			      << " ( " << mcp->getPDG() << " )" << std::endl
			      << "SimCaloHit " << simCaloHit << "  "
			      << "Pos SimCaloHit = " << "( " << simCaloHit->getPosition()[0] << "," << simCaloHit->getPosition()[1] << "," << simCaloHit->getPosition()[2] 
			      << " )" << "  " << "CaloHit " << caloHit << "  "
			      << "Pos CaloHit = " << "( " << caloHit->getPosition()[0] << "," << caloHit->getPosition()[1] << "," << caloHit->getPosition()[2] << " )" 
			      <<  std::endl
			      << "E CaloHit (full hit) = " << EHit << "  " << "CaloHit type = " << type << "  " << "E SimCaloHit = " << ESimHit << "  " << "EHitCalculated = " 
			      << EHitCalculated << "  " << "EHitContributed = " << EHitContribution << std::endl << std::endl;		    
			    */
			    
			    
			}
			  
			  delete navCalorimeter;
			  navCalorimeter = 0;
			  
			}
			
		      }
		      
		    }

		  }
		  
		}
		
	      }
	      
	    }
	    
	  
      	    // determine ESumCalorimeterHits and ESumSubCalorimeterHits
	  
	    //debug
	    /*
	      std::cout << "collectedCalorimeterHitsWithEnergies.size() = " << collectedCalorimeterHitsWithEnergies.size() << "  " 
	      << "collectedSubCalorimeterHitsWithEnergies.size() = " << collectedSubCalorimeterHitsWithEnergies.size() << std::endl;
	    */
	    
	    for(unsigned int i = 0; i < collectedCalorimeterHitsWithEnergies.size(); ++i) {
	      
	      // debug
	      // CalorimeterHit* caloHit = collectedCalorimeterHitsWithEnergies.at(i).first;
	      float EHit = collectedCalorimeterHitsWithEnergies.at(i).second;
	      
	      // debug
	      // CalorimeterHit* caloSubHit = collectedSubCalorimeterHitsWithEnergies.at(i).first;
	      float ESubHit = collectedSubCalorimeterHitsWithEnergies.at(i).second;
	      
	      
	      ESumCalorimeterHits += EHit;
	      ESumSubCalorimeterHits += ESubHit;
	    
	      // debug
	      // std::cout << "caloHit = " << caloHit << "  " << "EHit = " << EHit << "  " << "caloSubHit = " << caloSubHit << "  " << "ESubHit = " << ESubHit << std::endl;
	      
	    }       	    

	    /*
	    double pMCP = sqrt( (mcpOfTrack->getMomentum()[0])*(mcpOfTrack->getMomentum()[0]) + (mcpOfTrack->getMomentum()[1])*(mcpOfTrack->getMomentum()[1]) +
				(mcpOfTrack->getMomentum()[2])*(mcpOfTrack->getMomentum()[2]) );
	    */

	    double eMCP = mcpOfTrack->getEnergy();

	    
	    // debug
	    /*
	    std::cout << "MCP = " << MarlinUtil::getMCName(mcpOfTrack->getPDG()) << " ( " << mcpOfTrack->getPDG() << " )" << "  " << "pMCP = " << pMCP << "  " 
		      << "EMC = " << eMCP << "  " << std::endl 
		      << "E related to track (hit accuracy) = " << ESumCalorimeterHits << "  " 
		      << "E related to track perfectly (sub-hit accuracy) = " << ESumSubCalorimeterHits << std::endl << std::endl;	      
	    */

	  
	    delete navTracks;
	    navTracks = 0;
	  
	    
	    if ( ( eMCP > _eMin ) && ( ESumSubCalorimeterHits > (_k * eMCP)) ) {
	      
	      std::cout << std::endl
			<< "--------------------------------------------------------------------------------------------------------------------------------------" 
			<< std::endl << std::endl
			<< " ==> EVENT WITH GEANT4 BUG FOUND <=="
			<< std::endl << std::endl
			<< "     EVENT WILL BE DISCARDED ..."
			<< std::endl << std::endl
			<< "--------------------------------------------------------------------------------------------------------------------------------------" 
			<< std::endl << std::endl;

	      invalidTrackFound = true;
	      break;
	      
	    }
	  
	  }
	  catch(DataNotAvailableException &e){
	    std::cout << "Collection " << _colNameRelationTrackToMCP << " or " << _colNameRelationCaloHitToSimCaloHit  << " not available in event." << std::endl;
	  };
	
	} // end of loop over tracks
    
	if (invalidTrackFound) break; // break loop over LCCollection
      
      }
      
    } // end of loop over LCLollections in this LCEvent
  
  }

  catch(DataNotAvailableException &e) {std::cout << "no valid collection in event " << _nEvt << std::endl; };
    

  bool isEventValid = !invalidTrackFound;
  setReturnValue(isEventValid);

  ++_nEvt;

}


void CutOnGEANT4Bug::check( LCEvent * evt ) {
 
}


void CutOnGEANT4Bug::end() {

}

