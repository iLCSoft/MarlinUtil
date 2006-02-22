
#include "MarlinUtil.h"
 



//============================================================================
void MarlinUtil::getMC_Balance(LCEvent * evt, double* accumulatedEnergies){
//============================================================================
  //FIXME : boundery check for the array accumulatedEnergies is needed, or a different data exchange

  int idpdg;
  const double* mom;
  float enr;
  double mass;

  // FIXME : hard coded name of MC Collection (this is the convention)
  std::string colNameMC("MCParticle");
  
  try {

    LCCollection* mcpCol = evt->getCollection(colNameMC);
    
        //-----------------------------------------------------------------------
    // Calculate balance at IP taking into account everything
    //-----------------------------------------------------------------------
    double px,py,pz,pt,ttet;
    
    double e_to_tube  = 0.;
    double e_to_tubex = 0.;
    double e_to_tubey = 0.;
    double e_to_tubez = 0.;
    int    n_to_tube = 0; 
    
    double e_neutr = 0.;   
    double e_neutrx= 0.;   
    double e_neutry= 0.;   
    double e_neutrz= 0.;   
    int    n_neutr= 0;    
    
    double e_muon = 0.;  
    double e_muonx= 0.;  
    double e_muony= 0.;  
    double e_muonz= 0.;  
    int    n_muon= 0;   
    
    double e_elect = 0.; 
    double e_electx= 0.; 
    double e_electy= 0.; 
    double e_electz= 0.; 
    int    n_elect= 0;  
    
    double e_photon = 0.;
    double e_photonx= 0.;
    double e_photony= 0.;
    double e_photonz= 0.;
    int    n_photon= 0; 
    
    double e_pi0 = 0.;
    double e_pi0x= 0.;
    double e_pi0y= 0.;
    double e_pi0z= 0.;
    int    n_pi0= 0; 
    
    double e_llhadr = 0.;
    double e_llhadrx= 0.;
    double e_llhadry= 0.;
    double e_llhadrz= 0.;
    int    n_llhadr= 0; 
    
    double e_slhadr = 0.;
    double e_slhadrx= 0.;
    double e_slhadry= 0.;
    double e_slhadrz= 0.;
    int    n_slhadr= 0; 
    
    double e_chadr = 0.; 
    double e_chadrx= 0.; 
    double e_chadry= 0.; 
    double e_chadrz= 0.; 
    int    n_chadr= 0;  

    for(int i=0; i<mcpCol->getNumberOfElements() ; i++){
      MCParticle* imc = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
      idpdg = imc-> getPDG (); 
      mom = imc-> getMomentum (); 
      enr = imc-> getEnergy (); 
      mass = imc-> getMass (); 
      if( imc-> getGeneratorStatus() == 1) { // stable particles only   
	px = mom[0]; 
	py = mom[1]; 
	pz = mom[2];
	pt = hypot(px,py);
	ttet = atan2(pt,pz);
	if ((fabs(ttet) < 0.1) || (fabs(M_PI-ttet) < 0.1)) {
	  e_to_tube  += enr;
	  e_to_tubex += px;
	  e_to_tubey += py;
	  e_to_tubez += pz;
	  n_to_tube ++;
	  continue;
	} 
	if((abs(idpdg)==12)||(abs(idpdg)==14)||(abs(idpdg)==16)) {
	  e_neutr  += enr;
	  e_neutrx += px;
	  e_neutry += py;
	  e_neutrz += pz;
	  n_neutr ++;
	  continue;
	} 
	if(abs(idpdg)==13) { // mu+ mu- 
	  e_muon  += enr;
	  e_muonx += px;
	  e_muony += py;
	  e_muonz += pz;
	  n_muon ++;
	  continue;
	} 
	if(abs(idpdg)==11) { //  e+ e-
	  e_elect  += enr;
	  e_electx += px;
	  e_electy += py;
	  e_electz += pz;
	  n_elect ++;
	  continue;
	} 
	if(idpdg == 111) { // Pi0 as stable 
	  e_pi0  += enr;
	  e_pi0x += px;
	  e_pi0y += py;
	  e_pi0z += pz;
	  n_pi0 ++;
	  continue;
	} 
	if(idpdg == 22) { // photon
	  e_photon  += enr;
	  e_photonx += px;
	  e_photony += py;
	  e_photonz += pz;
	  n_photon ++;
	  continue;
	} 
	if(    // long lived neutral hadrons
	   (abs(idpdg)==2112)|| // neutron
	   (abs(idpdg)== 130)   // KoL
	   ) {                     
	  e_llhadr  += enr;
	  e_llhadrx += px;
	  e_llhadry += py;
	  e_llhadrz += pz;
	  n_llhadr ++;
	  continue;
	}
	if(  // short lived neutral hadrons
	   (abs(idpdg)== 310)|| // KoS
	   (abs(idpdg)==3122)|| // Lambda0
	   (abs(idpdg)==3212)|| // Sigma0
	   (abs(idpdg)==3322)   // Xi0
	   ) {
	  e_slhadr  += enr;
	  e_slhadrx += px;
	  e_slhadry += py;
	  e_slhadrz += pz;
	  n_slhadr ++;
	  continue;
	}
	if(!(abs(idpdg)==12) && !(abs(idpdg)==14) && !(abs(idpdg)==16) &&  // neutrinos   
	   !(abs(idpdg)==13) && // mu+ mu- 
	   !(abs(idpdg)==11) && //  e+ e-
	   !(idpdg == 111) && // Pi0
	   !(idpdg == 22) &&  // photon
	   !(abs(idpdg)==2112) && !(abs(idpdg)== 311) && // neutral hadrons
	   !(abs(idpdg)== 130) && !(abs(idpdg)== 310) && // neutral hadrons
	   !(abs(idpdg)==3122) && !(abs(idpdg)==3212)) { // neutral hadrons
	  e_chadr  += enr;
	  e_chadrx += px;
	  e_chadry += py;
	  e_chadrz += pz;
	  n_chadr ++;
	  continue;
	}
	std::cout <<" Unknow for this program  ID is " <<idpdg<< std::endl;
      }    // Stable particles only
    }      // End for for MCParticles 
    double e_sum = e_elect+e_muon+e_chadr+e_pi0+e_photon+e_llhadr+e_slhadr+e_neutr;
    int n_sum =  n_elect+n_muon+n_chadr+n_pi0+n_photon+n_llhadr+n_slhadr+n_neutr;
    double evt_energy = e_sum;
    int n_evt = n_sum;

    //   Minus muon loosed energy = 1.6 GeV in average
    double e_mu_lost = e_muon  - n_muon*1.6;
    double e_lost = e_neutr + e_mu_lost + e_to_tube;
    int n_lost = n_neutr + n_to_tube;

    double e_real = evt_energy - e_lost;
    int n_real = n_evt - n_lost;

    if (e_real< 0.0) e_real = 0.000001;
    if (n_real < 0)  n_real = 0;

    /**  
	 std::cout <<" =============================================================="<< std::endl;
	 std::cout << " ========   Record Balance  ======="<< std::endl ;
	 std::cout <<" =============================================================="<< std::endl;
	 std::cout <<" ==============  Possible lost  ==================="<< std::endl;
	 std::cout <<"  Neutrino energy      = "<<e_neutr<<",  in "<< n_neutr<<" neutrinos"<< std::endl;
	 std::cout <<"  Energy to beam tube  = "<<e_to_tube<<",  in "<<n_to_tube<<" particles"<< std::endl;
	 std::cout <<"  Muons energy lost    = "<<e_mu_lost<<"  in "<<n_muon<<" muons"<< std::endl;
	 std::cout <<"  --------------------------------------------------"<< std::endl;
	 std::cout <<"  Total Event energy at IP = "<<evt_energy<<" [GeV]"<< std::endl;
	 std::cout <<"  --------------------------------------------------"<< std::endl;
	 std::cout <<"  Whole lost Energy        = "<<e_lost<< std::endl;
	 std::cout <<"  Available Energy in calo.= "<<e_real<< std::endl;
	 std::cout <<" =============================================================="<< std::endl;
	 std::cout <<"  Muon energy               = "<<e_muon  <<'\t'<<"  in "<< n_muon  <<" muons"<< std::endl;
	 std::cout <<"  Electron energy           = "<<e_elect <<'\t'<<"  in "<< n_elect <<" electrons"<< std::endl;
	 std::cout <<"  Charged hadron energy     = "<<e_chadr <<'\t'<<"  in "<< n_chadr <<" hadrons"<< std::endl;
	 std::cout <<"  -------------------------------------------------------------"<< std::endl;
	 std::cout <<"  Pi0 energy (if stable)    = "<<e_pi0   <<'\t'<<"  in "<< n_pi0   <<" Pi zeros"<< std::endl;
	 std::cout <<"  Photon energy             = "<<e_photon<<'\t'<<"  in "<< n_photon<<" photons"<< std::endl;
	 std::cout <<"  -------------------------------------------------------------"<< std::endl;
	 std::cout <<"  Long lived hadron energy  = "<<e_llhadr<<'\t'<<"  in "<< n_llhadr<<" hadrons"<< std::endl;
	 std::cout <<"  Short lived hadron energy = "<<e_slhadr<<'\t'<<"  in "<< n_slhadr<<" hadrons"<< std::endl;
	 std::cout <<" =============================================================="<< std::endl;
    */
    

    accumulatedEnergies[0]  = e_real;
    accumulatedEnergies[1]  = e_to_tube;
    accumulatedEnergies[2]  = e_neutr;
    accumulatedEnergies[3]  = e_elect;
    accumulatedEnergies[4]  = e_muon;
    accumulatedEnergies[5]  = e_photon;
    accumulatedEnergies[6]  = e_pi0;
    accumulatedEnergies[7]  = e_llhadr;
    accumulatedEnergies[8]  = e_slhadr;
    accumulatedEnergies[9]  = e_chadr;
    accumulatedEnergies[10] = n_real;
    accumulatedEnergies[11] = n_to_tube;
    accumulatedEnergies[12] = n_neutr;
    accumulatedEnergies[13] = n_elect;
    accumulatedEnergies[14] = n_muon;
    accumulatedEnergies[15] = n_photon;
    accumulatedEnergies[16] = n_pi0;
    accumulatedEnergies[17] = n_llhadr;
    accumulatedEnergies[18] = n_slhadr;
    accumulatedEnergies[19] = n_chadr;
    
  }
  catch(DataNotAvailableException &e){
    std::cout << "Cannot find MC Particle Collection in event "
	      << evt->getEventNumber () << std::endl ;
  };
  


 } // End  MC_Balance
