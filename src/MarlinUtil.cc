#include "MarlinUtil.h"

//#ifdef USE_SEPARATE_HEPPDT
//#include <HepPDT/TableBuilder.hh>
//#include <HepPDT/ParticleDataTable.hh>
//#include <HepPDT/TempParticleData.hh>
//#else
//#include <CLHEP/HepPDT/TableBuilder.hh>
//#include <CLHEP/HepPDT/ParticleDataTable.hh>
//#include <CLHEP/HepPDT/TempParticleData.hh>
//#endif


void MarlinUtil::printMCParticle(MCParticle* MCP, bool printDaughters) {


  double pxMC = MCP->getMomentum()[0];
  double pyMC = MCP->getMomentum()[1];
  double pzMC = MCP->getMomentum()[2];
	  
  double absPMC = sqrt( pow((pxMC),2) + pow((pyMC),2) +  pow((pzMC),2) );
  int nOfDaughters = MCP->getDaughters().size();
  
  std::cout << "MC Particle:" << std::endl
	    << "----------------------------------------------------------------------------------------------------------" << std::endl 
	    << "MCP  : " << MCP->getPDG() << "  " << " m = " << MCP->getMass() << "  " << "|p| = " << absPMC << " " 
	    << "(" << pxMC << "," << pyMC << "," << pzMC << ")" << "  " << "E = " << MCP->getEnergy() << "  " << " q = " << MCP->getCharge() << std::endl
	    << "nDaughters : " << nOfDaughters << std::endl;


  if (printDaughters) { 
    
    for(int j=0; j<nOfDaughters; ++j) {
      MCParticle* MCPDaughter = MCP->getDaughters()[j];
      
      double pxMCDaughter = MCPDaughter->getMomentum()[0];
      double pyMCDaughter = MCPDaughter->getMomentum()[1];
      double pzMCDaughter = MCPDaughter->getMomentum()[2];
      
      double absPMCDaughter = sqrt( pow((pxMCDaughter),2) + pow((pyMCDaughter),2) +  pow((pzMCDaughter),2) );
      
      int nOfDaughtersD = MCPDaughter->getDaughters().size();
      
      std::cout << j << "  " << "MCP : " << MCPDaughter->getPDG() << "  " << " m = " << MCPDaughter->getMass() << "  " << "|p| = " << absPMCDaughter << " " 
		<< "(" << pxMCDaughter << "," << pyMCDaughter << "," << pzMCDaughter << ")" << "  " << "E = " << MCPDaughter->getEnergy() << "  " << " q = " 
		<< MCPDaughter->getCharge() << std::endl
		<< "nDaughtersD : " << nOfDaughtersD << std::endl;
      
      for(int k=0; k<nOfDaughtersD; ++k) {
	MCParticle* MCPDaughterD = MCPDaughter->getDaughters()[k];
	
	
	double pxMCDaughterD = MCPDaughterD->getMomentum()[0];
	double pyMCDaughterD = MCPDaughterD->getMomentum()[1];
	double pzMCDaughterD = MCPDaughterD->getMomentum()[2];
	
	double absPMCDaughterD = sqrt( pow((pxMCDaughterD),2) + pow((pyMCDaughterD),2) +  pow((pzMCDaughterD),2) );
	
	std::cout << "---> " << k << "  " << "MCP : " << MCPDaughterD->getPDG() << "  " << " m = " << MCPDaughterD->getMass() << "  " << "|p| = " << absPMCDaughterD << " " 
		  << "(" << pxMCDaughterD << "," << pyMCDaughterD << "," << pzMCDaughterD << ")" << "  " << "E = " << MCPDaughterD->getEnergy() << "  " << " q = " 
		  << MCPDaughterD->getCharge() << std::endl;
	
      }
      
    }

  }
  
  std::cout  << "----------------------------------------------------------------------------------------------------------" << std::endl;
  
}



// ____________________________________________________________________________________________________



std::string MarlinUtil::getMCName(int PDGCode) {


  double mass = 0.0;
  double errmassp = 0.0;
  double errmassn = 0.0;
  double width = 0.0;
  double errwidthp = 0.0;
  double errwidthn = 0.0;
  std::string I("");
  std::string G("");
  std::string J("");
  std::string P("");
  std::string C("");
  std::string A("");
  int PDGCodeRead = 0;
  std::string Charge("");
  int R = 0;
  std::string S("");
  std::string Name("");
  std::string Quarks("");



  std::string Line;
  
  std::ifstream FileStream;
  CSVParser ParseStream;
  // FIXME: do not always open and close textfile
  // FIXME: use variable for filename
  FileStream.open("mass_width_2006.csv");
  if (!FileStream) {
    std::cout << "Cannot open 'mass_width_2006.csv'" << std::endl;
  }

  bool endOfFile = false;
  while (!endOfFile) {
    
    endOfFile = FileStream.eof();
    getline(FileStream, Line); // Get a line
    if (Line == "") continue;
    
    ParseStream << Line; // Feed the line to the ParseStream
    
    ParseStream >> mass  >> errmassp  >> errmassn 
		>> width  >> errwidthp  >> errwidthn 
		>> I  >> G  >> J  >> P 
		>> C  >> A  >> PDGCodeRead  >> Charge 
		>> R  >> S  >> Name  >> Quarks;
    
    if (abs(PDGCode) == PDGCodeRead) break;	  
    
  };
  
  if (endOfFile) {
    
    std::cout << std::endl << "Cannot find particle with PDG code " << PDGCode 
	      << " in file 'mass_width_2004.csv'" << std::endl;
    
    mass = 0.0;
    errmassp = 0.0;
    errmassn = 0.0;
    width = 0.0;
    errwidthp = 0.0;
    errwidthn = 0.0;
    I = "";
    G = "";
    J = "";
    P = "";
    C = "";
    A = "";
    PDGCodeRead = 0;
    Charge = "";
    R = 0;
    S = "";
    Name = "unknown";
    Quarks ="";

  }
  
  FileStream.close();


  // debug
  // std::cout << "Name = " << Name << std::endl;
  
  if ( (Name == "") || (Name == " ") ) Name = "unknown";


  // remove the blanks at the end of Name
  int blankPosition = Name.find(" ",0);
  int nameLength = Name.length();

  // debug
  // std::cout << "name = " << Name << "  " << "blankPosition = " << blankPosition << "  " << "nameLength = " << nameLength << std::endl;

  if ( blankPosition > 0 )  Name.erase(blankPosition,nameLength-blankPosition);


  return Name;




  // FIXME: take information from CLHEP
  /*
  const char pdgfile[] = "mass_width_2006.mc";
  std::ifstream pdfile;
  pdfile.open( pdgfile );
  if( !pdfile ) { 
    std::cerr << "cannot open " << pdgfile << std::endl;
    return "ERROR opening file";
  }
  // construct empty PDT
  DefaultConfig::ParticleDataTable datacol( "PDG Table" );
  {
    // Construct table builder
    HepPDT::TableBuilder  tb(datacol);
    // read the input - put as many here as you want
    if( !addPDGParticles( pdfile, tb ) ) { std::cout << "error reading PDG file " << std::endl; }
  }   // the tb destructor fills datacol

  std::string name = datacol.particle( HepPDT::ParticleID(PDGCode) )->name();

  pdfile.close();
  
  return name;
  */

}



// ____________________________________________________________________________________________________



int MarlinUtil::getPDGCode(std::string name) {
  /*
  const char pdgfile[] = "mass_width_2004.mc";
  std::ifstream pdfile;
  pdfile.open( pdgfile );
  if( !pdfile ) { 
    std::cerr << "cannot open " << pdgfile << std::endl;
    return -1;
  }
  // construct empty PDT
  DefaultConfig::ParticleDataTable datacol( "PDG Table" );
  {
    // Construct table builder
    HepPDT::TableBuilder  tb(datacol);
    // read the input - put as many here as you want
    if( !addPDGParticles( pdfile, tb ) ) { std::cout << "error reading PDG file " << std::endl; }
  }   // the tb destructor fills datacol


  int PDGCode = 0;

  DefaultConfig::ParticleData* pp = datacol.particle(name); // FIXME: Not yet implemented in CLHEP

  //  int PDGCode = datacol.particle(name)->pid();

  pdfile.close();
  
  return PDGCode;
  */

  return 0;

}


// ____________________________________________________________________________________________________


MCParticleVec MarlinUtil::getAllMCParents( MCParticle* mcPart )  {

    const MCParticleVec& parents = mcPart->getParents();
    MCParticleVec result;
    
    bool isFirstParent = parents.size() == 0;
    bool isStableOnGeneratorLevel = mcPart->getGeneratorStatus() == 1;

    // get only the first parent of the MC tree with generator status 1 (stable during generation)
    if ( isFirstParent || isStableOnGeneratorLevel ) {

      // check for invalid MC tree (e.g. 'loop' in the tree)
      MCParticleVec::const_iterator position = find(result.begin(),result.end(),mcPart);
      if ( position != result.end() ) {
	
	std::cout << "Warning: invalid MC tree found" << std::endl;
	MCParticleVec empty;
	return empty;
	
      }


      // debug
      /*
      std::cout << "parent: " << mcPart << "  " << mcPart->getPDG() << "  " 
		<< "GenStat: " << mcPart->getGeneratorStatus() << "  " << "CrSim: " << mcPart->isCreatedInSimulation() << "  " 
		<< "BckSc: " << mcPart->isBackscatter() << "  " << "VTXNotEnd: " << mcPart->vertexIsNotEndpointOfParent() << "  " 
		<< "DecTr: " << mcPart->isDecayedInTracker() << "  " << "DecCa: " << mcPart->isDecayedInCalorimeter() << "  " 
		<< "LeftDet: " << mcPart->hasLeftDetector() << "  " << "Stop: " << mcPart->isStopped() << "  " 
		<< mcPart->getEndpoint()[0] << "|" << mcPart->getEndpoint()[1] << "|" << mcPart->getEndpoint()[2] << std::endl;
      */
      
	

      result.push_back( mcPart );

    }
    else {

    
      for ( MCParticleVec::const_iterator it = parents.begin(); it != parents.end(); ++it ) {

	// debug
	/*
	std::cout << "MCP:    " << (*it) << "  " << (*it)->getPDG() << "  " 
		  << "GenStat: " << (*it)->getGeneratorStatus() << "  " << "CrSim: " << (*it)->isCreatedInSimulation() << "  " 
		  << "BckSc: " << (*it)->isBackscatter() << "  " << "VTXNotEnd: " << (*it)->vertexIsNotEndpointOfParent() << "  " 
		  << "DecTr: " << (*it)->isDecayedInTracker() << "  " << "DecCa: " << (*it)->isDecayedInCalorimeter() << "  " 
		  << "LeftDet: " << (*it)->hasLeftDetector() << "  " << "Stop: " << (*it)->isStopped() << "  " 
		  << (*it)->getEndpoint()[0] << "|" << (*it)->getEndpoint()[1] << "|" << (*it)->getEndpoint()[2]<< std::endl;
	*/


	MCParticleVec pparents = getAllMCParents( *it );

	result.insert( result.end(), pparents.begin(), pparents.end() );

      }

    }
    
    return result;
   
}


// ____________________________________________________________________________________________________



MCParticleVec MarlinUtil::getAllMCDaughters( MCParticle* mcPart )  {

  const MCParticleVec& daughters = mcPart->getDaughters();
  MCParticleVec result;
  
  // use this if-else to get only the last daughters of the MC tree
  //     if ( daughters.size() == 0 )
  //       {
  
  // check for invalid MC tree (e.g. 'loop' in the tree)
  MCParticleVec::const_iterator position = find(result.begin(),result.end(),mcPart);
  if ( position != result.end() ) {
    
    std::cout << "Warning: invalid MC tree found" << std::endl;
    MCParticleVec empty;
    return empty;
    
  }
  
  result.push_back( mcPart );
  
  //       }
  //     else
  
  
  for ( MCParticleVec::const_iterator it = daughters.begin(); it != daughters.end(); ++it ) {
    
    MCParticleVec ddaughters = getAllMCDaughters( *it );
    result.insert( result.end(), ddaughters.begin(), ddaughters.end() );
  }
  
  return result;

}


// ____________________________________________________________________________________________________


bool MarlinUtil::isDaughterOf( MCParticle* daughter, MCParticle* parent )  {


  const MCParticleVec& daughters = parent->getDaughters();
  bool isDaughter = false;

  for ( MCParticleVec::const_iterator it = daughters.begin(); it != daughters.end(); ++it ) {
    
    if ( (*it) == daughter) {
      
      isDaughter = true;
      break;      

    }

  }

  
  if ( !isDaughter ) {
  
    for ( MCParticleVec::const_iterator it = daughters.begin(); it != daughters.end(); ++it ) {
      
      bool isDDaughter = isDaughterOf( daughter,(*it) );

      if (isDDaughter) {

	isDaughter = true;
	break; 
	
      }
     
    }
    
  }

  return isDaughter;
  
}


// ____________________________________________________________________________________________________


bool MarlinUtil::DecayChainInTree(std::vector<int> DecayChannel, LCEvent* evt)
{

  // FIXME rewrite the whole method!

  // fixed to Z0 to hadrons

  bool flag = false;

  LCCollection* MCcol = evt->getCollection("MCParticle") ;

  int m = MCcol->getNumberOfElements();
  for(int j=0; j<m; ++j){
    
    MCParticle* MCP = dynamic_cast<MCParticle*>(MCcol->getElementAt(j));
    
    int idpdg = MCP->getPDG();
    
    if (idpdg == 23) {
    
      MCParticleVec daughters = MCP->getDaughters();

      int nOfDaughters = daughters.size();
      
      for(int i=0; i<nOfDaughters; ++i){
	  
	int PDGDaughter = daughters.at(i)->getPDG();

	if ( (PDGDaughter == 1) || (PDGDaughter == 2) || (PDGDaughter == 3) ) {

	  int nu = 0;
	  int nd = 0;
	  int ns = 0;
	  
	  for(int k=0; k<nOfDaughters; ++k){

	    int PDGDaughterCheck = abs(daughters.at(k)->getPDG());

	    if ( (PDGDaughterCheck == 1) || (PDGDaughterCheck == 2) || (PDGDaughterCheck == 3) || (PDGDaughterCheck == 21) ) {

	      if (PDGDaughterCheck == 1) ++nu;
	      if (PDGDaughterCheck == 2) ++nd;
	      if (PDGDaughterCheck == 3) ++ns;
	      
	    }
	    else {
	      flag = false;
	      break;
	    }

	  }
	  
	  if ( ( (nu==2) && (nd==0) && (ns==0) ) || ( (nu==0) && (nd==2) && (ns==0) ) || ( (nu==0) && (nd==0) && (ns==2) ) ) {
	    flag = true;
	    break;
	  }
	  
	}
	else flag = false;
	
      }
    }
    
    if (flag) break;

  }

  return flag;





  /*

  // Simple search on the tree. Return value equals to 1, if DecayChannel was found, and 0, if not.
  // Does't check if decay chain apperars more than one time in the tree. That means, this  method returns 1, if
  // at least one time the specific decay chain apperars.

  //  ChainOfParticlesPdg[0]=0;      // set first particle in chain as initial particle

  // for (int k=0; k<DecayChannel.size(); k++)   cout <<  ChainOfParticlesPdg[k] << endl;  // debug

 bool InChain = false;
 bool flag = true;
 int ParentN  = 0;
 
 for (int k=0; k<MCcol.size(); k++)
   {
     ParentN = k;
     flag = true;
     //cout << endl << "ParentN before " << ParentN << endl; //debug
     
     for (int l=DecayChannel.size()-1; l>1; l--) // search for ALL parents
       {


	 // GO ON HERE !!!



	 flag &= TreeParentID(ParentN)==ChainOfParticlesPdg[l-1]; 
	 //cout << l << "  " << TreeParentID(ParentN) << "  " << ChainOfParticlesPdg[l-1] << "  " << flag 
	 //     << "  " << ParentN;
	 ParentN = TreeParentN(ParentN);
	 //cout << "  " << ParentN << endl;
       }
      flag &= tree_stables_id[k]==ChainOfParticlesPdg[DecayChannel.size()-1];  // compare with last
      flag &= tree_stables_id[k]!=TreeParentID(k);   // parent and particle should not be the same

      ParentN = k;
      for (int l=DecayChannel.size()-1; l>0; l--)  // search for initial 'particle' in chain
	{
	  ParentN = TreeParentN(ParentN);
	  if (l==1) flag &= ParentN==0; 
        }
      InChain |= flag;
    }
  return InChain;

  */

}


// ____________________________________________________________________________________________________


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
	//FIXME : Hard coded cut. Use GEAR instead
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
	   !(abs(idpdg)==3122) && !(abs(idpdg)==3212) && !(abs(idpdg)==3322) ) { // neutral hadrons
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

    //   Minus muon loosed energy = 1.6 GeV in average, i.e. Muons deposit 1.6 GeV an average in Calorimeter
    double e_mu_lost = e_muon  - n_muon*1.6;
    double e_lost = e_neutr + e_mu_lost + e_to_tube;
    int n_lost = n_neutr + n_to_tube;

    double e_real = evt_energy - e_lost;
    int n_real = n_evt - n_lost;

    if (e_real< 0.0) e_real = 0.000001;
    if (n_real < 0)  n_real = 0;

      
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
    accumulatedEnergies[20] = evt_energy - e_neutr - e_to_tube;






    
  }
  catch(DataNotAvailableException &e){
    std::cout << "Cannot find MC Particle Collection in event "
	      << evt->getEventNumber () << std::endl ;
  };
  


 } // End  MC_Balance



// ____________________________________________________________________________________________________



void MarlinUtil::printTrack(Track* track, double bField) {

  const double* p = getMomentum(track,bField);
  const double pAbs = getAbsMomentum(track,bField);
  const double pt = sqrt( p[0]*p[0] + p[1]*p[1] );
  double d0 = track->getD0();
  double z0 = track->getZ0();
  double omega = track->getOmega();
  double phi0  = track->getPhi();
  double tanlambda = track->getTanLambda();
  
  std::cout << "Track id() : " << track->id() << "  " << "|p| = " << pAbs << "  " << "pt = " << pt << "  " 
	    << "p = " << "(" << p[0] << "," << p[1] << "," << p[2] << ")" << std::endl
	    << "(d0,z0,phi0,omega,tanL) = " << "(" << d0 << "," << z0 << "," << phi0 << "," << omega << "," << tanlambda << ")" << std::endl 
	    << "Ri = " << track->getRadiusOfInnermostHit() << "  " << "# of SubTracks = " << track->getTracks().size() 
	    << "   " << "# of tracker hits = " << track->getTrackerHits().size() << std::endl << std::endl;
  

  delete[] p;

}



// ____________________________________________________________________________________________________



const double* MarlinUtil::getMomentum(Track* track, double bField) {

  // user need to care about deletion of the array
  
  double d0 = track->getD0();
  double z0 = track->getZ0();
  double omega = track->getOmega();
  double phi0  = track->getPhi();
  double tanlambda = track->getTanLambda();
  
  HelixClass* helix = new HelixClass();
  helix->Initialize_Canonical(phi0, d0, z0, omega, tanlambda, bField);

  double* p = new double[3];

  for (int k=0; k < 3; ++k) p[k] = helix->getMomentum()[k];
  
  delete helix;
  helix = 0;

  return p;

}



// ____________________________________________________________________________________________________



const double MarlinUtil::getAbsMomentum(Track* track, double bField) {
  
  const double* p = getMomentum(track,bField);
  
  double pAbs = 0.0;

  for (int i = 0; i < 3 ; ++i) pAbs += p[i]*p[i];
  pAbs = sqrt(pAbs);

  const double pAbsReturn = pAbs;

  return pAbsReturn;
  
}


// ____________________________________________________________________________________________________

void MarlinUtil::printCluster(Cluster* cluster) {
    
  double energy = cluster->getEnergy();
  unsigned int nHits = cluster->getCalorimeterHits().size();
  double x = cluster->getPosition()[0];
  double y = cluster->getPosition()[1];	 
  double z = cluster->getPosition()[2];
  double r = sqrt(x*x + y*y + z*z);
  
  std::cout << "Cluster id() : " << cluster->id() << "  " << "|r| = " << r << "  " << "(" << x << "," << y << "," << z << ")" << "  " << "E = " << energy << "  "
	    << "# of hits = " << nHits << "  " << "# of Particle IDs = " << cluster->getParticleIDs().size()
	    << "  " << "# of SubClusters = " << cluster->getClusters().size() << std::endl << std::endl;
  
}



// ____________________________________________________________________________________________________



void MarlinUtil::printRecoParticle(ReconstructedParticle* recoParticle, double bField) {


  double xRef = recoParticle->getReferencePoint()[0];
  double yRef = recoParticle->getReferencePoint()[1];
  double zRef = recoParticle->getReferencePoint()[2];
  
  double pxReco = recoParticle->getMomentum()[0];
  double pyReco = recoParticle->getMomentum()[1];
  double pzReco = recoParticle->getMomentum()[2];
  double absPReco = sqrt(pxReco*pxReco + pyReco*pyReco + pzReco*pzReco); 
  double Energy = recoParticle->getEnergy();
  double mass = recoParticle->getMass();
  
  double checkEnergy = sqrt(absPReco*absPReco + mass*mass);      

  int type = recoParticle->getType();
  std::string typeName("not assigned");
  
  switch (type) {
  case 0  : typeName = "no type";                        break;
  case 1  : typeName = "electron/positron";              break;
  case 2  : typeName = "charged hadron (Wolf: or muon)"; break;
  case 3  : typeName = "photon";                         break;
  case 4  : typeName = "neutral hadron";                 break;
  case 5  : typeName = "muon";                           break;
  case 15 : typeName = "compound object";                break;
  }
  
  int nParticleIDs = recoParticle->getParticleIDs().size();
  
  const TrackVec Tracks = recoParticle->getTracks();
  const ClusterVec Clusters = recoParticle->getClusters();
  int nOfTracks = Tracks.size();
  int nOfClusters = Clusters.size();
  
      
  std::cout << "--------------------------------------------------------------------------------------------------------------------------------------" << std::endl 
	    << typeName << "  " << "|p| = " << absPReco << "  " << "(" << pxReco << "," << pyReco << "," << pzReco << ")" 
	    << "  " << "E = " << Energy << " (" << checkEnergy << ")" << "  " << " m = " << recoParticle->getMass() << "  " 
	    << "q = " << recoParticle->getCharge() 	<< std::endl 
	    << "RefPoint = (" << xRef << "," << yRef << "," << zRef << ")" << "  " << "( LCObjectID = " << recoParticle->id() << " )" << "  " 
	    << "# of Particle IDs = " << nParticleIDs << "  ";
  
  for(int j=0; j<nParticleIDs; ++j) {
    std::cout << "(" << recoParticle->getParticleIDs()[j];
    if (j<nParticleIDs-1) std::cout << ",";
    else std::cout << ")" << "  ";
  }
  
  std::cout << "ID used = " << recoParticle->getParticleIDUsed() << std::endl
	    << "nTracks = " << nOfTracks << std::endl;
      
  for (int iOfTracks = 0; iOfTracks<nOfTracks; ++iOfTracks) {
    
    Track* track = recoParticle->getTracks()[iOfTracks];

    printTrack(track,bField);

  }

  // FIXME:  same loop twice, not so nice 060905 OW
  double energyOfClusters = 0.0;
  for (int iOfClusters = 0; iOfClusters<nOfClusters; ++iOfClusters) {
    
    Cluster* cluster = recoParticle->getClusters()[iOfClusters];
    energyOfClusters += cluster->getEnergy();
    
  }

  std::cout << "nClusters = " << nOfClusters << "  " << "EClusters = " << energyOfClusters << std::endl;
      
  for (int iOfClusters = 0; iOfClusters<nOfClusters; ++iOfClusters) {
    
    Cluster* cluster = recoParticle->getClusters()[iOfClusters];
    
    double energy = cluster->getEnergy();
    double x = cluster->getPosition()[0];
    double y = cluster->getPosition()[1];	 
    double z = cluster->getPosition()[2];
    double r = sqrt(x*x + y*y + z*z);
    
    std::cout << iOfClusters << "  " << "E = " << energy << "  " << "|r| = " << r << "  " << "(" << x << "," << y << "," << z << ")" << "  "
	      << "# of Particle IDs = " << cluster->getParticleIDs().size()
	      << "  " << "# of SubClusters = " << cluster->getClusters().size() << std::endl;
    
  }

  std::cout << "--------------------------------------------------------------------------------------------------------------------------------------" << std::endl 
	    << std::endl;
      

}



// ____________________________________________________________________________________________________



int MarlinUtil::countAllSimTrackerHits(LCEvent* evt,MCParticle* MCP) {

  int counter = 0;

  std::vector< std::string >::const_iterator iter;
  const std::vector< std::string >* ColNames = evt->getCollectionNames();
  
  for( iter = ColNames->begin() ; iter != ColNames->end() ; iter++) {
    
    LCCollection* col = evt->getCollection( *iter ) ;
    
    if ( col->getTypeName() == LCIO::SIMTRACKERHIT ) {
      
      int n = col->getNumberOfElements();

      for(int i=0; i<n; ++i){

	SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(col->getElementAt(i));
	
	if ( hit->getMCParticle() == MCP )  ++counter;
	
      }

    }
    
  }
  
  return counter;  
  
}



// ____________________________________________________________________________________________________



int MarlinUtil::countAllSimCalorimeterHits(LCEvent* evt,MCParticle* MCP,double& accumulatedSimCaloEnergy) {

  // initialise
  int counter = 0;
  accumulatedSimCaloEnergy = 0.0;

  std::vector< std::string >::const_iterator iter;
  const std::vector< std::string >* ColNames = evt->getCollectionNames();
  
  for( iter = ColNames->begin() ; iter != ColNames->end() ; iter++) {
    
    LCCollection* col = evt->getCollection( *iter ) ;
    
    if ( col->getTypeName() == LCIO::SIMCALORIMETERHIT ) {
      
      int n = col->getNumberOfElements();

      for(int i=0; i<n; ++i){

	SimCalorimeterHit* hit = dynamic_cast<SimCalorimeterHit*>(col->getElementAt(i));
	
	for(int j=0; j<hit->getNMCContributions (); ++j){
	 
	  if ( hit->getParticleCont(j) == MCP ) {

	    ++counter;
	    accumulatedSimCaloEnergy += hit->getEnergyCont(j);

	  }
	  
	}
	
      }
      
    }
    
  }
  
  return counter;  
  
}



// ____________________________________________________________________________________________________


double MarlinUtil::getEnergyDepositedInFullCalorimeter(LCEvent* evt) {

  // initialise
  double accumulatedCaloEnergy = 0.0;

  std::vector< std::string >::const_iterator iter;
  const std::vector< std::string >* ColNames = evt->getCollectionNames();
  
  for( iter = ColNames->begin() ; iter != ColNames->end() ; iter++) {
    
    LCCollection* col = evt->getCollection( *iter ) ;
    
    if ( col->getTypeName() == LCIO::CALORIMETERHIT ) {
      
      int n = col->getNumberOfElements();

      for(int i=0; i<n; ++i){

	CalorimeterHit* hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(i));
	
	accumulatedCaloEnergy += hit->getEnergy();
		       	
      }
      
    }
    
  }
  
  return accumulatedCaloEnergy;  
  
}



// ____________________________________________________________________________________________________
// ____________________________________________________________________________________________________





MCParticleHelper::MCParticleHelper() {


  _pdgCodesMCParticles.clear();
  _massMCParticles.clear();
  _nameMCParticles.clear();
  _chargeMCParticles.clear();


  double mass = 0.0;
  double errmassp = 0.0;
  double errmassn = 0.0;
  double width = 0.0;
  double errwidthp = 0.0;
  double errwidthn = 0.0;
  std::string I("");
  std::string G("");
  std::string J("");
  std::string P("");
  std::string C("");
  std::string A("");
  int PDGCodeRead = 0;
  std::string Charge("");
  int R = 0;
  std::string S("");
  std::string Name("");
  std::string Quarks("");
  
  

  std::string Line;
  
  std::ifstream FileStream;
  CSVParser ParseStream;
  // FIXME: do not always open and close textfile
  // FIXME: use variable for filename
  FileStream.open("mass_width_2006.csv");
  if (!FileStream) {
    std::cout << "Cannot open 'mass_width_2006.csv'" << std::endl;
  }

  bool endOfFile = false;
  unsigned int numberOfLines = 0;

  while (!endOfFile) {
    
    endOfFile = FileStream.eof();
    getline(FileStream, Line); // Get a line
    if (Line == "") continue;
    
    ParseStream << Line; // Feed the line to the ParseStream
    
    ParseStream >> mass  >> errmassp  >> errmassn 
		>> width  >> errwidthp  >> errwidthn 
		>> I  >> G  >> J  >> P 
		>> C  >> A  >> PDGCodeRead  >> Charge 
		>> R  >> S  >> Name  >> Quarks;

    if ( ( (P != "+") && (P != "-") && (P != "?") && (P != "") ) || (Charge.length() == 0) ) continue;

    ++numberOfLines;

     // debug
    /*
    std::cout << "m: " << mass << "  " << "Em: " << errmassp << "|" << errmassn << "  " << "width: " << width << "  " << "Ew: " << errwidthp << "|" << errwidthn << "  " 
	      << "I: " << I << "  " << "G: " << G << "  " << "J: " << J << "  " << "P: " << P << "  " << "C: " << C << "  " << "A: " << A << "  " << "PDG: " << PDGCodeRead 
	      << "  " << "Q: " << Charge << "  " << "R: " << R << "  " << "S: " << S << "  " << Name << "  " << Quarks << "  " << Charge.length() << std::endl;
    */


    if ( (Name == "") || (Name == " ") ) Name = "unknown";


    // remove the blanks at the end of Name
    int blankPosition = Name.find(" ",0);
    int nameLength = Name.length();

    // debug
    // std::cout << "name = " << Name << "  " << "blankPosition = " << blankPosition << "  " << "nameLength = " << nameLength << std::endl;
    
    if ( blankPosition > 0 )  Name.erase(blankPosition,nameLength-blankPosition);

    _pdgCodesMCParticles.push_back(PDGCodeRead);
    _massMCParticles.push_back(mass);
    _nameMCParticles.push_back(Name);
    _chargeMCParticles.push_back(Charge);
    
  };
  
  FileStream.close();

  // debug
  //  std::cout << "numberOfLines: " << numberOfLines << std::endl;





  // FIXME: take information from CLHEP
  /*
  const char pdgfile[] = "mass_width_2006.mc";
  std::ifstream pdfile;
  pdfile.open( pdgfile );
  if( !pdfile ) { 
    std::cerr << "cannot open " << pdgfile << std::endl;
    return "ERROR opening file";
  }
  // construct empty PDT
  DefaultConfig::ParticleDataTable datacol( "PDG Table" );
  {
    // Construct table builder
    HepPDT::TableBuilder  tb(datacol);
    // read the input - put as many here as you want
    if( !addPDGParticles( pdfile, tb ) ) { std::cout << "error reading PDG file " << std::endl; }
  }   // the tb destructor fills datacol

  std::string name = datacol.particle( HepPDT::ParticleID(PDGCode) )->name();

  pdfile.close();
  
  */


  
}


std::string MCParticleHelper::getMCCharge(int PDGCode) {  


  std::vector<int>::const_iterator position = find(_pdgCodesMCParticles.begin(),_pdgCodesMCParticles.end(),abs(PDGCode)); 


  if ( position == _pdgCodesMCParticles.end() ) {
    
    std::string returnValue("");
    return returnValue;

  }    
  else {

    int index = position - _pdgCodesMCParticles.begin();

    std::string chargeString = _chargeMCParticles.at(index);
    
    if ( PDGCode < 0 ) {
      
      if ( (chargeString.length() > 0) && (chargeString.at(0) == '-') ) chargeString.replace(0,1,"+");
      else if ( (chargeString.length() > 0) && (chargeString.at(0) == '+') ) chargeString.replace(0,1,"-");
      
      // special cases of ++ and -- particles
      else if ( (chargeString.length() > 1) && (chargeString.at(0) == '-') && (chargeString.at(1) == '-') ) chargeString.replace(0,2,"++"); 
      else if ( (chargeString.length() > 1) && (chargeString.at(0) == '+') && (chargeString.at(1) == '+') ) chargeString.replace(0,2,"--"); 

      else { } // charge is 0, do nothing
      

    }


    return chargeString;

  }
  
}


// ____________________________________________________________________________________________________

