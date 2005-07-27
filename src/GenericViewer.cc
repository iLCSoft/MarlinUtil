#include "GenericViewer.h"
#include <EVENT/LCObject.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/Cluster.h>
#include <EVENT/Track.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>
#include <iostream>
#include "ClusterShapes.h"
#include <ced_cli.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>
#include "HelixClass.h"
#include <math.h>

using namespace lcio ;
using namespace marlin ;


GenericViewer aGenericViewer ;


GenericViewer::GenericViewer() : Processor("GenericViewer") {

  _description = "Drawing Utility" ;
  
  std::vector<std::string> simCaloHitCollections;

  simCaloHitCollections.push_back(std::string("ecal02_EcalBarrel"));
  simCaloHitCollections.push_back(std::string("ecal02_EcalEndcap"));
  simCaloHitCollections.push_back(std::string("hcalFeRPC1_HcalBarrelEnd"));
  simCaloHitCollections.push_back(std::string("hcalFeRPC1_HcalBarrelReg"));
  simCaloHitCollections.push_back(std::string("hcalFeRPC1_HcalEndCaps"));


  registerProcessorParameter( "SimCaloHitCollections" , 
			      "Sim Calo Hit Collection Names" ,
			      _simCaloHitCollections ,
			       simCaloHitCollections);

  std::vector<std::string> caloHitCollections;
  caloHitCollections.push_back(std::string("ECAL"));
  caloHitCollections.push_back(std::string("HCAL"));

  registerProcessorParameter("CaloHitCollections" , 
			     "Calo Hit Collection Names" , 
			     _caloHitCollections , 
			     caloHitCollections);
  
  std::vector<std::string> simTrackerHitCollections;
  simTrackerHitCollections.push_back("tpc03_TPC");

  registerProcessorParameter("SimTrackerHitCollections" , 
			     "Sim Tracker Hit Collection Names" , 
			     _simTrackerHitCollections , 
			     simTrackerHitCollections);

  std::vector<std::string> trackerHitCollections;
  trackerHitCollections.push_back("TPCTrackerHits");

  registerProcessorParameter("TrackerHitCollections" , 
			     "Tracker Hit Collection Names" , 
			     _trackerHitCollections , 
			     trackerHitCollections);

  
  registerProcessorParameter("TrueClusterCollection" , 
			     "True Cluster Collection Name" , 
			     _trueClustersCollection , 
			     std::string("TrueClusters"));

  
  registerProcessorParameter("ClusterCollection" , 
			     "Cluster Collection Name" , 
			     _clustersCollection , 
			     std::string("ClustersAR"));

  
  registerProcessorParameter("TrueTrackCollection" , 
			     "True Track Collection Name" , 
			     _trueTracksCollection , 
			     std::string("TrueTracks"));
  
  registerProcessorParameter("TrackCollection" , 
			     "Track Collection Name" , 
			     _tracksCollection , 
			     std::string("TPC_Tracks"));

  registerProcessorParameter("ParticleCollection" , 
			     "Particle Collection Name" , 
			     _particleCollection , 
			     std::string("RecoParticles"));

  registerProcessorParameter("LayerCaloHit", 
			     "Layer for Calo Hits",
			     _layerCaloHit, 
			     (int)-1);  

  registerProcessorParameter("LayerTrackerHit", 
			     "Layer for Tracker Hits",
			     _layerTrackerHit, 
			     (int)-1);

  registerProcessorParameter("LayerTrueClusters", 
			     "Layer for True Clusters",
			     _layerTrueClusters, 
			     (int)-1);

  registerProcessorParameter("LayerClusters", 
			     "Layer for Reco Clusters",
			     _layerClusters, 
			     (int)-1);

  registerProcessorParameter("LayerTrueTracks", 
			     "Layer for True Tracks",
			     _layerTrueTracks, 
			     (int)-1);

  registerProcessorParameter("LayerTracks", 
			     "Layer for Tracks",
			     _layerTracks, 
			     (int)-1);

  registerProcessorParameter("LayerSimTrackerHit", 
			     "Layer for Sim Tracker Hits",
			     _layerSimTrackerHit, 
			     (int)-1);

  registerProcessorParameter("LayerSimCaloHit", 
			     "Layer for Sim Calo Hits",
			     _layerSimCaloHit, 
			     (int)-1);

  registerProcessorParameter("LayerReco",
			     "Layer for Reco Particles",
			     _layerReco,
			     (int)9);

  registerProcessorParameter("DetectorModel",
			     "Detector Model",
			     _detModel,
			     (int)0);

  registerProcessorParameter("MagneticField",
			     "Magnetic Field",
			     _bField,
			     (float)4.0);

}

void GenericViewer::init() {

    _nRun = -1;
    _nEvt = 0;
    ced_client_init("localhost",7286);
    ced_register_elements();

}


void GenericViewer::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
  _nEvt = 0;
} 

void GenericViewer::processEvent( LCEvent * evt ) { 

    std::cout << std::endl;
    std::cout << " +++++++++++++++++++++++++++" << std::endl;
    std::cout << "       Event Number : " << _nEvt << std::endl;
    std::cout << std::endl;

    _mcpList.clear();
    
    LCCollection * mcpcol ;
    try {
	mcpcol = evt->getCollection("MCParticle");
	int nelem = mcpcol->getNumberOfElements();
	std::cout << " # of MCParticles = " << nelem << std::endl;
	for (int ielem(0); ielem < nelem; ++ielem) {
	    MCParticle * mcp = 
		dynamic_cast<MCParticle*>(mcpcol->getElementAt(ielem));
	    _mcpList[mcp] = ielem;
	    /*	    std::cout << "Particle : " << ielem 
		      << " ; PDG : " << mcp->getPDG()
		      << " ; Energy : " << mcp->getEnergy() 
		      << " ; Charge : " << mcp->getCharge() 
		      << " ; Mass   : " << mcp->getMass() 
		      << std::endl; */
	    if (abs(mcp->getCharge()) > 0.) {
	      HelixClass * helix = new HelixClass();
	      float ver[3];
	      float mom[3];
	      float charg = (float)mcp->getCharge();
	      for (int ii=0; ii<3; ++ii) {
		ver[ii] = (float)mcp->getVertex()[ii];
		mom[ii] = (float)mcp->getMomentum()[ii];
	      }
	      helix->Initialize_VP(ver,mom,charg,_bField);
	      /* std::cout << "Track parameters ---> " << std::endl;
	      std::cout << "d0 = " << helix->getD0()
			<< " ; z0 = " << helix->getZ0()
			<< " ; omega = " << helix->getOmega()
			<< " ; tanlam = " << helix->getTanLambda()
			<< " ; tan(phi0)  = " << tan(helix->getPhi0())
			<< " ; pt = " << helix->getPXY()
			<< std::endl; */
	      delete helix;
	    } 
	}
    }
    catch(DataNotAvailableException &e) {
	std::cout << "No MCParticle collection found " << std::endl;
	exit(1);
    }

    std::cout << std::endl;
    std::cout << std::endl;
    

// Drawing Geometry

    ced_new_event();  // Reset drawing buffer and START drawing collection

    //          The Simplest geometry has Z axis as beam axis

    static CED_GeoCylinder geoCylindersLDC[] = {       // for TESLA Detector Geometry
      //      {    50.0,  6,  0.0, 5658.5, -5658.5, 0xff      }, // beam tube
      {   380.0, 24,  0.0, 2658.5, -2658.5, 0xff      }, // inner TPC
      {  1840.0,  8, 22.5, 2700.0, -2700.0, 0x7f7f1f  }, // inner ECAL
      {  2045.7,  8, 22.5, 2700.0, -2700.0, 0x7f7f1f  }, // outer ECAL
      {  2045.7,  8, 22.5, 101.00,  2820.0, 0x7f7f1f  }, // endcap ECAL
      {  2045.7,  8, 22.5, 101.00, -3022.0, 0x7f7f1f  }, // endcap ECAL
      {  3000.0, 16,  0.0, 2658.5, -2658.5, 0xcf00    }, // outer HCAL
      {  3000.0,  8, 22.5, 702.25,  2826.0, 0xcf00    }, // endcap HCAL
      {  3000.0,  8, 22.5, 702.25, -4230.5, 0xcf00    }, // endcap HCAL
	// radius, poligon order, angle degree, 1/2 length, shift in z, color
    }; 
    
    static CED_GeoCylinder geoCylindersSiD[] = {       // for SiD Detector Geometry
      //	{    12.0,  100,  0.0, 5658.5, -5658.5, 0xff      }, // beam tube
	{  186.35,   20,  0.0,  271.0,  -271.0, 0xff      }, // 1st SiD layer
	{  448.85,   20,  0.0,  621.0,  -621.0, 0xff      }, // 2nd SiD layer
	{  711.35,   20,  0.0,  971.0,  -971.0, 0xff      }, // 3d  SiD layer
	{  973.85,   20,  0.0, 1321.0, -1321.0, 0xff      }, // 4th SiD layer
	{ 1236.35,   20,  0.0, 1649.4, -1649.4, 0xff      }, // 5th SiD layer
	{ 1270.00,   30,  0.0, 1682.5, -1682.5, 0x7f7f1f  }, // inner ECAL
	{ 1385.00,   30,  0.0, 1682.5, -1682.5, 0x7f7f1f  }, // outer ECAL
	{ 1385.00,   30,  0.0,  562.5,  1682.5, 0x7f7f1f  }, // endcap ECAL
	{ 1385.00,   30,  0.0,  562.5, -2807.5, 0x7f7f1f  }, // endcap ECAL
	{ 2500.00,   30,  0.0, 1795.0, -1795.0, 0xcf00    }, // outer  HCAL  
	{ 2500.00,   30,  0.0, 487.25 , 1795.0, 0xcf00    }, // endcap HCAL
	{ 2500.00,   30,  0.0, 487.25 ,-2770.0, 0xcf00    }, // endcap HCAL      
    };


    static CED_GeoCylinder geoCylindersGLD[] = {       // for Huge Detector Geometry
      //      {    50.0,  6,  0.0, 5658.5, -5658.5, 0xff      }, // beam tube
      {   400.0, 24,  0.0, 2600.0, -2600.0, 0xff      }, // inner TPC
      {  2100.0, 30, 22.5, 2700.0, -2700.0, 0x7f7f1f  }, // inner CALO
      {  3405.0, 30, 22.5, 2700.0, -2700.0, 0x7f7f1f  }, // outer CAL0
      {  3405.0, 30, 22.5,  737.5,  2700.0, 0x7f7f1f  }, // endcap CAL0
      {  3405.0, 30, 22.5,  737.5, -4175.0, 0x7f7f1f  }, // endcap ECAL
      //      {  3000.0, 16,  0.0, 2658.5, -2658.5, 0xcf00    }, // outer HCAL
      //      {  3000.0,  8, 22.5, 702.25,  2826.0, 0xcf00    }, // endcap HCAL
      //      {  3000.0,  8, 22.5, 702.25, -4230.5, 0xcf00    }, // endcap HCAL
	// radius, poligon order, angle degree, 1/2 length, shift in z, color
    }; 

    /*
      static CED_GeoCylinder geoCylinders[] = {    // for Prototype
      {    180.0,  4,  45.0, 110.0, 0.0, 0xff },   // beam tube
      {    500.0,  4,  45.0, 250.0, 220., 0xff }   // inner TPC
      };
    */

    if (_detModel == 0) {
      ced_geocylinders(sizeof(geoCylindersLDC)/sizeof(CED_GeoCylinder),geoCylindersLDC);
    }
    else if (_detModel == 1) {
      ced_geocylinders(sizeof(geoCylindersSiD)/sizeof(CED_GeoCylinder),geoCylindersSiD);
    }
    else if (_detModel == 2) {
      ced_geocylinders(sizeof(geoCylindersGLD)/sizeof(CED_GeoCylinder),geoCylindersGLD);
    }


// Drawing Simulated Tracker Hits 
    if (_layerSimTrackerHit >= 0) {
	for ( unsigned int icol = 0; 
	      icol < _simTrackerHitCollections.size(); ++ icol) {
	    try {
		LCCollection * col = 
		  evt->getCollection(_simTrackerHitCollections[icol].c_str()) ;
		int nelem = col->getNumberOfElements();
		int counter = 0;
		for (int ielem = 0; ielem < nelem; ++ielem) {
		  SimTrackerHit * hit = 
		    dynamic_cast<SimTrackerHit*>( col->getElementAt(ielem));
		  const MCParticle * mcp = hit->getMCParticle();
		  if (mcp == NULL) {
		    counter = 0;
		  }
		  else {
		    std::map<MCParticle *, int >::iterator pos;
		    for (pos = _mcpList.begin(); pos != _mcpList.end(); ++pos) 
		      {
			MCParticle * part = pos->first; 
			counter = pos->second;
			if (part == mcp)
			  break;
		      }
		  }
		  
		  int color = returnColor(counter);
		  float x = (float)hit->getPosition()[0];
		  float y = (float)hit->getPosition()[1];
		  float z = (float)hit->getPosition()[2];
		  ced_hit(x,y,z, 
			  2 | (_layerSimTrackerHit<<CED_LAYER_SHIFT),30,color);
		}	    
	    }
	    catch(DataNotAvailableException &e) {}
	    ced_send_event();
	}
    }

// Drawing Simulated Calorimeter Hits
    if (_layerSimCaloHit >= 0) {
      float totenergy = 0.0;
      int tothits = 0;
      for ( unsigned int icol = 0; 
	    icol < _simCaloHitCollections.size(); ++ icol) {
	try {
	  LCCollection * col = 
	    evt->getCollection(_simCaloHitCollections[icol].c_str()) ;
	  int nelem = col->getNumberOfElements();
	  for (int ielem = 0; ielem < nelem; ++ielem) {
	    SimCalorimeterHit * hit = 
	      dynamic_cast<SimCalorimeterHit*>( col->getElementAt(ielem));
	    int counter = 0;
	    totenergy += hit->getEnergy();
	    if (hit->getNMCParticles() > 0) {
	      MCParticle * par = hit->getParticleCont(0);
	      counter = _mcpList[par];
	    }
	    else {
	      counter = 0;
	    }
	    
	    tothits++;
	    int color = returnColor(counter);
	    float x = (float)hit->getPosition()[0];
	    float y = (float)hit->getPosition()[1];
	    float z = (float)hit->getPosition()[2];
	    ced_hit(x,y,z, _layerSimTrackerHit<<CED_LAYER_SHIFT,3,color);
	  }
	}
	catch(DataNotAvailableException &e) {}	
	ced_send_event();
      }
    }
    

// Drawing Tracker Hits 
    if (_layerTrackerHit >= 0 ) {
	for ( unsigned int icol = 0; icol < _trackerHitCollections.size(); ++ icol) {
	    try {
		LCCollection * col = 
		  evt->getCollection(_trackerHitCollections[icol].c_str()) ;
		int nelem = col->getNumberOfElements();
		for (int ielem = 0; ielem < nelem; ++ielem) {
		    TrackerHit * hit = 
		      dynamic_cast<TrackerHit*>( col->getElementAt(ielem));
		    float x = (float)hit->getPosition()[0];
		    float y = (float)hit->getPosition()[1];
		    float z = (float)hit->getPosition()[2];
		    ced_hit(x,y,z, 
			    _layerTrackerHit<<CED_LAYER_SHIFT,2,0x7cf774);
		}
	    
	    }
	    catch(DataNotAvailableException &e) {}
	    ced_send_event();
	}
    }


// Drawing Calorimeter Hits
    if (_layerCaloHit >= 0) {
      int tothits = 0;
	for ( unsigned int icol = 0; 
	      icol < _caloHitCollections.size(); ++ icol) {
	    try {
		LCCollection * col = evt->getCollection(_caloHitCollections[icol].c_str()) ;
		int nelem = col->getNumberOfElements();
		for (int ielem = 0; ielem < nelem; ++ielem) {
		  CalorimeterHit * hit = 
		    dynamic_cast<CalorimeterHit*>( col->getElementAt(ielem));
		  float x = (float)hit->getPosition()[0];
		  float y = (float)hit->getPosition()[1];
		  float z = (float)hit->getPosition()[2];
		  ced_hit(x,y,z, _layerCaloHit<<CED_LAYER_SHIFT,2,0xFFFFFF);
		  tothits++;
		}	    
	    }
	    catch(DataNotAvailableException &e) {}
	    ced_send_event();
	}      
	std::cout << "Total number of calo hits : " << tothits << std::endl;
    }


// Drawing Reconstructed Particles
    if (_layerReco >= 0) {
      try {
	LCCollection * col = evt->getCollection(_particleCollection.c_str());
	int nelem = col->getNumberOfElements();

	float TotEn = 0.0;
	float TotPX = 0.0;
	float TotPY = 0.0;
	float TotPZ = 0.0;

	for (int ip(0); ip < nelem; ++ip) {
	  ReconstructedParticle * part = 
	    dynamic_cast<ReconstructedParticle*>
	    (col->getElementAt(ip));
	  
	  TrackVec trackVec = part->getTracks();
	  int nTracks =  (int)trackVec.size();
	  ClusterVec clusterVec = part->getClusters();
	  int nClusters = (int)clusterVec.size();

	  float ene = part->getEnergy();
	  float px  = (float)part->getMomentum()[0];
	  float py  = (float)part->getMomentum()[1];
	  float pz  = (float)part->getMomentum()[2];
	  int type = (int)part->getType();

	  TotEn += ene;
	  TotPX += px;
	  TotPY += py;
	  TotPZ += pz;

	  std::cout << "Particle : " << ip
		    << " type : " << type 
		    << " PX = " << px
		    << " PY = " << py
		    << " PZ = " << pz
		    << " E  = " << ene << std::endl;

	  int color = returnColor(ip);

	  if (nClusters > 0 ) {
	    Cluster * cluster = clusterVec[0];
	    CalorimeterHitVec hitvec = cluster->getCalorimeterHits();
	    int nHits = (int)hitvec.size();
	    for (int iHit = 0; iHit < nHits; ++iHit) {
	      CalorimeterHit * hit = hitvec[iHit];
	      float x = hit->getPosition()[0];
	      float y = hit->getPosition()[1];
	      float z = hit->getPosition()[2];	      
	      ced_hit(x,y,z,_layerReco<<CED_LAYER_SHIFT,2,color);
	    }
	  }

	  if (nTracks > 0 ) {
	    Track * track = trackVec[0];
	    TrackerHitVec hitvec = track->getTrackerHits();
	    std::cout << "track parameters --->" << std::endl;
	    std::cout << "d0 = " << track->getD0()
		      << " ; z0 = " << track->getZ0() 
		      << " ; omega = " << track->getOmega() 
		      << " ; tanlam = " << track->getTanLambda()
		      << " ; tan(phi0) = " << tan(track->getPhi()) 
		      << std::endl;
	    int nHits = (int)hitvec.size();
	    for (int iHit = 0; iHit < nHits; ++iHit) {
	      TrackerHit * hit = hitvec[iHit];
	      float x = (float)hit->getPosition()[0];
	      float y = (float)hit->getPosition()[1];
	      float z = (float)hit->getPosition()[2];	      
	      ced_hit(x,y,z,_layerReco<<CED_LAYER_SHIFT,2,color);
	    }
	  }
	  ced_send_event();
	}

	std::cout << std::endl;
	std::cout << "Total Energy and Momentum Balance of Event" << std::endl;
	std::cout << "Energy = " << TotEn
		  << " PX = " << TotPX
		  << " PY = " << TotPY
		  << " PZ = " << TotPZ << std::endl;
	std::cout << std::endl;

	  
      }
      catch( DataNotAvailableException &e){}
    }


// Drawing True Clusters 
    if (_layerTrueClusters >= 0) {
	try {
	    LCCollection * col = 
	      evt->getCollection(_trueClustersCollection.c_str());
	    int nelem = col->getNumberOfElements();
	    std::cout << std::endl;
	    std::cout << "++++++++++++++++++++++++++++++++++++++++" 
		      << std::endl;
	    std::cout << "      DISPLAYING TRUE CLUSTERS " << std::endl;
	    std::cout << "++++++++++++++++++++++++++++++++++++++++" 
		      << std::endl;
	    std::cout << "Number of True Clusters : " << nelem << std::endl;
	    for (int iclust(0); iclust < nelem; ++iclust) {
		Cluster * cluster = 
		  dynamic_cast<Cluster*>(col->getElementAt(iclust));
		const LCCollection * relcol;
		try {
		    relcol = evt->getCollection("TrueClusterToMCP");	
		}
		catch(DataNotAvailableException &e){
		    std::cout << 
		      "No relation exist between True Clusters and MCPs" 
			      << std::endl; 	
		    std::cout << "Quitting program" << std::endl; 
		    exit(1);
		}
		LCRelationNavigator navigate(relcol);
      		LCObjectVec objectVec = navigate.getRelatedToObjects(cluster);
		MCParticle * mcp = NULL;  
		int iPDG = 0;
		int color = 0;
		if (objectVec.size() > 0) {
		    mcp = dynamic_cast<MCParticle*> (objectVec[0]);
		    iPDG = mcp->getPDG();
		    color = _mcpList[mcp];
		}
		CalorimeterHitVec hitvec = cluster->getCalorimeterHits();
		int nhits = (int)hitvec.size();
		std::cout << "True Cluster : "<< iclust 
			  << " # hits : " << nhits 
			  << " PDG  : " << iPDG << std::endl;  		
		for (int ihit(0); ihit < nhits; ++ihit) {
		    CalorimeterHit * calhit = hitvec[ihit];
		    float x = calhit->getPosition()[0];
		    float y = calhit->getPosition()[1];
		    float z = calhit->getPosition()[2];
		    int kcol = returnColor(color);
		    ced_hit(x,y,z,_layerTrueClusters<<CED_LAYER_SHIFT,2,kcol);
		} 		    
	    }
	}
	catch(DataNotAvailableException &e){}	
	ced_send_event();
    }


// Drawing Reconstructed Clusters
    if (_layerClusters >= 0) {
	try {
	    LCCollection * col = 
	      evt->getCollection(_clustersCollection.c_str());
	    int nelem = col->getNumberOfElements();
	    for (int iclust(0); iclust < nelem; ++iclust) {
		Cluster * cluster = 
		  dynamic_cast<Cluster*>(col->getElementAt(iclust));
		CalorimeterHitVec hitvec = cluster->getCalorimeterHits();
		int nhits = (int)hitvec.size();
		for (int ihit(0); ihit < nhits; ++ihit) {
		    CalorimeterHit * calhit = hitvec[ihit];
		    float x = calhit->getPosition()[0];
		    float y = calhit->getPosition()[1];
		    float z = calhit->getPosition()[2];
		    int color = returnColor(iclust);
		    ced_hit(x,y,z,_layerClusters<<CED_LAYER_SHIFT,2,color);
		} 		    
	    }
	}
	catch(DataNotAvailableException &e){}	
	ced_send_event();
    }

// Drawing True Tracks
    if (_layerTrueTracks >= 0) {
	try {
	    LCCollection * col = 
	      evt->getCollection(_trueTracksCollection.c_str());
	    int nelem = col->getNumberOfElements();
	    std::cout << std::endl;
	    std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
	    std::cout << "      DISPLAYING TRUE TRACKS " << std::endl;
	    std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
	    std::cout << "Number of True Tracks : " << nelem << std::endl;
	    for (int iclust(0); iclust < nelem; ++iclust) {
		Track * track = 
		    dynamic_cast<Track*>(col->getElementAt(iclust));
		const LCCollection * relcol;
		try {
		    relcol = evt->getCollection("TrueTrackToMCP");	
		}
		catch(DataNotAvailableException &e){
		    std::cout << "No relation exist between Track and MCP" 
			      << std::endl; 	
		    std::cout << "Quitting program" << std::endl; 
		    exit(1);
		}
		LCRelationNavigator navigate(relcol);
      		LCObjectVec objectVec = navigate.getRelatedToObjects( track );
		MCParticle * mcp = NULL;  
		int iPDG = 0;
		int color = 0;
		if (objectVec.size() > 0) {
		    mcp = dynamic_cast<MCParticle*> (objectVec[0]);
		    iPDG = mcp->getPDG();
		    color = _mcpList[mcp];
		}
		TrackerHitVec hitvec = track->getTrackerHits();
		int nhits = (int)hitvec.size();
		std::cout << "True Track : "<< iclust 
			  << " # hits : " << nhits 
			  << " PDG  : " << iPDG << std::endl;  
		float * ah = new float[nhits];
		float * xh = new float[nhits];
		float * yh = new float[nhits];
		float * zh = new float[nhits];
		float zmin = 1.0E+10;
		float zmax = -1.0E+10;
		for (int ihit(0); ihit < nhits; ++ihit) {
		    TrackerHit * hit = hitvec[ihit];
		    float x = (float)hit->getPosition()[0];
		    float y = (float)hit->getPosition()[1];
		    float z = (float)hit->getPosition()[2];
		    int kcol = returnColor(color);
		    ced_hit(x,y,z, _layerTrueTracks<<CED_LAYER_SHIFT,2,kcol);
		    ah[ihit] = 1.0;
		    xh[ihit] = x;
		    yh[ihit] = y;
		    zh[ihit] = z;
		    if (z < zmin)
			zmin = z;
		    if (z > zmax)
			zmax = z;
		} 	
		ClusterShapes * shapes = new ClusterShapes(nhits,ah,xh,yh,zh);
		float dz = (zmax - zmin) / 500.;
		float par[5];
		float dpar[5];
		float chi2;
		float distmax;
		shapes->FitHelix(500, 0, 1, par, dpar, chi2, distmax);
		float x0 = par[0];
		float y0 = par[1];
		float r0 = par[2];
		float bz = par[3];
		float phi0 = par[4];
		if (chi2 > 0. && chi2 < 10.) {
		  for (int iz(0); iz < 500; ++iz) {
		    float z1 = zmin + iz*dz;
		    float z2 = z1 + dz;
		    float x1 = x0 + r0*cos(bz*z1+phi0);
		    float y1 = y0 + r0*sin(bz*z1+phi0);
		    float x2 = x0 + r0*cos(bz*z2+phi0);
		    float y2 = y0 + r0*sin(bz*z2+phi0);
		    ced_line(x1,y1,z1,x2,y2,z2,
			     _layerTrueTracks<<CED_LAYER_SHIFT,2,0xFFFFFF);
		  }
		}
		ced_send_event();
		delete shapes;
		delete[] xh;
		delete[] yh;
		delete[] zh;
		delete[] ah;
	    }
	}
	catch(DataNotAvailableException &e){}	
	ced_send_event();
    }


// Drawing Reconstructed Tracks
    if (_layerTracks >= 0) {
	try {
	    LCCollection * col = evt->getCollection(_tracksCollection.c_str());
	    int nelem = col->getNumberOfElements();
	    std::cout << "Number of Tracks : " << nelem << std::endl;
	    for (int iclust(0); iclust < nelem; ++iclust) {
		Track * track = 
		  dynamic_cast<Track*>(col->getElementAt(iclust));
		TrackerHitVec hitvec = track->getTrackerHits();
		int nhits = (int)hitvec.size();
		std::cout << iclust << " : # hits = " << nhits << std::endl;
		std::cout << " Track parameters ---> " << std::endl; 
		std::cout << "d0 = " << track->getD0()
			  << " ; z0 = " << track->getZ0() 
			  << " ; omega = " << track->getOmega() 
			  << " ; tanlam = " << track->getTanLambda()
			  << " ; tan(phi0) = " << tan(track->getPhi()) 
			  << std::endl;
		float * ah = new float[nhits];
		float * xh = new float[nhits];
		float * yh = new float[nhits];
		float * zh = new float[nhits];
		float zmin = 1.0E+10;
		float zmax = -1.0E+10;
		for (int ihit(0); ihit < nhits; ++ihit) {
		    TrackerHit * hit = hitvec[ihit];
		    float x = (float)hit->getPosition()[0];
		    float y = (float)hit->getPosition()[1];
		    float z = (float)hit->getPosition()[2];
		    int color = returnColor(iclust);
		    ced_hit(x,y,z,_layerTracks<<CED_LAYER_SHIFT,2,color);
		    ah[ihit] = 1.0;
		    xh[ihit] = x;
		    yh[ihit] = y;
		    zh[ihit] = z;
		    if (z < zmin)
			zmin = z;
		    if (z > zmax)
			zmax = z;
		} 	
		ClusterShapes * shapes = new ClusterShapes(nhits,ah,xh,yh,zh);
		float zBegin, zEnd;
		if (fabs(zmin)<fabs(zmax)) {
		  zBegin = zmin;
		  zEnd   = zmax;
		}
		else {
		  zBegin = zmax;
		  zEnd   = zmin;
		}
		float signPz = zEnd - zBegin;		
		float dz = (zmax - zmin) / 500.;
		float par[5];
		float dpar[5];
		float chi2;
		float distmax;
      		shapes->FitHelix(500, 0, 1, par, dpar, chi2, distmax);
		float x0 = par[0];
		float y0 = par[1];
		float r0 = par[2];
		float bz = par[3];
		float phi0 = par[4];
		std::cout << "Chi2 Of Fit " << chi2 << std::endl;
		HelixClass * helix = new HelixClass();
		helix->Initialize_BZ(x0, y0, r0, 
				     bz, phi0, _bField,signPz,
				     zBegin);
		std::cout << "d0 = " << helix->getD0()
			  << " ; z0 = " << helix->getZ0() 
			  << " ; omega = " << helix->getOmega() 
			  << " ; tanlam = " << helix->getTanLambda()
			  << " ; tan(phi0) = " << tan(helix->getPhi0()) 
			  << std::endl;		

		/*		for (int ihit = 0; ihit < nhits; ++ihit) {
		  float xPoint[3];
		  xPoint[0] = xh[ihit];
		  xPoint[1] = yh[ihit];
		  xPoint[2] = zh[ihit];
		  float Distance[3];
		  helix->getDistanceToPoint(xPoint,Distance);
		  std::cout << ihit 
			    << " " << Distance[0]
			    << " " << Distance[1]
			    << " " << Distance[2]
			    << std::endl;
		}*/

		ced_send_event();
		if (chi2 > 0. && chi2 < 10.) {
		  for (int iz(0); iz < 500; ++iz) {
		    float z1 = zmin + iz*dz;
		    float z2 = z1 + dz;
		    float x1 = x0 + r0*cos(bz*z1+phi0);
		    float y1 = y0 + r0*sin(bz*z1+phi0);
		    float x2 = x0 + r0*cos(bz*z2+phi0);
		    float y2 = y0 + r0*sin(bz*z2+phi0);			
		    ced_line(x1,y1,z1,x2,y2,z2,_layerTrueTracks<<CED_LAYER_SHIFT,1,0xFFFFFF);
		    
		  }		
		}
		delete shapes;
		delete helix;
		delete[] xh;
		delete[] yh;
		delete[] zh;
		delete[] ah;
	    }
	}
	catch(DataNotAvailableException &e){}	
	ced_send_event();
    }

    getchar();
    _nEvt++;

}


void GenericViewer::check( LCEvent * evt ) { }
  
void GenericViewer::end(){ } 

int GenericViewer::returnColor(int counter) {

	int icol =  counter % 16;
	int kcol;
	if (icol==0)  kcol = 0x00ff00;
	if (icol==1)  kcol = 0xAA00ff;
	if (icol==2)  kcol = 0xff0000;
	if (icol==3)  kcol = 0x00ffff;
	if (icol==4)  kcol = 0xffff00;
	if (icol==5)  kcol = 0xff00ff;
	if (icol==6)  kcol = 0xffffff;
	if (icol==7)  kcol = 0x0fffff;
	if (icol==8)  kcol = 0x000088;
	if (icol==9)  kcol = 0x008800;
	if (icol==10) kcol = 0x880000;
	if (icol==11) kcol = 0x008888;
	if (icol==12) kcol = 0x888800;
	if (icol==13) kcol = 0x880088;
	if (icol==14) kcol = 0x888888;
	if (icol==15) kcol = 0x00A888;

	return kcol;

}
