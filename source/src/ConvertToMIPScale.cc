#include "ConvertToMIPScale.h"


using namespace lcio;
using namespace marlin;



ConvertToMIPScale aConvertToMIPScale;



ConvertToMIPScale::ConvertToMIPScale() : Processor("ConvertToMIPScale") {
  
  _description = "transforms the energy amplitude of each calorimeter hit passing the cut into the MIP scale";
  

  registerInputCollection( LCIO::CALORIMETERHIT, 
			    "ECALInputCollectionName", 
			    "ECAL Input Collection Name",
			    _inputEcalCollection,
			    std::string("ECAL") );
  
  
  registerInputCollection( LCIO::CALORIMETERHIT, 
			    "HCALInputCollectionName", 
			    "HCAL Input Collection Name",
			    _inputHcalCollection,
			    std::string("HCAL") );
  

  registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "ECALOutputCollectionName",
			    "ECAL Output Collection Name, stores calorimeter hits with amplitudes in MIP energy scale",
			    _outputEcalCollection,
			    std::string("ECAL_MIP") ); 


  registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "HCALOutputCollectionName",
			    "HCAL Output Collection Name, stores calorimeter hits with amplitudes in MIP energy scale",
			    _outputHcalCollection,
			    std::string("HCAL_MIP") ); 


  registerProcessorParameter("CutECAL", 
			     "Cut in MIPs on the amplitudes in the ECAL collection",
			     _cutEcal,
			     (float)0.0);


  registerProcessorParameter("CutHCAL", 
			     "Cut in MIPs on the amplitudes in the HCAL collection",
			     _cutHcal,
			     (float)0.0);



  std::vector<float> mipCoeffEcal;
  mipCoeffEcal.push_back(1.0);
  mipCoeffEcal.push_back(1.0);

  registerProcessorParameter("MIPCoeffEcal", 
			     "Coefficients for the MIP calibration in the ECAL in GeV/MIP",
			     _mipCoeffEcal,
			     mipCoeffEcal);

  

  std::vector<float> mipCoeffHcal;
  mipCoeffHcal.push_back(1.0);

  registerProcessorParameter("MIPCoeffHcal",
			     "Coefficients for the MIP calibration in the HCAL in GeV/MIP",
			     _mipCoeffHcal,
			     mipCoeffHcal);


}


void ConvertToMIPScale::init() {
  
    _nRun = -1;
    _nEvt = 0;
    
}



void ConvertToMIPScale::processRunHeader( LCRunHeader* ) {

  ++_nRun;
  _nEvt = 0;

} 



void ConvertToMIPScale::processEvent( LCEvent * evt ) { 
    

  LCCollectionVec* ecalColMIPScale = new LCCollectionVec(LCIO::CALORIMETERHIT);
  LCCollectionVec* hcalColMIPScale = new LCCollectionVec(LCIO::CALORIMETERHIT);


  LCFlagImpl flag;

  flag.setBit(LCIO::CHBIT_LONG); // CHBIT_LONG = 31, stores position of each LCIO::CALORIMETERHIT
  
  ecalColMIPScale->setFlag(flag.getFlag());
  hcalColMIPScale->setFlag(flag.getFlag());



  // convert  ECAL hits to MIP scale
  
  try {

    LCCollection* col = evt->getCollection(_inputEcalCollection.c_str());
    
    unsigned int n = col->getNumberOfElements();

    for (unsigned int j = 0; j < n; ++j) {

      CalorimeterHit* hit = dynamic_cast<CalorimeterHit*>( col->getElementAt(j) );

      if ( hit->getEnergy() >= _cutEcal ) {
	
	
	// type 0 : ECAL sampling 1
	// type 1 : ECAL sampling 2
	// type 2 : HCAL
	if ( hit->getType() == 0 ) {

	  CalorimeterHitImpl* hitMIPScale = new CalorimeterHitImpl();

	  hitMIPScale->setCellID0(hit->getCellID0());
	  hitMIPScale->setCellID1(hit->getCellID1());

	  hitMIPScale->setTime(hit->getTime());
	  hitMIPScale->setPosition(hit->getPosition());
	  hitMIPScale->setType(hit->getType());

	  hitMIPScale->setEnergy( (hit->getEnergy())/(_mipCoeffEcal.at(0)) );

	  ecalColMIPScale->addElement(hitMIPScale);

	}
	else if ( hit->getType() == 1 ) {

	  CalorimeterHitImpl* hitMIPScale = new CalorimeterHitImpl();

	  hitMIPScale->setCellID0(hit->getCellID0());
	  hitMIPScale->setCellID1(hit->getCellID1());

	  hitMIPScale->setTime(hit->getTime());
	  hitMIPScale->setPosition(hit->getPosition());
	  hitMIPScale->setType(hit->getType());

	  hitMIPScale->setEnergy( (hit->getEnergy())/(_mipCoeffEcal.at(1)) );

	  ecalColMIPScale->addElement(hitMIPScale);

	}
	else {
	
	  std::cout << "Warning: CalorimeterHit type in ECAL collection " << _inputEcalCollection << "not set properly. No MIP conversion will be done for this hit." 
		    << _nEvt << std::endl;
	  
	}

      }
      
    }
    
  }
  catch(DataNotAvailableException &e){std::cout << "Collection " << _inputEcalCollection << " not found in LCEvent " << _nEvt << std::endl; };
  

  


  // convert  HCAL hits to MIP scale
  
  try {

    LCCollection* col = evt->getCollection(_inputHcalCollection.c_str());
    
    unsigned int n = col->getNumberOfElements();

    for (unsigned int j = 0; j < n; ++j) {

      CalorimeterHit* hit = dynamic_cast<CalorimeterHit*>( col->getElementAt(j) );

      if ( hit->getEnergy() >= _cutHcal ) {
	
	
	// type 0 : ECAL sampling 1
	// type 1 : ECAL sampling 2
	// type 2 : HCAL
	if ( hit->getType() == 2 ) {
	  
	  CalorimeterHitImpl* hitMIPScale = new CalorimeterHitImpl();

	  hitMIPScale->setCellID0(hit->getCellID0());
	  hitMIPScale->setCellID1(hit->getCellID1());

	  hitMIPScale->setTime(hit->getTime());
	  hitMIPScale->setPosition(hit->getPosition());
	  hitMIPScale->setType(hit->getType());

	  hitMIPScale->setEnergy( (hit->getEnergy())/(_mipCoeffHcal.at(0)) );

	  hcalColMIPScale->addElement(hitMIPScale);

	}
	else {
	
	  std::cout << "Warning: CalorimeterHit type in HCAL collection " << _inputHcalCollection << "not set properly. No MIP conversion will be done for this hit." 
		    << _nEvt << std::endl;
	  
	}

      }
      
    }
    
  }
  catch(DataNotAvailableException &e){std::cout << "Collection " << _inputHcalCollection << " not found in LCEvent " << _nEvt << std::endl; };





  

  evt->addCollection(ecalColMIPScale,_outputEcalCollection.c_str());
  evt->addCollection(hcalColMIPScale,_outputHcalCollection.c_str());

  ++_nEvt;

}



void ConvertToMIPScale::check( LCEvent* ) { }


  
void ConvertToMIPScale::end(){ } 
