/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "FPCCDPixelHit.h"
#include "FPCCDData.h"

#include <IMPL/LCGenericObjectImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/SimTrackerHit.h>
#include <UTIL/LCTOOLS.h>
#include <iostream>

// =====================================================================
FPCCDData::FPCCDData(int maxlayer, int maxladder): _maxlayer(maxlayer), _maxladder(maxladder)
{
  //std::cout << "***FPCCDData class: constructor!**** " << std::endl; 
  //std::cout << "_maxlayer = " << _maxlayer << std::endl; 
  //std::cout << "_maxladder = " << _maxladder << std::endl; 
  //The above check shows that Both _maxlayer and _maxladder are always same, 6 and 17.

  _pxHits.resize(_maxlayer); //--> "resize" method belongs to std::vector. This changes the number of elements. In this case, resized to 6 elements. (_maxlayer is 6 by default.)
  //_pxHits is defined in FPCCDData.h. _pxHits is PixelDataBuf_t, namely std::vector< std::vector<std::map<unsigned int, FPCCDPixelHit*>> > .
  //This structure is matrix whose elements are FPCCDPixelHit*.
  
  for(int i=0;i<_maxlayer;i++){ _pxHits[i].resize(_maxladder); }
  //--> 2nd elements are resized. By default, _maxlayer is 17 regardless that each layer has different number of ladders. 
  //Namely, when you do with 1st layer, 2nd elements has some value until 10th and has unknown values from 11th to 17th.
}


// ===================================================================
void FPCCDData::clear()
{
  for(int layer=0;layer<_maxlayer;layer++){
    for(int ladder=0;ladder<_maxladder;ladder++){
            PixelHitMap_t::iterator it = _pxHits[layer][ladder].begin();
      while( it != _pxHits[layer][ladder].end()){
        delete (*it).second;
        it++ ;
      }
      _pxHits[layer][ladder].clear();
    }
  }
}

// ===================================================================
void FPCCDData::dump()
{
  std::cout << "FPCCDData::dump() prints out all Pixel hits" << std::endl;
  for(int layer=0;layer<_maxlayer;layer++) {
    for(int ladder=0;ladder<_maxladder;ladder++) {
      if( _pxHits[layer][ladder].size() > 0 ) {
        PixelHitMap_t::iterator it=_pxHits[layer][ladder].begin();
        while( it != _pxHits[layer][ladder].end() ) {
          (*it).second->print();
          it++;
        }
      }
    }
  }
}

// ====================================================================
void FPCCDData::addPixelHit(FPCCDPixelHit &aHit, bool isSignal)
{
  int layer=aHit.getLayerID();
  int ladder=aHit.getLadderID();

  unsigned int hitid=aHit.encodeCellWord();
  FPCCDPixelHit::HitQuality_t addedQuality;
  { isSignal ? addedQuality = FPCCDPixelHit::kSingle : addedQuality = FPCCDPixelHit::kBKG ; }  

  PixelHitMap_t::iterator it=_pxHits[layer][ladder].find(hitid);//map (or multimap) returns end if not found.
  // Append ADC if same pixelID exist in the map
  if( it != _pxHits[layer][ladder].end() ) {
    FPCCDPixelHit *findhit=(*it).second;  
    //    addedQuality = aHit.getQuality();
    findhit->addPixelHit(aHit, addedQuality); //Don't regard this addPixelHit as a recursive method. This method comes from FPCCDPixelHit.
  }
  // Add new hit if not exist in the map
  else {
    aHit.setQuality(addedQuality);
    _pxHits[layer][ladder].insert(PixelHitMap_t::value_type(hitid, new FPCCDPixelHit(aHit)));
  }
}

// ===================================================================
void FPCCDData::packPixelHits(EVENT::LCCollection &colvec)
{
  // This function writes all CCD data into hitvec.
  // _pxHits are cleared after copying to save space
  for(int layer=0;layer<_maxlayer;layer++) {
    for(int ladder=0;ladder<_maxladder;ladder++) {
      PixelHitMap_t::iterator it=_pxHits[layer][ladder].begin();//PixelHitMap_t is "std::map<unsigned int, FPCCDPixelHit*>".
      //The type of _pxHits[some number][some number] is PixelHitMap_t. This is OK.
      if( it != _pxHits[layer][ladder].end() ) {
        IMPL::LCGenericObjectImpl *out=new IMPL::LCGenericObjectImpl();
        unsigned int index=0;//In C++, initialization is repeated. So, in 2nd loop, this discription also occuers.
        // store layer number, ladder number in the first word
        // Most significant 4 bits (left 4 bits) should be kept unused 
        // for future use as a data format type
        int word0=( (layer & 0x000000FF ) << 8 | ( ladder & 0x000000FF ) );
        //mori: This process is bit operation. The followings is from FPCCDData.h
        /*************from FPCCDData.h******
          Format of LCGenericObject is as follows:
          PixelHits in one ladder are packed in one element.
          First 32 bit word contains layerID and ladderID of the element.
          Counting bit position from left to right, layerID and ladderID
          information are stored as follows. 
          bit 31-16: reserved for future use --> Mori [simthits address which makes layerID and ladderID of hit points should be put here.]
          bit 15-8 : layerID
          bit 7-0  : ladderID
          - Last words are for pixel hits in a ladder. Two word is used 
            for each pixel hit. 
            -- First word
            bit 31: 0 
            bit 30-29 : hit quality bit
            0 = hit is created by a single SimTrackerHit of 
            a signal event.
            1 = hit is created by multiple SimTrackerHits of 
            a signal event
            2 = hit is created by a sigle SimTrackerHit of 
            a signal event and hits by background events
            3 = hit is created by multiple SimTrackerHits 
            of a signal event and background events. 
            bit 28-16 : xiID 
            bit 15-0  : zetaID
            -- Second word
            dEdx value given by SimTrackerHit object.
            Float value is stored with a help of union.
        **********************************/
        //Namely, word0 = 0x0000LLDD (LL is layerID and DD is ladderID.)
        //I think layerID doesn't need 8bits, only needs 4bits(ultimately speaking 3bits but it is no useful for reading code. ). 
        //The type of word0 is 4bite, namely 32bits. The bit digit of 0x00000000 is same as the one of "int".
        //Namely, word0 has unused 16bits (but can take only minus value.). 

        out->setIntVal(index, word0 );//virtual void  setIntVal (unsigned index, int value) --> Sets the integer value at the given index. 
        index++;
        while( it !=_pxHits[layer][ladder].end() ) {
          FPCCDPixelHit *aHit=(*it).second;
          // First word, from left to right,
          // MSB =0 to indicate hitID word --> MSB means Most Significant Bit or Bite. Most left bit or bite is MSB. ex) MSB of 00110101 = 0 ( in bit ) and MSB of b189ff77 is b1 ( in bite ) 
          //     next 2 bit for quality
          //     next 13 + 16 bit for hit id ( xi and zeta )
          unsigned int hitid=(unsigned int)aHit->encodeCellWord();//encodeCellWord() is from FPCCDPixelHit. 
          //If already _xiID and _zetaID is set in the instance pointed by *aHit,
          //then it returns unsigned int number which has 32bits data (effectively, 29bits) of _xiID and _zetaID.  
          //When packing data, this extra 3bit is deleted. Don't put additional value in these 3bits. 
          unsigned int quality=(unsigned int)aHit->getQuality();//from FPCCDPixelHit.h. 
          //(define) HitQuality_t getQuality(){ return _quality; } 
          //(define) typedef enum { kSingle=0, kSignalOverlap=0x01, kBKGOverlap=0x02, kBKG=0x03} HitQuality_t;
          int hitwd= ( ( (unsigned int)quality << 29 & 0x60000000 ) |
                       ( (hitid & 0x1FFFFFFF ) ) );
          //0x60000000 is 1100....000 (number of 0 is 29.)
          //0x1FFFFFFF is 00010000....0000 (number of 0 from 1 is 28.)
          out->setIntVal(index++, hitwd );
          // 2nd word is edep
          union intfloat { float edep; int iedep; } edepout; //union is like struct. usage: edepout.edep and edpeout.iedep.

          edepout.edep=aHit->getEdep();//If already _edep is set, then return _edep.
          out->setIntVal(index++, edepout.iedep);
          
          it++;
             
          //*************important change is here. 20121204 Mori**************//   
          //New information, which simthit was used, is added into LCGenericObject here.
          int osize = aHit->getSizeOfOrderID();
          out->setIntVal(index++, osize);
          for(int oi = 0; oi < osize; oi++){ out->setIntVal(index++, aHit->getOrderID(oi)); }
          //**************************end**************************************//

        } // moving data in aHit to LCGenericObjectImpl
        colvec.addElement(out);  // add one element --> Adds the given element to (end of) the collection.
        // maybe like std::vector::push_back()
      }
    }
  }
  clear();//Contents of _pxHits is cleared.
}

// =====================================================
int FPCCDData::unpackPixelHits(EVENT::LCCollection &col)
{
// Convert pixelhit data in col to _pxHits 
// returns total number of PixelHits converted.

//In FPCCDClustering.cc, The pointer of LCCollection whose name is LCGenericObject is replaced for &col.  

  int nhits=0;
  int nElements=col.getNumberOfElements();

//   UTIL::LCTOOLS::printLCGenericObjects(col);

  clear();
  for(int ie=0;ie<nElements;ie++){
    int ig=0;
    EVENT::LCGenericObject *obj=dynamic_cast<EVENT::LCGenericObject*>(col.getElementAt(ie));
    unsigned int iw0=obj->getIntVal(ig);//Returns the integer value for the given index. 
    int layer=( ( iw0>>8 ) & 0x000000FF );
    int ladder= ( iw0 & 0x000000FF ) ;
    ig++;
    while( ig < obj->getNInt() ) { //obj->getNInt()  -->  Number of integer values stored in this object. 
      unsigned int iw=obj->getIntVal(ig);
      FPCCDPixelHit *aHit=new FPCCDPixelHit(layer, ladder);
      // Converting first word
      aHit->decodeCellWord(iw);// xiID and zetaID are set here. 
      unsigned int qwd=(iw & 0x60000000 ) >> 29 ;
      switch ( qwd ) {
      case 0 : aHit->setQuality(FPCCDPixelHit::kSingle); break;
      case 1 : aHit->setQuality(FPCCDPixelHit::kSignalOverlap) ; break ;
      case 2 : aHit->setQuality(FPCCDPixelHit::kBKGOverlap) ; break ;
      case 3 : aHit->setQuality((FPCCDPixelHit::HitQuality_t)(FPCCDPixelHit::kSignalOverlap|FPCCDPixelHit::kBKGOverlap)) ; break ;
      //case 3 --> first a pair of blaces is cast. second pair is casted value.
      }
      ig++;
      // converting second word ; 
      union intfloat { float edep ; int iedep ; } edepin ;
      edepin.iedep=obj->getIntVal(ig++);
      aHit->setEdep(edepin.edep);
       
      //************************important change 20121204 mori********************************//
      int osize = obj->getIntVal(ig++);//osize shows multiplicity of original simthits in one pixel hit.
      for(int oi = 0; oi < osize; oi++ ){ 
         aHit->setOrderID( obj->getIntVal(ig++) );
      }
      //***************************************************************************************//

      unsigned int hitid=aHit->encodeCellWord();
      _pxHits[layer][ladder].insert(PixelHitMap_t::value_type(hitid, aHit));
      nhits++;
//      aHit->print();
    }
  }
  return nhits;
}

// =====================================================
void FPCCDData::Add(FPCCDData &bgHit)
{
  FPCCDPixelHit *aHit;
  
  for(int i=0; i<_maxlayer; i++){  
    for(int j=0; j<_maxladder; j++){
      
      
      PixelHitMap_t::iterator it= bgHit.itBegin(i, j);
      while( it != bgHit.itEnd(i, j)){
        
        aHit=(*it).second;
        //std::cout << "Add is starting!" << std::endl;
	   //std::cout << "getQuality = " << aHit->getQuality() << std::endl;
	   bool isSignal = (aHit->getQuality() == FPCCDPixelHit::kSingle);
	   //std::cout << "isSignal = " << isSignal  << std::endl;
        FPCCDData::addPixelHit( *aHit, isSignal);
        
	   //FPCCDData::addPixelHit( *aHit, false);
        
        it++;
      }
      
    }
  }
}
