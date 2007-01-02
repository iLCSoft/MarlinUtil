
#include "MarlinDrawUtil.h"
#include <cmath>
 


// ____________________________________________________________________________________________________


int MarlinDrawUtil::getColor(int pdgCode) {
  
  switch (abs(pdgCode)) {
    
    // CHARGED OBJECTS:

    // leptons
  case 11   : return 0x2bffb8; // electron : turquoise = 0x2bffb8, similar color for corresponding neutrinos
  case 13   : return 0x0cff2c; // muon     : green     = 0x0cff2c, similar color for corresponding neutrinos
  case 15   : return 0x122aff; // tauon    : blue      = 0x122aff, similar color for corresponding neutrinos


    // charged hadrons --> refine
    
    // baryons:
  case 2212 : return 0xff5765; // proton : redpink       = 0xff5765
  case 3222 : return 0xde9698; // sigma+ : light redpink = 0xde9698
  case 3112 : return 0xde9698; // sigma- : light redpink = 0xde9698


    // mesons:
  case 211  : return 0xff7402; // Pi+/-  : orange  = 0xff7402
  case 321  : return 0xd500ff; // K+/-   : violett = 0xd500ff

    
    // gauge bosons:





    
    // NEUTRAL OBJECTS:

    // leptons    
  case 12   : return 0x24d59a; // nue  : turquoise    = 0x24d59a, similar color for corresponding charged leptons
  case 14   : return 0x0acd24; // numu : green        = 0x0acd24, similar color for corresponding charged leptons
  case 16   : return 0x0c1db1; // nue  : blue         = 0x0c1db1, similar color for corresponding charged leptons


    // neutral hadrons --> refine

    // baryons:
  case 2112 : return 0xff0627; // neutron : red = 0xff0627
  case 3122 : return 0x8c84de; // lambda0 : dark light blue = 0x8c84de
  case 3212 : return 0x84c9de; // sigma0  : dark light turquoise = 0x84c9de
  case 3322 : return 0x95de88; // Xi0     : dark light green = 0x95de88 
    

    // mesons:
  case 111  : return 0xf6a8f7; // Pi0 : pink         = 0xf6a8f7
  case 130  : return 0xd500ff; // K0L : violett      = 0xd500ff
  case 310  : return 0xb000d3; // K0S : dark violett = 0xb000d3
  case 311  : return 0xb000d3; // K0  : dark violett = 0xb000d3


    // gauge bosons:
  case 22   : return 0xf9f920; // gamma : yellow = 0xf9f920
 

    
    // DEFAULT

  default   : return 0xffffff; // default : white  = 0xffffff


  }

}


// ____________________________________________________________________________________________________


int MarlinDrawUtil::getColorAmplitude(float amplitude, float max_amplitude,std::string mode, float limit) {

  float ratio = amplitude/max_amplitude;
  int r = 0;
  int g = 0;
  int b = 0;
  int color = 0;

  const int color_categories = 5;
  const float color_categories_limit = limit/(float)color_categories;


  if (mode == "rainbow") {
    if ( (0.0<=ratio)&&(ratio<=color_categories_limit) ) {
      r = 255;
      g = (int)floor(255*( ratio/0.02 )); //change floor
      b = 0;
      color = (r<<16) + (g<<8) + b;
      return color;
    }
    else if ( (color_categories_limit<ratio)&&(ratio<=2*color_categories_limit) ) {
      r = (int)(255 - floor(255*( (ratio - 0.02)/(0.04 - 0.02) ))); //change floor
      g = 255;
      b = 0;
      color = (r<<16) + (g<<8) + b;
      return color;
    }
    else if ( (2*color_categories_limit<ratio)&&(ratio<=3*color_categories_limit) ) {
      r = 0;
      g = 255;
      b = (int)floor(255*( (ratio - 0.04)/(0.06 - 0.04) )); //change floor 
      color = (r<<16) + (g<<8) + b;
      return color;
    }
    else if ( (3*color_categories_limit<ratio)&&(ratio<=4*color_categories_limit) ) {
      r = 0;
      g = (int)(255 - floor(255*( (ratio - 0.06)/(0.08 - 0.06) ))); //change floor 
      b = 255;
      color = (r<<16) + (g<<8) + b;
      return color;
    }
    else if ( (4*color_categories_limit<ratio)&&(ratio<=5*color_categories_limit) ) {
      r = (int)floor(255*( (ratio - 0.08)/(0.10 - 0.08) )); //change floor
      g = 0;
      b = 255;
      color = (r<<16) + (g<<8) + b;
      return color;
    }    
    else {
      r = 255;
      g = 255;
      b = 255;
      color = (r<<16) + (g<<8) + b;     
      return color;
    }
  }

  else return 0;

}

// ____________________________________________________________________________________________________

