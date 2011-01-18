
#include "MarlinDrawUtil.h"
#include <cmath>

//SJA:FIXED:added to make gcc4.3 compliant
#include <cstdlib>




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


unsigned int MarlinDrawUtil::getColorAmplitude(float amplitude, float max_amplitude,std::string mode, float limit) {


  const unsigned int color_categories = 4; // in rgb there are 4 steps to scan a rainbow spectrum from red to blue
  const unsigned int color_category_width = 256;

  float ratio = amplitude/max_amplitude;

  
  if (limit < 0.0) limit = 0.0;
  if (limit > 1.0) limit = 1.0;

  
  unsigned int ratioInt = (unsigned int)(color_categories*color_category_width*(ratio/limit));


  // debug  
  // std::cout << "ratio = " << ratio << "  " << ratioInt << std::endl;


  unsigned int r = 0;
  unsigned int g = 0;
  unsigned int b = 0;
  unsigned int color = 0;




  if (mode == "rainbow") {
    if ( (0 <= ratioInt) && (ratioInt < 1*color_category_width) ) {
      r = 255;
      g = 0 + ratioInt;
      b = 0;
      color = (r<<16) + (g<<8) + b;
      
      // debug
      // std::cout << r << "  " << g << "  " <<  b << std::endl;
      
      return color;
    }
    else if ( (1*color_category_width <= ratioInt) && (ratioInt < 2*color_category_width ) ) {
      r = 255 - (ratioInt - 1*color_category_width);
      g = 255;
      b = 0;
      color = (r<<16) + (g<<8) + b;

      // debug
      // std::cout << r << "  " << g << "  " <<  b << std::endl;

      return color;
    }
    else if ( (2*color_category_width <= ratioInt) && (ratioInt < 3*color_category_width) ) {
      r = 0;
      g = 255;
      b = 0 + (ratioInt - 2*color_category_width);
      color = (r<<16) + (g<<8) + b;

      // debug
      // std::cout << r << "  " << g << "  " <<  b << std::endl;

      return color;
    }
    else if ( (3*color_category_width <= ratioInt) && (ratioInt < 4*color_category_width) ) {
      r = 0;
      g = 255 - (ratioInt - 3*color_category_width);
      b = 255;
      color = (r<<16) + (g<<8) + b;

      // debug
      // std::cout << r << "  " << g << "  " <<  b << std::endl;

      return color;
    }
    else {
      r = 255;
      g = 255;
      b = 255;
      color = (r<<16) + (g<<8) + b;     
      
      // debug
      // std::cout << r << "  " << g << "  " <<  b << std::endl;
      
      return color;
    }
  }





   

  // old version calculating in floating point representation
  /*     
  const float color_categories_limit = limit/(float)color_categories;

  if (mode == "rainbowFloat") {
    if ( (0.0<=ratio)&&(ratio<=color_categories_limit) ) {
      r = 255;
      g = (unsigned int)floor(255.*color_categories*(ratio-0.0/color_categories)) + 0; //change floor
      b = 0;
      color = (r<<16) + (g<<8) + b;
      
      // debug
      std::cout << r << "  " << g << "  " <<  b << std::endl;
      
      return color;
    }
    else if ( (color_categories_limit<ratio)&&(ratio<=2*color_categories_limit) ) {
      r = (unsigned int)floor(-255.*color_categories*(ratio-1.0/color_categories)) + 255; //change floor
      g = 255;
      b = 0;
      color = (r<<16) + (g<<8) + b;

      // debug
      std::cout << r << "  " << g << "  " <<  b << std::endl;

      return color;
    }
    else if ( (2*color_categories_limit<ratio)&&(ratio<=3*color_categories_limit) ) {
      r = 0;
      g = 255;
      b = (unsigned int)floor(255.*color_categories*(ratio-2.0/color_categories)) + 0; //change floor 
      color = (r<<16) + (g<<8) + b;

      // debug
      std::cout << r << "  " << g << "  " <<  b << std::endl;

      return color;
    }
    else if ( (3*color_categories_limit<ratio)&&(ratio<=4*color_categories_limit) ) {
      r = 0;
      g = (int)floor(-255.*color_categories*(ratio-3.0/color_categories)) + 255; //change floor 
      b = 255;
      color = (r<<16) + (g<<8) + b;

      // debug
      std::cout << r << "  " << g << "  " <<  b << std::endl;

      return color;
    }
    else {
      r = 255;
      g = 255;
      b = 255;
      color = (r<<16) + (g<<8) + b;     
      
      // debug
      std::cout << r << "  " << g << "  " <<  b << std::endl;
      
      return color;
    }
  }
  */
  



  else return 0;

}

// ____________________________________________________________________________________________________

