#ifndef MarlinDrawUtil_h
#define MarlinDrawUtil_h 1

#include <iostream>
#include <string>

#include <lcio.h>

using namespace lcio;


class MarlinDrawUtil {

 public:

  static int getColor(int pdgCode);
  static unsigned int getColorAmplitude(float amplitude, float max_amplitude, std::string mode, float limit);

};

#endif
