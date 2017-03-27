#include "PseudoHistogram.h"

//SJA:FIXED:added to make gcc4.3 compliant
#include <cstdlib>



//=============================================================================

PseudoHistogram::PseudoHistogram(int NOfBins, double min, double max):

  _FullNumberOfBins(NOfBins + 2),
  _NumberOfBins(NOfBins),
  _MinValue(min),
  _MaxValue(max),
  _BinWidth((fabs(_MaxValue - _MinValue)) / _NumberOfBins),
  _NOfEntries(_FullNumberOfBins, 0),
  _Content(_FullNumberOfBins, 0.0),
  _UpperIntervalLimit(_FullNumberOfBins-1, 0.0)

{
  for (int i = 0; i < _FullNumberOfBins-1; ++i) {
    _UpperIntervalLimit[i] = _MinValue + i*_BinWidth;
  }

}

//=============================================================================

PseudoHistogram::~PseudoHistogram() {
}

//=============================================================================

void PseudoHistogram::clearContent() {

  for (int i = 0; i < _FullNumberOfBins; ++i) {
    _NOfEntries[i] = 0;
    _Content[i] = 0.0;
  }

}

//=============================================================================

void PseudoHistogram::fill(double x, double w) {

  if (x < _MinValue) {
    ++_NOfEntries[0];
    _Content[0] += w;  
  }
  else if (x > _MaxValue) {
    ++_NOfEntries[_FullNumberOfBins-1];
    _Content[_FullNumberOfBins-1] += w;  
  }
  else if (x == _MaxValue) {
    ++_NOfEntries[_FullNumberOfBins-2];
    _Content[_FullNumberOfBins-2] += w;  
  }
  else {
    int index = (int)floor((((double)_NumberOfBins)/
			    (_MaxValue - _MinValue))*(x - _MinValue)) + 1;
    ++_NOfEntries[index];
    _Content[index] += w;
  }

}

//=============================================================================

int PseudoHistogram::findBin(double x) {

  if (x < _MinValue) {
    return 0;
  }
  else if (x > _MaxValue) {
    return _FullNumberOfBins-1;
  }
  else if (x == _MaxValue) {
    return _FullNumberOfBins-2;
  }
  else {
  return (int)floor((((double)_NumberOfBins)/
            (_MaxValue - _MinValue))*(x - _MinValue)) + 1;
  }
}

//=============================================================================

double PseudoHistogram::getBinContent(int bin) {

  if ( isInRange(bin) ) {
    return _Content[bin];
  }
  else {
    std::cout << "Requested bin not in range of the 'PseudoHistogram'." << std::endl;
    return -1.0;
  }
  
}

//=============================================================================

int PseudoHistogram::getNumberOfEntries(int bin) {

  if ( isInRange(bin) ) { 
    return _NOfEntries[bin];
  }
  else {
    std::cout << "Requested bin not in range of the 'PseudoHistogram'." << std::endl;
    return -1;
  }

}

//=============================================================================

bool PseudoHistogram::isInRange(int bin) {

  if ( (bin >= 0) && (bin < _FullNumberOfBins) ) { 
    return true;
  }
  else return false;

}

//=============================================================================

double PseudoHistogram::integral(int startbin, int endbin) {

  if ( isInRange(startbin) && isInRange(endbin) ) {    
    double result = 0.0;
    for (int i = startbin; i < abs(endbin - startbin)+startbin+1; ++i) {
      result += getBinContent(i);
    }
    return result;
  }
  else {
    std::cout << "At least one requested bin is not in range of the 'PseudoHistogram'."
	      << std::endl;
    return -1.0;
  }
  
}

//=============================================================================

void PseudoHistogram::printContent() {

  std::cout << "bin" << "     " << "content" << "     " << "entries" << "     " 
	    << "interval" << std::endl;

  std::cout.precision(3);
  for (int i = 0; i < _FullNumberOfBins; ++i) {
    std::cout << i << "          " << _Content[i] << "          " << _NOfEntries[i] 
	      << "          ";
    if (i==0) {
      std::cout << "(-inf," << _UpperIntervalLimit[i] << "]" << std::endl;
    }
    else if (i==_FullNumberOfBins-1) {
      std::cout << "[" << _UpperIntervalLimit[i-1] << ",inf)" << std::endl;
    }
    else {
      std::cout << "[" << _UpperIntervalLimit[i-1] << "," <<  _UpperIntervalLimit[i]
		<< "]" << std::endl;
    }
  }
}

//=============================================================================
