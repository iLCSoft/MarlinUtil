#include <iostream>
#include <math.h>

/**
 *    Simple class with a histogram-like array without any display features
 *
 *    @authors O. Wendt (DESY)
 *    @version $Id: PseudoHistogram.h,v 1.1 2005-10-19 14:36:34 owendt Exp $
 *
 */
class PseudoHistogram {

 public:

  /**
   *    Constructor with number of bins and lower and upper boundaries
   */
  PseudoHistogram(int NOfBins, double min, double max);
  ~PseudoHistogram();
  
  /**
   *    Clear the content but leave the structure of the object, i.e. bins and boundaries
   */
  void clearContent();

  /**
   *    Fill a value x with the weight w to the pseudo-histogram
   */
  void fill(double x, double w);

  /**
   *    Returns content of bin
   */
  double getBinContent(int bin);
  
  /**
   *    Returns number of entries in bin
   */
  int getNumberOfEntries(int bin);

  /**
   *    Checks if bin is in range of the pseudo-histogram (over- and underflow bins are
   *    taken into account)
   */
  bool isInRange(int bin);

  /**
   *    Returns the weighted sum of the pseudo-histogram within startbin and
   *    endbin
   */
  double integral(int startbin, int endbin);

  /**
   *    Prints content of the pseudo-histogram on the standard output 
   */
  void printContent();


 private:
  int _FullNumberOfBins; // number of bins plus the over- and undeflow bin
  int _NumberOfBins; // number of bins without the over- and undeflow bin
  double _MinValue;
  double _MaxValue;
  double _BinWidth;
  int* _NOfEntries;
  double* _Content;
  double* _UpperIntervalLimit;

};
