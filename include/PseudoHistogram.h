#include <iostream>
#include <math.h>

/**
 *    Simple class with a histogram-like array without any display features
 *
 *    @author O. Wendt (DESY)
 *    @version $Id: PseudoHistogram.h,v 1.3 2005-11-04 16:54:59 owendt Exp $
 *
 */
class PseudoHistogram {

 public:

  /**
   * Constructor with number of bins and lower and upper boundaries
   *
   * @param NOfBins : number of bins in the histogram. Over- and underflow bin are
   *                  created additionally.
   * @param min     : smallest value in the histogram
   * @param max     : largest value in the histogram
   */
  PseudoHistogram(int NOfBins, double min, double max);

  /**
   * Destructor
   */
  ~PseudoHistogram();
  
  /**
   * Clears the content but leaves the structure of the object, i.\ e.\ bins and 
   * boundaries
   */
  void clearContent();

  /**
   * Fill a value x with the weight w to the pseudo-histogram
   *
   * @param x : value to fill in the histogram
   * @param w : weight of the value x
   */
  void fill(double x, double w);

  /**
   * Returns content of bin
   *
   * @param bin : number of the bin
   */
  double getBinContent(int bin);
  
  /**
   * Returns number of entries in bin
   *
   * @param bin : number of the bin
   */
  int getNumberOfEntries(int bin);

  /**
   * Checks if bin is in range of the pseudo-histogram (over- and underflow bins are
   * taken into account)
   *
   * @param bin : number of the bin
   */
  bool isInRange(int bin);

  /**
   * Returns the weighted sum of the pseudo-histogram within startbin and
   * endbin
   *
   * @param startbin : number of the start bin for the integral
   * @param endbin   : number of the end bin for the integral
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
