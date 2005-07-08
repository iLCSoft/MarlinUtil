/**
 * A header file which defines the basics of a tpc
 */
#ifndef tpc_h
#define tpc_h 1

#include <vector>
#include <map>
#include <lcio.h>
#include <EVENT/SimTrackerHit.h>

using namespace std;
using namespace lcio;

class Voxel_tpc;

class Tpc {

 public:
  Tpc();
  Tpc(float, float, float, unsigned int, float, float, float);
  ~Tpc();

  /**
   * an inline function which get the row hits 
   * being inline makes arow=the_tpc.getRowHits as fast as arow=row_hits[row]
   * to make it inline the implentation of the function must be declared here
   */
  vector <Voxel_tpc *> getRowHits(unsigned int row ){return row_hits[row];};

  void putRowHit(unsigned int , Voxel_tpc *);
  void putSimTrackerHit(Voxel_tpc *, SimTrackerHit *);
  SimTrackerHit * getSimTrackerHit(Voxel_tpc *);

  float getOuterRadius(){return outer_radius;};
  float getInnerRadius(){return inner_radius;};
  float getHalfLength(){return half_length;};
  int   getNumberOfRows(){return number_rows;};
  float getPadRPitch(){return pad_r_pitch;};
  float getPixelPhiWidth(){return pixel_phi_width;};
  float getTimeWidth(){return time_width;};
  int   getNumberOfTimeSlices(){return number_time_slices;};
  float getTpcRphiResMax(){return tpc_rphi_res_max;};
  float getTpcZRes(){return tpc_z_res;};
  int   getNumberOfPhiSegments(int);
  double getRowRadius(int);
  float gettpc_rphi_res(float);



 private:

  float outer_radius;
  float inner_radius;
  float half_length;
  int number_rows;
  float pad_r_pitch;
  float pixel_phi_width;
  float time_width;
  int number_time_slices;

  float tpc_rphi_res_max;
  float tpc_z_res;

  vector<int> number_of_phi_segments;
  vector<double> row_radius;

  vector< vector <Voxel_tpc *> > row_hits;

  map< Voxel_tpc *,SimTrackerHit *> tpc_hit_map;

};

// Global pointer to the Tk_Te_Bank structure which is defined in tktebank.cc 
extern Tpc * the_tpc;

#endif
