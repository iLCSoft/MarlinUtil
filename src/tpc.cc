#include <iostream>
#include <vector>
#include <stdexcept>
#include "tpc.h"
#include "constants.h"


using namespace std;
using namespace constants;

// Global pointer to the TPC
Tpc * the_tpc;

Tpc::Tpc(){
}

Tpc::Tpc(float o_radius, float i_radius, float h_length, unsigned int n_rows, float p_r_pitch, float p_phi_width, float t_width, float rphires, float zres){

  outer_radius = o_radius;
  inner_radius = i_radius;
  half_length = h_length;
  number_rows = n_rows;
  pad_r_pitch = p_r_pitch;
  pixel_phi_width = p_phi_width;
  time_width = t_width;
  number_time_slices = (int)(2.* half_length / time_width); 
  tpc_rphi_res_max = rphires;
  tpc_z_res = zres;

  // set size of row_hits to hold (n_rows) vectors
  row_hits.resize(n_rows);



  // set row radius and number of phi section per row
  
  for (int i=0; i<number_rows; i++){
    row_radius.push_back( (pad_r_pitch/2.) + ( (float)i*pad_r_pitch) + inner_radius);
    number_of_phi_segments.push_back(int((twopi*row_radius[i]) / pixel_phi_width));
  }

  cout << "trying: A tpc has been instaniated" << endl;  

  //  for (unsigned int i=0; i<n_rows; i++){
  //    vector <Voxel_tpc *> next_row;
  //    row_hits.push_back(next_row);
  //  }
  
}

Tpc::~Tpc(){
}

void Tpc::putRowHit(unsigned int row, Voxel_tpc * p_voxel){


  if (row_hits.size() > row ){
    row_hits[row].push_back(p_voxel);
  }

  else{ 
    cout << "row index out of range" << endl;
    exit(1);
  }

  //  cout << "row = " << row << endl; 
  //  cout << "row_hits.size[row] = " << row_hits[row].size() << endl;

}

void Tpc::putSimTrackerHit(Voxel_tpc * voxel, SimTrackerHit * hit){
  tpc_hit_map[voxel] = hit;
}


SimTrackerHit * Tpc::getSimTrackerHit(Voxel_tpc * voxel){
  if (SimTrackerHit * hit = tpc_hit_map[voxel]) return hit;
  else return NULL;
}

int Tpc::getNumberOfPhiSegments(int row){
  
  //  if(row>=0&&row<number_of_phi_segments.size()) return number_of_phi_segments.at(row);
  //  else throw out_of_range("row out of range");

  // use the .at(row) method instead of [row]. This method throws and exception 
  return number_of_phi_segments.at(row);
};

double Tpc::getRowRadius(int row){

  //  if(row>0&&row<row_radius.size()) return row_radius[row];
  //  else throw out_of_range("row out of range");
  return row_radius.at(row);



};

float Tpc::gettpc_rphi_res(float z){
return tpc_rphi_res_max-fabs(z)/half_length*0.10;
};


