#include "GroupTracks.h"

GroupTracks::GroupTracks() {
  _trackARVec.clear();
}

GroupTracks::GroupTracks( TrackExtended * track ) {
  _trackARVec.clear();
  _trackARVec.push_back( track );
}

GroupTracks::~GroupTracks() {}

void GroupTracks::addTrackExtended( TrackExtended * track ) {  
  _trackARVec.push_back( track );
}

void GroupTracks::ClearTrackExtendedVec() {
  _trackARVec.clear();
}

TrackExtendedVec & GroupTracks::getTrackExtendedVec() {
  return _trackARVec;
}
