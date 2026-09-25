#ifndef _SPACEPOINTLAYERINPUT_H_SEEN_
#define _SPACEPOINTLAYERINPUT_H_SEEN_

#include <vector>

#include "TTree.h"

#include "TMS_Geom.h"

// Point layers (TMS_SpacePoint::GetLayer()) for space points read back from a
// reco file. Files converted since 2026-09-25 carry them in SpacePointLayer;
// older files (BothNeighbors pairing) don't, but there every point sits
// exactly on its X-bar plane's z, so the nearest plane index is its layer.
class SpacePointLayerInput {
  public:
    SpacePointLayerInput(TTree *reco_tree, int max_points)
      : layers_(max_points, -1), from_file_(reco_tree->GetBranch("SpacePointLayer") != nullptr) {
      if (from_file_) reco_tree->SetBranchAddress("SpacePointLayer", layers_.data());
    }

    // Layer of point i (after the tree entry has been read), whose z is z.
    int Layer(int i, double z) const {
      return from_file_ ? layers_[i] : TMS_Geom::GetInstance().GetPlaneIndexNearestZ(z);
    }

    bool FromFile() const { return from_file_; }

  private:
    std::vector<int> layers_;
    bool from_file_;
};

#endif
