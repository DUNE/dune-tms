#ifndef _TMS_LIGHTCONTRIBCOLUMNS_H_SEEN_
#define _TMS_LIGHTCONTRIBCOLUMNS_H_SEEN_

#include <vector>

#include "Rtypes.h"

#include "TMS_TrueHit.h"

// Flat output columns for the per-particle light breakdown of the true hits. A true hit's entries are
// [offset, offset + count) of these vectors, one per particle that produced detected photons. Shared
// by the two tree writers (TMS_TreeWriter, TMS_ReadoutTreeWriter) so that the layout, which the
// validation scripts rely on, is defined in one place.
struct TMS_LightContribColumns {
  std::vector<int> PrimaryId;
  std::vector<Long64_t> VertexId;
  std::vector<int> Photons;
  std::vector<float> Share;

  void Clear() {
    PrimaryId.clear();
    VertexId.clear();
    Photons.clear();
    Share.clear();
  }

  // Appends the contributions of one true hit and returns where they start and how many there are
  void Append(const TMS_TrueHit& true_hit, int& offset, int& count) {
    offset = static_cast<int>(Photons.size());
    count = static_cast<int>(true_hit.GetNLightContributions());
    const int total = true_hit.GetNPhotons();
    for (int c = 0; c < count; c++) {
      PrimaryId.push_back(true_hit.GetLightContributionPrimaryId(c));
      VertexId.push_back(true_hit.GetLightContributionVertexGlobalId(c));
      Photons.push_back(true_hit.GetLightContributionPhotons(c));
      Share.push_back(total > 0 ? static_cast<float>(true_hit.GetLightContributionPhotons(c)) / total : 0.f);
    }
  }
};

#endif
