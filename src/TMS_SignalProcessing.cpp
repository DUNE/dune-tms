#include "TMS_SignalProcessing.h"
#include "TMS_Readout_Manager.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

void TMS_SignalProcessing::MergeCoincidentHits(TMS_Event &event) {
  std::vector<TMS_Hit> &TMS_Hits = event.GetHitsRawRef();

  std::sort(TMS_Hits.begin(), TMS_Hits.end(), TMS_Hit::SortByZThenT);

  const double readout_time = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_ReadoutTime();

  // Loop over the original hits
  for (std::vector<TMS_Hit>::iterator it = TMS_Hits.begin();
      it != TMS_Hits.end(); it++) {
    if ((*it).GetPedSup()) continue; // Skip hits which are already removed
    // Maybe this hit has already been counted
    double z = (*it).GetZ();
    double y = (*it).GetNotZ();
    //double e = hit.GetE();
    double t = (*it).GetT();

    // Look ahead to find duplicates, but stop when z != z2
    std::vector<std::vector<TMS_Hit>::iterator> duplicates;
    for (std::vector<TMS_Hit>::iterator jt = it+1; jt != TMS_Hits.end() && jt->GetZ() == z; ++jt) {

      TMS_Hit hit2 = *(jt);
      if (hit2.GetPedSup()) continue; // Skip hits which are already removed
      double z2 = hit2.GetZ();
      double y2 = hit2.GetNotZ();
      //double e2 = hit2.GetE();
      double t2 = hit2.GetT();

      // Merge
      if (z == z2 && y == y2 && fabs(t2-t) < readout_time) {
        (*it).MergeWith(hit2);
        // Phase III: merge the event-level truth side table by HitId alongside the reco-level
        // merge above, since TMS_TrueHit is no longer embedded in TMS_Hit.
        event.MergeTrueHit((*it).GetHitId(), hit2.GetHitId());
        // todo, we may want to store an array of true hits. One way would be to move the merging code within the hit class
        duplicates.push_back(jt);
      }
    }
    // Now flag the duplicates for removal
    for (auto rit = duplicates.rbegin(); rit != duplicates.rend(); ++rit) {
      auto hit_to_erase = *rit;
      hit_to_erase->SetPedSup(true);
    }
  }
  // Now erase all hits that are set as ped supped
  std::vector<TMS_Hit> remaining_hits;
  for (auto& hit : TMS_Hits) {
    if (!hit.GetPedSup()) {
      remaining_hits.push_back(hit);
      if (hit.GetE() > 10000)  std::cout << "Warning: Found hit higher than 10 GeV. Seems unlikely. Hit E = " << (hit.GetE() / 1000.0) << " GeV." << std::endl;
    } else {
      event.EraseTrueHit(hit.GetHitId());
    }
  }
  TMS_Hits = std::move(remaining_hits);
}

void TMS_SignalProcessing::SimulatePedestalSubtraction(TMS_Event &event) {
  // Don't actually remove the hits because they may be relevant for other processes
  std::vector<TMS_Hit> &TMS_Hits = event.GetHitsRawRef();

  for (auto& hit : TMS_Hits) {
     double hit_pe = hit.GetPE();
     if (hit_pe < TMS_Readout_Manager::GetInstance().Get_Sim_Readout_PedestalSubtractionThreshold()) {
       hit.SetPedSup(true);
     }
  }
}
