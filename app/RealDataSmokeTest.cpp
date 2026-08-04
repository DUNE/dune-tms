// Phase III Stage C: synthetic real-data smoke test.
//
// No real DAQ input format exists yet for TMS (confirmed against the repo/issue
// history), so there's no real file to run this against. Instead this builds a
// handful of TMS_Hits directly via the truthless constructor (no TG4HitSegment, no
// TMS_TrueHit at all) and pushes them through the same pipeline stages a real DAQ
// hit stream would need: TMS_SignalProcessing::MergeCoincidentHits() ->
// TMS_TimeSlicer::RunTimeSlicer() -> [per slice] TMS_TrackFinder::FindTracks() ->
// TMS_TreeWriter::Fill(). (The time slicer wasn't in the original Stage C plan's
// abbreviated pipeline list, but TMS_TreeWriter::Fill() throws without it -- it
// populates TimeSliceBounds, a hard prerequisite for any event reaching Fill(), real
// data included.) The point is to exercise the "reco path runs with the truth table
// entirely empty" claim end to end, not to produce physically meaningful
// reconstruction.
//
// Geometry still needs to come from somewhere -- an edep-sim file's embedded
// TGeoManager is reused purely as a source of bar/plane geometry (same as
// TrackLengthTester.cpp); none of that file's TG4Event/truth branches are read.
//
// Hit positions come directly from a walk of the real geometry tree (see
// FindRepresentativeBarPositions() below), not a guessed coordinate grid: an earlier
// version of this test scanned a coordinate line/grid hoping to land inside a bar and
// consistently found zero, because module-layer daughters ("Module" nodes) each only
// span their own module's x-slot rather than the full detector width -- there's no
// x value that's "generically" inside a bar independent of which module it falls in.
// Walking the tree and reading real placement matrices sidesteps that assumption
// entirely and works regardless of a geometry's specific module/bar partitioning.

#include <iostream>
#include <functional>
#include <set>
#include <array>
#include "TFile.h"
#include "TGeoManager.h"

// Geometry singleton
#include "TMS_Geom.h"
// Event class
#include "TMS_Event.h"
// Real-or-simulated signal processing (hit merging)
#include "TMS_SignalProcessing.h"
// Reconstructor
#include "TMS_Reco.h"
// Time slicer
#include "TMS_TimeSlicer.h"
// TTree writer
#include "TMS_TreeWriter.h"
// General manager
#include "TMS_Manager.h"

namespace {

// Walks the geometry tree from the top, and for each distinct ModuleLayer node number
// (i.e. each distinct plane) encountered for the first time, records the global center
// of that plane's first module's first scintillator-bar leaf. One real, known-good
// position per plane, up to maxPositions -- mirrors TMS_Geom::SurveyNodeRecursive()'s
// traversal pattern but only descends one representative bar per plane instead of
// visiting every bar.
std::vector<std::array<double, 3>> FindRepresentativeBarPositions(TGeoManager *g, int maxPositions) {
  TMS_Manager &manager = TMS_Manager::GetInstance();
  const std::string moduleLayerName = manager.Get_GEOMETRY_VOLUME_ModuleLayer();

  std::vector<std::array<double, 3>> positions;
  std::set<int> seenPlanes;

  std::function<void()> recurse = [&]() {
    if ((int)positions.size() >= maxPositions) return;
    TGeoNode *node = g->GetCurrentNode();
    std::string name(node->GetName());

    if (name.find(moduleLayerName) != std::string::npos) {
      int planeNum = node->GetNumber();
      if (seenPlanes.insert(planeNum).second && node->GetNdaughters() > 0) {
        g->CdDown(0); // first module
        TGeoNode *moduleNode = g->GetCurrentNode();
        if (moduleNode->GetNdaughters() > 0) {
          g->CdDown(0); // first bar
          const double local[3] = {0, 0, 0};
          double master[3];
          g->GetCurrentMatrix()->LocalToMaster(local, master);
          positions.push_back({master[0], master[1], master[2]});
          g->CdUp();
        }
        g->CdUp();
      }
      return; // don't descend further into this module layer's subtree
    }

    int nDaughters = node->GetNdaughters();
    for (int i = 0; i < nDaughters && (int)positions.size() < maxPositions; i++) {
      g->CdDown(i);
      recurse();
      g->CdUp();
    }
  };

  g->CdTop();
  recurse();
  g->CdTop();
  return positions;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Need one argument: an edep-sim file to source detector geometry from "
                 "(its TG4Event truth is not read)" << std::endl;
    return -1;
  }

  std::string filename = std::string(argv[1]);
  TMS_Manager::GetInstance().SetFileName(filename);

  TFile *input = new TFile(filename.c_str(), "open");
  TGeoManager *geom = (TGeoManager*)input->Get("EDepSimGeometry");
  if (!geom) {
    std::cerr << "FAIL: could not find EDepSimGeometry in " << filename << std::endl;
    return -1;
  }
  TMS_Geom::GetInstance().SetGeometry(geom);

  TMS_Event event;
  event.SetSpillNumber(0);

  const int maxPositions = 60;
  std::vector<std::array<double, 3>> positions =
      FindRepresentativeBarPositions(TMS_Geom::GetInstance().GetGeometry(), maxPositions);
  std::cout << "Found " << positions.size() << " representative real bar positions "
               "(one per distinct plane, from the geometry tree)" << std::endl;

  int nHitsAdded = 0;
  for (const auto &pos : positions) {
    TMS_Hit candidate(pos[0], pos[1], pos[2], /*energy=*/2.0, /*time=*/nHitsAdded * 1.0, /*pe=*/50.0);
    if (candidate.GetPlaneNumber() < 0 || candidate.GetBarNumber() < 0) {
      std::cerr << "WARNING: real bar-leaf position (" << pos[0] << ", " << pos[1] << ", " << pos[2]
                << ") did not resolve to a valid bar via TMS_Bar(x,y,z) -- skipping" << std::endl;
      continue;
    }
    candidate.SetHitId(event.NextHitId());
    event.GetHitsRawRef().push_back(candidate);
    nHitsAdded++;
  }

  std::cout << "Built " << nHitsAdded << " truthless hits from real bar positions" << std::endl;
  if (nHitsAdded == 0) {
    std::cerr << "FAIL: none of the real bar-leaf positions resolved via TMS_Bar(x,y,z) -- "
                 "the truthless constructor's geometry lookup disagrees with the tree walk "
                 "above, which found real bar leaf nodes directly" << std::endl;
    return -1;
  }

  // Confirm the truth table really is empty for every hit -- this is the actual
  // "real data" claim under test, not just the absence of a crash.
  for (const auto &hit : event.GetHitsRawRef()) {
    if (event.GetTrueHit(hit.GetHitId()) != nullptr) {
      std::cerr << "FAIL: hit " << hit.GetHitId()
                << " unexpectedly has an entry in TrueHitByHitId" << std::endl;
      return -1;
    }
  }
  std::cout << "Confirmed: TrueHitByHitId has no entry for any of the " << nHitsAdded
            << " hits" << std::endl;

  // Run the same real-or-simulated pipeline stages a real DAQ hit stream would need.
  TMS_SignalProcessing::GetInstance().MergeCoincidentHits(event);
  std::cout << "After MergeCoincidentHits: " << event.GetNHits() << " hits" << std::endl;

  // TMS_TreeWriter::Fill() unconditionally reads event.GetTimeSliceBounds() (see
  // TMS_TreeWriter.cpp's TimeSliceStartTime/EndTime fill), which throws if
  // TimeSliceBounds is empty -- only TMS_TimeSlicer::RunTimeSlicer() populates it. This
  // is a hard prerequisite for any event reaching Fill(), real data included, not a
  // det-sim/MC-only step, so it belongs here even though the original Stage C plan
  // text didn't name it explicitly.
  int nslices = TMS_TimeSlicer::GetSlicer().RunTimeSlicer(event);
  std::cout << "RunTimeSlicer produced " << nslices << " slice(s)" << std::endl;

  for (int slice = 0; slice < nslices; slice++) {
    TMS_Event event_slice = TMS_Event(event, slice);

    TMS_TrackFinder::GetFinder().FindTracks(event_slice);
    std::cout << "Slice " << slice << ": FindTracks completed, "
              << event_slice.GetTracks().size() << " reconstructed track(s)" << std::endl;

    TMS_TreeWriter::GetWriter().Fill(event_slice);
  }
  TMS_TreeWriter::GetWriter().Write();

  std::cout << "PASS: synthetic truthless-hit smoke test completed with no crash" << std::endl;

  input->Close();
  return 0;
}
