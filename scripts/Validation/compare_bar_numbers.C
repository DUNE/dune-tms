// Diagnostic for the geometry-restructuring work on feature/usher_geometry-from-gdml:
// dumps every branch whose value comes from TMS_Bar::GetBarNumber() (the value that
// changed when TMS_Bar.cpp's BarNumber formula was switched from the hardcoded
// TMS_Const::TMS_Start_Exact to the geometry-survey-derived GetXStartOfTMS()/GetYStartOfTMS())
// from two ConvertToTMSTree output files to text, so they can be diffed with plain `diff`.
//
// Usage (inside the SL7/Apptainer build environment, so ROOT resolves the same way the
// build does):
//   root -b -q 'compare_bar_numbers.C("before.root", "after.root")'
// Then:
//   diff before_scan/*.txt after_scan/*.txt
//
// Deliberately excludes RecoTrackKalmanPlaneBarView (per-hit-per-track, up to
// __TMS_MAX_TRACKS__ * __TMS_MAX_LINE_HITS__ * 3 entries) from the Scan dump -- it's the
// same BarNumber value as the others, just per-hit-in-track instead of per-track-endpoint,
// and dumping it in full would be unreadably large for a first pass. Extend the branch
// list below if that finer-grained check is ever needed.

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

void DumpTree(const char* filename, const char* treename, const char* branches, const char* outfile) {
  TFile *f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Could not open " << filename << std::endl;
    return;
  }
  TTree *t = (TTree*)f->Get(treename);
  if (!t) {
    std::cerr << "Tree " << treename << " not found in " << filename << std::endl;
    f->Close();
    return;
  }
  t->SetScanField(0); // no row limit

  gSystem->RedirectOutput(outfile, "w");
  std::cout << "# " << filename << " : " << treename << " : " << t->GetEntries() << " entries" << std::endl;
  t->Scan(branches, "", "colsize=12");
  gSystem->RedirectOutput(0);

  f->Close();
}

void compare_bar_numbers(const char* file_before, const char* file_after) {
  gSystem->mkdir("before_scan", true);
  gSystem->mkdir("after_scan", true);

  // RecoHitBar: per-reconstructed-hit bar index, bounded by nHits.
  // Tree is registered as "Line_Candidates" -- that's the ROOT tree name; the C++
  // member variable that owns it is called Branch_Lines (see TMS_TreeWriter.cpp:41).
  DumpTree(file_before, "Line_Candidates", "EventNo:SliceNo:SpillNo:nHits:RecoHitBar", "before_scan/RecoHitBar.txt");
  DumpTree(file_after,  "Line_Candidates", "EventNo:SliceNo:SpillNo:nHits:RecoHitBar", "after_scan/RecoHitBar.txt");

  // Kalman track-endpoint bar views: bounded by nTracks, tree Reco_Tree
  DumpTree(file_before, "Reco_Tree",
      "EventNo:SliceNo:SpillNo:nTracks:RecoTrackKalmanFirstPlaneBarView:RecoTrackKalmanLastPlaneBarView:"
      "RecoTrackKalmanFirstPlaneBarViewTrue:RecoTrackKalmanLastPlaneBarViewTrue",
      "before_scan/RecoTrackKalmanBarView.txt");
  DumpTree(file_after, "Reco_Tree",
      "EventNo:SliceNo:SpillNo:nTracks:RecoTrackKalmanFirstPlaneBarView:RecoTrackKalmanLastPlaneBarView:"
      "RecoTrackKalmanFirstPlaneBarViewTrue:RecoTrackKalmanLastPlaneBarViewTrue",
      "after_scan/RecoTrackKalmanBarView.txt");

  // TrueHitBar: per-true-hit bar index, bounded by NTrueHits, tree Truth_Info
  DumpTree(file_before, "Truth_Info", "EventNo:SpillNo:NTrueHits:TrueHitBar", "before_scan/TrueHitBar.txt");
  DumpTree(file_after,  "Truth_Info", "EventNo:SpillNo:NTrueHits:TrueHitBar", "after_scan/TrueHitBar.txt");

  std::cout << "Dumped before_scan/*.txt and after_scan/*.txt -- run:" << std::endl;
  std::cout << "  diff before_scan/RecoHitBar.txt after_scan/RecoHitBar.txt" << std::endl;
  std::cout << "  diff before_scan/RecoTrackKalmanBarView.txt after_scan/RecoTrackKalmanBarView.txt" << std::endl;
  std::cout << "  diff before_scan/TrueHitBar.txt after_scan/TrueHitBar.txt" << std::endl;
}
