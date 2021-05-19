// python code was very slow...
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TGeoManager.h"
#include "TStopwatch.h"

// For now write the efficiency here...
#include "TH1D.h"

// EDepSim includes
#include "EDepSim/TG4Event.h"
#include "EDepSim/TG4PrimaryVertex.h"

// TMS includes
// Geometry singleton
#include "TMS_Geom.h"
// Event class
#include "TMS_Event.h"
// Event viewer singleton
#include "TMS_EventViewer.h"
// Reconstructor
#include "TMS_Reco.h"

bool ConvertToTMSTree(std::string filename, std::string output_filename) {
  std::cout << "Got " << filename << ", writing to " << output_filename << std::endl;

  // The input file
  TFile *input = new TFile(filename.c_str());
  // The EDepSim events
  TTree *events = (TTree*)input->Get("EDepSimEvents");
  // The generator pass-through information
  TTree *gRoo = (TTree*)input->Get("DetSimPassThru/gRooTracker");
  // Get the detector geometry
  TGeoManager *geom = (TGeoManager*)input->Get("EDepSimGeometry");

  // Get the event
  TG4Event *event = NULL;
  events->SetBranchAddress("Event", &event);

  // Load up the geometry
  TMS_Geom::GetInstance().SetGeometry(geom);
  TMS_Geom::GetInstance().SetFileName(filename);

  int N_entries = events->GetEntries();

  std::cout << "Starting loop over " << N_entries << " entries..." << std::endl;
  TStopwatch Timer;
  Timer.Start();

  int i = 0;
  for (; i < N_entries; ++i) {
    if (i > 1000) break;
    events->GetEntry(i);
    gRoo->GetEntry(i);

    if (i % (N_entries/10) == 0) {
      std::cout << "Processed " << i << "/" << N_entries << " (" << double(i)*100./N_entries << "%)" << std::endl;
    }

    // Make a TMS event
    TMS_Event tms_event = TMS_Event(*event);
    // Dump information
    //tms_event.Print();

    // Try finding some tracks
    TMS_TrackFinder::GetFinder().FindTracks(tms_event);

    // View it
    TMS_EventViewer::GetViewer().Draw(tms_event);
  } // End loop over all the events

  /*
  TH1D* eff = TMS_TrackFinder::GetFinder().GetEfficiencyHist();
  TH1D* total = TMS_TrackFinder::GetFinder().GetTotalHist();
  TH1D* eff_ratio = TMS_TrackFinder::GetFinder().GetEfficiency();
  TFile *file = new TFile("test_eff.root", "recreate");
  file->cd();
  eff->Write();
  total->Write();
  eff_ratio->Write();
  file->Close();
  */

  Timer.Stop();
  std::cout << "Event loop took " << Timer.RealTime() << "s for " << i << " entries (" << Timer.RealTime()/N_entries << " s/entries)" << std::endl;

  return true;
}

int main(int argc, char **argv) {
  if (argc != 2 && argc != 3) {
    std::cerr << "Need one or two arguments: [EDepSim output file] [Output filename]" << std::endl;
    return -1;
  }

  std::string EDepSimFile = std::string(argv[1]);
  std::string OutputFile;

  // If two arguments are given
  if (argc == 2) {
    std::string filename = std::string(argv[1]);
    OutputFile = filename.substr(0, filename.find(".root"));
    OutputFile += "_output.root";
  } else {
    OutputFile = std::string(argv[2]);
  }

  bool ok = ConvertToTMSTree(EDepSimFile, OutputFile);
  if (ok) return 0;
  else return -1;
}
