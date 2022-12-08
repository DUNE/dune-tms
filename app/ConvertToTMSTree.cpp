// python code was very slow...
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TGeoManager.h"
#include "TStopwatch.h"

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
// TTree writer
#include "TMS_TreeWriter.h"
// General manager
#include "TMS_Manager.h"

bool ConvertToTMSTree(std::string filename, std::string output_filename) {
  std::cout << "Got " << filename << ", writing to " << output_filename << std::endl;

  // The input file
  TFile *input = new TFile(filename.c_str(), "open");
  // The EDepSim events
  TTree *events = (TTree*)(input->Get("EDepSimEvents")->Clone("events"));
  // The generator pass-through information
  TTree *gRoo = (TTree*) (input->Get("DetSimPassThru/gRooTracker"));
  if (gRoo)
    gRoo = (TTree*) (gRoo->Clone("gRoo"));
  // Get the detector geometry
  TGeoManager *geom = (TGeoManager*)input->Get("EDepSimGeometry");

  // Get the event
  TG4Event *event = NULL;
  events->SetBranchAddress("Event", &event);
  // Get the true neutrino vector from the gRooTracker object
  int StdHepPdg[__EDEP_SIM_MAX_PART__];
  double StdHepP4[__EDEP_SIM_MAX_PART__][4];
  if (gRoo){
    gRoo->SetBranchStatus("*", false);
    gRoo->SetBranchStatus("StdHepPdg", true);
    gRoo->SetBranchStatus("StdHepP4", true);
    gRoo->SetBranchAddress("StdHepPdg", StdHepPdg);
    gRoo->SetBranchAddress("StdHepP4", StdHepP4);
  }
  // The global manager
  TMS_Manager::GetInstance().SetFileName(filename);

  // Load up the geometry
  TMS_Geom::GetInstance().SetGeometry(geom);

  int N_entries = events->GetEntries();
  
  int event_counter = 0;

  bool DrawPDF = TMS_Manager::GetInstance().Get_DrawPDF();

  std::cout << "Starting loop over " << N_entries << " entries..." << std::endl;
  TStopwatch Timer;
  Timer.Start();

  int i = 0;
  //N_entries = 500;

  // Do we overlay events
  bool Overlay = false;
  // How many events do we want to overlay?
  int nOverlays = 3;
  // The vector carrying our events that we want to overlay
  std::vector<TMS_Event> overlay_events;

  for (; i < N_entries; ++i) {
    events->GetEntry(i);
    if (gRoo)
      gRoo->GetEntry(i);

#ifndef DEBUG
    if (i % (N_entries/10) == 0) {
#endif
      std::cout << "Processed " << i << "/" << N_entries << " (" << double(i)*100./N_entries << "%)" << std::endl;
    }

    // Make a TMS event
    TMS_Event tms_event = TMS_Event(*event);
    // Fill up truth information from the GRooTracker object
    if (gRoo){
      tms_event.FillTruthFromGRooTracker(StdHepPdg, StdHepP4);
    }

    // Keep filling up the vector and move on to the next event
    if (Overlay && i % nOverlays != 0) {
      overlay_events.push_back(tms_event);
      continue;
    }

    // Add event information and truth from another event
    if (Overlay && i % nOverlays == 0) {
      // Now loop over previous events
      for (auto &event : overlay_events) tms_event.AddEvent(event);
      overlay_events.clear();
    }

    // Dump information
    //tms_event.Print();

    // Try finding some tracks

    int nslices = TMS_TimeSlicer::GetSlicer().RunTimeSlicer(tms_event);
    for (int slice = 0; slice < nslices; slice++) {
    
      auto slice_hits = tms_event.GetHits(slice);
      /*double min_time = 1e9;
      double max_time = -1e9;
      for (auto hit : slice_hits) {
        if (hit.GetPedSup()) std::cout<<"Found a ped supped hit in slice"<<std::endl;
        min_time = std::min(min_time, hit.GetT());
        max_time = std::max(max_time, hit.GetT());
      }
      if (max_time - min_time > (10000.0 / 52)) std::cout<<"Found a time range larger than expected: "<<(max_time - min_time)<<", min="<<min_time<<", max="<<max_time<<", slice="<<slice<<std::endl; */
      TMS_Event tms_event_slice = TMS_Event();
      tms_event_slice.SetHitsRaw(slice_hits);
      tms_event_slice.SetSliceNumber(slice);
      tms_event_slice.SetEventNumber(event_counter);
      event_counter += 1;
      int spill_number = i;
      tms_event_slice.SetSpillNumber(spill_number);
      
      TMS_TrackFinder::GetFinder().FindTracks(tms_event_slice);

#ifdef DUNEANAOBJ_ENABLED
      caf::SRTMS srtms = TMS_Utils::ConvertEvent();
#endif

      // View it
      if (DrawPDF) TMS_EventViewer::GetViewer().Draw(tms_event_slice);
      // Write it
      TMS_TreeWriter::GetWriter().Fill(tms_event_slice);
    }
  } // End loop over all the events

  Timer.Stop();
  std::cout << "Event loop took " << Timer.RealTime() << "s for " << i << " entries (" << Timer.RealTime()/N_entries << " s/entries)" << std::endl;

  TMS_TreeWriter::GetWriter().Write();

  delete events;
  if (gRoo)
    delete gRoo;
  input->Close();

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
