// python code was very slow...
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TGeoManager.h"
#include "TStopwatch.h"
#include "TParameter.h"

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
// TTree writer for det sim
#include "TMS_ReadoutTreeWriter.h"
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
  // Nice list of vars to avoid future confusion (internet search for genie groottracker)
  // https://internal.dunescience.org/doxygen/read__t2k__rootracker_8C.html
  int StdHepPdg[__EDEP_SIM_MAX_PART__];
  double StdHepP4[__EDEP_SIM_MAX_PART__][4];
  double EvtVtx[__EDEP_SIM_MAX_PART__][4];
  if (gRoo){
    gRoo->SetBranchStatus("*", false);
    gRoo->SetBranchStatus("StdHepPdg", true);
    gRoo->SetBranchStatus("StdHepP4", true);
    gRoo->SetBranchStatus("EvtVtx", true);
    gRoo->SetBranchAddress("StdHepPdg", StdHepPdg);
    gRoo->SetBranchAddress("StdHepP4", StdHepP4);
    gRoo->SetBranchAddress("EvtVtx", EvtVtx);
  }
  // The global manager
  TMS_Manager::GetInstance().SetFileName(filename);

  // Load up the geometry
  TMS_Geom::GetInstance().SetGeometry(geom);

  int N_entries = events->GetEntries();
  const int max_n_events_config = TMS_Manager::GetInstance().Get_MaximumNEvents();
  if (max_n_events_config >= 0 && max_n_events_config < N_entries) N_entries = max_n_events_config;
  
  int event_counter = 0;

  bool DrawPDF = TMS_Manager::GetInstance().Get_DrawPDF();

  std::cout << "Starting loop over " << N_entries << " entries..." << std::endl;
  TStopwatch Timer;
  Timer.Start();

  int i = 0;
  int truth_info_entry_number = 0;
  //N_entries = 500;

  // Do we overlay events
  bool Overlay = false;
  // How many events do we want to overlay?
  int nOverlays = 50;
  // The vector carrying our events that we want to overlay
  std::vector<TMS_Event> overlay_events;
  
  bool NerscOverlay = false;
  TParameter<double>* spillPeriod_s = (TParameter<double>*)input->Get("spillPeriod_s");
  if (spillPeriod_s != NULL) NerscOverlay = true;
  double SpillPeriod = 0;
  if (NerscOverlay) {
    std::cout<<"Combining spills"<<std::endl;
    SpillPeriod = spillPeriod_s->GetVal() * 1e9; // convert to ns
    std::cout<<"Found spillSeriod_s of "<<SpillPeriod<<"ns"<<std::endl;
    if (SpillPeriod < 1e7 || SpillPeriod > 1e11) {
      std::cout<<"Fatal: Found spillSeriod_s that is unusually high or low. Expecting something like ~1.2e9 ns"<<std::endl;
      exit(1);
    }
    TMS_Manager::GetInstance().Set_Nersc_Spill_Period(SpillPeriod);
  }
  int current_spill_number = 0;

  for (; i < N_entries; ++i) {
    if (N_entries <= 10 || i % (N_entries/10) == 0) {
      std::cout << "Processed " << i << "/" << N_entries << " (" << double(i)*100./N_entries << "%)" << std::endl;
    }
    
    events->GetEntry(i);
    // todo, gRoo has a different indexing than events with overlay
    if (gRoo)
      gRoo->GetEntry(i);
      
    // Todo: This should no longer be needed when this bug is fixed in the spill builder
    // https://github.com/DUNE/2x2_sim/issues/54
    event->EventId = i;
    //if (event->Primaries.size() > 0)
    //  std::cout<<"Entry "<<i<<", interaction number of vtx 0: "<<event->Primaries[0].GetInteractionNumber()<<", vs event.EventId "<<event->EventId<<std::endl;

    // Make a TMS event
    TMS_Event tms_event = TMS_Event(*event);
    tms_event.SetSpillNumber(i);
    // Fill up truth information from the GRooTracker object
    if (gRoo){
      tms_event.FillTruthFromGRooTracker(StdHepPdg, StdHepP4, EvtVtx);
    }

    // Keep filling up the vector and move on to the next event
    if (Overlay && i % nOverlays != 0) {
      overlay_events.push_back(tms_event);
      continue;
    }
    
    // Keep filling up the vector if within spill
    if (NerscOverlay) {
      double next_spill_time = (current_spill_number + 0.5) * TMS_Manager::GetInstance().Get_Nersc_Spill_Period();
      double current_spill_time = event->Primaries.begin()->Position.T();
      // Check that this neutrino is within spill, but not last event
      if (current_spill_time < next_spill_time && i != N_entries - 1) {
        overlay_events.push_back(tms_event);
        continue;
      }
    }

    // Add event information and truth from another event
    if (overlay_events.size() > 0) {
      // Now loop over previous events
      std::cout<<"Overlaying "<<overlay_events.size()<<" events"<<std::endl;
      // The first event should be the starting point so reverse it
      std::reverse(overlay_events.begin(), overlay_events.end());
      TMS_Event last_event = overlay_events.back();
      overlay_events.pop_back();
      for (auto &event : overlay_events) last_event.AddEvent(event);
      // Make sure to set the spill number correctly
      last_event.SetSpillNumber(current_spill_number);
      overlay_events.clear();
      // Now add this event as the first event in the next set
      overlay_events.push_back(tms_event);
      if (NerscOverlay) current_spill_number += 1;
      // ... and make this event the combined spill called "last_event"
      tms_event = last_event;
    }
    
    // Apply the det sim now, after overlaying events
    // This doesn't work right now
    tms_event.ApplyReconstructionEffects();

    // Dump information
    //tms_event.Print();

    // Calculate the mapping between vertex ID and visible energy for the primary event
    tms_event.GetVertexIdOfMostVisibleEnergy();
    
    // Save det sim information
    TMS_ReadoutTreeWriter::GetWriter().Fill(tms_event);
    
    int nslices = TMS_TimeSlicer::GetSlicer().RunTimeSlicer(tms_event);
    std::cout<<"Sliced event "<<i<<" into "<<nslices<<" slices"<<std::endl;
    
    // Check if this is not pileup
    if (gRoo && event->Primaries.size() == 1 && tms_event.GetNVertices() == 1) {
      // Fill the info of the one and only true vertex in the spill
      auto primary_vertex = event->Primaries[0];
      int interaction_number = primary_vertex.GetInteractionNumber();
      gRoo->GetEntry(interaction_number);
      tms_event.FillTruthFromGRooTracker(StdHepPdg, StdHepP4, EvtVtx);
    }

    TMS_TreeWriter::GetWriter().FillSpill(tms_event, truth_info_entry_number, nslices);
    truth_info_entry_number += nslices;
    
    // Could save per spill info here
    
    //std::cout<<"Ran time slicer"<<std::endl;
    for (int slice = 0; slice < nslices; slice++) {
      // First make an event based on the slice
      TMS_Event tms_event_slice;
      // If the time slicer is off, use the entire old TMS_Event. That way muon KE branch is copied.
      if (!TMS_Manager::GetInstance().Get_Reco_TIME_RunTimeSlicer()) tms_event_slice = tms_event;
      else tms_event_slice = TMS_Event(tms_event, slice);

      // Fill truth info, but only for slice != 0 (but with no time slicer, all slices = 1 so do it anyway.
      if (gRoo && (slice != 0 || nslices == 1)) {
        // First find the vertex which contributed most to the event
        // This also sets the slice specific variables
        int primary_vertex_id = tms_event_slice.GetVertexIdOfMostVisibleEnergy();
        if (primary_vertex_id >= 0) {
          // Now find out how much that true vertex contributed in general
          auto map = tms_event.GetTrueVisibleEnergyPerVertex();
          
          if (map.find(primary_vertex_id) == map.end()) 
              std::cout<<"Warning: Didn't find primary_vertex_id "<<primary_vertex_id<<" inside map of size "<<map.size()<<std::endl;
          double visible_energy_from_vertex = map[primary_vertex_id];
          tms_event_slice.SetTotalVisibleEnergyFromVertex(visible_energy_from_vertex);
          
          gRoo->GetEntry(primary_vertex_id);
          tms_event_slice.FillTruthFromGRooTracker(StdHepPdg, StdHepP4, EvtVtx);
          tms_event_slice.SetLeptonInfoUsingVertexID(primary_vertex_id);
        }
      }
      
      event_counter += 1;
      
      // Try finding some tracks
      TMS_TrackFinder::GetFinder().FindTracks(tms_event_slice);

#ifdef DUNEANAOBJ_ENABLED
      caf::SRTMS srtms = TMS_Utils::ConvertEvent();
#endif
      //tms_event_slice.Print();

      // View it
      if (DrawPDF) TMS_EventViewer::GetViewer().Draw(tms_event_slice);
      // Write it
      TMS_TreeWriter::GetWriter().Fill(tms_event_slice);
    }
  } // End loop over all the events

  Timer.Stop();
  std::cout << "Event loop took " << Timer.RealTime() << "s for " << i << " entries (" << Timer.RealTime()/N_entries << " s/entries)" << std::endl;

  TMS_TreeWriter::GetWriter().Write();
  TMS_ReadoutTreeWriter::GetWriter().Write();

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
