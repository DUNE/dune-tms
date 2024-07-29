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
  int nOverlays = 3;
  // The vector carrying our events that we want to overlay
  std::vector<TMS_Event> overlay_events;

  for (; i < N_entries; ++i) {
    events->GetEntry(i);
    // todo, gRoo has a different indexing than events with overlay
    if (gRoo)
      gRoo->GetEntry(i);

#ifndef DEBUG
    if (N_entries <= 10 || i % (N_entries/10) == 0) {
#endif
      std::cout << "Processed " << i << "/" << N_entries << " (" << double(i)*100./N_entries << "%)" << std::endl;
    }

    // Make a TMS event
    TMS_Event tms_event = TMS_Event(*event);

    // Fill up truth information from the GRooTracker object
    if (gRoo){
      tms_event.FillTruthFromGRooTracker(StdHepPdg, StdHepP4, EvtVtx);
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

    // Calculate the mapping between vertex ID and visible energy for the primary event
    tms_event.GetVertexIdOfMostVisibleEnergy();
    
    // Save det sim information
    TMS_ReadoutTreeWriter::GetWriter().Fill(tms_event);

    int nslices = TMS_TimeSlicer::GetSlicer().RunTimeSlicer(tms_event);
    
    // Check if this is not pileup
    if (event->Primaries.size() == 1 && tms_event.GetNVertices() == 1) {
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
      if (slice == 0 || slice == nslices - 1) std::cout<<"Processing slice "<<slice<<" of event number "<<i<<" / "<<N_entries<<std::endl;
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
          
          if ((unsigned long) primary_vertex_id >= event->Primaries.size())
              std::cout<<"Warning: primary_vertex_id "<<primary_vertex_id<<" is above event->Primaries.size "<<event->Primaries.size()<<std::endl;
          else {
            // Now set the remaining information from the gRoo tchain.
            auto primary_vertex = event->Primaries[primary_vertex_id];
            int interaction_number = primary_vertex.GetInteractionNumber();
            gRoo->GetEntry(interaction_number);
            tms_event_slice.FillTruthFromGRooTracker(StdHepPdg, StdHepP4, EvtVtx);
            // And the lepton info
            int lepton_index = -1;
            int current_index = 0;
            for (auto particle : primary_vertex.Particles) {
              int pdg = std::abs(particle.GetPDGCode());
              if (pdg >= 11 && pdg <= 16) {
                lepton_index = current_index;
                break;
              }
              current_index += 1;
            }
            if (lepton_index >= 0) {
              auto lepton = primary_vertex.Particles[lepton_index];
              int lepton_pdg = lepton.GetPDGCode();
              auto lepton_position = primary_vertex.GetPosition();
              auto lepton_momentum = lepton.GetMomentum();
              tms_event_slice.FillTrueLeptonInfo(lepton_pdg, lepton_position, lepton_momentum);
            }
          }
        }
      }
      
      event_counter += 1;
      int spill_number = i;
      tms_event_slice.SetSpillNumber(spill_number);
      std::cout << "Slice number: " << slice << std::endl;
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
