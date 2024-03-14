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
  int StdHepPdg[__EDEP_SIM_MAX_PART__];
  double StdHepP4[__EDEP_SIM_MAX_PART__][4];
  double StdHepX4[__EDEP_SIM_MAX_PART__][4];
  if (gRoo){
    gRoo->SetBranchStatus("*", false);
    gRoo->SetBranchStatus("StdHepPdg", true);
    gRoo->SetBranchStatus("StdHepP4", true);
    gRoo->SetBranchStatus("StdHepX4", true);
    gRoo->SetBranchAddress("StdHepPdg", StdHepPdg);
    gRoo->SetBranchAddress("StdHepP4", StdHepP4);
    gRoo->SetBranchAddress("StdHepX4", StdHepX4);
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

  int SpillNumber = 0;
  // Do we have a spill that we need to combine?
  bool Spill = false;
  TParameter<float>* spillPeriod_s = (TParameter<float>*)input->Get("spillPeriod_s");
  if (spillPeriod_s != NULL) Spill = true;
  double SpillPeriod = 0;
  if (Spill) {
    std::cout<<"Combining spills"<<std::endl;
    SpillPeriod = spillPeriod_s->GetVal() * 1e9; // convert to ns
    std::cout<<"Found spillSeriod_s of "<<SpillPeriod<<std::endl;
  }
  // The vector carrying our events that we want to overlay
  std::vector<TMS_Event> spill_events;

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
    // todo, gRoo has a different indexing than events with overlay
    if (gRoo)
      gRoo->GetEntry(i);

//#ifndef DEBUG
    if (N_entries <= 10 || i % (N_entries/10) == 0) {
      std::cout << "Processed " << i << "/" << N_entries << " (" << double(i)*100./N_entries << "%)" << std::endl;
    }
//#endif
    
    
    double primary_event_time = -999;
    if (event->Primaries.size() > 0)
      primary_event_time = event->Primaries.begin()->Position.T(); // ns

    // Remove spill time from truth info
    // Need to subtract spill time from hit time
    // Or else floats in tree writer will lose precision for spills > 0
    // ie. 50 + 1e9 = 1e9 for floats but not doubles, but using doubles crashes
    // But also we only care about hit time within spill anyway since we match by spill number
    double timeoffset = SpillNumber * SpillPeriod;
    // Check if end is beyond current spill
    while (primary_event_time > timeoffset + SpillPeriod) { 
      primary_event_time += SpillPeriod;
      SpillNumber += 1;
    }
    if (Spill) {
      std::cout<<"Adjusting hit times with time offset of "<<timeoffset<<std::endl;
      // ... interaction vertex
      for (std::vector<TG4PrimaryVertex>::iterator v = event->Primaries.begin(); v != event->Primaries.end(); ++v) {
        //v->Position.T() = event_time;
        double old_event_time = v->Position.T();
        double event_time = old_event_time - timeoffset;
        v->Position.SetT(event_time);
      }

      // ... trajectories
      std::cout<<"2Adjusting hit times with time offset of "<<timeoffset<<std::endl;
      for (std::vector<TG4Trajectory>::iterator t = event->Trajectories.begin(); t != event->Trajectories.end(); ++t) {
        // loop over all points in the trajectory
        for (std::vector<TG4TrajectoryPoint>::iterator p = t->Points.begin(); p != t->Points.end(); ++p) {
          double new_time = p->Position.T() - timeoffset;
          p->Position.SetT(new_time);
        }
      }

      // ... and, finally, energy depositions
      std::cout<<"3Adjusting hit times with time offset of "<<timeoffset<<std::endl;
      for (auto d = event->SegmentDetectors.begin(); d != event->SegmentDetectors.end(); ++d) {
        for (std::vector<TG4HitSegment>::iterator h = d->second.begin(); h != d->second.end(); ++h) {
          double start_time = h->Start.T() - timeoffset;
          double stop_time = h->Stop.T() - timeoffset;
          //std::cout<<"start time original: "<<h->Start.T()<<", start_time_new: "<<start_time<<std::endl;
          h->Start.SetT(start_time);
          h->Stop.SetT(stop_time);
        }
      }
    }
    
    // Make a TMS event
    TMS_Event tms_event = TMS_Event(*event, true, timeoffset, SpillNumber);
    // Fill up truth information from the GRooTracker object
    if (gRoo){
      tms_event.FillTruthFromGRooTracker(StdHepPdg, StdHepP4);
      tms_event.FillAdditionalTruthFromGRooTracker(StdHepX4);
    }

    // Keep filling up the vector and move on to the next event
    if (Overlay && i % nOverlays != 0) {
      overlay_events.push_back(tms_event);
      continue;
    }

    // Add event information and truth from another event
    if (Overlay && i % nOverlays == 0) {
      // Now loop over previous events
      std::cout << "Overlaying" << std::endl;
      for (auto &event : overlay_events) tms_event.AddEvent(event);
      overlay_events.clear();
    }

    // Check if this event is part of the same spill, or the start of the next spill
    if (Spill) {
      // Need some truth info for the next part    
      //tms_event.FillTruthFromGRooTracker(StdHepPdg, StdHepP4);
      //tms_event.FillAdditionalTruthFromGRooTracker(StdHepX4);
      //std::cout<<"StdHepX4: "<<StdHepX4[0][0]<<","<<StdHepX4[0][1]<<","<<StdHepX4[0][2]<<","<<StdHepX4[0][3]<<std::endl;
    
      if (event->Primaries.size() > 0) {
        //double event_time = tms_event.GetNeutrinoX4().T(); // This is always 0
        //double alt_event_time = tms_event.GetHits().size() == 0 ? -9999 : tms_event.GetHits()[0].GetT();
        double end_of_spill = (SpillNumber + 0.5) * SpillPeriod; 
        //std::cout<<"i="<<i<<", event_time = "<<event_time<<", primary_event_time = "<<primary_event_time<<", alt_event_time = "<<alt_event_time<<",end of spill = "<<end_of_spill<<std::endl;
        //std::cout<<"primary_event_time="<<primary_event_time<<", end_of_spill="<<end_of_spill<<std::endl;
        if (primary_event_time <= end_of_spill) {
          // This event is within the window of the current spill,
          // so add to list and continue iterating
          spill_events.push_back(tms_event);
          // Continue to next event unless this is the last event in the file
          if (i < N_entries - 1) continue;
        }
      }
      else { 
        std::cout<<"Warning: Found g4 event without a primary. Skipping"<<std::endl;
      }
    }

    // Now combine the spills into a single tms_event
    if (Spill) {
      if (spill_events.size() == 0) { 
        std::cout<<"Found empty spill_events "<<spill_events.size()<<" for spill number "<<SpillNumber<<std::endl;
        continue; 
      }
      std::cout<<"Combining "<<spill_events.size()<<" events into a single spill"<<std::endl;
      // Get the first event in the list
      //auto temp_event = spill_events[0];
      // Now add all the events pass the first event
      //for (auto it = std::begin(spill_events) + 1; it != std::end(spill_events); ++it) temp_event.AddEvent(*it);
      TMS_Event temp_event;
      std::cout << "spill2 in" << std::endl << std::flush;
      for (auto &event : spill_events) temp_event.AddEvent(event);
      spill_events.clear();
      std::cout << "spill2 oot" << std::endl << std::flush;
      // Now save the current event as the first event of the next spill
      spill_events.push_back(tms_event);
      // And finally run the next bit of code on the combined event
      tms_event = temp_event;
      SpillNumber += 1;
    }

    // Dump information
    //tms_event.Print();
    
    // Now do the detector sim after overlay/spill code combines multiple events
    tms_event.ApplyReconstructionEffects();

    // Calculate the mapping between vertex ID and visible energy for the primary event
    tms_event.GetVertexIdOfMostVisibleEnergy();
    
    // Save det sim information
    TMS_ReadoutTreeWriter::GetWriter().Fill(tms_event);

    int nslices = TMS_TimeSlicer::GetSlicer().RunTimeSlicer(tms_event);
    
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
          
          double visible_energy_from_vertex = map[primary_vertex_id];
          tms_event_slice.SetTotalVisibleEnergyFromVertex(visible_energy_from_vertex);
          
          // Now set the remaining information from the gRoo tchain.
          auto primary_vertex = event->Primaries[primary_vertex_id];
          int interaction_number = primary_vertex.GetInteractionNumber();
          gRoo->GetEntry(interaction_number);
          tms_event_slice.FillTruthFromGRooTracker(StdHepPdg, StdHepP4);
          tms_event_slice.FillAdditionalTruthFromGRooTracker(StdHepX4);
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
      
      event_counter += 1;
      int spill_number = i;
      tms_event_slice.SetSpillNumber(spill_number);

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
