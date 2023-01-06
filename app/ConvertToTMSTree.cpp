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
  double StdHepX4[__EDEP_SIM_MAX_PART__][4];
  //gRoo->Show(0);
  //gRoo->Show(1);
  //events->Show(0);
  //events->Show(1);
  //exit(0);
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
      //tms_event.FillTruthInformation(*event);
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

    int nslices = TMS_TimeSlicer::GetSlicer().RunTimeSlicer(tms_event);
    //std::cout<<"Ran time slicer"<<std::endl;
    for (int slice = 0; slice < nslices; slice++) {
      std::cout<<"Processing slice "<<slice<<std::endl;
      // First make an event based on the slice
      TMS_Event tms_event_slice = TMS_Event(tms_event, slice);
      //tms_event_slice.SetEventNumber(event_counter);
    
      //auto slice_hits = tms_event.GetHits(slice);
      //auto all_hits = tms_event.GetHits(-1);
      //std::cout<<"Slice: "<<slice<<", N hits in slice: "<<slice_hits.size()<<", n hits total: "<<all_hits.size()<<std::endl;
      /*double min_time = 1e9;
      double max_time = -1e9;
      for (auto hit : slice_hits) {
        if (hit.GetPedSup()) std::cout<<"Found a ped supped hit in slice"<<std::endl;
        min_time = std::min(min_time, hit.GetT());
        max_time = std::max(max_time, hit.GetT());
      }
      if (max_time - min_time > (10000.0 / 52)) std::cout<<"Found a time range larger than expected: "<<(max_time - min_time)<<", min="<<min_time<<", max="<<max_time<<", slice="<<slice<<std::endl; */
      //tms_event_slice.SetHitsRaw(slice_hits);
      //tms_event_slice.SetSliceNumber(slice);
      //std::cout<<"Made tms_event_slice"<<std::endl;
      
      // Fill truth info, but only for slice != 0 (but with no time slicer, all slices = 1 so do it anyway.
      // TODO add back
      if (false && gRoo && (slice != 0 || nslices == 1)) {
        // First find the vertex which contributed most to the event
        // This also sets the slice specific variables
        int primary_vertex_id = tms_event_slice.GetVertexIdOfMostVisibleEnergy();
        if (primary_vertex_id >= 0) {
          //std::cout<<"Got primary vertex id: "<<primary_vertex_id<<std::endl;
          // Now find out how much that true vertex contributed in general
          auto map = tms_event.GetTrueVisibleEnergyPerVertex();
          
          double visible_energy_from_vertex = map[primary_vertex_id];
          tms_event_slice.SetTotalVisibleEnergyFromVertex(visible_energy_from_vertex);
          //std::cout<<"For slice "<<slice<<", found total visible energy of "<<tms_event_slice.GetTotalVisibleEnergyFromVertex()<<" for id "<<primary_vertex_id<<", compare to "<<tms_event_slice.GetVisibleEnergyFromVertexInSlice()<<std::endl;
          //std::cout<<"N hits full: "<<tms_event.GetNHits()<<", N hits slice: "<<tms_event_slice.GetNHits()<<std::endl;
          //std::cout<<"Set total visible energy"<<std::endl;
          
          // Now set the remaining information from the gRoo tchain.
          auto primary_vertex = event->Primaries[primary_vertex_id];
          int interaction_number = primary_vertex.GetInteractionNumber();
          gRoo->GetEntry(interaction_number);
          tms_event_slice.FillTruthFromGRooTracker(StdHepPdg, StdHepP4);
          tms_event_slice.FillAdditionalTruthFromGRooTracker(StdHepX4);
          //std::cout<<"Set truth info"<<std::endl;
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
          //std::cout<<"Found particle"<<std::endl;
          if (lepton_index >= 0) {
            auto lepton = primary_vertex.Particles[lepton_index];
            int lepton_pdg = lepton.GetPDGCode();
            auto lepton_position = primary_vertex.GetPosition();
            auto lepton_momentum = lepton.GetMomentum();
            tms_event_slice.FillTrueLeptonInfo(lepton_pdg, lepton_position, lepton_momentum);
          }
          //std::cout<<"Filled lepton truth info"<<std::endl;
        }
        //std::cout<<"E vtx: "<<tms_event_slice.GetVisibleEnergyFromVertexInSlice()<<", E tot: "<<tms_event_slice.GetTotalVisibleEnergyFromVertex()<<", E other: "<<tms_event_slice.GetVisibleEnergyFromOtherVerticesInSlice()<<std::endl;
      
        /*
        bool found_correct_primary = false;
        auto time_range = tms_event_slice.GetEventTimeRange();
        double min_time = time_range.first;
        double max_time = time_range.second;
        const double buffer = 5; // ns because reco and true time aren't always the same
        int best_candidate = -999;
        double best_time_score = 1e9;
        for (TG4PrimaryVertexContainer::iterator it = event->Primaries.begin(); it != event->Primaries.end(); ++it) {
          double time = it->GetPosition().T();
          // Find neutrinos in slice
          if (min_time - buffer <= time && time <= max_time + buffer) {
            // Find event nearest to the front of the slice
            double time_score = std::abs(time - min_time);
            // todo, exclude NC events?
            if (best_candidate == -999 || time_score < best_time_score) {
              best_candidate = it - event->Primaries.begin();
              best_time_score = time_score;
            }
          }
        }
        if (best_candidate != -999) {
          int interaction_number = event->Primaries[best_candidate].GetInteractionNumber();
          gRoo->GetEntry(interaction_number);
          found_correct_primary = true;
          //std::cout<<"For slice "<<slice<<", found best candidate interaction_number="<<interaction_number<<", score="<<best_time_score<<std::endl;
        }
        if (found_correct_primary) { 
          tms_event_slice.FillTruthFromGRooTracker(StdHepPdg, StdHepP4);
          tms_event_slice.FillAdditionalTruthFromGRooTracker(StdHepX4);
          auto primary_vertex = event->Primaries[best_candidate];
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
            // Find all the energy associated with this vertex
            double total_energy_in_slice = 0;
            double total_energy_of_primary_neutrino_in_slice = 0;
            double total_other_energy_in_slice = 0;
            double total_energy_of_primary_neutrino = 0;
            for (auto hit : all_hits) {
              int slice_number = hit.GetSlice();
              //int primary_id = hit.GetTrueHit().GetPrimaryId();
              bool is_primary_neutrino = false;
              int vertex_id = hit.GetTrueHit().GetVertexId();
              if (vertex_id < 0) { std::cout<<"Found vertex id < 0: "<<vertex_id<<std::endl; exit(0); }
              if (vertex_id == best_candidate) {
                is_primary_neutrino = true;
              }
              if (slice_number == slice) { 
                total_energy_in_slice += hit.GetE();
                if (is_primary_neutrino) total_energy_of_primary_neutrino_in_slice += hit.GetE();
                else total_other_energy_in_slice += hit.GetE();
              }
              if (is_primary_neutrino) total_energy_of_primary_neutrino += hit.GetE();
              //hit.Print();
            }
            std::cout<<"Total energy in slice: "<<total_energy_in_slice<<std::endl;
            std::cout<<"Total primary energy in slice: "<<total_energy_of_primary_neutrino_in_slice<<std::endl;
            std::cout<<"Total other energy in slice: "<<total_other_energy_in_slice<<std::endl;
            std::cout<<"Total energy of primary: "<<total_energy_of_primary_neutrino<<std::endl;
            //if (total_energy_of_primary_neutrino_in_slice > 0) exit(0);
          }
          
        } */
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
