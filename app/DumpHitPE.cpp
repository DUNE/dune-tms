// Dumps per-hit true energy deposit (E), simulated photoelectrons (PE), and the pedestal-
// suppression flag (PedSup) for every TMS_Hit in a spill file, straight after
// TMS_Event::FinalizeEvent() runs the full detector-response chain (optical model, merging,
// readout noise, pedestal subtraction -- see TMS_Event::ApplyReconstructionEffects()).
//
// Deliberately skips time-slicing/track-finding/space-point-building -- none of it changes a
// hit's E/PE/PedSup, which are already final at the whole-spill TMS_Event level -- so this is
// much cheaper than a full ConvertToTMSTree pass. Mirrors ConvertToTMSTree's event-overlay
// logic (Overlay/NerscOverlay) so the hit population -- and therefore what merges/pedestal-
// suppresses -- matches what the real pipeline sees.

#include <fstream>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TGeoManager.h"
#include "TParameter.h"

#include "EDepSim/TG4Event.h"

#include "TMS_Event.h"
#include "TMS_Geom.h"
#include "TMS_Manager.h"
#include "TMS_TrueHit.h"
#include "TMS_Utils.h"

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <input_edep_sim_spills.root> <output_hits.csv>" << std::endl;
    return -1;
  }
  const std::string input_filename = argv[1];
  const std::string output_csv_path = argv[2];

  TFile *input = new TFile(input_filename.c_str(), "open");
  if (!input || input->IsZombie()) {
    std::cerr << "Failed to open input file: " << input_filename << std::endl;
    return -1;
  }

  TTree *events_raw = (TTree *)input->Get("EDepSimEvents");
  if (!events_raw) {
    std::cerr << "Input file is missing the required 'EDepSimEvents' tree: " << input_filename << std::endl;
    return -1;
  }
  TTree *events = (TTree *)events_raw->Clone("events");

  TGeoManager *geom = (TGeoManager *)input->Get("EDepSimGeometry");
  if (!geom) {
    std::cerr << "Input file is missing 'EDepSimGeometry': " << input_filename << std::endl;
    return -1;
  }
  TMS_Geom::GetInstance().SetGeometry(geom);
  TMS_Manager::GetInstance().SetFileName(input_filename);

  TG4Event *event = NULL;
  events->SetBranchAddress("Event", &event);

  bool NerscOverlay = false;
  TParameter<double> *spillPeriod_s = (TParameter<double> *)input->Get("spillPeriod_s");
  double SpillPeriod = 0;
  if (spillPeriod_s != NULL) {
    NerscOverlay = true;
    SpillPeriod = spillPeriod_s->GetVal() * 1e9;  // s -> ns
    std::cout << "Combining spills, spillPeriod_s = " << SpillPeriod << " ns" << std::endl;
    TMS_Manager::GetInstance().Set_Nersc_Spill_Period(SpillPeriod);
  }
  int current_spill_number = 0;

  std::ofstream out(output_csv_path);
  // recoE = hit.GetE() -- NOT the true deposit. TMS_DetectorSimulation::SimulateOpticalModel()
  // overwrites it at the end (TMS_DetectorSimulation.cpp ~168-171) to
  // `PE * Get_RECO_CALIBRATION_EnergyCalibration()`, a deterministic rescaling of PE by a fixed
  // MeV/PE constant -- so recoE/PE is a constant, not a physics correlation. trueE comes from the
  // separate TMS_TrueHit object (TMS_Event::GetTrueHit()), which SimulateOpticalModel reads from
  // but never overwrites, and is the real energy deposit.
  out << "recoE,trueE,PE,pedsup,pdg\n";

  int N_entries = events->GetEntries();
  std::vector<TMS_Event> overlay_events;
  long long n_hits_written = 0;

  for (int i = 0; i < N_entries; ++i) {
    if (N_entries <= 10 || i % (N_entries / 10) == 0) {
      std::cout << "Processed " << i << "/" << N_entries << std::endl;
    }
    events->GetEntry(i);
    event->EventId = i;

    TMS_Event tms_event = TMS_Event(*event);
    tms_event.SetSpillNumber(i);

    if (NerscOverlay) {
      double next_spill_time = (current_spill_number + 0.5) * TMS_Manager::GetInstance().Get_Nersc_Spill_Period();
      double current_spill_time = event->Primaries.begin()->Position.T();
      if (current_spill_time < next_spill_time && i != N_entries - 1) {
        overlay_events.push_back(tms_event);
        continue;
      }
    }

    if (overlay_events.size() > 0) {
      std::reverse(overlay_events.begin(), overlay_events.end());
      TMS_Event last_event = overlay_events.back();
      overlay_events.pop_back();
      last_event.OverlayEvents(overlay_events);
      last_event.SetSpillNumber(current_spill_number);
      overlay_events.clear();
      overlay_events.push_back(tms_event);
      if (NerscOverlay) current_spill_number += 1;
      tms_event = last_event;
    }

    // Runs the full detector-response chain (optical model, merging, readout noise, pedestal
    // subtraction) on the whole-spill event -- see TMS_Event::ApplyReconstructionEffects().
    tms_event.FinalizeEvent();

    for (const auto &hit : tms_event.GetHits(-1, /*include_ped_sup=*/true)) {
      // Dominant true contributor to this hit, same lookup TMS_TreeWriter uses for
      // RecoHitPrimary*Energy (see TMS_TreeWriter.cpp:1542-1553).
      int pdg = 0;  // 0 = no truth match (shouldn't normally happen on MC input)
      const auto particle_info = TMS_Utils::GetPrimaryIdsByEnergy({hit}, tms_event);
      if (!particle_info.energies.empty()) {
        const long long vertex_global_id = particle_info.vertexglobalids[0];
        const int track_id = particle_info.trackids[0];
        const int particle_index = tms_event.GetTrueParticleIndex(vertex_global_id, track_id);
        if (particle_index >= 0) {
          pdg = tms_event.GetTrueParticles()[particle_index].GetPDG();
        }
      }
      const TMS_TrueHit *true_hit = tms_event.GetTrueHit(hit.GetHitId());
      const double true_e = true_hit ? true_hit->GetE() : -1.0;  // -1 = no truth (shouldn't happen on MC)

      out << hit.GetE() << "," << true_e << "," << hit.GetPE() << ","
          << (hit.GetPedSup() ? 1 : 0) << "," << pdg << "\n";
      ++n_hits_written;
    }
  }

  out.close();
  std::cout << "Wrote " << n_hits_written << " hits to " << output_csv_path << std::endl;

  delete events;
  input->Close();
  return 0;
}
