// Dumps every charged Geant4 trajectory that reaches the TMS in one spill,
// read straight from the raw TG4Event (no TMS detector response), for a
// whole-spill truth event display (2026-09-30, TMS studies talk).
//
// Usage: DumpSpillTrajectories <spill_file.root> <spillNumber> <out.json>
//
// Spills are grouped exactly as ConvertToTMSTree does with a NERSC-overlay
// file: interactions are stored in time order, and an interaction belongs to
// spill n when its first primary vertex time lies in
// [(n - 0.5) * period, (n + 0.5) * period), period = spillPeriod_s. So
// spillNumber here matches the converter's SpillNumber.
//
// For each trajectory with at least one point inside the TMS box
// (TMS_Geom::IsInsideTMS), ALL of its points are written, so the plotting
// side can clip segments that cross the box boundary. Neutral particles and
// nuclear fragments are skipped: they leave no track of their own.
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "TDatabasePDG.h"
#include "TFile.h"
#include "TParameter.h"
#include "TParticlePDG.h"
#include "TTree.h"

#include "EDepSim/TG4Event.h"

#include "TMS_Geom.h"
#include "TMS_Manager.h"

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <spill_file.root> <spillNumber> <out.json>\n";
    return 1;
  }
  const std::string filename = argv[1];
  const int spill_number = std::stoi(argv[2]);
  const std::string out_path = argv[3];

  TFile *input = new TFile(filename.c_str(), "open");
  if (!input || input->IsZombie()) {
    std::cerr << "Could not open " << filename << "\n";
    return 1;
  }
  TGeoManager *geom = (TGeoManager *)input->Get("EDepSimGeometry");
  TMS_Manager::GetInstance().SetFileName(filename);
  TMS_Geom::GetInstance().SetGeometry(geom);

  TParameter<double> *period_s = (TParameter<double> *)input->Get("spillPeriod_s");
  if (!period_s) {
    std::cerr << "No spillPeriod_s in " << filename << ": not a spill file\n";
    return 1;
  }
  const double period_ns = period_s->GetVal() * 1e9;
  const double spill_t0 = spill_number * period_ns;

  TTree *events = (TTree *)input->Get("EDepSimEvents");
  TG4Event *event = nullptr;
  events->SetBranchAddress("Event", &event);

  TDatabasePDG *pdg_db = TDatabasePDG::Instance();
  const TMS_Geom &tms = TMS_Geom::GetInstance();
  const TVector3 lo = tms.GetStartOfTMS(), hi = tms.GetEndOfTMS();

  std::ofstream json(out_path);
  json << "{\"spill\":" << spill_number << ",\"spill_t0_ns\":" << std::fixed << spill_t0
       << ",\"tms_box\":[[" << lo.X() << "," << lo.Y() << "," << lo.Z() << "],[" << hi.X() << "," << hi.Y() << ","
       << hi.Z() << "]],\"trajectories\":[";
  bool first = true;
  int n_interactions = 0, n_written = 0;
  for (Long64_t i = 0; i < events->GetEntries(); ++i) {
    events->GetEntry(i);
    if (event->Primaries.empty()) continue;
    const double t = event->Primaries.begin()->Position.T();
    const int spill_of_entry = (int)std::floor(t / period_ns + 0.5);
    if (spill_of_entry < spill_number) continue;
    if (spill_of_entry > spill_number) break;  // entries are in time order
    ++n_interactions;
    for (const TG4Trajectory &traj : event->Trajectories) {
      const int pdg = traj.GetPDGCode();
      if (std::abs(pdg) > 1000000000) continue;  // nuclei and fragments
      const TParticlePDG *part = pdg_db->GetParticle(pdg);
      if (!part || part->Charge() == 0) continue;
      bool reaches_tms = false;
      for (const TG4TrajectoryPoint &p : traj.Points)
        if (tms.IsInsideTMS(p.GetPosition().Vect())) { reaches_tms = true; break; }
      if (!reaches_tms) continue;
      if (!first) json << ",";
      first = false;
      ++n_written;
      const TLorentzVector &p4 = traj.GetInitialMomentum();
      json << "{\"entry\":" << i << ",\"track\":" << traj.GetTrackId() << ",\"parent\":" << traj.GetParentId()
           << ",\"pdg\":" << pdg << ",\"p_mev\":" << p4.P() << ",\"pts\":[";
      for (size_t k = 0; k < traj.Points.size(); ++k) {
        const TLorentzVector &x = traj.Points[k].GetPosition();
        if (k) json << ",";
        json << "[" << x.X() << "," << x.Y() << "," << x.Z() << "," << (x.T() - spill_t0) << "]";
      }
      json << "]}";
    }
  }
  json << "],\"n_interactions\":" << n_interactions << "}";
  json.close();
  std::cout << "Spill " << spill_number << ": " << n_interactions << " interactions, " << n_written
            << " charged trajectories reaching the TMS. Wrote " << out_path << "\n";
  return 0;
}
