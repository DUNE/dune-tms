// Dumps a true particle's own Geant4 trajectory points (TMS_TrueParticle's
// stepping points, read straight from the raw TG4Event -- no TMS detector
// response involved), restricted to inside the TMS. This is a genuinely
// different thing from what the Kalman event displays call "true_trajectory"
// there: that field is built from RECONSTRUCTED space points that are
// truth-matched on both views, which inherits the space-point builder's own
// combinatorial ghosting (one real hit can pair with several partners in an
// adjacent plane -- see TMS_SpacePointBuilder.cpp's own header comment) and
// so can zigzag between several real-but-ghosted candidates at one nominal
// z. This tool exists to give a real, continuous G4 path to compare against
// instead (2026-09-22, "seven clean muons" display).
//
// Usage: DumpTrueTrajectory <geom.root> <out.json> <vertexId1> [vertexId2 ...]
//
// Each Truth_Spill "VertexGlobalID" ends in a local VertexID (e.g.
// 115002871 -> 2871) -- confirmed empirically (2026-09-22, by matching
// Truth_Spill's own BirthPosition against every raw event's primary-vertex
// position) that this local VertexID equals the EDepSimEvents TTree entry
// number directly: one raw G4 event (one neutrino interaction) per entry,
// entry N holding the vertex with local id N. So no separate spill lookup is
// needed -- just read entry `vertexId` for each target.
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

#include "EDepSim/TG4Event.h"

#include "TMS_Event.h"
#include "TMS_Geom.h"
#include "TMS_Manager.h"
#include "TMS_TrueParticle.h"

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <geom.root> <out.json> <vertexId1> [vertexId2 ...]\n";
    return 1;
  }
  const std::string geom_filename = argv[1];
  const std::string out_path = argv[2];
  std::vector<int> vertex_ids;
  for (int i = 3; i < argc; ++i) vertex_ids.push_back(std::stoi(argv[i]));

  TFile *input = new TFile(geom_filename.c_str(), "open");
  if (!input || input->IsZombie()) {
    std::cerr << "Could not open " << geom_filename << "\n";
    return 1;
  }
  TGeoManager *geom = (TGeoManager *)input->Get("EDepSimGeometry");
  TMS_Manager::GetInstance().SetFileName(geom_filename);
  TMS_Geom::GetInstance().SetGeometry(geom);

  TTree *events = (TTree *)input->Get("EDepSimEvents");
  if (!events) {
    std::cerr << "No EDepSimEvents tree in " << geom_filename << "\n";
    return 1;
  }
  TG4Event *event = nullptr;
  events->SetBranchAddress("Event", &event);

  std::ofstream json(out_path);
  json << std::fixed << "{";
  bool first_muon = true;
  int found = 0;
  for (int vertex_id : vertex_ids) {
    if (vertex_id < 0 || vertex_id >= events->GetEntries()) {
      std::cerr << "vertexId " << vertex_id << " out of range (" << events->GetEntries() << " entries)\n";
      continue;
    }
    events->GetEntry(vertex_id);
    TMS_Event tms_event(*event, true);
    // GetPositionPoints(z_start,z_end,onlyInsideTMS) isn't const, so work on
    // a local copy rather than the event's own const-returned vector.
    std::vector<TMS_TrueParticle> true_particles = tms_event.GetTrueParticles();
    bool found_this_one = false;
    for (TMS_TrueParticle &part : true_particles) {
      if (part.GetTrackId() != 0) continue;
      if (std::abs(part.GetPDG()) != 13) continue;
      found_this_one = true;
      ++found;
      if (!first_muon) json << ",";
      first_muon = false;
      // Wide open z range (the whole detector plus margin) -- onlyInsideTMS
      // already does the real clipping via TMS_Geom::StaticIsInsideTMS.
      std::vector<TVector3> pts = part.GetPositionPoints(-1.0e6, 1.0e6, true);
      json << "\"" << vertex_id << "\":[";
      for (size_t i = 0; i < pts.size(); ++i) {
        if (i) json << ",";
        json << "[" << pts[i].X() << "," << pts[i].Y() << "," << pts[i].Z() << "]";
      }
      json << "]";
      std::cout << "vertexId=" << vertex_id << " pdg=" << part.GetPDG()
                << " raw G4 trajectory points=" << part.GetPositionPoints().size()
                << " inside-TMS points=" << pts.size() << "\n";
      break;  // exactly one primary muon per vertex expected here
    }
    if (!found_this_one) std::cerr << "vertexId " << vertex_id << ": no primary muon found\n";
  }
  json << "}";
  json.close();
  std::cout << "Found " << found << "/" << vertex_ids.size() << " target muons. Wrote " << out_path << "\n";
  return 0;
}
