// Lists every electronics channel of the detector geometry as "plane,bar,view" (the TMS_ChannelId
// the simulation and the bad channel list use), one per scintillator bar node. Used to draw bad
// channel lists. Usage: DumpChannels <input edep-sim file> [channels.csv]
// (or --gdml-file <file.gdml> in place of the input file). Needs TMS_DIR set.
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"

#include "TMS_Bar.h"
#include "TMS_Geom.h"
#include "TMS_Hit.h"
#include "TMS_Manager.h"

namespace {

struct Center { double x, y, z; };

// Collects the master-frame center of every scintillator bar node, as TMS_Bar::FindModules() uses it.
// Bars are only constructed after the walk, since building one moves the geometry navigator.
void Walk(TGeoManager* geom, const std::vector<std::string>& bar_names, std::vector<Center>& centers) {
  TGeoNode* node = geom->GetCurrentNode();
  const std::string name(node->GetName());
  for (const std::string& bar_name : bar_names) {
    if (name.find(bar_name) == std::string::npos) continue;
    TGeoBBox* box = dynamic_cast<TGeoBBox*>(node->GetVolume()->GetShape());
    if (!box) return;
    const double* origin = box->GetOrigin();
    double master[3];
    geom->GetCurrentMatrix()->LocalToMaster(origin, master);
    TMS_Geom::GetInstance().Scale(master[0], master[1], master[2]);
    centers.push_back({master[0], master[1], master[2]});
    return;
  }
  for (int i = 0; i < node->GetNdaughters(); ++i) {
    geom->CdDown(i);
    Walk(geom, bar_names, centers);
    geom->CdUp();
  }
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: DumpChannels <input edep-sim file | --gdml-file file.gdml> [channels.csv]\n";
    return 1;
  }
  TGeoManager* geom = nullptr;
  int next = 2;
  if (std::string(argv[1]) == "--gdml-file") {
    if (argc < 3) return 1;
    geom = TGeoManager::Import(argv[2]);
    next = 3;
  } else {
    TFile* input = TFile::Open(argv[1]);
    if (input && !input->IsZombie()) geom = (TGeoManager*)input->Get("EDepSimGeometry");
  }
  if (!geom) { std::cerr << "No geometry\n"; return 1; }
  const std::string out_name = argc > next ? argv[next] : "channels.csv";

  TMS_Geom::GetInstance().SetGeometry(geom);
  TMS_Manager& manager = TMS_Manager::GetInstance();
  const std::vector<std::string> bar_names = {manager.Get_GEOMETRY_VOLUME_ScintLayer(),
      manager.Get_GEOMETRY_VOLUME_ScintLayerOrtho(), manager.Get_GEOMETRY_VOLUME_ScintLayerParallel()};

  std::vector<Center> centers;
  geom->CdTop();
  Walk(geom, bar_names, centers);
  geom->CdTop();

  std::map<TMS_ChannelId, int> channels;
  int unplaced = 0;
  for (const Center& c : centers) {
    TMS_Bar bar(c.x, c.y, c.z);
    if (bar.GetPlaneNumber() < 0 || bar.GetBarNumber() < 0 || bar.GetBarType() == TMS_Bar::kError) { ++unplaced; continue; }
    channels[{bar.GetPlaneNumber(), bar.GetBarNumber(), bar.GetBarTypeNumber()}]++;
  }
  int duplicates = 0;
  for (const auto& kv : channels) if (kv.second > 1) ++duplicates;

  std::ofstream out(out_name);
  out << "plane,bar,view\n";
  for (const auto& kv : channels) out << kv.first.Plane << "," << kv.first.Bar << "," << kv.first.View << "\n";

  std::cout << centers.size() << " bar nodes, " << channels.size() << " distinct channels written to " << out_name
            << "; " << unplaced << " could not be placed, " << duplicates << " channel ids shared by more than one node\n";
  return (unplaced || duplicates) ? 2 : 0;
}
