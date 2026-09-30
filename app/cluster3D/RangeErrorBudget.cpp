// Error budget for the range momentum of muons that stop in the TMS: how much of the
// reconstructed range momentum's spread is physics (range straggling, and our energy-loss model),
// how much is sampling the track only at the scintillator planes, and how much is reconstruction.
//
// For each muon of a Cluster3DRecoTruth muons.csv (ND-LAr box vertex, stops in the TMS; the row of
// the slice holding most of its hits), the momentum needed to traverse, and stop in, the material
// along its TRUE Geant4 trajectory (read from the raw TG4Event), walked backward from the 20 MeV/c
// floor with the same Bethe-Bloch-per-material calculation the Kalman follower uses, over:
//   p_true_path      the whole trajectory inside the TMS (entry to stopping point)
//   p_from_first_hit from the muon's first hit plane to its stopping point
//   p_hits           from its first hit plane to its last hit plane (what planes can sample)
//   p_straight       a straight line from the TMS entry point to the stopping point
// next to the true momentum entering the TMS and the reconstructed range momentum of its best
// Cluster3D track (both from muons.csv). The muon's hits come from the reco file's per-hit table
// (hits whose largest true contributor is the muon).
//
// Usage: RangeErrorBudget <spill.EDEPSIM_SPILLS.root> <reco.root> <muons.csv> <out.csv>
// Only muons.csv rows whose sourcefile is <reco.root> are used.
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TGeoManager.h"
#include "TTree.h"
#include "TVector3.h"

#include "BetheBloch.h"
#include "EDepSim/TG4Event.h"
#include "Material.h"
#include "TMS_Event.h"
#include "TMS_Geom.h"
#include "TMS_Manager.h"
#include "TMS_TrueParticle.h"
#include "TMS_VertexId.h"

namespace {

constexpr double kMinMomentumMeV = 20.0;  // the follower's floor (TMS_KalmanFollower.cpp)

// Longest energy-loss sub-step (mm); RANGE_BUDGET_STEP_MM overrides (1 mm by default).
double StepMM() {
  static const double step = std::getenv("RANGE_BUDGET_STEP_MM") ? std::atof(std::getenv("RANGE_BUDGET_STEP_MM")) : 1.0;
  return step;
}

// Scale on the Bethe-Bloch stopping power (RANGE_BUDGET_DEDX_SCALE, 1 by default): used to measure
// what scale would bring the true-path range momentum onto Geant4's (the floor read +1.7%, 2026-09-29).
double DEdxScale() {
  static const double s = std::getenv("RANGE_BUDGET_DEDX_SCALE") ? std::atof(std::getenv("RANGE_BUDGET_DEDX_SCALE")) : 1.0;
  return s;
}

// Momentum (MeV/c) at the start of a polyline that stops exactly at its end: energy loss through
// the geometry's materials, walked backward from the floor, as TMS_KalmanFollower's
// RangeMomentumMeV() does for a straight segment.
double RangeMomentum(const std::vector<TVector3> &points) {
  std::vector<std::pair<TGeoMaterial *, double> > steps;
  for (std::size_t i = 1; i < points.size(); ++i) {
    if ((points[i] - points[i - 1]).Mag() < 1e-6) continue;
    const auto seg = TMS_Geom::GetInstance().GetMaterials(points[i - 1], points[i]);
    steps.insert(steps.end(), seg.begin(), seg.end());
  }
  BetheBloch_Calculator bethe(Material::kPolyStyrene);
  double energy = std::sqrt(kMinMomentumMeV * kMinMomentumMeV + BetheBloch_Utils::Mm * BetheBloch_Utils::Mm);
  const double scale = TMS_Geom::GetInstance().Scale(1.0);
  for (auto it = steps.rbegin(); it != steps.rend(); ++it) {
    double density = it->first->GetDensity() / (CLHEP::g / CLHEP::cm3) / std::pow(scale, 3);
    const double thickness = TMS_Geom::GetInstance().Scale(it->second / 10.0);  // mm -> cm
    try {
      Material matter(density);
      bethe.fMaterial = matter;
    } catch (const std::invalid_argument &) {
      continue;
    }
    // Sub-steps of at most StepMM(): dE/dx is steep near the stopping point, so evaluating a
    // whole material step at the energy of its downstream end overestimates the loss.
    const int n = std::max(1, static_cast<int>(std::ceil(thickness * 10.0 / StepMM())));
    for (int k = 0; k < n; ++k) {
      const double loss = DEdxScale() * bethe.Calc_dEdx(energy) * density * thickness / n;
      if (std::isfinite(loss)) energy += loss;
    }
  }
  return BetheBloch_Utils::EnergyToMomentum(BetheBloch_Utils::Mm, energy);
}

// The part of a trajectory (in its own order) between z0 and z1: from where it first reaches z0 to
// where it last is at or below z1, with interpolated end points.
std::vector<TVector3> Clip(const std::vector<TVector3> &pts, double z0, double z1) {
  std::vector<TVector3> out;
  if (pts.size() < 2) return out;
  std::size_t first = pts.size(), last = 0;
  for (std::size_t i = 0; i < pts.size(); ++i)
    if (pts[i].Z() >= z0) { first = i; break; }
  for (std::size_t i = pts.size(); i-- > 0;)
    if (pts[i].Z() <= z1) { last = i; break; }
  if (first >= pts.size() || last < first) return out;
  auto interp = [](const TVector3 &a, const TVector3 &b, double z) {
    const double t = (b.Z() != a.Z()) ? (z - a.Z()) / (b.Z() - a.Z()) : 0.0;
    return a + t * (b - a);
  };
  if (first > 0) out.push_back(interp(pts[first - 1], pts[first], z0));
  for (std::size_t i = first; i <= last; ++i) out.push_back(pts[i]);
  if (last + 1 < pts.size() && pts[last + 1].Z() > z1) out.push_back(interp(pts[last], pts[last + 1], z1));
  return out;
}

std::vector<std::string> Split(const std::string &line) {
  std::vector<std::string> v;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) v.push_back(item);
  return v;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 5) {
    std::cerr << "Usage: " << argv[0] << " <spill.root> <reco.root> <muons.csv> <out.csv>\n";
    return 1;
  }
  const std::string spill_name = argv[1], reco_name = argv[2];
  TFile spill(spill_name.c_str());
  TGeoManager *geom = (TGeoManager *)spill.Get("EDepSimGeometry");
  if (!geom) { std::cerr << "No EDepSimGeometry in " << spill_name << "\n"; return 1; }
  TMS_Manager::GetInstance().SetFileName(spill_name);
  TMS_Geom::GetInstance().SetGeometry(geom);
  TTree *events = (TTree *)spill.Get("EDepSimEvents");
  TG4Event *event = nullptr;
  events->SetBranchAddress("Event", &event);

  // Muons: the main-slice row per muon, ND-LAr box, stopping in the TMS.
  std::ifstream in(argv[3]);
  std::string line;
  std::getline(in, line);
  std::map<std::string, int> col;
  {
    const auto head = Split(line);
    for (std::size_t i = 0; i < head.size(); ++i) col[head[i]] = i;
  }
  const std::string reco_base = reco_name.substr(reco_name.find_last_of('/') + 1);
  std::map<std::pair<long long, int>, std::vector<std::string> > muons;
  while (std::getline(in, line)) {
    const auto v = Split(line);
    const std::string src = v[col["sourcefile"]];
    if (src.substr(src.find_last_of('/') + 1) != reco_base) continue;
    if (v[col["vertex_in_lar_box"]] != "1" || v[col["stops_in_tms"]] != "1") continue;
    const auto key = std::make_pair(std::stoll(v[col["vertexglobalid"]]), std::stoi(v[col["trackid"]]));
    auto it = muons.find(key);
    if (it == muons.end() || std::stoi(v[col["true_hits_in_slice"]]) > std::stoi(it->second[col["true_hits_in_slice"]])) muons[key] = v;
  }
  std::cout << muons.size() << " stopping ND-LAr muons in " << reco_base << "\n";

  // Optional (RANGE_BUDGET_RECO_ENDS=<csv: sourcefile,vgid,trackid,reco_first_z,reco_last_z>): the
  // first and last hit z of each muon's reconstructed track, to also walk the TRUE path between the
  // RECONSTRUCTED ends -- separating endpoint errors from the reconstruction's walk itself.
  std::map<std::pair<long long, int>, std::pair<double, double> > recoEnds;
  if (const char *path = std::getenv("RANGE_BUDGET_RECO_ENDS")) {
    std::ifstream ends(path);
    std::string l;
    std::getline(ends, l);
    while (std::getline(ends, l)) {
      const auto v = Split(l);
      if (v.size() < 5 || v[0] != reco_base) continue;
      recoEnds[std::make_pair(std::stoll(v[1]), std::stoi(v[2]))] = std::make_pair(std::stod(v[3]), std::stod(v[4]));
    }
  }

  // Their hits' z range, from the reco file's per-hit table.
  TFile reco(reco_name.c_str());
  TTree *rt = (TTree *)reco.Get("Reco_Tree");
  const int kMax = 20000;
  static float h_z[kMax];
  static long long h_vg[kMax];
  static int h_tk[kMax], h_ped[kMax];
  int nh = 0;
  rt->SetBranchAddress("nSpacePointHits", &nh);
  rt->SetBranchAddress("SpacePointHitZ", h_z);
  rt->SetBranchAddress("SpacePointHitTrueVertexGlobalId", h_vg);
  rt->SetBranchAddress("SpacePointHitTrueTrackId", h_tk);
  rt->SetBranchAddress("SpacePointHitPedSup", h_ped);

  std::ofstream out(argv[4]);
  out << "sourcefile,vgid,trackid,p_true_enter,p_reco_range,found,true_last_hit_z,best_last_hit_z,"
         "z_first_hit,z_last_hit,z_entry,z_stop,p_true_path,p_from_first_hit,p_hits,p_straight,"
         "p_true_first_hit,p_path_first_hit,p_hits_first_hit,reco_first_z,reco_last_z,p_true_path_reco_ends\n";
  for (const auto &kv : muons) {
    const auto &v = kv.second;
    rt->GetEntry(std::stoll(v[col["entry"]]));
    double zmin = 1e30, zmax = -1e30;
    for (int h = 0; h < nh && h < kMax; ++h)
      if (!h_ped[h] && h_vg[h] == kv.first.first && h_tk[h] == kv.first.second) {
        zmin = std::min(zmin, (double)h_z[h]);
        zmax = std::max(zmax, (double)h_z[h]);
      }
    const int local_vertex = static_cast<int>(kv.first.first % TMS_VertexIdScale);
    if (local_vertex < 0 || local_vertex >= events->GetEntries() || zmin > zmax) continue;
    events->GetEntry(local_vertex);
    TMS_Event tms_event(*event, true);
    std::vector<TMS_TrueParticle> parts = tms_event.GetTrueParticles();
    std::vector<TVector3> pts, allPts, allMom;
    for (TMS_TrueParticle &part : parts)
      if (part.GetTrackId() == kv.first.second && std::abs(part.GetPDG()) == 13) {
        pts = part.GetPositionPoints(-1.0e6, 1.0e6, true);
        for (const TLorentzVector &x : part.GetPositionPoints()) allPts.push_back(x.Vect());
        allMom = part.GetMomentumPoints();
        break;
      }
    if (pts.size() < 2) continue;
    const double zentry = pts.front().Z(), zstop = pts.back().Z();
    const double p_path = RangeMomentum(pts);
    const double p_first = RangeMomentum(Clip(pts, zmin, 1e9));
    const double p_hits = RangeMomentum(Clip(pts, zmin, zmax));
    const double p_straight = RangeMomentum({pts.front(), pts.back()});
    // The reconstruction reports its range momentum at its first hit, but the truth's "momentum
    // entering the TMS" is at the first stored trajectory point inside the TMS box -- a median
    // 138 mm (84%: 268 mm) further in, since Geant4 stores points sparsely. So also: the true
    // momentum interpolated at the first hit's z along the full trajectory, and the range walks
    // from there.
    double p_true_first_hit = -1.0;
    if (allPts.size() == allMom.size())
      for (std::size_t i = 1; i < allPts.size(); ++i)
        if (allPts[i - 1].Z() <= zmin && allPts[i].Z() > zmin) {
          const double t = (zmin - allPts[i - 1].Z()) / (allPts[i].Z() - allPts[i - 1].Z());
          p_true_first_hit = allMom[i - 1].Mag() + t * (allMom[i].Mag() - allMom[i - 1].Mag());
          break;
        }
    const std::vector<TVector3> fromFirst = Clip(allPts, zmin, pts.back().Z());
    const double p_path_first_hit = fromFirst.size() >= 2 ? RangeMomentum(fromFirst) : -1.0;
    const std::vector<TVector3> hitsOnly = Clip(allPts, zmin, zmax);
    const double p_hits_first_hit = hitsOnly.size() >= 2 ? RangeMomentum(hitsOnly) : -1.0;
    double reco_first = -1.0, reco_last = -1.0, p_reco_ends = -1.0;
    {
      auto it = recoEnds.find(kv.first);
      if (it != recoEnds.end()) {
        reco_first = it->second.first;
        reco_last = it->second.second;
        const std::vector<TVector3> between = Clip(allPts, reco_first, reco_last);
        if (between.size() >= 2) p_reco_ends = RangeMomentum(between);
      }
    }
    if (std::getenv("RANGE_BUDGET_DEBUG")) {
      // Path length, and areal density per material, along the trajectory and along the chord.
      auto budget = [](const std::vector<TVector3> &p, const char *what) {
        double len = 0.0;
        std::map<std::string, double> gcm2;
        for (std::size_t i = 1; i < p.size(); ++i) {
          len += (p[i] - p[i - 1]).Mag();
          for (const auto &st : TMS_Geom::GetInstance().GetMaterials(p[i - 1], p[i]))
            gcm2[st.first->GetName()] += st.first->GetDensity() / (CLHEP::g / CLHEP::cm3) * st.second / 10.0;
        }
        std::cout << "   " << what << ": " << p.size() << " points, length " << len << " mm;";
        for (const auto &m : gcm2) std::cout << " " << m.first << " " << m.second << " g/cm2";
        std::cout << "\n";
      };
      std::cout << "muon " << kv.first.first << ":" << kv.first.second << " p_true_enter " << v[col["true_momentum_tms_mev"]]
                << " path " << p_path << " straight " << p_straight << "\n";
      budget(pts, "trajectory");
      budget({pts.front(), pts.back()}, "chord");
      int nonmono = 0;
      for (std::size_t i = 1; i < pts.size(); ++i) if (pts[i].Z() < pts[i - 1].Z()) ++nonmono;
      std::cout << "   z decreasing steps: " << nonmono << "; first " << pts.front().X() << "," << pts.front().Y() << "," << pts.front().Z()
                << " last " << pts.back().X() << "," << pts.back().Y() << "," << pts.back().Z() << "\n";
    }
    out << v[col["sourcefile"]] << "," << kv.first.first << "," << kv.first.second << "," << v[col["true_momentum_tms_mev"]] << ","
        << v[col["best_range_momentum_mev"]] << "," << v[col["found"]] << "," << v[col["true_last_hit_z"]] << ","
        << v[col["best_last_hit_z"]] << "," << zmin << "," << zmax << "," << zentry << "," << zstop << "," << p_path << ","
        << p_first << "," << p_hits << "," << p_straight << "," << p_true_first_hit << "," << p_path_first_hit << ","
        << p_hits_first_hit << "," << reco_first << "," << reco_last << "," << p_reco_ends << "\n";
  }
  return 0;
}
