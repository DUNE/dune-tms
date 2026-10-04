// Single-stage scans of the detector simulation, for the stage-by-stage validation plots
// (scripts/Validation/DetSim/stages/). Like ArtificialResegmentationTest, synthetic muon crossings
// of one known bar are pushed through the real TMS_Event::FinalizeEvent(); the readout
// configuration (TMS_READOUT_TOML) selects the pipeline, so run once per configuration.
//
// Modes:
//   path  Perpendicular crossings of the reference bar at minimum-ionizing dE/dx, scanning the
//         path length through the bar, up to its full 16 mm thickness (and, for a 9 mm crossing, the
//         deposited energy). One row
//         per throw: the light before and after the threshold, the hit time, ToT and photon count.
//   position  Full-thickness (16 mm) MIP crossings of the reference bar at a series of distances from its
//         readout end, along the bar. One row per throw: the light and hit time, and the distance
//         from the readout as the simulation itself computes it (for the PE-vs-position check
//         against the PDR light-yield expectation, and the fiber propagation delay).
//   pileup  Several particles (tracks 1..N) crossing the same bar at staggered times and with
//         different energies, to check how one channel combines them: readouts, summed light and
//         energy, hit time, and the light provenance (which particle made most of the light, and
//         which made the first photon).
//   pair  Two crossings of the reference bar a time dt apart, scanning dt, plus the two halves of
//         a split X-bar crossed at once. One row per throw: the number of readouts.
//
// Usage: DetSimStageScan <path|position|pileup|pair> <edep_sim_file_for_geometry> <output_csv> [nthrows]
//                        [distance_mm]   (path mode only: distance of the crossing from the readout end;
//                                         default is the reference position, 2.7 m)

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TGeoManager.h"
#include "TLorentzVector.h"

#include "EDepSim/TG4Event.h"
#include "EDepSim/TG4PrimaryVertex.h"
#include "EDepSim/TG4HitSegment.h"

#include "TMS_Bar.h"
#include "TMS_Event.h"
#include "TMS_Geom.h"
#include "TMS_Manager.h"

namespace {

constexpr int kTrackId = 1;
constexpr int kPDGCodeMuMinus = 13;
constexpr double kSpeedOfLight_mm_per_ns = 299.792458;
// Minimum-ionizing dE/dx of the reference crossing in ArtificialResegmentationTest
constexpr double kMipdEdx_MeV_per_mm = 1.8207009 / 10.164372;

struct Crossing {
  double x, y, z;  // mm, midpoint
  double t;        // ns
  double path_mm;  // along z
  double edep_MeV;
  int track_id = 1;
};

// Same reference bar as ArtificialResegmentationTest (plane 28, bar 28, view 1), but centered on the
// bar in z. The scintillator is 16 mm thick (z = 12924.5 to 12940.5 mm, center 12932.5, from the
// geometry), so a 16 mm path is a full perpendicular crossing; the 9 mm default path, kept for the
// scans that compare against ArtificialResegmentationTest, is 56% of it.
constexpr double kBarThickness_mm = 16.0;
const Crossing kRefBar = {-2434.968, -2302.8132, 12932.5, 0.0, 9.0, 9.0 * kMipdEdx_MeV_per_mm};
// The two halves of split X-bar row y=-744.23 in plane 81 (bar 128 at x<0, bar 129 at x>0)
const Crossing kXBarNeg = {-1000.0, -744.23, 18502.5, 0.0, 9.0, 9.0 * kMipdEdx_MeV_per_mm};
const Crossing kXBarPos = {1000.0, -744.23, 18502.5, 0.0, 9.0, 9.0 * kMipdEdx_MeV_per_mm};

TG4Event BuildEvent(const std::vector<Crossing>& crossings, int event_id) {
  TG4Event event;
  event.RunId = 1;
  event.EventId = event_id;
  TG4PrimaryVertex vtx;
  vtx.Position = TLorentzVector(crossings[0].x, crossings[0].y, crossings[0].z - 100, 0);
  vtx.GeneratorName = "DetSimStageScan";
  vtx.Reaction = "synthetic";
  // One muon per distinct track id, each with a trajectory from before its first crossing to
  // after its last
  std::map<int, std::pair<const Crossing*, const Crossing*>> tracks;
  for (const auto& c : crossings) {
    auto it = tracks.find(c.track_id);
    if (it == tracks.end()) tracks[c.track_id] = {&c, &c};
    else it->second.second = &c;
  }
  for (const auto& kv : tracks) {
    TG4PrimaryParticle primary;
    primary.TrackId = kv.first;
    primary.PDGCode = kPDGCodeMuMinus;
    primary.Momentum = TLorentzVector(0, 0, 1, 1);
    vtx.Particles.push_back(primary);
  }
  event.Primaries.push_back(vtx);
  for (const auto& kv : tracks) {
    TG4Trajectory traj;
    traj.TrackId = kv.first;
    traj.ParentId = -1;
    traj.PDGCode = kPDGCodeMuMinus;
    traj.Name = "mu-";
    TG4TrajectoryPoint p0, p1;
    p0.Position = vtx.Position;
    p1.Position = TLorentzVector(kv.second.second->x, kv.second.second->y, kv.second.second->z + 100, kv.second.second->t);
    traj.Points.push_back(p0);
    traj.Points.push_back(p1);
    event.Trajectories.push_back(traj);
  }
  TG4HitSegmentContainer segs;
  for (const auto& c : crossings) {
    TG4HitSegment seg;
    seg.PrimaryId = c.track_id;
    seg.Contrib.push_back(c.track_id);
    seg.EnergyDeposit = c.edep_MeV;
    seg.SecondaryDeposit = 0;
    seg.Start = TLorentzVector(c.x, c.y, c.z - c.path_mm / 2, c.t);
    seg.Stop = TLorentzVector(c.x, c.y, c.z + c.path_mm / 2, c.t + c.path_mm / kSpeedOfLight_mm_per_ns);
    seg.TrackLength = c.path_mm;
    segs.push_back(seg);
  }
  event.SegmentDetectors["volTMS"] = segs;
  return event;
}

// The reference bar's type, length and readout end, found from a probe crossing, and the
// coordinate along the bar that puts a crossing at a given distance from the readout end (the
// same readout-end definition as TMS_DetectorSimulation)
struct BarAxis {
  bool is_xbar = false;
  double length = 0, center = 0;
  bool Probe();
  void MoveTo(Crossing& c, double distance_mm) const;
};

bool BarAxis::Probe() {
  TG4Event event = BuildEvent({kRefBar}, 0);
  TMS_Event tms_event(event);
  tms_event.FinalizeEvent();
  auto hits = tms_event.GetHits(-1, /*include_ped_sup=*/true);
  if (hits.empty()) { std::cerr << "probe crossing made no hit\n"; return false; }
  const TMS_Bar& bar = hits.front().GetBar();
  is_xbar = bar.GetBarType() == TMS_Bar::kXBar;
  length = bar.GetBarLength();
  center = bar.GetAxisReadoutCenter();
  return true;
}

void BarAxis::MoveTo(Crossing& c, double d) const {
  const TMS_Geom& geo = TMS_Geom::GetInstance();
  if (is_xbar) c.x = (kRefBar.x < 0) ? geo.XBarNegReadoutLocation(center, length) + d
                                      : geo.XBarPosReadoutLocation(center, length) - d;
  else c.y = geo.YBarReadoutLocation(center, length) - d;
}

void RunPathScan(std::ofstream& out, int nthrows, double distance_mm) {
  out << "path_mm,edep_MeV,throw_index,n_hits_total,n_hits_surviving,pe_all,pe_surviving,reco_energy,"
         "hit_time,true_time,tot,n_photons\n";
  // (path, energy scale): the path-length scan at MIP dE/dx, then brighter 9 mm crossings for ToT
  std::vector<std::pair<double, double>> cells;
  for (double p : {0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.5, 9.0, 12.0, kBarThickness_mm}) cells.push_back({p, 1.0});
  for (double s : {2.0, 3.0, 4.0, 6.0, 8.0}) cells.push_back({9.0, s});
  BarAxis axis;
  if (distance_mm > 0 && !axis.Probe()) return;
  int event_id = 0;
  for (const auto& cell : cells) {
    Crossing c = kRefBar;
    if (distance_mm > 0) axis.MoveTo(c, distance_mm);
    c.path_mm = cell.first;
    c.edep_MeV = cell.first * kMipdEdx_MeV_per_mm * cell.second;
    for (int t = 0; t < nthrows; ++t) {
      TG4Event event = BuildEvent({c}, event_id++);
      TMS_Event tms_event(event);
      tms_event.FinalizeEvent();
      int n_total = 0, n_surviving = 0;
      double pe_all = 0, pe_surviving = 0, reco_e = 0;
      double hit_time = 0, true_time = -999, tot = -999;
      int n_photons = -999;
      double best_pe = -1;
      for (const auto& hit : tms_event.GetHits(-1, /*include_ped_sup=*/true)) {
        if (hit.GetBarNumber() < 0) continue;
        n_total++;
        pe_all += hit.GetPE();
        const TMS_TrueHit* true_hit = tms_event.GetTrueHit(hit.GetHitId());
        if (true_hit) {
          true_time = true_hit->GetT();
          // -999 with the legacy pipeline, which draws no individual photons
          if (true_hit->GetNPhotons() >= 0) n_photons = std::max(n_photons, 0) + true_hit->GetNPhotons();
        }
        if (hit.GetPedSup()) continue;
        n_surviving++;
        pe_surviving += hit.GetPE();
        reco_e += hit.GetE();
        if (hit.GetPE() > best_pe) {
          best_pe = hit.GetPE();
          hit_time = hit.GetT();
          tot = hit.GetToT();
        }
      }
      out << c.path_mm << "," << c.edep_MeV << "," << t << "," << n_total << "," << n_surviving << "," << pe_all
          << "," << pe_surviving << "," << reco_e << "," << hit_time << "," << true_time << "," << tot << ","
          << n_photons << "\n";
    }
  }
  std::cout << "path scan: " << cells.size() << " cells x " << nthrows << " throws\n";
}

void RunPositionScan(std::ofstream& out, int nthrows) {
  out << "distance_mm,x,y,bar_type,bar_number,bar_length_mm,throw_index,n_hits_surviving,pe_all,pe_surviving,"
         "hit_time,true_time\n";
  BarAxis axis;
  if (!axis.Probe()) return;
  const bool is_xbar = axis.is_xbar;
  const double length = axis.length;
  std::vector<double> distances;
  for (double f : {0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98}) distances.push_back(f * length);
  int event_id = 0;
  for (double d : distances) {
    Crossing c = kRefBar;
    c.path_mm = kBarThickness_mm;
    c.edep_MeV = kBarThickness_mm * kMipdEdx_MeV_per_mm;
    axis.MoveTo(c, d);
    for (int t = 0; t < nthrows; ++t) {
      TG4Event event = BuildEvent({c}, event_id++);
      TMS_Event tms_event(event);
      tms_event.FinalizeEvent();
      int n_surviving = 0, bar_number = -1;
      double pe_all = 0, pe_surviving = 0, hit_time = 0, true_time = -999;
      for (const auto& hit : tms_event.GetHits(-1, /*include_ped_sup=*/true)) {
        if (hit.GetBarNumber() < 0) continue;
        bar_number = hit.GetBarNumber();
        pe_all += hit.GetPE();
        const TMS_TrueHit* true_hit = tms_event.GetTrueHit(hit.GetHitId());
        if (true_hit) true_time = true_hit->GetT();
        if (hit.GetPedSup()) continue;
        n_surviving++;
        pe_surviving += hit.GetPE();
        hit_time = hit.GetT();
      }
      out << d << "," << c.x << "," << c.y << "," << (is_xbar ? "X" : "Y") << "," << bar_number << "," << length << ","
          << t << "," << n_surviving << "," << pe_all << "," << pe_surviving << "," << hit_time << "," << true_time << "\n";
    }
  }
  std::cout << "position scan: " << distances.size() << " positions x " << nthrows << " throws, bar type "
            << (is_xbar ? "X" : "Y") << ", length " << length << " mm\n";
}

void RunPileupScan(std::ofstream& out, int nthrows) {
  out << "pattern,n_particles,dt,throw_index,n_readouts,pe_surviving,reco_energy,hit_time,first_true_time,"
         "n_photons,top_light_id,top_light_share,first_photon_id,n_contrib,min_contrib_id,max_energy_id\n";
  // Energy of each particle's crossing, relative to a full MIP crossing, for each pattern
  struct Pattern { const char* name; std::vector<double> scale; };
  const std::vector<Pattern> patterns = {
      {"equal", {1, 1, 1, 1}},
      {"bright_first", {1, 0.5, 0.25, 0.125}},
      {"bright_last", {0.125, 0.25, 0.5, 1}}};
  const std::vector<double> delays = {0, 5, 20, 50, 100, 150, 250};
  int event_id = 0;
  auto run = [&](const Pattern& pat, int n, double dt) {
    std::vector<Crossing> crossings;
    int max_energy_id = 1;
    for (int i = 0; i < n; ++i) {
      Crossing c = kRefBar;
      c.path_mm = kBarThickness_mm;
      c.edep_MeV = kBarThickness_mm * kMipdEdx_MeV_per_mm * pat.scale[i];
      c.t = i * dt;
      c.track_id = i + 1;
      if (pat.scale[i] > pat.scale[max_energy_id - 1]) max_energy_id = i + 1;
      crossings.push_back(c);
    }
    for (int t = 0; t < nthrows; ++t) {
      TG4Event event = BuildEvent(crossings, event_id++);
      TMS_Event tms_event(event);
      tms_event.FinalizeEvent();
      int n_readouts = 0, n_photons = 0, top_id = -999, first_id = -999, n_contrib = 0, min_contrib_id = -999;
      double pe = 0, reco_e = 0, hit_time = 1e9, first_true = 1e9, top_share = -999, best_pe = -1;
      for (const auto& hit : tms_event.GetHits(-1, /*include_ped_sup=*/false)) {
        if (hit.GetBarNumber() < 0) continue;
        n_readouts++;
        pe += hit.GetPE();
        reco_e += hit.GetE();
        const TMS_TrueHit* true_hit = tms_event.GetTrueHit(hit.GetHitId());
        if (true_hit && true_hit->GetNPhotons() > 0) n_photons += true_hit->GetNPhotons();
        // Time and provenance of the brightest readout
        if (hit.GetPE() > best_pe) {
          best_pe = hit.GetPE();
          if (true_hit) {
            top_id = true_hit->GetPrimaryIdByLight();
            top_share = true_hit->GetLightShare();
            first_id = true_hit->GetFirstPhotonPrimaryId();
            n_contrib = static_cast<int>(true_hit->GetNLightContributions());
            for (size_t c = 0; c < true_hit->GetNLightContributions(); ++c) {
              const int id = true_hit->GetLightContributionPrimaryId(c);
              if (min_contrib_id == -999 || id < min_contrib_id) min_contrib_id = id;
            }
          }
        }
        if (hit.GetT() < hit_time) hit_time = hit.GetT();
        if (true_hit && true_hit->GetT() < first_true) first_true = true_hit->GetT();
      }
      if (n_readouts == 0) { hit_time = -999; first_true = -999; }
      out << pat.name << "," << n << "," << dt << "," << t << "," << n_readouts << "," << pe << "," << reco_e << ","
          << hit_time << "," << first_true << "," << n_photons << "," << top_id << "," << top_share << "," << first_id
          << "," << n_contrib << "," << min_contrib_id << "," << max_energy_id << "\n";
    }
  };
  int cells = 0;
  run(patterns[0], 1, 0); cells++;
  for (const auto& pat : patterns)
    for (int n : {2, 3, 4})
      for (double dt : delays) { run(pat, n, dt); cells++; }
  std::cout << "pileup scan: " << cells << " cells x " << nthrows << " throws\n";
}

void RunPairScan(std::ofstream& out, int nthrows) {
  out << "case,dt,throw_index,n_readouts,n_channels,pe_surviving\n";
  std::vector<double> delays;
  for (double dt = 0; dt <= 1000; dt += 10) delays.push_back(dt);
  int event_id = 0;
  auto run = [&](const std::string& label, double dt, const std::vector<Crossing>& crossings) {
    for (int t = 0; t < nthrows; ++t) {
      TG4Event event = BuildEvent(crossings, event_id++);
      TMS_Event tms_event(event);
      tms_event.FinalizeEvent();
      int n = 0;
      double pe = 0;
      std::map<TMS_ChannelId, int> channels;
      for (const auto& hit : tms_event.GetHits(-1, /*include_ped_sup=*/false)) {
        n++;
        pe += hit.GetPE();
        channels[hit.GetChannelId()]++;
      }
      out << label << "," << dt << "," << t << "," << n << "," << channels.size() << "," << pe << "\n";
    }
  };
  for (double dt : delays) {
    Crossing second = kRefBar;
    second.t = dt;
    run("same_bar", dt, {kRefBar, second});
  }
  run("xbar_halves", 0, {kXBarNeg, kXBarPos});
  run("single", 0, {kRefBar});
  std::cout << "pair scan: " << delays.size() + 2 << " cells x " << nthrows << " throws\n";
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <path|position|pileup|pair> <edep_sim_file_for_geometry> <output_csv> [nthrows] [distance_mm]\n";
    return 1;
  }
  const std::string mode = argv[1];
  const std::string geom_file = argv[2];
  const std::string out_csv = argv[3];
  const int nthrows = (argc >= 5) ? std::atoi(argv[4]) : (mode == "pair" ? 200 : mode == "position" ? 1000 : mode == "pileup" ? 200 : 2000);
  const double distance_mm = (argc >= 6) ? std::atof(argv[5]) : -1;
  if (mode != "path" && mode != "position" && mode != "pileup" && mode != "pair") {
    std::cerr << "Unknown mode " << mode << "\n";
    return 1;
  }

  TMS_Manager::GetInstance().SetFileName(geom_file);
  TFile* input = new TFile(geom_file.c_str(), "open");
  TGeoManager* geom = dynamic_cast<TGeoManager*>(input->Get("EDepSimGeometry"));
  if (!geom) {
    std::cerr << "Could not load EDepSimGeometry from " << geom_file << "\n";
    return 1;
  }
  TMS_Geom::GetInstance().SetGeometry(geom);

  std::ofstream out(out_csv);
  if (mode == "path") RunPathScan(out, nthrows, distance_mm);
  else if (mode == "position") RunPositionScan(out, nthrows);
  else if (mode == "pileup") RunPileupScan(out, nthrows);
  else RunPairScan(out, nthrows);
  out.close();
  std::cout << "Wrote " << out_csv << "\n";
  return 0;
}
