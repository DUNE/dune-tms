// Single-stage scans of the detector simulation, for the stage-by-stage validation plots
// (scripts/Validation/DetSim/stages/). Like ArtificialResegmentationTest, synthetic muon crossings
// of one known bar are pushed through the real TMS_Event::FinalizeEvent(); the readout
// configuration (TMS_READOUT_TOML) selects the pipeline, so run once per configuration.
//
// Modes:
//   path  Perpendicular crossings of the reference bar at minimum-ionizing dE/dx, scanning the
//         path length through the bar (and, for a full crossing, the deposited energy). One row
//         per throw: the light before and after the threshold, the hit time, ToT and photon count.
//   position  Full-thickness MIP crossings of the reference bar at a series of distances from its
//         readout end, along the bar. One row per throw: the light and hit time, and the distance
//         from the readout as the simulation itself computes it (for the PE-vs-position check
//         against the PDR light-yield expectation, and the fiber propagation delay).
//   pair  Two crossings of the reference bar a time dt apart, scanning dt, plus the two halves of
//         a split X-bar crossed at once. One row per throw: the number of readouts.
//
// Usage: DetSimStageScan <path|position|pair> <edep_sim_file_for_geometry> <output_csv> [nthrows]

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
};

// Same reference bar as ArtificialResegmentationTest (plane 28, bar 28, view 1); a 9 mm path
// centered on it stays inside the 10 mm bar
const Crossing kRefBar = {-2434.968, -2302.8132, 12935.819, 0.0, 9.0, 9.0 * kMipdEdx_MeV_per_mm};
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
  TG4PrimaryParticle primary;
  primary.TrackId = kTrackId;
  primary.PDGCode = kPDGCodeMuMinus;
  primary.Momentum = TLorentzVector(0, 0, 1, 1);
  vtx.Particles.push_back(primary);
  event.Primaries.push_back(vtx);
  TG4Trajectory traj;
  traj.TrackId = kTrackId;
  traj.ParentId = -1;
  traj.PDGCode = kPDGCodeMuMinus;
  traj.Name = "mu-";
  TG4TrajectoryPoint p0, p1;
  p0.Position = vtx.Position;
  p1.Position = TLorentzVector(crossings.back().x, crossings.back().y, crossings.back().z + 100, crossings.back().t);
  traj.Points.push_back(p0);
  traj.Points.push_back(p1);
  event.Trajectories.push_back(traj);
  TG4HitSegmentContainer segs;
  for (const auto& c : crossings) {
    TG4HitSegment seg;
    seg.PrimaryId = kTrackId;
    seg.Contrib.push_back(kTrackId);
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

void RunPathScan(std::ofstream& out, int nthrows) {
  out << "path_mm,edep_MeV,throw_index,n_hits_total,n_hits_surviving,pe_all,pe_surviving,reco_energy,"
         "hit_time,true_time,tot,n_photons\n";
  // (path, energy scale): the path-length scan at MIP dE/dx, then brighter full crossings for ToT
  std::vector<std::pair<double, double>> cells;
  for (double p : {0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.5, 9.0}) cells.push_back({p, 1.0});
  for (double s : {2.0, 3.0, 4.0, 6.0, 8.0}) cells.push_back({9.0, s});
  int event_id = 0;
  for (const auto& cell : cells) {
    Crossing c = kRefBar;
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
  // Find the reference bar's type, length and readout end from a probe event, then walk the
  // crossing along the bar axis using the same readout-end definition as TMS_DetectorSimulation
  bool is_xbar = false;
  double length = 0, center = 0;
  {
    TG4Event event = BuildEvent({kRefBar}, 0);
    TMS_Event tms_event(event);
    tms_event.FinalizeEvent();
    auto hits = tms_event.GetHits(-1, /*include_ped_sup=*/true);
    if (hits.empty()) { std::cerr << "probe crossing made no hit\n"; return; }
    const TMS_Bar& bar = hits.front().GetBar();
    is_xbar = bar.GetBarType() == TMS_Bar::kXBar;
    length = bar.GetBarLength();
    center = bar.GetAxisReadoutCenter();
  }
  const TMS_Geom& geo = TMS_Geom::GetInstance();
  // Distance from the readout end d -> coordinate along the bar. The reference X-bar sits at x < 0.
  auto coordinate = [&](double d) {
    if (is_xbar) return (kRefBar.x < 0) ? geo.XBarNegReadoutLocation(center, length) + d
                                        : geo.XBarPosReadoutLocation(center, length) - d;
    return geo.YBarReadoutLocation(center, length) - d;
  };
  std::vector<double> distances;
  for (double f : {0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98}) distances.push_back(f * length);
  int event_id = 0;
  for (double d : distances) {
    Crossing c = kRefBar;
    (is_xbar ? c.x : c.y) = coordinate(d);
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
    std::cerr << "Usage: " << argv[0] << " <path|position|pair> <edep_sim_file_for_geometry> <output_csv> [nthrows]\n";
    return 1;
  }
  const std::string mode = argv[1];
  const std::string geom_file = argv[2];
  const std::string out_csv = argv[3];
  const int nthrows = (argc >= 5) ? std::atoi(argv[4]) : (mode == "pair" ? 200 : mode == "position" ? 1000 : 2000);
  if (mode != "path" && mode != "position" && mode != "pair") {
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
  if (mode == "path") RunPathScan(out, nthrows);
  else if (mode == "position") RunPositionScan(out, nthrows);
  else RunPairScan(out, nthrows);
  out.close();
  std::cout << "Wrote " << out_csv << "\n";
  return 0;
}
