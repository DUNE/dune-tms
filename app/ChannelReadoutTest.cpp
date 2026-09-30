// Targeted check of the per-channel readout model (TMS_DetectorSimulation::SimulateChannelReadout,
// used with Sim.DetSim.UseResponseElements=true): two muon crossings of known bars at a known time
// difference are pushed through the real TMS_Event::FinalizeEvent(), and the number of final
// readouts in the channels involved is compared with what the configured readout window, deadtime
// and zombie time imply. The readout configuration is read once per process (TMS_READOUT_TOML),
// so run this once per configuration to be tested.
//
// Cases:
// - same bar, second crossing inside / after the readout window (and, with a deadtime
//   configured, inside the deadtime, in the zombie time, and after it);
// - the two halves of a split X-bar crossed at the same time: separate channels (bars 128 and 129
//   of plane 81 share z and NotZ, so a (z, NotZ)-keyed merge wrongly combines them).
//
// Usage: ChannelReadoutTest <edep_sim_file_for_geometry> [nthrows]
// Exit status is the number of failing cases.

#include <cmath>
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

#include "TMS_Event.h"
#include "TMS_Geom.h"
#include "TMS_Manager.h"
#include "TMS_Readout_Manager.h"

namespace {

constexpr int kTrackId = 1;
constexpr int kPDGCodeMuMinus = 13;
constexpr double kPathLengthMM = 9.0;  // along z, inside a 10 mm bar
constexpr double kEdepMeV = 5.0;       // well above the pedestal threshold after attenuation

struct Crossing {
  double x, y, z;  // mm, midpoint
  double t;        // ns
};

// Same reference bar as ArtificialResegmentationTest (plane 28, bar 28, view 1)
const Crossing kRefBar = {-2434.968, -2302.8132, 12935.819, 0.0};
// The two halves of split X-bar row y=-744.23 in plane 81 (bar 128 at x<0, bar 129 at x>0)
const Crossing kXBarNeg = {-1000.0, -744.23, 18502.5, 0.0};
const Crossing kXBarPos = {1000.0, -744.23, 18502.5, 0.0};

TG4HitSegment MakeSegment(const Crossing& c) {
  TG4HitSegment seg;
  seg.PrimaryId = kTrackId;
  seg.Contrib.push_back(kTrackId);
  seg.EnergyDeposit = kEdepMeV;
  seg.SecondaryDeposit = 0;
  seg.Start = TLorentzVector(c.x, c.y, c.z - kPathLengthMM / 2, c.t);
  seg.Stop = TLorentzVector(c.x, c.y, c.z + kPathLengthMM / 2, c.t + 0.03);
  seg.TrackLength = kPathLengthMM;
  return seg;
}

TG4Event BuildEvent(const std::vector<Crossing>& crossings, int event_id) {
  TG4Event event;
  event.RunId = 1;
  event.EventId = event_id;
  TG4PrimaryVertex vtx;
  vtx.Position = TLorentzVector(crossings[0].x, crossings[0].y, crossings[0].z - 100, 0);
  vtx.GeneratorName = "ChannelReadoutTest";
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
  for (const auto& c : crossings) segs.push_back(MakeSegment(c));
  event.SegmentDetectors["volTMS"] = segs;
  return event;
}

// Readouts expected in one channel for two crossings dt apart, per the configured model
int ExpectedSameChannelReadouts(double dt) {
  const double readout = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_ReadoutTime();
  const double dead = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_Deadtime();
  const double zombie = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_ZombieTime();
  if (dt < readout) return 1;
  if (dead > 0 && dt < readout + dead) {
    if (zombie > 0 && dt >= readout + dead - zombie) return 2;
    return 1;
  }
  return 2;
}

struct Case {
  std::string label;
  std::vector<Crossing> crossings;
  int expected;  // readouts summed over the channels touched
};

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <edep_sim_file_for_geometry> [nthrows]\n";
    return 1;
  }
  const std::string geom_file = argv[1];
  const int nthrows = (argc >= 3) ? std::atoi(argv[2]) : 200;
  TMS_Manager::GetInstance().SetFileName(geom_file);
  TFile* input = new TFile(geom_file.c_str(), "open");
  TGeoManager* geom = dynamic_cast<TGeoManager*>(input->Get("EDepSimGeometry"));
  if (!geom) {
    std::cerr << "Could not load EDepSimGeometry from " << geom_file << "\n";
    return 1;
  }
  TMS_Geom::GetInstance().SetGeometry(geom);

  const bool response_elements = TMS_Readout_Manager::GetInstance().Get_Sim_DetSim_UseResponseElements();
  const double readout = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_ReadoutTime();
  const double dead = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_Deadtime();
  const double zombie = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_ZombieTime();
  std::cout << "UseResponseElements=" << response_elements << " ReadoutTime=" << readout
            << " Deadtime=" << dead << " ZombieTime=" << zombie << "\n";

  // Second-crossing delays, kept >= 30 ns from every window edge so photon-arrival and
  // electronic time jitter cannot move a crossing across one
  std::vector<double> delays = {50, readout + 80};
  if (dead > 0) {
    delays.push_back(readout + 0.5 * (dead - std::max(zombie, 0.0)));
    if (zombie > 0) delays.push_back(readout + dead - 0.5 * zombie);
    delays.push_back(readout + dead + 80);
  }
  std::vector<Case> cases;
  for (double dt : delays) {
    Crossing second = kRefBar;
    second.t = dt;
    cases.push_back({"same_bar_dt" + std::to_string(static_cast<int>(dt)), {kRefBar, second}, ExpectedSameChannelReadouts(dt)});
  }
  // Split X-bar halves at the same time: two channels. The default pipeline's (z, NotZ)-keyed
  // merge combines them (reported, not required, with the flag off).
  cases.push_back({"xbar_halves_dt0", {kXBarNeg, kXBarPos}, 2});

  int n_fail = 0;
  int event_id = 0;
  for (const Case& c : cases) {
    std::map<int, int> count_histogram;
    double pe_sum = 0;
    int n_channels_min = 1000, n_channels_max = 0;
    for (int t = 0; t < nthrows; ++t) {
      TG4Event event = BuildEvent(c.crossings, event_id++);
      TMS_Event tms_event(event);
      tms_event.FinalizeEvent();
      int n = 0;
      std::map<TMS_ChannelId, int> channels;
      for (const auto& hit : tms_event.GetHits(-1, /*include_ped_sup=*/false)) {
        n++;
        pe_sum += hit.GetPE();
        channels[hit.GetChannelId()]++;
      }
      count_histogram[n]++;
      n_channels_min = std::min(n_channels_min, static_cast<int>(channels.size()));
      n_channels_max = std::max(n_channels_max, static_cast<int>(channels.size()));
    }
    const bool pass = count_histogram.size() == 1 && count_histogram.begin()->first == c.expected;
    const bool required = response_elements || c.label.rfind("xbar", 0) != 0;
    if (!pass && required) n_fail++;
    std::cout << (pass ? "PASS  " : (required ? "FAIL  " : "info  ")) << c.label << ": expected " << c.expected
              << " readout(s); got";
    for (const auto& h : count_histogram) std::cout << " " << h.first << "x" << h.second;
    std::cout << " throws; channels " << n_channels_min << "-" << n_channels_max
              << "; mean total PE " << pe_sum / nthrows << "\n";
  }
  std::cout << (n_fail == 0 ? "ALL PASSED" : "FAILURES") << ": " << n_fail << " failing case(s)\n";
  return n_fail;
}
