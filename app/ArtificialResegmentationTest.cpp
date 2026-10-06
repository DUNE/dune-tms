// Controlled artificial-resegmentation benchmark (Test A: synthetic uniform MIP).
//
// A fixed physical muon traversal (real bar, real geometry, uniform dE/dx, no delta
// ray) is represented as Nseg artificial TG4HitSegments for Nseg in {1,2,5,10,20}, and
// each representation is pushed through the exact same production detector-response
// code (TMS_Event::FinalizeEvent()) many times. Because only the bookkeeping
// segmentation changes -- physical truth is held fixed -- any resulting difference
// across Nseg is unambiguously a segmentation artifact, not physics variation. See
// reports/2026-09-13_segmentation_timing_benchmark/ for the observational precursor to
// this test and the study plan this implements.
//
// Generalized from the original single-scenario version (2026-09-13) to run several
// scenarios (varying crossing angle, path length, and deposited energy around the same
// validated reference point) and to emit truth-provenance columns
// (TrueHitPrimaryId/VertexId-by-energy, see TMS_TrueHit::GetPrimaryIdByEnergy()) so
// segmentation-invariance can be checked at the truth-attribution level too, not just
// PE/timing. This is intended as the Phase 2 "capture a baseline before the redesign
// lands" step of the detector-sim response-element redesign plan -- run it now against
// the untouched pipeline, then re-run after each later phase behind the
// Sim.DetSim.UseResponseElements config flag and diff against this baseline.
//
// NOT YET DONE: an optical-segment-length convergence scan (0.5/1/2/5mm, per the
// review's own suggestion) -- there is no such parameter in the simulation yet
// (Sim.DetSim.OpticalSegmentLength doesn't exist until the later re-segmentation phase
// lands), so that axis can't be exercised until then. Don't add a fake/no-op CLI knob
// for it now; add it when the parameter it's supposed to sweep actually exists.
//
// Usage: ArtificialResegmentationTest <edep_sim_file_for_geometry> <output_csv> [nthrows]

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
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

namespace {

// Reference physical traversal, taken from a real, already-validated single-segment
// (TrueNTrueParticles==1) ~10mm muon-dominated hit in the completed observational
// benchmark (before/0000001_RecoCandidates.root, Truth_Info entry index=3, hit=92):
// x=-2434.968 y=-2302.8132 z=12935.819 mm, dx=10.164372 mm, TrueHitE=1.8207009 MeV,
// bar=28 plane=28 view=1. All scenarios below share this exact (x,y,z) crossing
// midpoint, since it's confirmed to resolve to a real scintillator bar -- only the
// crossing direction, path length, and deposited energy vary per scenario.
constexpr double kRefX = -2434.968;
constexpr double kRefY = -2302.8132;
constexpr double kRefZCenter = 12935.819;
// The real hit this was taken from had dx=10.164372 mm, but centered on its z the 20-step
// version puts an end fragment's midpoint outside the bar, so that fragment is dropped before
// any detector response and ~5% of the energy silently vanishes (found 2026-09-21 via the
// true_energy_in_hits column; the ~5% PE deficit reported for this test on 2026-09-13 was this
// artifact, not a detector effect). Use a 9 mm crossing at the same dE/dx instead, which keeps
// every fragment midpoint inside the bar at all tested step counts.
// The scintillator bar is 16 mm thick in z (12924.5 to 12940.5 mm), so this is 9/16 = 56% of a
// perpendicular crossing; "full_thickness" in the scenario label is historical.
constexpr double kNominalPathLengthMM = 9.0;
constexpr double kNominalEdepMeV = 1.8207009 * (9.0 / 10.164372);
constexpr int kTrackId = 1;
constexpr int kPDGCodeMuMinus = 13;
constexpr double kSpeedOfLight_mm_per_ns = 299.792458;  // relativistic muon, beta~1

std::vector<int> NsegScan() { return {1, 2, 5, 10, 20}; }

struct Direction3 {
  double x, y, z;
};

Direction3 Normalized(Direction3 d) {
  const double n = std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
  return {d.x / n, d.y / n, d.z / n};
}

// One physical scenario: a crossing of kPathLengthMM through the reference point along
// kDirection, depositing kTotalEdepMeV total. All share the same validated midpoint.
struct Scenario {
  std::string label;
  Direction3 direction;
  double path_length_mm;
  double total_edep_MeV;
};

std::vector<Scenario> ScenarioList() {
  return {
    // Baseline: perpendicular full-thickness crossing, same as the original single-
    // scenario version of this test.
    {"perpendicular_full_thickness", Normalized({0, 0, 1}), kNominalPathLengthMM, kNominalEdepMeV},
    // Tilted crossing: same path length, angled ~25 degrees off the plane normal in x.
    // Probes whether segmentation artifacts depend on crossing angle, not just
    // thickness-direction crossings.
    {"tilted_25deg", Normalized({0.42, 0, 1}), kNominalPathLengthMM, kNominalEdepMeV},
    // Short/grazing crossing: targets the sub-2mm true-path-length regime flagged in
    // the Phase 0 benchmark as a persistent ~20% survival floor unchanged by prior
    // fixes -- this scenario lets that regime be probed under controlled resegmentation
    // too, once later phases land.
    {"grazing_short_path", Normalized({0, 0, 1}), 1.5, kNominalEdepMeV * (1.5 / kNominalPathLengthMM)},
    // Higher local dE/dx: same path length, more deposited energy, to check the
    // segmentation artifacts aren't specific to one energy scale.
    {"higher_dedx", Normalized({0, 0, 1}), kNominalPathLengthMM, kNominalEdepMeV * 2.0},
  };
}

TG4Event BuildSyntheticEvent(const Scenario& scenario, int nseg, int run_id, int event_id) {
  const Direction3 dir = scenario.direction;
  const double half = scenario.path_length_mm / 2.0;
  const TLorentzVector start(kRefX - half * dir.x, kRefY - half * dir.y, kRefZCenter - half * dir.z, 0.0);
  const TLorentzVector stop(kRefX + half * dir.x, kRefY + half * dir.y, kRefZCenter + half * dir.z,
                             scenario.path_length_mm / kSpeedOfLight_mm_per_ns);

  TG4Event event;
  event.RunId = run_id;
  event.EventId = event_id;

  // One primary vertex, one primary muon -- required so TrackId=kTrackId resolves
  // cleanly in TMS_Event::ProcessTG4Event()'s vertex-mapping (avoids the "track id not
  // found" fallback path).
  TG4PrimaryVertex vtx;
  vtx.Position = start;
  vtx.GeneratorName = "ArtificialResegmentationTest";
  vtx.Reaction = "synthetic";
  TG4PrimaryParticle primary;
  primary.TrackId = kTrackId;
  primary.PDGCode = kPDGCodeMuMinus;
  primary.Momentum = TLorentzVector(dir.x, dir.y, dir.z, 1);  // direction only matters qualitatively here
  vtx.Particles.push_back(primary);
  event.Primaries.push_back(vtx);

  // Matching trajectory so TMS_TrueParticle bookkeeping resolves.
  TG4Trajectory traj;
  traj.TrackId = kTrackId;
  traj.ParentId = -1;
  traj.PDGCode = kPDGCodeMuMinus;
  traj.Name = "mu-";
  TG4TrajectoryPoint p0, p1;
  p0.Position = start;
  p1.Position = stop;
  traj.Points.push_back(p0);
  traj.Points.push_back(p1);
  event.Trajectories.push_back(traj);

  // Nseg artificial, equal-length segments, each carrying Edep/Nseg (uniform dE/dx),
  // linearly interpolated position and constant-velocity time -- only the bookkeeping
  // boundaries change with Nseg, not the total physical content.
  TG4HitSegmentContainer segs;
  for (int i = 0; i < nseg; ++i) {
    const double f0 = static_cast<double>(i) / nseg;
    const double f1 = static_cast<double>(i + 1) / nseg;
    TG4HitSegment seg;
    seg.PrimaryId = kTrackId;
    seg.Contrib.push_back(kTrackId);
    seg.EnergyDeposit = scenario.total_edep_MeV / nseg;
    seg.SecondaryDeposit = 0;
    seg.Start = TLorentzVector(start.X() + f0 * (stop.X() - start.X()),
                                start.Y() + f0 * (stop.Y() - start.Y()),
                                start.Z() + f0 * (stop.Z() - start.Z()),
                                start.T() + f0 * (stop.T() - start.T()));
    seg.Stop = TLorentzVector(start.X() + f1 * (stop.X() - start.X()),
                               start.Y() + f1 * (stop.Y() - start.Y()),
                               start.Z() + f1 * (stop.Z() - start.Z()),
                               start.T() + f1 * (stop.T() - start.T()));
    seg.TrackLength = scenario.path_length_mm / nseg;
    segs.push_back(seg);
  }
  event.SegmentDetectors["volTMS"] = segs;

  return event;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <edep_sim_file_for_geometry> <output_csv> [nthrows]\n";
    return 1;
  }
  const std::string geom_file = argv[1];
  const std::string out_csv = argv[2];
  const int nthrows = (argc >= 4) ? std::atoi(argv[3]) : 5000;

  TMS_Manager::GetInstance().SetFileName(geom_file);

  TFile* input = new TFile(geom_file.c_str(), "open");
  TGeoManager* geom = dynamic_cast<TGeoManager*>(input->Get("EDepSimGeometry"));
  if (!geom) {
    std::cerr << "Could not load EDepSimGeometry from " << geom_file << "\n";
    return 1;
  }
  TMS_Geom::GetInstance().SetGeometry(geom);

  std::ofstream out(out_csv);
  out << "scenario,nseg,throw_index,n_hits_total,n_hits_surviving,total_pe,total_reco_energy,"
         "min_hit_time,max_hit_time,true_time,"
         "primary_id_by_energy,vertex_id_by_energy,n_contributors,min_contributor_energy_share,"
         "true_energy_in_hits,injected_energy\n";

  int global_event_id = 0;
  int n_diagnostic_failures = 0;
  for (const Scenario& scenario : ScenarioList()) {
    for (int nseg : NsegScan()) {
      for (int t = 0; t < nthrows; ++t) {
        TG4Event event = BuildSyntheticEvent(scenario, nseg, /*run_id=*/1, /*event_id=*/global_event_id++);
        TMS_Event tms_event(event);
        tms_event.FinalizeEvent();

        const auto hits = tms_event.GetHits(-1, /*include_ped_sup=*/true);
        int n_total = 0, n_surviving = 0;
        double total_pe = 0, total_e = 0;
        double min_t = 1e18, max_t = -1e18;
        double true_time = -999;
        // Total true energy that actually reached a hit (any bar, pedestal-suppressed or
        // not); compared with the injected energy this exposes fragments dropped before
        // the detector response, e.g. ones whose midpoint falls outside a bar.
        double true_energy_in_hits = 0;
        // Truth-provenance columns are read from the highest-PE surviving hit (the
        // "main" hit for this throw) -- with nseg>1 there can be more than one final
        // hit if resegmentation/merging fails, which n_hits_surviving already flags;
        // these columns track whether the *attributed* contributor is invariant too.
        int primary_id_by_energy = -1;
        long long vertex_id_by_energy = -1;
        int n_contributors = -1;
        double min_contributor_energy_share = -1;
        double best_pe = -1;
        for (const auto& hit : hits) {
          if (hit.GetBarNumber() < 0) continue;  // shouldn't happen given the reference position
          n_total++;
          if (!hit.GetPedSup()) {
            n_surviving++;
            total_pe += hit.GetPE();
            total_e += hit.GetE();
            min_t = std::min(min_t, hit.GetT());
            max_t = std::max(max_t, hit.GetT());
            if (hit.GetPE() > best_pe) {
              best_pe = hit.GetPE();
              const TMS_TrueHit* true_hit_for_main = tms_event.GetTrueHit(hit.GetHitId());
              if (true_hit_for_main) {
                primary_id_by_energy = true_hit_for_main->GetPrimaryIdByEnergy();
                vertex_id_by_energy = true_hit_for_main->GetVertexGlobalIdByEnergy();
                n_contributors = static_cast<int>(true_hit_for_main->GetNTrueParticles());
                min_contributor_energy_share = true_hit_for_main->GetEnergyShare(0);
                for (size_t i = 1; i < true_hit_for_main->GetNTrueParticles(); i++) {
                  min_contributor_energy_share = std::min(min_contributor_energy_share,
                                                            true_hit_for_main->GetEnergyShare(i));
                }
              }
            }
          }
          const TMS_TrueHit* true_hit = tms_event.GetTrueHit(hit.GetHitId());
          if (true_hit) {
            true_time = true_hit->GetT();
            true_energy_in_hits += true_hit->GetE();
          }
        }
        if (scenario.label == "perpendicular_full_thickness" && nseg == 1 && t == 0 && n_total == 0) {
          n_diagnostic_failures++;
          std::cerr << "WARNING: reference traversal produced zero hits at all -- check "
                       "that the reference (x,y,z) resolves to a real scintillator bar.\n";
        }
        out << scenario.label << "," << nseg << "," << t << "," << n_total << "," << n_surviving << ","
            << total_pe << "," << total_e << ","
            << (n_surviving ? min_t : 0.0) << "," << (n_surviving ? max_t : 0.0) << ","
            << true_time << ","
            << primary_id_by_energy << "," << vertex_id_by_energy << ","
            << n_contributors << "," << min_contributor_energy_share << ","
            << true_energy_in_hits << "," << scenario.total_edep_MeV << "\n";
      }
    }
  }
  out.close();

  const size_t total_throws = ScenarioList().size() * NsegScan().size() * nthrows;
  std::cout << "Wrote " << out_csv << " (" << total_throws << " throws)\n";
  if (n_diagnostic_failures > 0) {
    std::cerr << "ERROR: reference traversal never produced a hit -- results are invalid.\n";
    return 1;
  }
  return 0;
}
