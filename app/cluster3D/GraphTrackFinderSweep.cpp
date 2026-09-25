// Runs TMS_GraphTrackFinder against every known DBSCAN-failure muon from the
// 2026-09-07 15-file clustering benchmark (26 ND-LAr-fiducial muons where
// DBSCAN+PCA failed to recover a track-like cluster) and reports how many of
// them Graph Track Finder can find a clean partial path for -- the real "Phase 0"
// characterization the Link-and-Tree proposal calls for before committing to
// full Kalman-follower development.
//
// Shares its truth-matching convention (Parent-chain collapse, X-hit-first-
// then-Y point labeling) with ClusterTruthEfficiency.cpp and
// GraphTrackFinderSliceTest.cpp so "muon-owned" means the same thing everywhere.
//
// Input: a CSV of file,vgid,trackid rows (one per known failure -- see
// /media/usher/Drive2/DUNE/TMS/reports/2026-09-07_clustering_benchmark/
// benchmark_summary.txt for where these come from). Output: a per-case CSV
// plus an aggregate summary on stdout.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "TFile.h"
#include "TTree.h"

#include "TMS_GraphTrackFinder.h"
#include "TMS_LayerGrouping.h"
#include "SpacePointLayerInput.h"
#include "TMS_SpacePoint.h"

namespace {

const int kMaxSpacePoints = 10000;
const int kMaxTrueParticles = 20000;

struct TrueLabel {
  long long vgid = -1;
  int trackid = -999;
  bool Valid() const { return vgid >= 0; }
  bool operator==(const TrueLabel &o) const { return vgid == o.vgid && trackid == o.trackid; }
};

struct LabelHash {
  size_t operator()(const TrueLabel &l) const {
    return std::hash<long long>()(l.vgid) ^ (std::hash<int>()(l.trackid) << 1);
  }
};

struct SpillParticles {
  int n = 0;
  std::vector<long long> vgid;
  std::vector<int> trackid;
  std::vector<int> pdg;
  std::vector<int> parent_trackid;
  std::unordered_map<TrueLabel, int, LabelHash> index_of;
  std::vector<int> collapsed_trackid;
};

// See GraphTrackFinderSliceTest.cpp for the rationale (a muon's own delta-ray
// hits should count as the muon's, not "a different particle").
int CollapseTrackId(const SpillParticles &sp, int start_idx) {
  int idx = start_idx;
  int fallback_top_primary_trackid = sp.trackid[start_idx];
  int guard = 0;
  while (idx >= 0 && guard++ < 10000) {
    if (std::abs(sp.pdg[idx]) == 13) return sp.trackid[idx];
    fallback_top_primary_trackid = sp.trackid[idx];
    const int parent_tid = sp.parent_trackid[idx];
    if (parent_tid < 0) break;
    auto it = sp.index_of.find(TrueLabel{sp.vgid[idx], parent_tid});
    if (it == sp.index_of.end()) break;
    idx = it->second;
  }
  return fallback_top_primary_trackid;
}

// Everything this sweep needs to know about how one known failure fared.
struct CaseResult {
  bool ok = false;               // false if the file/target couldn't be loaded at all
  long long vgid = 0;
  int trackid = 0;
  int total_points = 0;
  int target_points_raw = 0;     // ghost-inflated raw space-point count
  int target_planes_total = 0;   // real distinct planes the target touches
  int best_path_points = 0;
  int best_path_planes = 0;
  double best_path_purity = 0.0; // purity (%) of whichever path has the most target-matched planes
  bool resource_limit_reached = false;
};

// Opens one RecoCandidates file, finds the slice where `target` has the most
// space points, builds real TMS_SpacePoints from it (native hit indices and
// truth labels straight from the already-written Reco_Tree branches -- no
// re-simulation needed), and runs the Graph Track Finder finder on it.
CaseResult RunOneCase(const std::string &input_filename, const TrueLabel &target,
                      const TMS_GraphTrackFinder::Config &config) {
  CaseResult res;
  res.vgid = target.vgid;
  res.trackid = target.trackid;

  TFile input(input_filename.c_str());
  if (input.IsZombie()) {
    std::cerr << "  [skip] failed to open " << input_filename << std::endl;
    return res;
  }
  TTree *reco_tree = (TTree *)input.Get("Reco_Tree");
  TTree *truth_spill = (TTree *)input.Get("Truth_Spill");
  if (!reco_tree || !truth_spill) {
    std::cerr << "  [skip] " << input_filename << " missing Reco_Tree/Truth_Spill" << std::endl;
    return res;
  }

  int spill_no_ts = 0, n_tp_ts = 0;
  static std::vector<long long> vgid_ts(kMaxTrueParticles);
  static std::vector<int> trackid_ts(kMaxTrueParticles);
  static std::vector<int> pdg_ts(kMaxTrueParticles);
  static std::vector<int> parent_ts(kMaxTrueParticles);
  truth_spill->SetBranchAddress("SpillNo", &spill_no_ts);
  truth_spill->SetBranchAddress("nTrueParticles", &n_tp_ts);
  truth_spill->SetBranchAddress("VertexGlobalID", vgid_ts.data());
  truth_spill->SetBranchAddress("TrackId", trackid_ts.data());
  truth_spill->SetBranchAddress("PDG", pdg_ts.data());
  truth_spill->SetBranchAddress("Parent", parent_ts.data());

  std::map<int, SpillParticles> spills;
  for (Long64_t e = 0; e < truth_spill->GetEntries(); ++e) {
    truth_spill->GetEntry(e);
    SpillParticles sp;
    sp.n = n_tp_ts;
    sp.vgid.assign(vgid_ts.begin(), vgid_ts.begin() + n_tp_ts);
    sp.trackid.assign(trackid_ts.begin(), trackid_ts.begin() + n_tp_ts);
    sp.pdg.assign(pdg_ts.begin(), pdg_ts.begin() + n_tp_ts);
    sp.parent_trackid.assign(parent_ts.begin(), parent_ts.begin() + n_tp_ts);
    for (int i = 0; i < n_tp_ts; ++i) sp.index_of[{sp.vgid[i], sp.trackid[i]}] = i;
    sp.collapsed_trackid.resize(n_tp_ts);
    for (int i = 0; i < n_tp_ts; ++i) sp.collapsed_trackid[i] = CollapseTrackId(sp, i);
    spills[spill_no_ts] = std::move(sp);
  }

  int n_space_points = 0, spill_no = 0, slice_no = 0;
  static std::vector<float> sp_x(kMaxSpacePoints), sp_y(kMaxSpacePoints), sp_z(kMaxSpacePoints);
  static std::vector<float> sp_time(kMaxSpacePoints);
  static std::vector<int> sp_x_hitidx(kMaxSpacePoints), sp_y_hitidx(kMaxSpacePoints);
  static std::vector<long long> sp_x_vgid(kMaxSpacePoints), sp_y_vgid(kMaxSpacePoints);
  static std::vector<int> sp_x_trackid(kMaxSpacePoints), sp_y_trackid(kMaxSpacePoints);
  reco_tree->SetBranchAddress("nSpacePoints", &n_space_points);
  reco_tree->SetBranchAddress("SpacePointX", sp_x.data());
  reco_tree->SetBranchAddress("SpacePointY", sp_y.data());
  reco_tree->SetBranchAddress("SpacePointZ", sp_z.data());
  const SpacePointLayerInput sp_layer(reco_tree, kMaxSpacePoints);
  reco_tree->SetBranchAddress("SpacePointTime", sp_time.data());
  reco_tree->SetBranchAddress("SpacePointXHitIndex", sp_x_hitidx.data());
  reco_tree->SetBranchAddress("SpacePointYHitIndex", sp_y_hitidx.data());
  reco_tree->SetBranchAddress("SpacePointXTrueVertexGlobalId", sp_x_vgid.data());
  reco_tree->SetBranchAddress("SpacePointXTrueTrackId", sp_x_trackid.data());
  reco_tree->SetBranchAddress("SpacePointYTrueVertexGlobalId", sp_y_vgid.data());
  reco_tree->SetBranchAddress("SpacePointYTrueTrackId", sp_y_trackid.data());
  reco_tree->SetBranchAddress("SpillNo", &spill_no);
  reco_tree->SetBranchAddress("SliceNo", &slice_no);

  Long64_t n_entries = reco_tree->GetEntries();
  int best_count = 0, best_total = 0;
  std::vector<TMS_SpacePoint> best_points;
  std::vector<TrueLabel> best_point_label;

  for (Long64_t entry = 0; entry < n_entries; ++entry) {
    reco_tree->GetEntry(entry);
    auto spill_it = spills.find(spill_no);
    if (spill_it == spills.end() || n_space_points <= 0) continue;
    const SpillParticles &sp = spill_it->second;

    auto collapse = [&](const TrueLabel &raw) -> TrueLabel {
      if (!raw.Valid()) return raw;
      auto it = sp.index_of.find(raw);
      if (it == sp.index_of.end()) return raw;
      return TrueLabel{raw.vgid, sp.collapsed_trackid[it->second]};
    };

    std::vector<TrueLabel> point_label(n_space_points);
    int count = 0;
    for (int i = 0; i < n_space_points; ++i) {
      TrueLabel x_label = collapse(TrueLabel{sp_x_vgid[i], sp_x_trackid[i]});
      TrueLabel y_label = collapse(TrueLabel{sp_y_vgid[i], sp_y_trackid[i]});
      point_label[i] = x_label.Valid() ? x_label : y_label;
      if (point_label[i] == target) ++count;
    }

    if (count > best_count) {
      best_count = count;
      best_total = n_space_points;
      best_points.clear();
      for (int i = 0; i < n_space_points; ++i) {
        best_points.push_back(TMS_SpacePoint(sp_x[i], sp_y[i], sp_z[i],
                                              sp_x_hitidx[i], sp_y_hitidx[i], sp_time[i],
                                              sp_layer.Layer(i, sp_z[i])));
      }
      best_point_label = point_label;
    }
  }

  if (best_count == 0) {
    std::cerr << "  [skip] target vgid=" << target.vgid << " trackid=" << target.trackid
              << " not found in " << input_filename << std::endl;
    return res;
  }

  res.total_points = best_total;
  res.target_points_raw = best_count;

  const TMS_GraphTrackFinder::Result result = TMS_GraphTrackFinder::Finder(config).Find(best_points);
  res.resource_limit_reached = result.Stats.ResourceLimitReached;

  const std::vector<int> z_layer = TMS_LayerGrouping::GroupIndexOfEachPoint(best_points, config.LayerZTolerance);
  std::set<int> target_layers_in_slice;
  for (std::size_t i = 0; i < best_points.size(); ++i) {
    if (best_point_label[i] == target) target_layers_in_slice.insert(z_layer[i]);
  }
  res.target_planes_total = static_cast<int>(target_layers_in_slice.size());

  // Report whichever path covers the most of the target's own planes (not
  // necessarily the top-scoring path -- a short, very pure path and a
  // longer, more complete one may both be useful outputs; here we want to
  // know the best the finder achieved at all, for this initial survey).
  int best_planes = 0, best_points_at_best = 0, best_matched_at_best = 0;
  for (const TMS_GraphTrackFinder::Path &path : result.Paths) {
    int matched = 0;
    std::set<int> matched_layers;
    for (std::size_t idx : path.SpacePointIndices) {
      if (best_point_label[idx] == target) {
        ++matched;
        matched_layers.insert(z_layer[idx]);
      }
    }
    if (static_cast<int>(matched_layers.size()) > best_planes) {
      best_planes = static_cast<int>(matched_layers.size());
      best_points_at_best = static_cast<int>(path.SpacePointIndices.size());
      best_matched_at_best = matched;
    }
  }
  res.best_path_planes = best_planes;
  res.best_path_points = best_points_at_best;
  res.best_path_purity = best_points_at_best > 0
      ? 100.0 * best_matched_at_best / best_points_at_best : 0.0;
  res.ok = true;
  return res;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc != 3 && argc != 5 && argc != 6 && argc != 8) {
    std::cerr << "Usage: " << argv[0]
              << " <targets.csv (file,vgid,trackid)> <output_results.csv>"
                 " [max_seed_layer_occupancy] [max_seed_hit_multiplicity]"
                 " [use_curvature_projection 0|1]"
                 " [occupancy_penalty] [hit_multiplicity_penalty]"
              << std::endl;
    return 1;
  }
  const std::string targets_path = argv[1];
  const std::string output_path = argv[2];
  // Real slices run to 100+ points/layer, far above what the seed gates were
  // implicitly tuned against (the ~9 points/layer synthetic test) -- see
  // linkandtree_validation.md: at the default 18/8, one known-viable case
  // (14 real planes, 204 space points) formed *zero* seeds anywhere in the
  // whole slice and so recovered nothing; raising these to 150/50 alone
  // turned it into an 8/14-plane, 100%-purity recovery. Exposed as CLI
  // overrides rather than changed defaults since retuning for real
  // occupancy scales is exactly the open question this sweep exists to
  // answer, not a settled decision.
  const bool has_seed_overrides = argc == 5 || argc == 6 || argc == 8;
  const std::size_t seed_occ_override = has_seed_overrides
      ? static_cast<std::size_t>(std::stoul(argv[3])) : 0;
  const std::size_t seed_mult_override = has_seed_overrides
      ? static_cast<std::size_t>(std::stoul(argv[4])) : 0;
  // Config::UseCurvatureProjection defaults to true -- without this override
  // this sweep would silently follow whatever that default currently is,
  // rather than deliberately choosing it, so make it explicit.
  const bool has_curvature_override = argc == 6 || argc == 8;
  const bool curvature_override = has_curvature_override ? (std::stoi(argv[5]) != 0) : true;
  const bool has_penalty_overrides = argc == 8;
  const double occ_penalty_override = has_penalty_overrides ? std::stod(argv[6]) : 0.0;
  const double mult_penalty_override = has_penalty_overrides ? std::stod(argv[7]) : 0.0;

  std::ifstream targets_file(targets_path);
  if (!targets_file) {
    std::cerr << "Could not open " << targets_path << std::endl;
    return 1;
  }
  std::string header;
  std::getline(targets_file, header);  // discard "file,vgid,trackid"

  TMS_GraphTrackFinder::Config config;  // default -- same as the synthetic/single-slice tests
  if (has_seed_overrides) {
    config.MaxSeedLayerOccupancy = seed_occ_override;
    config.MaxSeedHitMultiplicity = seed_mult_override;
  }
  if (has_curvature_override) config.UseCurvatureProjection = curvature_override;
  if (has_penalty_overrides) {
    config.OccupancyPenalty = occ_penalty_override;
    config.HitMultiplicityPenalty = mult_penalty_override;
  }

  std::ofstream out(output_path);
  out << "file,vgid,trackid,total_points,target_points_raw,target_planes_total,"
         "best_path_planes,best_path_points,best_path_purity,resource_limit_reached\n";

  int n_cases = 0, n_ok = 0, n_any_recovery = 0, n_majority_recovery = 0;
  int n_resource_limited = 0;
  double sum_plane_fraction = 0.0, sum_purity = 0.0;

  std::string line;
  while (std::getline(targets_file, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string file, vgid_str, trackid_str;
    std::getline(ss, file, ',');
    std::getline(ss, vgid_str, ',');
    std::getline(ss, trackid_str, ',');
    const TrueLabel target{std::stoll(vgid_str), std::stoi(trackid_str)};

    ++n_cases;
    std::cout << "[" << n_cases << "] " << file << " vgid=" << target.vgid
              << " trackid=" << target.trackid << " ... " << std::flush;
    const CaseResult res = RunOneCase(file, target, config);
    if (!res.ok) {
      std::cout << "FAILED TO LOAD" << std::endl;
      continue;
    }
    ++n_ok;
    if (res.resource_limit_reached) ++n_resource_limited;
    const double plane_fraction = res.target_planes_total > 0
        ? 100.0 * res.best_path_planes / res.target_planes_total : 0.0;
    sum_plane_fraction += plane_fraction;
    sum_purity += res.best_path_purity;
    if (res.best_path_planes > 0) ++n_any_recovery;
    if (plane_fraction >= 50.0) ++n_majority_recovery;

    std::cout << res.best_path_planes << "/" << res.target_planes_total << " planes ("
              << plane_fraction << "%), " << res.best_path_purity << "% purity"
              << (res.resource_limit_reached ? " [RESOURCE LIMITED]" : "") << std::endl;

    out << file << "," << target.vgid << "," << target.trackid << "," << res.total_points << ","
        << res.target_points_raw << "," << res.target_planes_total << "," << res.best_path_planes
        << "," << res.best_path_points << "," << res.best_path_purity << ","
        << (res.resource_limit_reached ? 1 : 0) << "\n";
  }

  std::cout << "\n=== Summary over " << n_ok << "/" << n_cases << " loaded cases ===\n"
            << "Any partial path recovered (>=1 plane): " << n_any_recovery << "/" << n_ok << "\n"
            << "Majority of target's planes recovered (>=50%): " << n_majority_recovery << "/"
            << n_ok << "\n"
            << "Mean plane-coverage fraction: " << (n_ok ? sum_plane_fraction / n_ok : 0.0) << "%\n"
            << "Mean best-path purity: " << (n_ok ? sum_purity / n_ok : 0.0) << "%\n"
            << "Cases that hit the hypothesis resource limit: " << n_resource_limited << "/" << n_ok
            << "\n"
            << "(For all " << n_ok << " cases, DBSCAN+PCA recovered 0 planes -- that's exactly why"
               " they're in this list: they're the 26 known ND-LAr-fiducial muon failures from the"
               " 2026-09-07 clustering benchmark.)\n";

  return 0;
}
