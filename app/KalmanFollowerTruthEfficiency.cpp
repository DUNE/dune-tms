// Full-population truth-based validation of the Kalman follower stage, on
// top of the exact same three-stage DBSCAN+PCA -> merged-cluster PCA ->
// Graph Track Finder pipeline GraphTrackFinderTruthEfficiency.cpp validates
// (same tolerances, same plurality-vote matching, same merge logic -- code
// duplicated rather than shared because the two tools' downstream needs
// diverge enough, per-muon, to make a shared loop harder to follow than two
// parallel ones). Whichever stage finds a muon's track-like object, that
// object's own points are handed to the Kalman follower: DBSCAN-direct and
// merged-cluster-PCA objects (unordered, no directed search behind them) go
// through TMS_KalmanFollower::Follower::RunBestSeed()'s multi-hypothesis
// first-z-layer seeding; the Graph Track Finder's own already-ordered best
// path goes through plain Run(), since its directed search already resolved
// this same first-point ambiguity.
//
// Scope, deliberately: efficiency / purity / completeness / ambiguity-
// resolution accuracy only. Momentum and charge resolution are real open
// questions (chi2 gate threshold, initial covariance realism -- Phase 2 of
// the project plan) but are NOT computed here; this tool's job is to answer
// "does the follower reliably find and correctly resolve the right hits,"
// not "how precise is the fitted momentum."

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "TFile.h"
#include "TGeoManager.h"
#include "TTree.h"

#include "TMS_FieldModel.h"
#include "TMS_Geom.h"
#include "TMS_GraphTrackFinder.h"
#include "TMS_KalmanFollower.h"
#include "TMS_SpacePoint.h"
#include "TMS_SpacePointCluster.h"
#include "TMS_SpacePointDBScan.h"

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
  std::vector<bool> tms_fiducial_start;
  std::vector<bool> lar_fiducial_start;
  std::unordered_map<TrueLabel, int, LabelHash> index_of;
  std::vector<int> collapsed_trackid;
};

// Same convention as ClusterTruthEfficiency.cpp / GraphTrackFinderTruthEfficiency.cpp.
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

// Same z-layer grouping TMS_GraphTrackFinder::Finder / TMS_KalmanFollower use internally.
std::vector<int> AssignZLayers(const std::vector<TMS_SpacePoint> &points, double tolerance) {
  std::vector<std::size_t> order(points.size());
  for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&points](std::size_t a, std::size_t b) {
    return points[a].GetZ() < points[b].GetZ();
  });
  std::vector<int> layer(points.size());
  int current_layer = -1;
  double layer_start_z = 0.0;
  for (std::size_t idx : order) {
    if (current_layer < 0 || points[idx].GetZ() - layer_start_z > tolerance) {
      ++current_layer;
      layer_start_z = points[idx].GetZ();
    }
    layer[idx] = current_layer;
  }
  return layer;
}

// Per-muon Kalman cross-check result, filled by RunFollowerAndScore() below --
// bundles the CSV columns so the three call sites (DBSCAN-direct, merged-PCA,
// Graph Track Finder) all fill them identically.
struct KalmanScore {
  bool ran = false;
  std::string seed_source;
  bool converged = false;
  int nodes_total = 0;
  int nodes_with_hit = 0;
  int gaps = 0;
  int correct_chosen = 0;
  int wrong_chosen = 0;
  int planes_covered = 0;
  int ambiguous_layers = 0;
  int ambiguous_truth_present = 0;
  int ambiguous_correct = 0;
  // Breaks down every target plane NOT in planes_covered by why: outside the
  // walked range entirely (before the seed's own first point -- the walk is
  // forward-only, see TMS_KalmanFollower.cpp's Run() -- or past
  // Config::MaxLayersBeyondSeed) vs. genuinely missed while walking through
  // the range (chi2-gate rejection, wrong pick, or a real gap). Added to
  // check whether the completeness ceiling is structural (out-of-range) or
  // a fit-quality problem (missed-in-range) -- see kalman_follower memory,
  // "investigate the completeness ceiling" (2026-09-15).
  int planes_before_walk = 0;
  int planes_after_walk = 0;
  int planes_missed_in_range = 0;
  std::string stop_reason = "not_started";
  // Of the gap layers (HasHit==false), how many had the truth-matched point
  // sitting right there among CandidateIndices but rejected by the chi2
  // gate (a tuning problem) vs. truth genuinely absent from that layer's
  // whole-slice candidate pool (a real hit-finding/reconstruction gap,
  // unrecoverable by any Follower::Config change).
  int gaps_truth_available = 0;
  int gaps_truth_absent = 0;
};

std::string StopReasonName(TMS_KalmanFollower::FitResult::StopReason r) {
  switch (r) {
    case TMS_KalmanFollower::FitResult::StopReason::ReachedRangeEnd: return "reached_range_end";
    case TMS_KalmanFollower::FitResult::StopReason::GapLimitExceeded: return "gap_limit_exceeded";
    case TMS_KalmanFollower::FitResult::StopReason::Diverged: return "diverged";
    default: return "not_started";
  }
}

KalmanScore ScoreFit(const TMS_KalmanFollower::FitResult &fit, const std::vector<TrueLabel> &point_label,
                     const TrueLabel &target, const std::vector<int> &z_layer_whole_slice,
                     const std::string &seed_source, const std::set<int> &target_layers_in_slice) {
  KalmanScore score;
  score.ran = true;
  score.seed_source = seed_source;
  score.converged = fit.Converged;
  score.nodes_total = (int)fit.Nodes.size();

  std::set<int> covered_layers;
  for (const TMS_KalmanFollower::FollowedNode &node : fit.Nodes) {
    if (!node.HasHit) {
      ++score.gaps;
      bool truth_present = false;
      for (std::size_t idx : node.CandidateIndices)
        if (point_label[idx] == target) truth_present = true;
      if (truth_present)
        ++score.gaps_truth_available;
      else
        ++score.gaps_truth_absent;
      continue;
    }
    ++score.nodes_with_hit;
    const bool correct = point_label[node.ChosenSpacePointIndex] == target;
    if (correct) {
      ++score.correct_chosen;
      covered_layers.insert(z_layer_whole_slice[node.ChosenSpacePointIndex]);
    } else {
      ++score.wrong_chosen;
    }
    if (node.CandidateIndices.size() > 1) {
      ++score.ambiguous_layers;
      bool truth_present = false;
      for (std::size_t idx : node.CandidateIndices)
        if (point_label[idx] == target) truth_present = true;
      if (truth_present) {
        ++score.ambiguous_truth_present;
        if (correct) ++score.ambiguous_correct;
      }
    }
  }
  score.planes_covered = (int)covered_layers.size();
  score.stop_reason = StopReasonName(fit.Stop);

  const int walk_start = fit.Nodes.empty() ? -1 : fit.Nodes.front().Layer;
  const int walk_end = fit.Nodes.empty() ? -1 : fit.Nodes.back().Layer;
  for (int layer : target_layers_in_slice) {
    if (fit.Nodes.empty() || layer < walk_start)
      ++score.planes_before_walk;
    else if (layer > walk_end)
      ++score.planes_after_walk;
  }
  score.planes_missed_in_range =
      (int)target_layers_in_slice.size() - score.planes_covered - score.planes_before_walk - score.planes_after_walk;
  return score;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 4 || argc > 5) {
    std::cerr << "Usage: " << argv[0]
              << " <edep_sim_geom_file> <input_reco_tree.root> <muons_output.csv> [append 0|1]"
              << std::endl;
    return -1;
  }

  const std::string geom_filename = argv[1];
  const std::string input_filename = argv[2];
  const std::string muons_csv_path = argv[3];
  const bool append = argc == 5 && std::stoi(argv[4]) != 0;

  // DBSCAN+PCA params: identical to ClusterTruthEfficiency.cpp / GraphTrackFinderTruthEfficiency.cpp.
  const int base_transverse_bars = 1;
  const int transverse_bars_per_plane_gap = 1;
  const int max_plane_gap = 3;
  const unsigned int min_points = 5;
  const double kLinearityThreshold = 0.8;
  const size_t kMinClusterSizeForTrack = 5;

  // Graph Track Finder fallback config: same validated real-data config as
  // GraphTrackFinderTruthEfficiency.cpp / KalmanFollowerSliceTest.cpp.
  TMS_GraphTrackFinder::Config lt_config;
  lt_config.MaxSeedLayerOccupancy = 150;
  lt_config.MaxSeedHitMultiplicity = 50;
  lt_config.OccupancyPenalty = 0.0;
  lt_config.HitMultiplicityPenalty = 0.0;
  lt_config.UseCurvatureProjection = false;

  TFile geom_input(geom_filename.c_str());
  if (geom_input.IsZombie()) {
    std::cerr << "Failed to open geometry source file: " << geom_filename << std::endl;
    return -1;
  }
  TGeoManager *geom = (TGeoManager *)geom_input.Get("EDepSimGeometry");
  if (!geom) {
    std::cerr << "Geometry source file is missing 'EDepSimGeometry': " << geom_filename << std::endl;
    return -1;
  }
  TMS_Geom::GetInstance().SetGeometry(geom);
  const double max_plane_pitch = TMS_Geom::GetInstance().GetMaxPlanePitch();
  const double bar_pitch = TMS_Geom::GetInstance().GetMaxBarPitch();
  if (max_plane_pitch <= 0 || bar_pitch <= 0) {
    std::cerr << "TMS_Geom found fewer than 2 surveyed planes or bars -- cannot derive a clustering tolerance."
              << std::endl;
    return -1;
  }
  const double worst_case_transverse = (base_transverse_bars + max_plane_gap * transverse_bars_per_plane_gap) * bar_pitch;
  const double broad_phase_radius =
      std::sqrt(worst_case_transverse * worst_case_transverse +
                std::pow(max_plane_pitch * (max_plane_gap + 1), 2));

  // Kalman follower: default Config, real (GDML-confirmed) 1.0T region field.
  // Stateless -- constructed once, reused (as a const&) for every muon.
  const RegionFieldModel field;
  const TMS_KalmanFollower::Config follower_config;
  const TMS_KalmanFollower::Follower follower(follower_config, field);

  TFile input(input_filename.c_str());
  if (input.IsZombie()) {
    std::cerr << "Failed to open input file: " << input_filename << std::endl;
    return -1;
  }
  TTree *reco_tree = (TTree *)input.Get("Reco_Tree");
  TTree *truth_info = (TTree *)input.Get("Truth_Info");
  TTree *truth_spill = (TTree *)input.Get("Truth_Spill");
  if (!reco_tree || !truth_info || !truth_spill) {
    std::cerr << "Input file is missing Reco_Tree/Truth_Info/Truth_Spill" << std::endl;
    return -1;
  }

  int spill_no_ts = 0, n_tp_ts = 0;
  static std::vector<long long> vgid_ts(kMaxTrueParticles);
  static std::vector<int> trackid_ts(kMaxTrueParticles);
  static std::vector<int> pdg_ts(kMaxTrueParticles);
  static std::vector<int> parent_ts(kMaxTrueParticles);
  static bool tms_fid_start_ts[kMaxTrueParticles];
  static bool lar_fid_start_ts[kMaxTrueParticles];
  truth_spill->SetBranchAddress("SpillNo", &spill_no_ts);
  truth_spill->SetBranchAddress("nTrueParticles", &n_tp_ts);
  truth_spill->SetBranchAddress("VertexGlobalID", vgid_ts.data());
  truth_spill->SetBranchAddress("TrackId", trackid_ts.data());
  truth_spill->SetBranchAddress("PDG", pdg_ts.data());
  truth_spill->SetBranchAddress("Parent", parent_ts.data());
  truth_spill->SetBranchAddress("TMSFiducialStart", tms_fid_start_ts);
  truth_spill->SetBranchAddress("LArFiducialStart", lar_fid_start_ts);

  std::map<int, SpillParticles> spills;
  for (Long64_t e = 0; e < truth_spill->GetEntries(); ++e) {
    truth_spill->GetEntry(e);
    SpillParticles sp;
    sp.n = n_tp_ts;
    sp.vgid.assign(vgid_ts.begin(), vgid_ts.begin() + n_tp_ts);
    sp.trackid.assign(trackid_ts.begin(), trackid_ts.begin() + n_tp_ts);
    sp.pdg.assign(pdg_ts.begin(), pdg_ts.begin() + n_tp_ts);
    sp.parent_trackid.assign(parent_ts.begin(), parent_ts.begin() + n_tp_ts);
    sp.tms_fiducial_start.assign(tms_fid_start_ts, tms_fid_start_ts + n_tp_ts);
    sp.lar_fiducial_start.assign(lar_fid_start_ts, lar_fid_start_ts + n_tp_ts);
    for (int i = 0; i < n_tp_ts; ++i) sp.index_of[{sp.vgid[i], sp.trackid[i]}] = i;
    sp.collapsed_trackid.resize(n_tp_ts);
    for (int i = 0; i < n_tp_ts; ++i) sp.collapsed_trackid[i] = CollapseTrackId(sp, i);
    spills[spill_no_ts] = std::move(sp);
  }
  std::cout << "Loaded Truth_Spill: " << spills.size() << " spills" << std::endl;

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
  reco_tree->SetBranchAddress("SpacePointTime", sp_time.data());
  reco_tree->SetBranchAddress("SpacePointXHitIndex", sp_x_hitidx.data());
  reco_tree->SetBranchAddress("SpacePointYHitIndex", sp_y_hitidx.data());
  reco_tree->SetBranchAddress("SpacePointXTrueVertexGlobalId", sp_x_vgid.data());
  reco_tree->SetBranchAddress("SpacePointXTrueTrackId", sp_x_trackid.data());
  reco_tree->SetBranchAddress("SpacePointYTrueVertexGlobalId", sp_y_vgid.data());
  reco_tree->SetBranchAddress("SpacePointYTrueTrackId", sp_y_trackid.data());
  reco_tree->SetBranchAddress("SpillNo", &spill_no);
  reco_tree->SetBranchAddress("SliceNo", &slice_no);

  int n_tp_ti = 0;
  static std::vector<int> true_nhits_slice(kMaxTrueParticles);
  truth_info->SetBranchAddress("nTrueParticles", &n_tp_ti);
  truth_info->SetBranchAddress("TrueNHitsInSlice", true_nhits_slice.data());

  std::ofstream muons_csv;
  if (append) {
    muons_csv.open(muons_csv_path, std::ios::app);
  } else {
    muons_csv.open(muons_csv_path);
    muons_csv << "sourcefile,entry,slice,spill,vertexglobalid,trackid,vertex_in_tms,vertex_in_lar_fiducial,"
                 "true_hits_in_slice,n_muons_in_slice,n_spacepoints_total,target_planes_total,"
                 "found_combined,"
                 "kalman_ran,kalman_seed_source,kalman_converged,"
                 "kalman_nodes_total,kalman_nodes_with_hit,kalman_gaps,"
                 "kalman_correct_chosen,kalman_wrong_chosen,kalman_purity_pct,"
                 "kalman_planes_covered,kalman_completeness_pct,"
                 "kalman_planes_before_walk,kalman_planes_after_walk,kalman_planes_missed_in_range,"
                 "kalman_stop_reason,kalman_gaps_truth_available,kalman_gaps_truth_absent,"
                 "kalman_ambiguous_layers,kalman_ambiguous_truth_present,kalman_ambiguous_correct,"
                 "probe_ran,probe_merged_size,probe_best_planes_covered,probe_best_purity_pct\n";
  }

  Long64_t n_entries = reco_tree->GetEntries();
  long n_slices_skipped_mismatch = 0, n_slices_seen = 0, n_slices_skipped_no_muon = 0;
  long n_muons_total = 0, n_found_combined = 0, n_kalman_ran = 0, n_kalman_converged = 0;

  for (Long64_t entry = 0; entry < n_entries; ++entry) {
    reco_tree->GetEntry(entry);
    truth_info->GetEntry(entry);

    auto spill_it = spills.find(spill_no);
    if (spill_it == spills.end()) continue;
    const SpillParticles &sp = spill_it->second;
    if (sp.n != n_tp_ti) {
      ++n_slices_skipped_mismatch;
      continue;
    }
    if (n_space_points <= 0) continue;

    std::vector<int> muon_particle_idx;
    for (int i = 0; i < sp.n; ++i) {
      if (std::abs(sp.pdg[i]) == 13 && true_nhits_slice[i] >= (int)min_points) {
        muon_particle_idx.push_back(i);
      }
    }
    if (muon_particle_idx.empty()) {
      ++n_slices_skipped_no_muon;
      continue;
    }
    ++n_slices_seen;

    auto collapse = [&](const TrueLabel &raw) -> TrueLabel {
      if (!raw.Valid()) return raw;
      auto it = sp.index_of.find(raw);
      if (it == sp.index_of.end()) return raw;
      return TrueLabel{raw.vgid, sp.collapsed_trackid[it->second]};
    };
    std::vector<TrueLabel> point_label(n_space_points);
    for (int i = 0; i < n_space_points; ++i) {
      const TrueLabel x_label = collapse(TrueLabel{sp_x_vgid[i], sp_x_trackid[i]});
      const TrueLabel y_label = collapse(TrueLabel{sp_y_vgid[i], sp_y_trackid[i]});
      point_label[i] = x_label.Valid() ? x_label : y_label;
    }

    const int n_muons_in_slice = (int)muon_particle_idx.size();

    // --- Pass A: DBSCAN+PCA, exactly as GraphTrackFinderTruthEfficiency.cpp. ---
    std::vector<TMS_SpacePoint> space_points;
    space_points.reserve(n_space_points);
    for (int i = 0; i < n_space_points; ++i) {
      space_points.emplace_back(sp_x[i], sp_y[i], sp_z[i], sp_x_hitidx[i], sp_y_hitidx[i], sp_time[i]);
    }
    std::vector<int> plane_index;
    plane_index.reserve(n_space_points);
    for (int i = 0; i < n_space_points; ++i) {
      plane_index.push_back(TMS_Geom::GetInstance().GetPlaneIndexNearestZ(sp_z[i]));
    }
    TMS_SpacePointDBScan dbscan(space_points, plane_index, min_points, bar_pitch, base_transverse_bars,
                                 transverse_bars_per_plane_gap, max_plane_gap, broad_phase_radius);
    std::vector<std::vector<int>> cluster_indices = dbscan.RunAndGetClusterIndices();
    std::vector<TMS_SpacePointCluster> clusters;
    clusters.reserve(cluster_indices.size());
    for (auto &indices : cluster_indices) clusters.emplace_back(space_points, indices);

    std::unordered_map<TrueLabel, std::vector<int>, LabelHash> muon_matches;
    for (size_t c = 0; c < clusters.size(); ++c) {
      const auto &cl = clusters[c];
      if (!cl.IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack)) continue;
      std::unordered_map<TrueLabel, int, LabelHash> votes;
      for (int idx : cluster_indices[c])
        if (point_label[idx].Valid()) votes[point_label[idx]]++;
      TrueLabel owner;
      int owner_count = 0;
      for (auto &kv : votes)
        if (kv.second > owner_count) { owner = kv.first; owner_count = kv.second; }
      int owner_pdg = 0;
      auto owner_idx_it = sp.index_of.find(owner);
      if (owner.Valid() && owner_idx_it != sp.index_of.end()) owner_pdg = sp.pdg[owner_idx_it->second];
      if (owner.Valid() && std::abs(owner_pdg) == 13) muon_matches[owner].push_back((int)c);
    }

    std::vector<int> point_cluster_id(n_space_points, 0);
    for (size_t c = 0; c < cluster_indices.size(); ++c)
      for (int idx : cluster_indices[c]) point_cluster_id[idx] = (int)c + 1;

    const std::vector<int> z_layer_whole_slice = AssignZLayers(space_points, lt_config.LayerZTolerance);

    // --- Pass B: merge-touching-clusters-plus-own-noise-and-re-PCA fallback,
    // identical logic to GraphTrackFinderTruthEfficiency.cpp. ---
    struct FallbackRun {
      TMS_GraphTrackFinder::Result result;
      std::vector<int> local_to_global;
      std::vector<int> z_layer_local;
    };

    std::unordered_map<int, std::vector<int>> merged_clusters;
    for (int pidx : muon_particle_idx) {
      TrueLabel label{sp.vgid[pidx], sp.trackid[pidx]};
      if (muon_matches.find(label) != muon_matches.end()) continue;

      std::set<int> touched_cluster_ids;
      std::vector<int> own_noise_points;
      for (int i = 0; i < n_space_points; ++i) {
        if (!(point_label[i] == label)) continue;
        const int cid = point_cluster_id[i];
        if (cid == 0) own_noise_points.push_back(i);
        else touched_cluster_ids.insert(cid);
      }
      std::vector<int> merged_indices = own_noise_points;
      for (int cid : touched_cluster_ids)
        for (int idx : cluster_indices[cid - 1]) merged_indices.push_back(idx);
      std::sort(merged_indices.begin(), merged_indices.end());
      merged_indices.erase(std::unique(merged_indices.begin(), merged_indices.end()), merged_indices.end());

      merged_clusters[pidx] = std::move(merged_indices);
    }

    std::unordered_map<int, FallbackRun> fallback_runs;
    std::unordered_map<int, bool> found_merged_by_pidx;
    for (auto &kv : merged_clusters) {
      const int pidx = kv.first;
      const std::vector<int> &merged_indices = kv.second;
      TrueLabel label{sp.vgid[pidx], sp.trackid[pidx]};
      found_merged_by_pidx[pidx] = false;
      if (merged_indices.size() < kMinClusterSizeForTrack) continue;

      TMS_SpacePointCluster merged_cluster(space_points, merged_indices);
      if (merged_cluster.IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack)) {
        std::unordered_map<TrueLabel, int, LabelHash> votes;
        for (int idx : merged_indices) if (point_label[idx].Valid()) votes[point_label[idx]]++;
        TrueLabel owner; int owner_count = 0;
        for (auto &v : votes) if (v.second > owner_count) { owner = v.first; owner_count = v.second; }
        if (owner == label) { found_merged_by_pidx[pidx] = true; continue; }
      }

      if (merged_indices.size() < lt_config.SeedLength) continue;
      std::vector<TMS_SpacePoint> local_points;
      local_points.reserve(merged_indices.size());
      for (int gi : merged_indices) local_points.push_back(space_points[gi]);
      FallbackRun run;
      run.z_layer_local = AssignZLayers(local_points, lt_config.LayerZTolerance);
      run.result = TMS_GraphTrackFinder::Finder(lt_config).Find(local_points);
      run.local_to_global = merged_indices;
      fallback_runs[pidx] = std::move(run);
    }

    // --- Probe (2026-09-16): what would GraphTrackFinder recover if run on
    // every cluster touching a muon's own points, REGARDLESS of whether
    // Stage 1/2's whole-cluster ownership check already succeeded? A real
    // case (vgid=1000000092001580, file 9) showed Stage 1 succeeding
    // trivially on a small isolated 4-plane tail fragment (25% completeness)
    // while the SAME GraphTrackFinder, run on the full touching-cluster set,
    // recovered 14/16 planes (87.5%) at 87.5% purity from a companion-
    // dominated cluster Stage 1's ownership gate had excluded -- because
    // Stage 2/3 are both skipped outright the instant Stage 1 succeeds
    // (`if (seedPath.empty())`), GraphTrackFinder never gets the chance to
    // even try decomposing that bigger cluster. Deliberately a SEPARATE,
    // parallel computation (not reusing/modifying merged_clusters/
    // fallback_runs above) so this shadow metric cannot alter the existing
    // validated found_dbscan/found_dbscan_merged/found_lt/kalman_* results
    // -- pure measurement, zero risk to the numbers this tool already
    // reports elsewhere.
    struct ProbeResult {
      bool ran = false;
      int merged_size = 0;
      int best_planes_covered = 0;
      double best_purity_pct = 0.0;
      std::vector<int> best_path_global_indices;  // ordered (GraphTrackFinder's own path order)
    };
    std::unordered_map<int, ProbeResult> probe_results;
    for (int pidx : muon_particle_idx) {
      TrueLabel label{sp.vgid[pidx], sp.trackid[pidx]};
      ProbeResult probe;

      std::set<int> touched_cluster_ids;
      std::vector<int> own_noise_points;
      for (int i = 0; i < n_space_points; ++i) {
        if (!(point_label[i] == label)) continue;
        const int cid = point_cluster_id[i];
        if (cid == 0) own_noise_points.push_back(i);
        else touched_cluster_ids.insert(cid);
      }
      std::vector<int> probe_indices = own_noise_points;
      for (int cid : touched_cluster_ids)
        for (int idx : cluster_indices[cid - 1]) probe_indices.push_back(idx);
      std::sort(probe_indices.begin(), probe_indices.end());
      probe_indices.erase(std::unique(probe_indices.begin(), probe_indices.end()), probe_indices.end());
      probe.merged_size = (int)probe_indices.size();

      if (probe_indices.size() >= lt_config.SeedLength) {
        std::vector<TMS_SpacePoint> probe_points;
        probe_points.reserve(probe_indices.size());
        for (int gi : probe_indices) probe_points.push_back(space_points[gi]);
        const TMS_GraphTrackFinder::Result probeResult = TMS_GraphTrackFinder::Finder(lt_config).Find(probe_points);
        probe.ran = true;
        for (const TMS_GraphTrackFinder::Path &path : probeResult.Paths) {
          std::size_t matched = 0;
          std::set<int> planesCovered;
          for (std::size_t li : path.SpacePointIndices) {
            const int gi = probe_indices[li];
            if (point_label[gi] == label) {
              ++matched;
              planesCovered.insert(z_layer_whole_slice[gi]);
            }
          }
          if ((int)planesCovered.size() > probe.best_planes_covered) {
            probe.best_planes_covered = (int)planesCovered.size();
            probe.best_purity_pct = path.SpacePointIndices.empty() ? 0.0 : 100.0 * matched / path.SpacePointIndices.size();
            probe.best_path_global_indices.clear();
            probe.best_path_global_indices.reserve(path.SpacePointIndices.size());
            for (std::size_t li : path.SpacePointIndices) probe.best_path_global_indices.push_back(probe_indices[li]);
          }
        }
      }
      probe_results[pidx] = probe;
    }

    for (int pidx : muon_particle_idx) {
      ++n_muons_total;
      TrueLabel label{sp.vgid[pidx], sp.trackid[pidx]};

      std::set<int> target_layers_in_slice;
      for (int i = 0; i < n_space_points; ++i)
        if (point_label[i] == label) target_layers_in_slice.insert(z_layer_whole_slice[i]);

      auto match_it = muon_matches.find(label);
      const bool found_dbscan = match_it != muon_matches.end();
      int best_dbscan_cluster = -1;
      if (found_dbscan) {
        size_t best_size = 0;
        for (int c : match_it->second)
          if (best_dbscan_cluster == -1 || clusters[c].GetSize() > best_size) {
            best_size = clusters[c].GetSize();
            best_dbscan_cluster = c;
          }
      }

      bool found_dbscan_merged = false;
      if (!found_dbscan) {
        auto it = found_merged_by_pidx.find(pidx);
        if (it != found_merged_by_pidx.end()) found_dbscan_merged = it->second;
      }

      bool found_lt = false;
      const TMS_GraphTrackFinder::Path *best_lt_path = nullptr;
      const FallbackRun *lt_run = nullptr;
      int best_lt_planes = 0;
      if (!found_dbscan && !found_dbscan_merged) {
        auto run_it = fallback_runs.find(pidx);
        if (run_it != fallback_runs.end()) {
          lt_run = &run_it->second;
          for (const TMS_GraphTrackFinder::Path &path : lt_run->result.Paths) {
            std::set<int> matched_layers;
            for (std::size_t local_idx : path.SpacePointIndices) {
              const int global_idx = lt_run->local_to_global[local_idx];
              if (point_label[global_idx] == label) matched_layers.insert(lt_run->z_layer_local[local_idx]);
            }
            if ((int)matched_layers.size() > best_lt_planes) {
              best_lt_planes = (int)matched_layers.size();
              best_lt_path = &path;
            }
          }
          found_lt = best_lt_planes > 0;
        }
      }

      // --- Which of the three original stages (if any) found something,
      // and how many of the target's own planes does ITS candidate object
      // actually touch? Needed to compare against the probe below on equal
      // footing (plane coverage of the raw candidate, before any Kalman
      // fit) -- see the 2026-09-16 "GraphTrackFinder never gets a chance to
      // run once Stage 1 succeeds" finding: Stage 1/2 succeeding at all,
      // even on a small isolated fragment, used to end the search here.
      auto CountOwnPlanes = [&](const std::vector<int> &idxs) {
        std::set<int> planes;
        for (int gi : idxs)
          if (point_label[gi] == label) planes.insert(z_layer_whole_slice[gi]);
        return (int)planes.size();
      };

      std::vector<int> stageCandidateIndices;
      std::string stageSeedSource;
      bool haveStageCandidate = false;
      if (found_dbscan) {
        stageCandidateIndices.assign(cluster_indices[best_dbscan_cluster].begin(), cluster_indices[best_dbscan_cluster].end());
        stageSeedSource = "dbscan_direct";
        haveStageCandidate = true;
      } else if (found_dbscan_merged) {
        stageCandidateIndices = merged_clusters[pidx];
        stageSeedSource = "merged_pca";
        haveStageCandidate = true;
      } else if (found_lt && best_lt_path != nullptr && best_lt_path->SpacePointIndices.size() >= 2) {
        stageCandidateIndices.reserve(best_lt_path->SpacePointIndices.size());
        for (std::size_t local_idx : best_lt_path->SpacePointIndices)
          stageCandidateIndices.push_back(lt_run->local_to_global[local_idx]);
        stageSeedSource = "graphtrack";
        haveStageCandidate = true;
      }
      const int stageCandidatePlanes = haveStageCandidate ? CountOwnPlanes(stageCandidateIndices) : 0;

      // --- Fix (2026-09-16): always compare against GraphTrackFinder run on
      // every cluster touching the target's own points, regardless of
      // whether an earlier stage already "succeeded" -- validated at full
      // population scale to recover 68.5%->94.1% completeness (ND-LAr-
      // fiducial: 70.5%->97.4%), at 92.6% mean purity on the paths that win,
      // for a real ~10x runtime cost (GraphTrackFinder now runs for
      // essentially every muon, not just the ~4% that used to reach it). ---
      const ProbeResult &probe = probe_results[pidx];
      const bool probeWins = probe.ran && probe.best_planes_covered > stageCandidatePlanes;

      std::vector<int> finalIndices;
      std::string finalSeedSource;
      bool haveFinal = false;
      if (probeWins) {
        finalIndices = probe.best_path_global_indices;
        finalSeedSource = "graphtrack_probe";
        haveFinal = true;
      } else if (haveStageCandidate) {
        finalIndices = stageCandidateIndices;
        finalSeedSource = stageSeedSource;
        haveFinal = true;
      }

      const bool found_combined = haveFinal;
      if (found_combined) ++n_found_combined;

      // --- Kalman follower: whichever candidate won above, fit it, always
      // via RunBestSeed()'s multi-hypothesis first-layer seeding -- even for
      // GraphTrackFinder-sourced (ordered) paths. Previously ordered paths
      // used plain Run(), trusting GraphTrackFinder's own directed search to
      // have already resolved the first-point ambiguity; that assumption
      // held when Stage 3 only ran as a last resort with nothing better
      // available, but broke down once it started winning competitively
      // against Stage 1/2 (the new graphtrack_probe pathway above) -- found
      // 2026-09-16 on a real case where GraphTrackFinder's own path anchored
      // its first TWO points on a companion particle (not the target),
      // seeding the whole fit's initial direction wrong and making it
      // confidently track the wrong particle for several layers before
      // losing the thread. RunBestSeed() tries every candidate at the
      // object's own first z-layer as an alternate seed hypothesis instead
      // of trusting a single one -- same machinery already validated for
      // DBSCAN-direct/merged-PCA seeds, just no longer withheld from
      // graphtrack-sourced ones. ---
      KalmanScore kscore;
      if (haveFinal) {
        const std::vector<std::size_t> objectIndices(finalIndices.begin(), finalIndices.end());
        const TMS_KalmanFollower::FitResult fit = follower.RunBestSeed(space_points, objectIndices);
        kscore = ScoreFit(fit, point_label, label, z_layer_whole_slice, finalSeedSource, target_layers_in_slice);
      }
      if (kscore.ran) {
        ++n_kalman_ran;
        if (kscore.converged) ++n_kalman_converged;
      }

      const double kalman_purity_pct =
          (kscore.correct_chosen + kscore.wrong_chosen) > 0
              ? 100.0 * kscore.correct_chosen / (kscore.correct_chosen + kscore.wrong_chosen)
              : 0.0;
      const double kalman_completeness_pct =
          !target_layers_in_slice.empty() ? 100.0 * kscore.planes_covered / target_layers_in_slice.size() : 0.0;

      muons_csv << input_filename << "," << entry << "," << slice_no << "," << spill_no << ","
                << label.vgid << "," << label.trackid << ","
                << (sp.tms_fiducial_start[pidx] ? 1 : 0) << "," << (sp.lar_fiducial_start[pidx] ? 1 : 0) << ","
                << true_nhits_slice[pidx] << "," << n_muons_in_slice << "," << n_space_points << ","
                << (int)target_layers_in_slice.size() << ","
                << (found_combined ? 1 : 0) << ","
                << (kscore.ran ? 1 : 0) << "," << kscore.seed_source << "," << (kscore.converged ? 1 : 0) << ","
                << kscore.nodes_total << "," << kscore.nodes_with_hit << "," << kscore.gaps << ","
                << kscore.correct_chosen << "," << kscore.wrong_chosen << "," << kalman_purity_pct << ","
                << kscore.planes_covered << "," << kalman_completeness_pct << ","
                << kscore.planes_before_walk << "," << kscore.planes_after_walk << "," << kscore.planes_missed_in_range << ","
                << kscore.stop_reason << "," << kscore.gaps_truth_available << "," << kscore.gaps_truth_absent << ","
                << kscore.ambiguous_layers << "," << kscore.ambiguous_truth_present << "," << kscore.ambiguous_correct << ","
                << (probe_results[pidx].ran ? 1 : 0) << "," << probe_results[pidx].merged_size << ","
                << probe_results[pidx].best_planes_covered << "," << probe_results[pidx].best_purity_pct
                << "\n";
    }

    if (n_slices_seen % 25 == 0) {
      std::cout << "  entry=" << entry << "/" << n_entries << " spill=" << spill_no << " slice=" << slice_no
                << " nSP=" << n_space_points << " muons_in_slice=" << n_muons_in_slice
                << " slices_seen=" << n_slices_seen << " kalman_ran=" << n_kalman_ran << std::endl;
    }
  }

  muons_csv.close();

  std::cout << "Done. " << n_muons_total << " muon candidates." << std::endl;
  std::cout << "Found (combined pipeline): " << n_found_combined << " ("
            << (n_muons_total > 0 ? 100.0 * n_found_combined / n_muons_total : 0.0) << "%)" << std::endl;
  std::cout << "Kalman follower ran: " << n_kalman_ran << ", converged: " << n_kalman_converged << " ("
            << (n_kalman_ran > 0 ? 100.0 * n_kalman_converged / n_kalman_ran : 0.0) << "%)" << std::endl;
  std::cout << "Slices seen (>=1 findable muon): " << n_slices_seen
            << ", skipped (no findable muon): " << n_slices_skipped_no_muon
            << ", skipped (nTrueParticles mismatch): " << n_slices_skipped_mismatch << std::endl;
  std::cout << "Wrote " << muons_csv_path << (append ? " (appended)" : "") << std::endl;

  return 0;
}
