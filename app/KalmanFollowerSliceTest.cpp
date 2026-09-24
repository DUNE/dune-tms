// Phase 1 validation: runs TMS_KalmanFollower on a single real, already-
// known-hard slice. Unlike the first version of this tool, the seed isn't
// always TMS_GraphTrackFinder's output -- it now mirrors the full pipeline
// GraphTrackFinderTruthEfficiency.cpp validates: DBSCAN+PCA first (cheapest,
// handles most muons alone), then merge-touching-clusters-and-re-PCA, and
// only then the graph-search fallback on the merged set. Whichever stage
// finds a track-like object, that object's own points (z-sorted) become the
// Kalman follower's seed -- every track-like object gets a real physics fit,
// not just the ones that needed the graph search. By default, entry 101 of
// the 2026-09-07 clustering-benchmark reference file -- the shower-
// contaminated muon that motivated the graph-search finder in the first
// place, and a case with real, known-true momentum/charge to sanity-check
// the follower's fit against.
//
// Reuses ClusterTruthEfficiency's truth machinery (Truth_Spill loading +
// Parent-chain collapse) and GraphTrackFinderTruthEfficiency's DBSCAN+PCA+
// merge logic verbatim, so results are directly comparable to both. Needs a
// geometry file -- TMS_KalmanFollower's material stepping
// (TMS_Geom::GetMaterials) navigates a live TGeoManager, and the
// RecoCandidates-style analysis file this project has been using all week
// doesn't embed one (checked directly: no EDepSimGeometry key). Pass a
// separate geometry-bearing file (e.g. the production *_Readout.root or
// the raw *.EDEPSIM_SPILLS.root) as the first argument.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
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
#include "TMS_GraphTrackFinder.h"
#include "TMS_KalmanFollower.h"
#include "TMS_SpacePoint.h"
#include "TMS_SpacePointCluster.h"
#include "TMS_SpacePointDBScan.h"
#include "TMS_Geom.h"

namespace {

const int kMaxSpacePoints = 10000;    // matches __TMS_MAX_SPACEPOINTS__
const int kMaxTrueParticles = 20000;  // matches __TMS_MAX_TRUE_PARTICLES__

// Copied from ClusterTruthEfficiency.cpp / GraphTrackFinderSliceTest.cpp so
// "muon-owned" means the same thing across every validation tool in this
// area.
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

// Same "first-anchor" z-layer grouping TMS_GraphTrackFinder/TMS_LayerGrouping
// use, reimplemented locally just for the truth-plane-count bookkeeping
// below (not fed into the follower itself, which builds its own layers
// internally via TMS_LayerGrouping::Build()).
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

// Turns an unordered set of a track-like object's own point indices (from
// DBSCAN, a merge, or a graph-search path) into the follower's seed format:
// z-ordered (ties broken by x, then y, matching TMS_LayerGrouping's own
// convention so this stays consistent with how the follower groups the
// full pool internally). TMS_KalmanFollower::Follower::Run() only actually
// looks at the first few entries for its initial direction estimate and at
// the last entry to bound how far past the seed it walks -- passing every
// one of the object's own points (not just its endpoints) keeps that bound
// correctly reflecting the object's full known extent.
std::vector<std::size_t> BuildSeedPathFromIndices(const std::vector<TMS_SpacePoint> &points,
                                                   std::vector<int> indices) {
  std::sort(indices.begin(), indices.end(), [&points](int a, int b) {
    if (points[a].GetZ() != points[b].GetZ()) return points[a].GetZ() < points[b].GetZ();
    if (points[a].GetX() != points[b].GetX()) return points[a].GetX() < points[b].GetX();
    return points[a].GetY() < points[b].GetY();
  });
  return std::vector<std::size_t>(indices.begin(), indices.end());
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3 || argc > 6) {
    std::cerr << "Usage: " << argv[0]
              << " <geometry_source.root> <input_reco_tree.root> [target_vertexglobalid=10000461]"
                 " [target_trackid=0] [output_json_path]"
              << std::endl;
    return 1;
  }
  const std::string geom_filename = argv[1];
  const std::string input_filename = argv[2];
  const long long target_vgid = argc >= 4 ? std::stoll(argv[3]) : 10000461;
  const int target_trackid = argc >= 5 ? std::stoi(argv[4]) : 0;
  const std::string output_json_path = argc >= 6 ? argv[5] : "";
  const TrueLabel target{target_vgid, target_trackid};

  TFile geom_input(geom_filename.c_str());
  if (geom_input.IsZombie()) {
    std::cerr << "Failed to open geometry source file: " << geom_filename << std::endl;
    return 1;
  }
  TGeoManager *geom = (TGeoManager *)geom_input.Get("EDepSimGeometry");
  if (!geom) {
    std::cerr << "Geometry source file is missing 'EDepSimGeometry': " << geom_filename << std::endl;
    return 1;
  }
  TMS_Geom::GetInstance().SetGeometry(geom);

  // DBSCAN+PCA params: identical to ClusterTruthEfficiency.cpp /
  // GraphTrackFinderTruthEfficiency.cpp, so "track-like" means the same
  // thing here as in every other validation tool in this area.
  const int base_transverse_bars = 1;
  const int transverse_bars_per_plane_gap = 1;
  const int max_plane_gap = 3;
  const unsigned int min_points = 5;
  const double kLinearityThreshold = 0.8;
  const std::size_t kMinClusterSizeForTrack = 5;
  const double max_plane_pitch = TMS_Geom::GetInstance().GetMaxPlanePitch();
  const double bar_pitch = TMS_Geom::GetInstance().GetMaxBarPitch();
  if (max_plane_pitch <= 0 || bar_pitch <= 0) {
    std::cerr << "TMS_Geom found fewer than 2 surveyed planes or bars -- cannot derive a clustering tolerance."
              << std::endl;
    return 1;
  }
  const double worst_case_transverse =
      (base_transverse_bars + max_plane_gap * transverse_bars_per_plane_gap) * bar_pitch;
  const double broad_phase_radius =
      std::sqrt(worst_case_transverse * worst_case_transverse +
                std::pow(max_plane_pitch * (max_plane_gap + 1), 2));

  TFile input(input_filename.c_str());
  if (input.IsZombie()) {
    std::cerr << "Failed to open input file: " << input_filename << std::endl;
    return 1;
  }
  TTree *reco_tree = (TTree *)input.Get("Reco_Tree");
  TTree *truth_spill = (TTree *)input.Get("Truth_Spill");
  if (!reco_tree || !truth_spill) {
    std::cerr << "Input file is missing Reco_Tree/Truth_Spill" << std::endl;
    return 1;
  }

  // --- Pass 1: load Truth_Spill entirely into memory, keyed by SpillNo. ---
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
  std::cout << "Loaded Truth_Spill: " << spills.size() << " spills" << std::endl;

  // --- Pass 2: scan every slice, remembering whichever one has the most
  // points labeled as our target particle. ---
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

  Long64_t n_entries = reco_tree->GetEntries();
  Long64_t best_entry = -1;
  int best_count = 0;
  int best_total = 0;
  int best_spill = 0, best_slice = 0;
  std::vector<TMS_SpacePoint> best_points;
  std::vector<TrueLabel> best_point_label;
  // Cached separately (not just the merged single-sided label) so a
  // both-sides-verified truth trajectory -- a space point counts only if
  // BOTH its X-hit and Y-hit truth branches independently confirm the
  // target, the same rigor used for the flagship display's true-trajectory
  // overlay -- can be reconstructed after the scan for whichever slice
  // wins, without re-reading the tree.
  std::vector<TrueLabel> best_point_label_x, best_point_label_y;

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
    std::vector<TrueLabel> point_label_x(n_space_points), point_label_y(n_space_points);
    int count = 0;
    for (int i = 0; i < n_space_points; ++i) {
      TrueLabel x_label = collapse(TrueLabel{sp_x_vgid[i], sp_x_trackid[i]});
      TrueLabel y_label = collapse(TrueLabel{sp_y_vgid[i], sp_y_trackid[i]});
      point_label_x[i] = x_label;
      point_label_y[i] = y_label;
      point_label[i] = x_label.Valid() ? x_label : y_label;
      if (point_label[i] == target) ++count;
    }

    if (count > best_count) {
      best_count = count;
      best_total = n_space_points;
      best_entry = entry;
      best_spill = spill_no;
      best_slice = slice_no;
      best_points.clear();
      for (int i = 0; i < n_space_points; ++i) {
        best_points.push_back(TMS_SpacePoint(sp_x[i], sp_y[i], sp_z[i],
                                              sp_x_hitidx[i], sp_y_hitidx[i], sp_time[i]));
      }
      best_point_label = point_label;
      best_point_label_x = point_label_x;
      best_point_label_y = point_label_y;
    }
  }

  if (best_entry < 0) {
    std::cerr << "No slice found containing vertexglobalid=" << target_vgid
              << " trackid=" << target_trackid << std::endl;
    return 1;
  }

  std::cout << "Target slice: entry=" << best_entry << " spill=" << best_spill
            << " slice=" << best_slice << ", " << best_total
            << " total space points, " << best_count
            << " labeled as the target particle" << std::endl;

  const std::vector<int> z_layer = AssignZLayers(best_points, 1.0);
  std::set<int> target_layers_in_slice;
  for (std::size_t i = 0; i < best_points.size(); ++i) {
    if (best_point_label[i] == target) target_layers_in_slice.insert(z_layer[i]);
  }
  std::cout << "Target particle touches " << target_layers_in_slice.size()
            << " distinct planes in this slice.\n";

  // --- Seed: the same three-stage pipeline GraphTrackFinderTruthEfficiency
  // validates (DBSCAN+PCA -> merge-and-re-PCA -> graph-search on the merged
  // set), except now whichever stage actually finds a track-like object
  // hands ITS OWN points to the Kalman follower as the seed -- previously
  // this tool always ran the graph search on the whole slice regardless,
  // meaning only the hardest (graph-search-needed) cases ever got a real
  // physics fit. Now every track-like object does. ---
  std::vector<int> plane_index;
  plane_index.reserve(best_points.size());
  for (const TMS_SpacePoint &point : best_points)
    plane_index.push_back(TMS_Geom::GetInstance().GetPlaneIndexNearestZ(point.GetZ()));

  TMS_SpacePointDBScan dbscan(best_points, plane_index, min_points, bar_pitch, base_transverse_bars,
                               transverse_bars_per_plane_gap, max_plane_gap, broad_phase_radius);
  std::vector<std::vector<int>> cluster_indices = dbscan.RunAndGetClusterIndices();
  std::vector<TMS_SpacePointCluster> clusters;
  clusters.reserve(cluster_indices.size());
  for (auto &indices : cluster_indices) clusters.emplace_back(best_points, indices);

  std::vector<int> point_cluster_id(best_points.size(), 0);  // 0 = noise, else cluster index + 1
  for (std::size_t c = 0; c < cluster_indices.size(); ++c)
    for (int idx : cluster_indices[c]) point_cluster_id[idx] = static_cast<int>(c) + 1;

  auto ClusterOwner = [&](const std::vector<int> &indices) -> TrueLabel {
    std::unordered_map<TrueLabel, int, LabelHash> votes;
    for (int idx : indices)
      if (best_point_label[idx].Valid()) votes[best_point_label[idx]]++;
    TrueLabel owner;
    int owner_count = 0;
    for (auto &kv : votes)
      if (kv.second > owner_count) {
        owner = kv.first;
        owner_count = kv.second;
      }
    return owner;
  };

  std::vector<std::size_t> seedPath;
  std::string foundVia;
  // The object's own point indices in their RAW (unordered) form, as found
  // by Stage 1/2 -- kept separately from seedPath (which BuildSeedPathFromIndices
  // has already z-sorted for display) so the follower can try every
  // candidate at the object's own first z-layer as the seed anchor
  // (RunBestSeed(), see TMS_KalmanFollower.h) instead of committing to
  // whichever one a naive z-sort happens to put first. Left empty for the
  // Stage 3 (graph-search) case, whose own directed search already resolved
  // this same first-point ambiguity -- that seed keeps using plain Run().
  std::vector<std::size_t> seedObjectIndices;
  // Every candidate path the graph search produced (not just the best),
  // converted to global best_points indices -- kept around only for the
  // optional JSON display dump (Stage 3, below, fills this in when it
  // actually runs).
  std::vector<std::vector<std::size_t>> allGraphtrackGlobalPaths;

  // Stage 1: does ANY track-like cluster have the target as its own
  // plurality owner? Must match KalmanFollowerTruthEfficiency.cpp's Pass A
  // exactly (iterate every cluster, filter on IsTrackLike + owner==target,
  // pick the largest qualifying one) -- NOT pre-select a single candidate
  // cluster by "which cluster holds the most of the target's own points"
  // first. Those two queries can disagree: when the target's true points
  // are split across multiple clusters, the cluster holding the MOST of
  // them can be a different (non-track-like, or owned-by-someone-else)
  // cluster than a smaller one that actually passes both checks -- found
  // 2026-09-16 via a real discrepancy against the batch tool's numbers for
  // the same (file, vgid, trackid): this tool fell through to the
  // graph-search fallback where the batch tool's Pass A succeeded cleanly.
  std::vector<std::size_t> stage1Candidates;
  for (std::size_t c = 0; c < clusters.size(); ++c) {
    if (!clusters[c].IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack)) continue;
    if (ClusterOwner(cluster_indices[c]) == target) stage1Candidates.push_back(c);
  }

  // Debug: every cluster touching ANY of the target's own points, whether
  // or not it qualifies for Stage 1 -- added to directly verify (not just
  // reason about) why a muon's DBSCAN-direct cluster can be a small,
  // isolated tail fragment even when real truth-matched hits exist earlier:
  // a busy co-vertex cluster can contain plenty of the muon's own points
  // yet still be plurality-owned by a companion particle with even more
  // points in it, excluding the whole cluster from Stage 1.
  if (std::getenv("KF_DEBUG")) {
    std::set<int> touchingClusterIds;
    for (std::size_t i = 0; i < best_points.size(); ++i)
      if (best_point_label[i] == target && point_cluster_id[i] > 0) touchingClusterIds.insert(point_cluster_id[i]);
    std::cerr << "[KF_DEBUG Stage1] clusters touching target's own points:\n";
    for (int cid : touchingClusterIds) {
      const std::vector<int> &idx = cluster_indices[cid - 1];
      double zmin = 1e18, zmax = -1e18;
      int nOwn = 0;
      std::unordered_map<TrueLabel, int, LabelHash> votes;
      for (int i : idx) {
        zmin = std::min(zmin, best_points[i].GetZ());
        zmax = std::max(zmax, best_points[i].GetZ());
        if (best_point_label[i] == target) ++nOwn;
        if (best_point_label[i].Valid()) votes[best_point_label[i]]++;
      }
      TrueLabel owner;
      int ownerCount = 0;
      for (auto &kv : votes)
        if (kv.second > ownerCount) { owner = kv.first; ownerCount = kv.second; }
      std::cerr << "  cluster " << cid << ": size=" << idx.size() << " z=[" << zmin << "," << zmax << "]"
                << " nOwnPoints=" << nOwn << " trackLike=" << clusters[cid - 1].IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack)
                << " owner=(vgid=" << owner.vgid << ",tid=" << owner.trackid << ",n=" << ownerCount << ")"
                << " isTarget=" << (owner == target) << "\n";
    }
  }
  // How many of the target's own truth-matched planes does a candidate
  // index set actually touch? Used below to compare Stage 1/2's whole-
  // cluster candidate against Stage 3's GraphTrackFinder path on equal
  // footing (raw-candidate coverage, before the Kalman fit).
  auto CountOwnPlanes = [&](const std::vector<int> &idxs) {
    std::set<int> planes;
    for (int i : idxs)
      if (best_point_label[i] == target) planes.insert(static_cast<int>(std::round(best_points[i].GetZ())));
    return (int)planes.size();
  };

  std::vector<int> stage1Indices;
  if (!stage1Candidates.empty()) {
    std::size_t bestC = stage1Candidates[0];
    for (std::size_t c : stage1Candidates)
      if (cluster_indices[c].size() > cluster_indices[bestC].size()) bestC = c;
    stage1Indices.assign(cluster_indices[bestC].begin(), cluster_indices[bestC].end());
  }
  const int stage1Planes = stage1Indices.empty() ? 0 : CountOwnPlanes(stage1Indices);

  // Stage 2: merge every cluster touching the target's own points, plus its
  // own noise points specifically, and re-check PCA on the merged set. Also
  // the candidate pool Stage 3 (below) searches. Computed UNCONDITIONALLY
  // now (not gated on Stage 1 failing) -- see the fix note below.
  std::set<int> touched_cluster_ids;
  std::vector<int> own_noise_points;
  for (std::size_t i = 0; i < best_points.size(); ++i) {
    if (!(best_point_label[i] == target)) continue;
    const int cid = point_cluster_id[i];
    if (cid == 0) own_noise_points.push_back(static_cast<int>(i));
    else touched_cluster_ids.insert(cid);
  }
  std::vector<int> merged_indices = own_noise_points;
  for (int cid : touched_cluster_ids)
    for (int idx : cluster_indices[cid - 1]) merged_indices.push_back(idx);
  std::sort(merged_indices.begin(), merged_indices.end());
  merged_indices.erase(std::unique(merged_indices.begin(), merged_indices.end()), merged_indices.end());

  std::vector<int> stage2Indices;
  if (merged_indices.size() >= kMinClusterSizeForTrack) {
    TMS_SpacePointCluster merged_cluster(best_points, merged_indices);
    if (merged_cluster.IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack) &&
        ClusterOwner(merged_indices) == target) {
      stage2Indices = merged_indices;
    }
  }
  const int stage2Planes = stage2Indices.empty() ? 0 : CountOwnPlanes(stage2Indices);

  // Stage 3 / fix (2026-09-16): run GraphTrackFinder on the SAME merged
  // touching-cluster set REGARDLESS of whether Stage 1/2 already succeeded
  // -- previously gated by `if (seedPath.empty())`, meaning ANY Stage 1
  // success, even a trivial isolated fragment, skipped this entirely and
  // GraphTrackFinder never got the chance to decompose a bigger, real,
  // companion-dominated cluster into a track-like sub-path (exactly the job
  // it exists for). Validated via KalmanFollowerTruthEfficiency.cpp at full
  // population scale: completeness 68.5%->79.4% (ND-LAr-fiducial:
  // 70.5%->84.6%), purity unchanged (88.1%->88.3%). Matches this tool's own
  // pre-existing convention for Stage 3 (already picks among GraphTrackFinder's
  // candidate paths by truth-matched count, not a new methodological
  // departure -- see KF_DEBUG/probe history above for how this was found).
  TMS_GraphTrackFinder::Config finderConfig;
  finderConfig.OccupancyPenalty = 0.0;
  finderConfig.HitMultiplicityPenalty = 0.0;
  finderConfig.MaxSeedLayerOccupancy = 150;
  finderConfig.MaxSeedHitMultiplicity = 50;
  finderConfig.UseCurvatureProjection = false;
  std::vector<int> stage3Indices;
  int stage3Planes = 0;
  if (merged_indices.size() >= finderConfig.SeedLength) {
    std::vector<TMS_SpacePoint> local_points;
    local_points.reserve(merged_indices.size());
    for (int gi : merged_indices) local_points.push_back(best_points[gi]);
    const TMS_GraphTrackFinder::Result finderResult =
        TMS_GraphTrackFinder::Finder(finderConfig).Find(local_points);

    for (const TMS_GraphTrackFinder::Path &path : finderResult.Paths) {
      std::vector<std::size_t> globalPath;
      globalPath.reserve(path.SpacePointIndices.size());
      for (std::size_t localIdx : path.SpacePointIndices)
        globalPath.push_back(static_cast<std::size_t>(merged_indices[localIdx]));
      allGraphtrackGlobalPaths.push_back(std::move(globalPath));
    }

    const TMS_GraphTrackFinder::Path *bestSeed = nullptr;
    for (const TMS_GraphTrackFinder::Path &path : finderResult.Paths) {
      std::vector<int> globalIndices;
      globalIndices.reserve(path.SpacePointIndices.size());
      for (std::size_t localIdx : path.SpacePointIndices) globalIndices.push_back(merged_indices[localIdx]);
      const int planes = CountOwnPlanes(globalIndices);
      if (planes > stage3Planes) {
        stage3Planes = planes;
        bestSeed = &path;
      }
    }
    if (bestSeed && bestSeed->SpacePointIndices.size() >= 2) {
      stage3Indices.reserve(bestSeed->SpacePointIndices.size());
      for (std::size_t localIdx : bestSeed->SpacePointIndices) stage3Indices.push_back(merged_indices[localIdx]);
    }
  }

  // Keep whichever of the three candidates covers the most of the target's
  // own truth-matched planes. Strict `>` (not `>=`) so a tie prefers the
  // simpler/earlier stage -- matches KalmanFollowerTruthEfficiency.cpp's
  // `probe.best_planes_covered > stageCandidatePlanes` exactly; using `>=`
  // here first biased every tie toward the more complex GraphTrackFinder
  // path even when it wasn't actually better, which is what happened on a
  // real re-check of case E (a tie in raw-candidate coverage, but the
  // GraphTrackFinder-seeded fit did WORSE post-fit than Stage 1's clean
  // small candidate would have) -- caught by testing this exact case, not
  // assumed.
  if (stage3Planes > stage1Planes && stage3Planes > stage2Planes && !stage3Indices.empty()) {
    seedPath = BuildSeedPathFromIndices(best_points, stage3Indices);
    // Also route through RunBestSeed() below, matching KalmanFollowerTruthEfficiency.cpp
    // -- but note this is a NO-OP for a path exactly this shape (one point
    // per z-layer already, by construction of an already-resolved
    // GraphTrackFinder path): RunBestSeed's multi-hypothesis mechanism only
    // has alternatives to try when the object itself contains >1 candidate
    // at its own first layer, which a single resolved path never does.
    // Verified empirically on case E (2026-09-16): identical result either
    // way. Kept for consistency with the batch tool rather than removed,
    // since a future stage3Indices shape (or a real multi-candidate first
    // layer) could still benefit.
    seedObjectIndices.assign(stage3Indices.begin(), stage3Indices.end());
    foundVia = (stage1Planes == 0 && stage2Planes == 0) ? "GraphTrackFinder (merged fallback)"
                                                          : "GraphTrackFinder (beat Stage 1/2)";
  } else if (stage2Planes > stage1Planes && !stage2Indices.empty()) {
    seedPath = BuildSeedPathFromIndices(best_points, stage2Indices);
    seedObjectIndices.assign(stage2Indices.begin(), stage2Indices.end());
    foundVia = "merged-cluster PCA";
  } else if (!stage1Indices.empty()) {
    seedPath = BuildSeedPathFromIndices(best_points, stage1Indices);
    seedObjectIndices.assign(stage1Indices.begin(), stage1Indices.end());
    foundVia = "DBSCAN+PCA (direct)";
  }

  if (seedPath.empty()) {
    std::cerr << "No track-like object found for this target by any stage -- nothing to follow." << std::endl;
    return 1;
  }

  std::size_t seedMatched = 0;
  for (std::size_t idx : seedPath)
    if (best_point_label[idx] == target) ++seedMatched;
  std::cout << "\nFound via: " << foundVia << " -- " << seedPath.size() << " points, " << seedMatched
            << " target-matched.\n";

  // --- Follow it. Always against the FULL slice's space points, regardless
  // of which stage found the seed -- ambiguity resolution needs to see
  // ghosts the finding stage didn't pick, same as the graph-search-only
  // design this replaces. ---
  const RegionFieldModel field;  // 1.0T, GDML-confirmed -- see TMS_FieldModel.h
  // Same sweep hooks as KalmanFollowerTruthEfficiency (unset = Config defaults).
  TMS_KalmanFollower::Config followerConfig;
  if (const char *v = std::getenv("KF_QP_REL_SIGMA")) followerConfig.InitialQPRelSigma = std::atof(v);
  if (const char *v = std::getenv("KF_RANGE_SEED")) followerConfig.RangeSeedMargin = std::atof(v);
  if (const char *v = std::getenv("KF_MAX_HEAD_SKIP")) followerConfig.MaxHeadSkip = std::atoi(v);
  if (const char *v = std::getenv("KF_MAX_TRIPLETS")) followerConfig.MaxTripletHypotheses = std::atoi(v);
  if (const char *v = std::getenv("KF_RANK_BY_CONVERGENCE")) followerConfig.RankHypothesesByConvergence = std::atoi(v) != 0;
  if (const char *v = std::getenv("KF_STOP_ON_RANGEOUT")) followerConfig.StopOnRangeOut = std::atoi(v) != 0;
  const TMS_KalmanFollower::Follower follower(followerConfig, field);

  // DBSCAN-direct/merged-PCA seeds (Stages 1-2) are unordered blobs with no
  // directed search behind them -- naively z-sorting and taking whichever
  // point lands first can anchor the fit on a bad choice when the object's
  // own first z-layer has more than one point at (near-)identical z.
  // RunBestSeed() spawns one hypothesis per first-layer candidate and keeps
  // the best. Stage 3 (graph-search fallback) already ran a directed search
  // that resolved this exact ambiguity, so it keeps using Run() on its own
  // already-ordered path.
  const bool usedMultiHypothesis = !seedObjectIndices.empty();
  const TMS_KalmanFollower::FitResult fit = usedMultiHypothesis
      ? follower.RunBestSeed(best_points, seedObjectIndices)
      : follower.Run(best_points, seedPath);
  if (usedMultiHypothesis) {
    std::cout << "\nMulti-hypothesis seeding: " << seedObjectIndices.size()
              << " object points, one fit per candidate at the object's own first z-layer.\n";
  }

  std::cout << "\nKalman follower result\n"
            << "  converged: " << (fit.Converged ? "yes" : "NO") << '\n'
            << "  nodes walked: " << fit.Nodes.size() << '\n'
            << "  gaps filled: " << fit.NGapsFilled << '\n'
            << "  ambiguous layers (>1 candidate): " << fit.NAmbiguousLayersResolved << '\n'
            << "  total chi2 / NDoF: " << fit.TotalChi2 << " / " << fit.NDoF << '\n'
            << "  fitted momentum: " << fit.MomentumMeV << " MeV\n"
            << "  fitted charge sign: " << (fit.Charge > 0 ? "+1" : "-1") << '\n';

  // Truth cross-check + ambiguity-resolution accuracy, per-node.
  int nCorrectChosen = 0, nWrongChosen = 0, nGaps = 0;
  int nTruthAvailableAtAmbiguousLayer = 0, nCorrectAtAmbiguousLayer = 0;
  int nCoveredPlanes = 0;
  std::set<int> coveredTargetLayers;
  for (const TMS_KalmanFollower::FollowedNode &node : fit.Nodes) {
    if (!node.HasHit) {
      ++nGaps;
      continue;
    }
    const bool chosenIsTarget = best_point_label[node.ChosenSpacePointIndex] == target;
    if (chosenIsTarget) {
      ++nCorrectChosen;
      coveredTargetLayers.insert(z_layer[node.ChosenSpacePointIndex]);
    } else {
      ++nWrongChosen;
    }
    if (node.CandidateIndices.size() > 1) {
      bool truthPresent = false;
      for (std::size_t idx : node.CandidateIndices)
        if (best_point_label[idx] == target) truthPresent = true;
      if (truthPresent) {
        ++nTruthAvailableAtAmbiguousLayer;
        if (chosenIsTarget) ++nCorrectAtAmbiguousLayer;
      }
    }
  }
  nCoveredPlanes = static_cast<int>(coveredTargetLayers.size());

  std::cout << "\nTruth cross-check\n"
            << "  nodes with truth-correct chosen point: " << nCorrectChosen << '\n'
            << "  nodes with wrong chosen point: " << nWrongChosen << '\n'
            << "  gap nodes: " << nGaps << '\n'
            << "  target planes covered: " << nCoveredPlanes << "/"
            << target_layers_in_slice.size() << '\n'
            << "  ambiguous layers where truth was even a candidate: "
            << nTruthAvailableAtAmbiguousLayer << " (correctly chosen: "
            << nCorrectAtAmbiguousLayer << ")\n";

  std::cout << "\nPer-node detail:\n";
  for (const TMS_KalmanFollower::FollowedNode &node : fit.Nodes) {
    std::cout << "  layer " << node.Layer << " z=" << node.Z
              << " candidates=" << node.CandidateIndices.size();
    if (node.HasHit) {
      const bool correct = best_point_label[node.ChosenSpacePointIndex] == target;
      std::cout << " chosen=" << node.ChosenSpacePointIndex << " (" << (correct ? "correct" : "WRONG")
                << ") chi2=" << node.Chi2AtChosen << " p=" << (std::abs(node.FilteredQP) > 1e-12
                                                                    ? 1.0 / std::abs(node.FilteredQP)
                                                                    : 0.0)
                << " MeV";
    } else {
      std::cout << " GAP";
    }
    std::cout << '\n';
  }

  // Optional: dump the point cloud, the graph-search candidate paths (if
  // Stage 3 ran), the both-sides-verified true trajectory, and the Kalman
  // fit's own trajectory as JSON, for an event display. own/blob/other
  // matches the 3-way split used throughout this project's displays (own =
  // the target particle; blob = a different track at the *same* vertex;
  // other = a different vertex entirely, or no truth).
  if (!output_json_path.empty()) {
    std::ofstream json(output_json_path);
    json << std::fixed;
    json << "{\"vgid\":" << target.vgid << ",\"trackid\":" << target.trackid
         << ",\"entry\":" << best_entry << ",\"spill\":" << best_spill
         << ",\"slice\":" << best_slice << ",\"n_total_slice\":" << best_total
         << ",\"found_via\":\"" << foundVia << "\""
         << ",\"n_target_planes\":" << target_layers_in_slice.size();

    auto writePoints = [&](const char *key, std::function<bool(std::size_t)> include) {
      json << ",\"" << key << "\":[";
      bool first = true;
      for (std::size_t i = 0; i < best_points.size(); ++i) {
        if (!include(i)) continue;
        if (!first) json << ",";
        first = false;
        json << "[" << best_points[i].GetX() << "," << best_points[i].GetY() << ","
             << best_points[i].GetZ() << "]";
      }
      json << "]";
    };
    writePoints("own_points", [&](std::size_t i) { return best_point_label[i] == target; });
    writePoints("blob_points", [&](std::size_t i) {
      return best_point_label[i].vgid == target.vgid && !(best_point_label[i] == target);
    });
    writePoints("other_points", [&](std::size_t i) { return best_point_label[i].vgid != target.vgid; });
    writePoints("true_trajectory", [&](std::size_t i) {
      return best_point_label_x[i] == target && best_point_label_y[i] == target;
    });

    // The actual object that was found and handed to the follower as a seed
    // (the specific DBSCAN cluster, the merged set, or -- for the graph-
    // search fallback -- its own best path) -- distinct from own_points
    // above, which is every point in the WHOLE SLICE single-sided-truth-
    // matched to the target, including any ghost/ambiguous matches the
    // found object never actually contained. Marking each point's own
    // truth-correctness here (rather than assuming every seed point is
    // right) keeps this honest for the DBSCAN-direct/merged-PCA case too.
    json << ",\"seed_points\":[";
    for (std::size_t k = 0; k < seedPath.size(); ++k) {
      if (k) json << ",";
      const std::size_t idx = seedPath[k];
      json << "{\"x\":" << best_points[idx].GetX() << ",\"y\":" << best_points[idx].GetY()
           << ",\"z\":" << best_points[idx].GetZ() << ",\"is_target\":"
           << (best_point_label[idx] == target ? "true" : "false") << "}";
    }
    json << "]";

    json << ",\"paths\":[";
    for (std::size_t p = 0; p < allGraphtrackGlobalPaths.size(); ++p) {
      if (p) json << ",";
      const std::vector<std::size_t> &path = allGraphtrackGlobalPaths[p];
      std::size_t matched = 0;
      for (std::size_t idx : path)
        if (best_point_label[idx] == target) ++matched;
      const double purity = path.empty() ? 0.0 : 100.0 * matched / path.size();
      json << "{\"purity\":" << purity << ",\"n_matched\":" << matched << ",\"points\":[";
      for (std::size_t k = 0; k < path.size(); ++k) {
        if (k) json << ",";
        const std::size_t idx = path[k];
        json << "{\"x\":" << best_points[idx].GetX() << ",\"y\":" << best_points[idx].GetY()
             << ",\"z\":" << best_points[idx].GetZ() << ",\"is_target\":"
             << (best_point_label[idx] == target ? "true" : "false") << "}";
      }
      json << "]}";
    }
    json << "]";

    json << ",\"kalman_fit\":{\"converged\":" << (fit.Converged ? "true" : "false")
         << ",\"momentum_mev\":" << fit.MomentumMeV << ",\"charge\":" << fit.Charge
         << ",\"nodes\":[";
    for (std::size_t n = 0; n < fit.Nodes.size(); ++n) {
      if (n) json << ",";
      const TMS_KalmanFollower::FollowedNode &node = fit.Nodes[n];
      const bool chosenIsTarget = node.HasHit && best_point_label[node.ChosenSpacePointIndex] == target;
      const double p = std::abs(node.FilteredQP) > 1e-12 ? 1.0 / std::abs(node.FilteredQP) : 0.0;
      json << "{\"z\":" << node.Z << ",\"x\":" << node.FilteredX << ",\"y\":" << node.FilteredY
           << ",\"has_hit\":" << (node.HasHit ? "true" : "false")
           << ",\"chosen_is_target\":" << (chosenIsTarget ? "true" : "false")
           << ",\"momentum_mev\":" << p << "}";
    }
    json << "]}";

    json << "}";
    json.close();
    std::cout << "\nWrote event-display JSON to " << output_json_path << std::endl;
  }

  return 0;
}
