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
  if (argc < 3 || argc > 5) {
    std::cerr << "Usage: " << argv[0]
              << " <geometry_source.root> <input_reco_tree.root> [target_vertexglobalid=10000461]"
                 " [target_trackid=0]"
              << std::endl;
    return 1;
  }
  const std::string geom_filename = argv[1];
  const std::string input_filename = argv[2];
  const long long target_vgid = argc >= 4 ? std::stoll(argv[3]) : 10000461;
  const int target_trackid = argc >= 5 ? std::stoi(argv[4]) : 0;
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
      best_entry = entry;
      best_spill = spill_no;
      best_slice = slice_no;
      best_points.clear();
      for (int i = 0; i < n_space_points; ++i) {
        best_points.push_back(TMS_SpacePoint(sp_x[i], sp_y[i], sp_z[i],
                                              sp_x_hitidx[i], sp_y_hitidx[i], sp_time[i]));
      }
      best_point_label = point_label;
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

  // Stage 1: does the cluster that plurality-owns the target's own points
  // (if any) already pass the PCA linearity check on its own?
  std::unordered_map<int, int> ownClusterVotes;
  for (std::size_t i = 0; i < best_points.size(); ++i)
    if (best_point_label[i] == target) ownClusterVotes[point_cluster_id[i]]++;
  int bestOwnClusterId = 0, bestOwnClusterVotes = 0;
  for (auto &kv : ownClusterVotes)
    if (kv.second > bestOwnClusterVotes) {
      bestOwnClusterId = kv.first;
      bestOwnClusterVotes = kv.second;
    }
  if (bestOwnClusterId > 0) {
    const TMS_SpacePointCluster &cl = clusters[bestOwnClusterId - 1];
    if (cl.IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack) &&
        ClusterOwner(cluster_indices[bestOwnClusterId - 1]) == target) {
      seedPath = BuildSeedPathFromIndices(best_points, cluster_indices[bestOwnClusterId - 1]);
      foundVia = "DBSCAN+PCA (direct)";
    }
  }

  // Stage 2: merge every cluster touching the target's own points, plus its
  // own noise points specifically, and re-check PCA on the merged set.
  std::vector<int> merged_indices;
  if (seedPath.empty()) {
    std::set<int> touched_cluster_ids;
    std::vector<int> own_noise_points;
    for (std::size_t i = 0; i < best_points.size(); ++i) {
      if (!(best_point_label[i] == target)) continue;
      const int cid = point_cluster_id[i];
      if (cid == 0) own_noise_points.push_back(static_cast<int>(i));
      else touched_cluster_ids.insert(cid);
    }
    merged_indices = own_noise_points;
    for (int cid : touched_cluster_ids)
      for (int idx : cluster_indices[cid - 1]) merged_indices.push_back(idx);
    std::sort(merged_indices.begin(), merged_indices.end());
    merged_indices.erase(std::unique(merged_indices.begin(), merged_indices.end()), merged_indices.end());

    if (merged_indices.size() >= kMinClusterSizeForTrack) {
      TMS_SpacePointCluster merged_cluster(best_points, merged_indices);
      if (merged_cluster.IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack) &&
          ClusterOwner(merged_indices) == target) {
        seedPath = BuildSeedPathFromIndices(best_points, merged_indices);
        foundVia = "merged-cluster PCA";
      }
    }
  }

  // Stage 3: still not track-like -- run the graph search on the merged
  // set (never the whole slice, matching GraphTrackFinderTruthEfficiency's
  // validated design), using the config already validated as needed for
  // real dense-slice occupancy (TMS_GraphTrackFinder::Config's stock
  // defaults find no usable seed at all on cases like this one).
  TMS_GraphTrackFinder::Config finderConfig;
  finderConfig.OccupancyPenalty = 0.0;
  finderConfig.HitMultiplicityPenalty = 0.0;
  finderConfig.MaxSeedLayerOccupancy = 150;
  finderConfig.MaxSeedHitMultiplicity = 50;
  finderConfig.UseCurvatureProjection = false;
  if (seedPath.empty() && merged_indices.size() >= finderConfig.SeedLength) {
    std::vector<TMS_SpacePoint> local_points;
    local_points.reserve(merged_indices.size());
    for (int gi : merged_indices) local_points.push_back(best_points[gi]);
    const TMS_GraphTrackFinder::Result finderResult =
        TMS_GraphTrackFinder::Finder(finderConfig).Find(local_points);

    std::size_t bestSeedMatched = 0;
    const TMS_GraphTrackFinder::Path *bestSeed = nullptr;
    for (const TMS_GraphTrackFinder::Path &path : finderResult.Paths) {
      std::size_t matched = 0;
      for (std::size_t localIdx : path.SpacePointIndices)
        if (best_point_label[merged_indices[localIdx]] == target) ++matched;
      if (matched > bestSeedMatched) {
        bestSeedMatched = matched;
        bestSeed = &path;
      }
    }
    if (bestSeed && bestSeed->SpacePointIndices.size() >= 2) {
      std::vector<int> globalIndices;
      globalIndices.reserve(bestSeed->SpacePointIndices.size());
      for (std::size_t localIdx : bestSeed->SpacePointIndices) globalIndices.push_back(merged_indices[localIdx]);
      seedPath = BuildSeedPathFromIndices(best_points, globalIndices);
      foundVia = "GraphTrackFinder (merged fallback)";
    }
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
  const TMS_KalmanFollower::Config followerConfig;
  const TMS_KalmanFollower::Follower follower(followerConfig, field);
  const TMS_KalmanFollower::FitResult fit = follower.Run(best_points, seedPath);

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

  return 0;
}
