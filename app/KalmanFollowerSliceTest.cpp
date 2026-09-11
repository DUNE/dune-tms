// Phase 1 validation: runs TMS_KalmanFollower on a single real, already-
// known-hard slice, seeded from TMS_GraphTrackFinder's own output (the same
// machinery GraphTrackFinderSliceTest.cpp validates). By default, entry 101
// of the 2026-09-07 clustering-benchmark reference file -- the shower-
// contaminated muon that motivated the graph-search finder in the first
// place, and a case with real, known-true momentum/charge to sanity-check
// the follower's fit against (see the README-style comment near main()).
//
// Reuses ClusterTruthEfficiency's truth machinery (Truth_Spill loading +
// Parent-chain collapse) and GraphTrackFinderSliceTest's slice-loading
// pattern, so results are directly comparable to both. Unlike those tools,
// this one DOES need a geometry file -- TMS_KalmanFollower's material
// stepping (TMS_Geom::GetMaterials) navigates a live TGeoManager, and the
// RecoCandidates-style analysis file this project has been using all week
// doesn't embed one (checked directly: no EDepSimGeometry key). Pass a
// separate geometry-bearing file (e.g. the production *_Readout.root or
// the raw *.EDEPSIM_SPILLS.root) as the first argument, same convention
// GraphTrackFinderTruthEfficiency.cpp already uses.

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

  // --- Seed: run TMS_GraphTrackFinder on the slice to get a topological
  // path, the same way GraphTrackFinderTruthEfficiency's fallback does.
  // This tool cares about the follower, not re-tuning the finder, so it
  // uses that same already-validated real-data config (relaxed seed gates,
  // zeroed occupancy/multiplicity penalty, quantization deadband on) rather
  // than TMS_GraphTrackFinder::Config's stock defaults, which don't find
  // any usable seed at all on this specific dense case. ---
  TMS_GraphTrackFinder::Config finderConfig;
  finderConfig.OccupancyPenalty = 0.0;
  finderConfig.HitMultiplicityPenalty = 0.0;
  finderConfig.MaxSeedLayerOccupancy = 150;
  finderConfig.MaxSeedHitMultiplicity = 50;
  finderConfig.UseCurvatureProjection = false;
  const TMS_GraphTrackFinder::Result finderResult =
      TMS_GraphTrackFinder::Finder(finderConfig).Find(best_points);

  std::size_t bestSeedMatched = 0;
  const TMS_GraphTrackFinder::Path *bestSeed = nullptr;
  for (const TMS_GraphTrackFinder::Path &path : finderResult.Paths) {
    std::size_t matched = 0;
    for (std::size_t idx : path.SpacePointIndices)
      if (best_point_label[idx] == target) ++matched;
    if (matched > bestSeedMatched) {
      bestSeedMatched = matched;
      bestSeed = &path;
    }
  }
  if (!bestSeed || bestSeed->SpacePointIndices.size() < 2) {
    std::cerr << "Graph Track Finder found no usable seed path for this target -- nothing to follow."
              << std::endl;
    return 1;
  }
  std::cout << "\nSeed path (Graph Track Finder, default config): "
            << bestSeed->SpacePointIndices.size() << " points, "
            << bestSeedMatched << " target-matched.\n";

  // --- Follow it. ---
  const RegionFieldModel field;  // placeholder Tesla magnitude -- see TMS_FieldModel.h
  const TMS_KalmanFollower::Config followerConfig;
  const TMS_KalmanFollower::Follower follower(followerConfig, field);
  const TMS_KalmanFollower::FitResult fit = follower.Run(best_points, bestSeed->SpacePointIndices);

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
