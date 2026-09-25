// Full-population truth-based validation of a three-stage DBSCAN+PCA ->
// merged-cluster PCA -> Graph Track Finder pipeline: run the existing, cheap
// DBSCAN+PCA clustering (identical to ClusterTruthEfficiency.cpp -- same
// tolerances, same plurality-vote matching) on every slice; for any muon
// candidate it misses, merge every DBSCAN cluster (track-like or not) that
// contains any of that muon's own points plus its own individually-
// unclustered ("noise") points, and re-run PCA on the merged set (a real
// track can legitimately get split across DBSCAN's own density boundaries,
// e.g. a sparse tail end going its own way -- merging first undoes exactly
// that, without pulling in the rest of the slice); only if the merged set
// is still not track-like does the much more expensive Graph Track Finder graph
// search run, on that same merged set. This mirrors how DBSCAN and Link-
// and-Tree are actually meant to be used together (a fallback for the
// cases DBSCAN structurally can't cluster, not a universal replacement),
// scoped to the muon's own split pieces rather than the whole slice or an
// arbitrary single cluster.

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

#include "TMS_Geom.h"
#include "TMS_GraphTrackFinder.h"
#include "TMS_SpacePoint.h"
#include "TMS_SpacePointCluster.h"
#include "TMS_SpacePointDBScan.h"
#include "TMS_LayerGrouping.h"
#include "SpacePointLayerInput.h"

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
  std::vector<float> momentum;
  std::vector<float> path_length_tms;
  std::vector<bool> tms_fiducial_start;
  std::vector<bool> lar_fiducial_start;
  std::unordered_map<TrueLabel, int, LabelHash> index_of;
  std::vector<int> collapsed_trackid;
};

// Same convention as ClusterTruthEfficiency.cpp.
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

  // DBSCAN+PCA params: identical to ClusterTruthEfficiency.cpp.
  const unsigned int min_points = 5;
  const double kLinearityThreshold = 0.8;
  const size_t kMinClusterSizeForTrack = 5;

  // Graph Track Finder fallback config: today's validated fix (seed gates at
  // real occupancy scale, occupancy/multiplicity growth-time scoring
  // penalty zeroed, quantization deadband on, curvature off). Deliberately
  // not the dense-layer search -- validated as a net regression.
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
  const double bar_pitch = TMS_Geom::GetInstance().GetMaxBarPitch();
  if (bar_pitch <= 0) {
    std::cerr << "TMS_Geom found fewer than 2 surveyed bars -- cannot derive a clustering tolerance." << std::endl;
    return -1;
  }
  // DBSCAN: TMS_SpacePointDBScan::DefaultTolerance() for this bar pitch.
  // Sweep hooks (Phase 1 baselines): DBSCAN_MAX_DZ_MM overrides the z window,
  // DBSCAN_MIN_POINTS the core-point threshold (which does NOT change the
  // muon population, still defined by min_points true hits).
  TMS_SpacePointDBScan::Tolerance dbscan_tolerance = TMS_SpacePointDBScan::DefaultTolerance(bar_pitch);
  if (const char *v = std::getenv("DBSCAN_MAX_DZ_MM")) dbscan_tolerance.MaxDzMM = std::atof(v);
  unsigned int dbscan_min_points = min_points;
  if (const char *v = std::getenv("DBSCAN_MIN_POINTS")) dbscan_min_points = std::atoi(v);

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
  static std::vector<float> mom_ts(kMaxTrueParticles * 4);
  static std::vector<float> path_ts(kMaxTrueParticles);
  static bool tms_fid_start_ts[kMaxTrueParticles];
  static bool lar_fid_start_ts[kMaxTrueParticles];
  truth_spill->SetBranchAddress("SpillNo", &spill_no_ts);
  truth_spill->SetBranchAddress("nTrueParticles", &n_tp_ts);
  truth_spill->SetBranchAddress("VertexGlobalID", vgid_ts.data());
  truth_spill->SetBranchAddress("TrackId", trackid_ts.data());
  truth_spill->SetBranchAddress("PDG", pdg_ts.data());
  truth_spill->SetBranchAddress("Parent", parent_ts.data());
  truth_spill->SetBranchAddress("MomentumTMSStart", mom_ts.data());
  truth_spill->SetBranchAddress("TruePathLengthInTMS", path_ts.data());
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
    sp.momentum.assign(mom_ts.begin(), mom_ts.begin() + n_tp_ts * 4);
    sp.path_length_tms.assign(path_ts.begin(), path_ts.begin() + n_tp_ts);
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
                 "momentum_mag,angle_deg,path_length_mm,true_hits_in_slice,n_muons_in_slice,"
                 "n_spacepoints_total,target_planes_total,"
                 "found_dbscan,dbscan_cluster_size,dbscan_matched_points,dbscan_purity,"
                 "found_dbscan_merged,merged_cluster_size,merged_linearity,"
                 "ran_linktree_fallback,found_linktree,"
                 "lt_path_points,lt_matched_points,lt_purity,"
                 "lt_planes,lt_resource_limit,found_combined\n";
  }

  Long64_t n_entries = reco_tree->GetEntries();
  long n_slices_skipped_mismatch = 0, n_slices_seen = 0, n_slices_skipped_no_muon = 0;
  long n_muons_total = 0, n_found_dbscan = 0, n_found_dbscan_merged = 0, n_ran_fallback = 0, n_found_via_fallback = 0;

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

    // --- Pass A: DBSCAN+PCA, exactly as ClusterTruthEfficiency.cpp. ---
    std::vector<TMS_SpacePoint> space_points;
    space_points.reserve(n_space_points);
    for (int i = 0; i < n_space_points; ++i) {
      space_points.emplace_back(sp_x[i], sp_y[i], sp_z[i], sp_x_hitidx[i], sp_y_hitidx[i], sp_time[i], sp_layer.Layer(i, sp_z[i]));
    }
    TMS_SpacePointDBScan dbscan(space_points, dbscan_min_points, dbscan_tolerance);
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

    // point_cluster_id[i]: 0 = unclustered ("noise"), else cluster index+1.
    std::vector<int> point_cluster_id(n_space_points, 0);
    for (size_t c = 0; c < cluster_indices.size(); ++c)
      for (int idx : cluster_indices[c]) point_cluster_id[idx] = (int)c + 1;

    // target_layers_in_slice (the completeness denominator) is independent
    // of which pass looked at the muon, so compute it once up front with a
    // single shared z-layer assignment over the whole slice.
    const std::vector<int> z_layer_whole_slice = TMS_LayerGrouping::GroupIndexOfEachPoint(space_points, lt_config.LayerZTolerance);

    const char *debug_vgid_env = std::getenv("LT_DEBUG_VGID");

    // --- Pass B: for each muon DBSCAN's per-cluster pass missed, build a
    // *merged* point set -- the union of every real cluster (track-like or
    // not) that contains any of the muon's own points, plus that muon's
    // own individually-unclustered ("noise") points specifically (not all
    // noise in the slice, just the muon's own orphaned hits). A single
    // real track can legitimately get split across DBSCAN's own density
    // boundaries (e.g. its sparse tail end going its own way, isolated
    // enough to form or fall into a separate group) -- merging first
    // undoes exactly that split, still without pulling in the rest of the
    // slice's unrelated activity. Re-run PCA on the merged set: if it's
    // now track-like on its own, no Graph Track Finder needed at all (this
    // muon was really findable by DBSCAN+PCA all along, just needed its
    // own split pieces reunited first). Only if it's still not track-like
    // does Graph Track Finder run, on the merged set.
    struct FallbackRun {
      TMS_GraphTrackFinder::Result result;
      std::vector<int> local_to_global;
      std::vector<int> z_layer_local;
    };

    std::unordered_map<int, std::vector<int>> merged_clusters;  // pidx -> merged global point indices
    for (int pidx : muon_particle_idx) {
      TrueLabel label{sp.vgid[pidx], sp.trackid[pidx]};
      if (muon_matches.find(label) != muon_matches.end()) continue;  // DBSCAN's own per-cluster pass already found it

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

      if (debug_vgid_env && label.vgid == std::stoll(debug_vgid_env)) {
        std::cerr << "[LT_DEBUG] vgid=" << label.vgid << " merged " << touched_cluster_ids.size()
                  << " clusters + " << own_noise_points.size() << " own noise points -> "
                  << merged_indices.size() << " total points (vs. whole slice " << n_space_points << ")"
                  << std::endl;
      }

      merged_clusters[pidx] = std::move(merged_indices);
    }

    // Re-run PCA/linearity on each muon's merged set; only queue the
    // Graph Track Finder fallback for the ones still not track-like on their own.
    std::unordered_map<int, FallbackRun> fallback_runs;  // keyed by pidx directly (per-muon, not shared)
    std::unordered_map<int, double> merged_linearity_by_pidx;
    std::unordered_map<int, int> merged_size_by_pidx;
    std::unordered_map<int, bool> found_merged_by_pidx;
    for (auto &kv : merged_clusters) {
      const int pidx = kv.first;
      const std::vector<int> &merged_indices = kv.second;
      TrueLabel label{sp.vgid[pidx], sp.trackid[pidx]};
      merged_size_by_pidx[pidx] = (int)merged_indices.size();
      found_merged_by_pidx[pidx] = false;
      if (merged_indices.size() < kMinClusterSizeForTrack) { merged_linearity_by_pidx[pidx] = 0.0; continue; }

      TMS_SpacePointCluster merged_cluster(space_points, merged_indices);
      merged_linearity_by_pidx[pidx] = merged_cluster.GetLinearity();
      if (merged_cluster.IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack)) {
        std::unordered_map<TrueLabel, int, LabelHash> votes;
        for (int idx : merged_indices) if (point_label[idx].Valid()) votes[point_label[idx]]++;
        TrueLabel owner; int owner_count = 0;
        for (auto &v : votes) if (v.second > owner_count) { owner = v.first; owner_count = v.second; }
        if (owner == label) { found_merged_by_pidx[pidx] = true; continue; }  // no LT needed
      }

      if (merged_indices.size() < lt_config.SeedLength) continue;  // can't possibly seed
      std::vector<TMS_SpacePoint> local_points;
      local_points.reserve(merged_indices.size());
      for (int gi : merged_indices) local_points.push_back(space_points[gi]);
      FallbackRun run;
      run.z_layer_local = TMS_LayerGrouping::GroupIndexOfEachPoint(local_points, lt_config.LayerZTolerance);
      run.result = TMS_GraphTrackFinder::Finder(lt_config).Find(local_points);
      run.local_to_global = merged_indices;
      fallback_runs[pidx] = std::move(run);
      ++n_ran_fallback;
    }

    for (int pidx : muon_particle_idx) {
      ++n_muons_total;
      TrueLabel label{sp.vgid[pidx], sp.trackid[pidx]};

      std::set<int> target_layers_in_slice;
      for (int i = 0; i < n_space_points; ++i)
        if (point_label[i] == label) target_layers_in_slice.insert(z_layer_whole_slice[i]);

      auto match_it = muon_matches.find(label);
      const bool found_dbscan = match_it != muon_matches.end();
      int dbscan_size = 0, dbscan_matched = 0;
      double dbscan_purity = 0.0;
      if (found_dbscan) {
        ++n_found_dbscan;
        int best_c = -1;
        size_t best_size = 0;
        for (int c : match_it->second)
          if (best_c == -1 || clusters[c].GetSize() > best_size) { best_size = clusters[c].GetSize(); best_c = c; }
        dbscan_size = clusters[best_c].GetSize();
        int owner_count = 0;
        for (int idx : cluster_indices[best_c]) if (point_label[idx] == label) ++owner_count;
        dbscan_matched = owner_count;
        dbscan_purity = dbscan_size > 0 ? 100.0 * owner_count / dbscan_size : 0.0;
      }

      bool found_dbscan_merged = false;
      int merged_size = 0;
      double merged_linearity = 0.0;
      if (!found_dbscan) {
        auto it = found_merged_by_pidx.find(pidx);
        if (it != found_merged_by_pidx.end()) found_dbscan_merged = it->second;
        auto sz_it = merged_size_by_pidx.find(pidx);
        if (sz_it != merged_size_by_pidx.end()) merged_size = sz_it->second;
        auto lin_it = merged_linearity_by_pidx.find(pidx);
        if (lin_it != merged_linearity_by_pidx.end()) merged_linearity = lin_it->second;
        if (found_dbscan_merged) ++n_found_dbscan_merged;
      }

      bool ran_fallback_this_muon = false;
      bool found_lt = false;
      int lt_points = 0, lt_matched = 0, lt_planes = 0;
      double lt_purity = 0.0;
      bool lt_resource_limit = false;
      if (!found_dbscan && !found_dbscan_merged) {
        auto run_it = fallback_runs.find(pidx);
        if (run_it != fallback_runs.end()) {
          ran_fallback_this_muon = true;
          const FallbackRun &run = run_it->second;
          lt_resource_limit = run.result.Stats.ResourceLimitReached;
          for (const TMS_GraphTrackFinder::Path &path : run.result.Paths) {
            int matched = 0;
            std::set<int> matched_layers;
            for (std::size_t local_idx : path.SpacePointIndices) {
              const int global_idx = run.local_to_global[local_idx];
              if (point_label[global_idx] == label) {
                ++matched;
                matched_layers.insert(run.z_layer_local[local_idx]);
              }
            }
            if ((int)matched_layers.size() > lt_planes) {
              lt_planes = (int)matched_layers.size();
              lt_points = (int)path.SpacePointIndices.size();
              lt_matched = matched;
            }
          }
          found_lt = lt_planes > 0;
          lt_purity = lt_points > 0 ? 100.0 * lt_matched / lt_points : 0.0;
          if (found_lt) ++n_found_via_fallback;
        }
      }

      const bool found_combined = found_dbscan || found_dbscan_merged || found_lt;

      const float *p4 = &sp.momentum[pidx * 4];
      const double momentum_mag = std::sqrt(p4[0] * p4[0] + p4[1] * p4[1] + p4[2] * p4[2]);
      const double angle_deg = momentum_mag > 0 ? std::acos(p4[2] / momentum_mag) * 180.0 / M_PI : 0.0;

      muons_csv << input_filename << "," << entry << "," << slice_no << "," << spill_no << ","
                << label.vgid << "," << label.trackid << ","
                << (sp.tms_fiducial_start[pidx] ? 1 : 0) << "," << (sp.lar_fiducial_start[pidx] ? 1 : 0) << ","
                << momentum_mag << "," << angle_deg << "," << sp.path_length_tms[pidx] << ","
                << true_nhits_slice[pidx] << "," << n_muons_in_slice << "," << n_space_points << ","
                << (int)target_layers_in_slice.size() << ","
                << (found_dbscan ? 1 : 0) << "," << dbscan_size << "," << dbscan_matched << "," << dbscan_purity << ","
                << (found_dbscan_merged ? 1 : 0) << "," << merged_size << "," << merged_linearity << ","
                << (ran_fallback_this_muon ? 1 : 0) << ","
                << (found_lt ? 1 : 0) << ","
                << lt_points << "," << lt_matched << "," << lt_purity << "," << lt_planes << ","
                << (lt_resource_limit ? 1 : 0) << "," << (found_combined ? 1 : 0) << "\n";
    }

    if (n_slices_seen % 25 == 0) {
      std::cout << "  entry=" << entry << "/" << n_entries << " spill=" << spill_no << " slice=" << slice_no
                << " nSP=" << n_space_points << " muons_in_slice=" << n_muons_in_slice
                << " slices_seen=" << n_slices_seen << " fallback_runs=" << n_ran_fallback << std::endl;
    }
  }

  muons_csv.close();

  std::cout << "Done. " << n_muons_total << " muon candidates." << std::endl;
  std::cout << "Found by DBSCAN+PCA (per-cluster): " << n_found_dbscan << " ("
            << (n_muons_total > 0 ? 100.0 * n_found_dbscan / n_muons_total : 0.0) << "%)" << std::endl;
  std::cout << "Additionally found by merged-cluster PCA (no Graph Track Finder needed): "
            << n_found_dbscan_merged << std::endl;
  std::cout << "Graph Track Finder fallback ran " << n_ran_fallback << " times (once per still-missing"
               " muon, on its merged cluster set), recovered " << n_found_via_fallback
            << " additional muons" << std::endl;
  const long n_total_found = n_found_dbscan + n_found_dbscan_merged + n_found_via_fallback;
  std::cout << "Combined found: " << n_total_found << " ("
            << (n_muons_total > 0 ? 100.0 * n_total_found / n_muons_total : 0.0)
            << "%)" << std::endl;
  std::cout << "Slices seen (>=1 findable muon): " << n_slices_seen
            << ", skipped (no findable muon): " << n_slices_skipped_no_muon
            << ", skipped (nTrueParticles mismatch): " << n_slices_skipped_mismatch << std::endl;
  std::cout << "Wrote " << muons_csv_path << (append ? " (appended)" : "") << std::endl;

  return 0;
}
