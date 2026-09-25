// Object-first truth validation of track finding + fitting: every track the
// reconstruction produces is classified against truth, instead of starting
// from a known true muon and asking whether it was found (the approach of
// ClusterTruthEfficiency / KalmanFollowerTruthEfficiency, which is blind to
// fake and duplicate tracks by construction).
//
// Reconstruction under test (reco-only, no truth used): DBSCAN over each
// slice's space points; every track-like cluster, largest first, goes to
// TMS_IterativeTrackFit::FitCluster(), which fits it with the Kalman
// follower and -- depending on the split mode -- keeps fitting the
// unclaimed remainder of clusters that merge more than one particle.
//
// Truth is used only afterwards, to score:
//   tracks CSV: one row per fitted track -- its truth owner (plurality of
//     "strict" labels: X-hit and Y-hit truth agree), strict purity, whether
//     the owner is a muon, and whether it duplicates an earlier track's owner.
//   muons CSV: one row per true muon with >= 5 true hits in the slice (the
//     same population as the muon-first tools) -- how many tracks it owns,
//     and the best one's strict plane coverage.
//
// Split mode (environment ITF_SPLIT): 0 = off (one fit per cluster, today's
// behavior), 1 = flagged clusters only (default), 2 = every cluster.
// Kalman follower environment hooks as in KalmanFollowerTruthEfficiency.
//
// Scope, deliberately: only track-like DBSCAN clusters are fitted. Muons
// buried in a non-track-like cluster (the muon-first tools' merge-and-re-PCA
// and graph-search stages, which pick clusters using the target's own truth
// labels) are not attempted here, so absolute efficiency is lower than the
// muon-first tools report. Use this tool to compare reconstruction variants
// on equal footing, and for fake/duplicate rates.

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
#include "TMS_IterativeTrackFit.h"
#include "TMS_KalmanFollower.h"
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
  bool operator<(const TrueLabel &o) const { return vgid != o.vgid ? vgid < o.vgid : trackid < o.trackid; }
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

// Same convention as ClusterTruthEfficiency.cpp / KalmanFollowerTruthEfficiency.cpp.
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
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0]
              << " <edep_sim_geom_file> <input_reco_tree.root> <tracks_output.csv> <muons_output.csv>" << std::endl;
    return -1;
  }
  const std::string geom_filename = argv[1];
  const std::string input_filename = argv[2];

  // DBSCAN+PCA params: identical to the muon-first tools.
  const unsigned int min_points = 5;
  const double kLinearityThreshold = 0.8;
  const size_t kMinClusterSizeForTrack = 5;

  TFile geom_input(geom_filename.c_str());
  TGeoManager *geom = geom_input.IsZombie() ? nullptr : (TGeoManager *)geom_input.Get("EDepSimGeometry");
  if (!geom) {
    std::cerr << "No EDepSimGeometry in " << geom_filename << std::endl;
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

  const RegionFieldModel field;
  TMS_KalmanFollower::Config follower_config;
  if (const char *v = std::getenv("KF_MAX_HEAD_SKIP")) follower_config.MaxHeadSkip = std::atoi(v);
  if (const char *v = std::getenv("KF_MAX_TRIPLETS")) follower_config.MaxTripletHypotheses = std::atoi(v);
  if (const char *v = std::getenv("KF_USE_TIME")) follower_config.UseTimeInSelection = std::atoi(v) != 0;
  if (const char *v = std::getenv("KF_TIME_SIGMA")) follower_config.TimeSigmaNs = std::atof(v);
  if (const char *v = std::getenv("KF_TIME_GATE")) follower_config.TimeGateNSigma = std::atof(v);
  const TMS_KalmanFollower::Follower follower(follower_config, field);

  TMS_IterativeTrackFit::Config itf;
  itf.DBScanMinPoints = dbscan_min_points;
  itf.DBScanTolerance = dbscan_tolerance;
  itf.LinearityThreshold = kLinearityThreshold;
  itf.MinClusterSizeForTrack = kMinClusterSizeForTrack;
  if (const char *v = std::getenv("ITF_SPLIT")) {
    const int mode = std::atoi(v);
    itf.Mode = mode == 0 ? TMS_IterativeTrackFit::Config::SplitMode::Off
             : mode == 2 ? TMS_IterativeTrackFit::Config::SplitMode::Always
                         : TMS_IterativeTrackFit::Config::SplitMode::Flagged;
  }
  if (const char *v = std::getenv("ITF_MAX_TRACKS")) itf.MaxTracksPerCluster = std::atoi(v);
  if (const char *v = std::getenv("ITF_MIN_HITS")) itf.MinHitsPerTrack = std::atoi(v);
  if (const char *v = std::getenv("ITF_MIN_SPLIT_HITS")) itf.MinHitsPerSplitTrack = std::atoi(v);
  if (const char *v = std::getenv("ITF_RESTRICT_SPLIT")) itf.RestrictSplitFitToCluster = std::atoi(v) != 0;

  TFile input(input_filename.c_str());
  TTree *reco_tree = input.IsZombie() ? nullptr : (TTree *)input.Get("Reco_Tree");
  TTree *truth_info = input.IsZombie() ? nullptr : (TTree *)input.Get("Truth_Info");
  TTree *truth_spill = input.IsZombie() ? nullptr : (TTree *)input.Get("Truth_Spill");
  if (!reco_tree || !truth_info || !truth_spill) {
    std::cerr << "Input file is missing Reco_Tree/Truth_Info/Truth_Spill" << std::endl;
    return -1;
  }

  // --- Truth_Spill into memory, keyed by SpillNo. ---
  int spill_no_ts = 0, n_tp_ts = 0;
  static std::vector<long long> vgid_ts(kMaxTrueParticles);
  static std::vector<int> trackid_ts(kMaxTrueParticles), pdg_ts(kMaxTrueParticles), parent_ts(kMaxTrueParticles);
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

  int n_tp_ti = 0;
  static std::vector<int> true_nhits_slice(kMaxTrueParticles);
  truth_info->SetBranchAddress("nTrueParticles", &n_tp_ti);
  truth_info->SetBranchAddress("TrueNHitsInSlice", true_nhits_slice.data());

  int n_space_points = 0, spill_no = 0, slice_no = 0;
  static std::vector<float> sp_x(kMaxSpacePoints), sp_y(kMaxSpacePoints), sp_z(kMaxSpacePoints), sp_time(kMaxSpacePoints);
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

  std::ofstream tracks_csv(argv[3]), muons_csv(argv[4]);
  tracks_csv << "sourcefile,entry,slice,cluster_id,cluster_size,cluster_flagged,iteration,n_hits,n_nodes,converged,"
                "momentum_mev,t0_ns,first_z,last_z,owner_vertexglobalid,owner_trackid,owner_pdg,owner_is_muon,"
                "owner_strict_hits,strict_purity,n_strict_labeled,duplicate_of_earlier\n";
  muons_csv << "sourcefile,entry,slice,vertexglobalid,trackid,vertex_in_tms,vertex_in_lar_fiducial,true_hits_in_slice,"
               "target_planes_strict,n_tracks_owned,best_planes_covered,best_completeness_pct,best_purity,"
               "best_iteration\n";

  long n_tracks = 0, n_split_tracks = 0;
  for (Long64_t entry = 0; entry < reco_tree->GetEntries(); ++entry) {
    reco_tree->GetEntry(entry);
    truth_info->GetEntry(entry);
    auto spill_it = spills.find(spill_no);
    if (spill_it == spills.end() || n_space_points <= 0) continue;
    const SpillParticles &sp = spill_it->second;
    if (sp.n != n_tp_ti) continue;

    auto collapse = [&](const TrueLabel &raw) -> TrueLabel {
      if (!raw.Valid()) return raw;
      auto it = sp.index_of.find(raw);
      if (it == sp.index_of.end()) return raw;
      return TrueLabel{raw.vgid, sp.collapsed_trackid[it->second]};
    };
    // Strict label: valid only where X-hit and Y-hit truth agree.
    std::vector<TrueLabel> strict(n_space_points);
    std::vector<TMS_SpacePoint> points;
    points.reserve(n_space_points);
    for (int i = 0; i < n_space_points; ++i) {
      const TrueLabel lx = collapse(TrueLabel{sp_x_vgid[i], sp_x_trackid[i]});
      const TrueLabel ly = collapse(TrueLabel{sp_y_vgid[i], sp_y_trackid[i]});
      if (lx.Valid() && lx == ly) strict[i] = lx;
      points.emplace_back(sp_x[i], sp_y[i], sp_z[i], sp_x_hitidx[i], sp_y_hitidx[i], sp_time[i], sp_layer.Layer(i, sp_z[i]));
    }
    const std::vector<int> z_layer = TMS_LayerGrouping::GroupIndexOfEachPoint(points, 1.0);

    // --- Reconstruction (no truth). ---
    TMS_SpacePointDBScan dbscan(points, dbscan_min_points, dbscan_tolerance);
    std::vector<std::vector<int> > clusters = dbscan.RunAndGetClusterIndices();
    std::vector<std::size_t> order;
    for (std::size_t c = 0; c < clusters.size(); ++c) {
      TMS_SpacePointCluster cl(points, clusters[c]);
      if (cl.IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack)) order.push_back(c);
    }
    // Largest first, so a big merged cluster claims its hits before any
    // small fragment beside it can.
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return clusters[a].size() > clusters[b].size(); });

    TMS_IterativeTrackFit::ClaimedHits claimed;
    struct Scored {
      TrueLabel owner;
      double purity = 0;
      int planes = 0;
      int iteration = 0;
    };
    std::vector<Scored> scored;
    std::set<TrueLabel> owners_seen;
    for (std::size_t c : order) {
      const std::vector<TMS_IterativeTrackFit::Track> tracks =
          TMS_IterativeTrackFit::FitCluster(points, clusters[c], follower, itf, claimed);
      for (const TMS_IterativeTrackFit::Track &t : tracks) {
        // --- Scoring (truth). ---
        std::unordered_map<TrueLabel, int, LabelHash> votes;
        std::unordered_map<TrueLabel, std::set<int>, LabelHash> planes;
        int n_hits = 0, n_labeled = 0;
        for (const TMS_KalmanFollower::FollowedNode &node : t.Fit.Nodes) {
          if (!node.HasHit) continue;
          ++n_hits;
          const TrueLabel &l = strict[node.ChosenSpacePointIndex];
          if (!l.Valid()) continue;
          ++n_labeled;
          ++votes[l];
          planes[l].insert(z_layer[node.ChosenSpacePointIndex]);
        }
        TrueLabel owner;
        int owner_n = 0;
        for (const auto &kv : votes)
          if (kv.second > owner_n || (kv.second == owner_n && kv.first < owner)) {
            owner = kv.first;
            owner_n = kv.second;
          }
        int owner_pdg = 0;
        if (owner.Valid()) {
          auto it = sp.index_of.find(owner);
          if (it != sp.index_of.end()) owner_pdg = sp.pdg[it->second];
        }
        const bool duplicate = owner.Valid() && owners_seen.count(owner) > 0;
        if (owner.Valid()) owners_seen.insert(owner);
        const double purity = n_hits > 0 ? static_cast<double>(owner_n) / n_hits : 0.0;
        scored.push_back({owner, purity, owner.Valid() ? (int)planes[owner].size() : 0, t.Iteration});
        ++n_tracks;
        if (t.Iteration > 0) ++n_split_tracks;

        tracks_csv << input_filename << "," << entry << "," << slice_no << "," << c + 1 << "," << clusters[c].size()
                   << "," << (t.ClusterFlagged ? 1 : 0) << "," << t.Iteration << "," << n_hits << ","
                   << t.Fit.Nodes.size() << "," << (t.Fit.Converged ? 1 : 0) << "," << t.Fit.MomentumMeV << ","
                   << t.Fit.TrackT0Ns << "," << (t.Fit.Nodes.empty() ? 0.0 : t.Fit.Nodes.front().Z) << ","
                   << (t.Fit.Nodes.empty() ? 0.0 : t.Fit.Nodes.back().Z) << "," << owner.vgid << "," << owner.trackid
                   << "," << owner_pdg << "," << (std::abs(owner_pdg) == 13 ? 1 : 0) << "," << owner_n << ","
                   << purity << "," << n_labeled << "," << (duplicate ? 1 : 0) << "\n";
      }
    }

    // --- Muon side: every true muon with >= min_points true hits here. ---
    for (int i = 0; i < sp.n; ++i) {
      if (std::abs(sp.pdg[i]) != 13 || true_nhits_slice[i] < (int)min_points) continue;
      const TrueLabel label{sp.vgid[i], sp.trackid[i]};
      std::set<int> target_planes;
      for (int p = 0; p < n_space_points; ++p)
        if (strict[p] == label) target_planes.insert(z_layer[p]);
      int n_owned = 0, best_planes = 0, best_iteration = -1;
      double best_purity = 0.0;
      for (const Scored &s : scored) {
        if (!(s.owner == label)) continue;
        ++n_owned;
        if (s.planes > best_planes) {
          best_planes = s.planes;
          best_purity = s.purity;
          best_iteration = s.iteration;
        }
      }
      muons_csv << input_filename << "," << entry << "," << slice_no << "," << label.vgid << "," << label.trackid
                << "," << (sp.tms_fiducial_start[i] ? 1 : 0) << "," << (sp.lar_fiducial_start[i] ? 1 : 0) << ","
                << true_nhits_slice[i] << "," << target_planes.size() << "," << n_owned << "," << best_planes << ","
                << (target_planes.empty() ? 0.0 : 100.0 * best_planes / target_planes.size()) << ","
                << best_purity << "," << best_iteration << "\n";
    }
  }
  std::cout << "Tracks fitted: " << n_tracks << " (from split remainders: " << n_split_tracks << ")" << std::endl;
  return 0;
}
