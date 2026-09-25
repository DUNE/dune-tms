// Truth-based validation of the space-point clustering (Stage 2.5 pipeline):
// track-finding efficiency, hit-finding efficiency and hit purity for muons,
// matched to real DBSCAN+PCA clusters by plurality vote against the space
// points' own truth linkage (SpacePointX/YTrueVertexGlobalId/TrueTrackId,
// already written per Reco_Tree entry -- no new branches needed).
//
// Truth_Info (one row per slice, same entry order as Reco_Tree) carries
// TrueNHitsInSlice[nTrueParticles], the number of ped-surviving reconstructed
// hits in *this slice* linked to each true particle. Truth_Spill (one row per
// spill) carries that same particle list's identity/kinematics (VertexGlobalID,
// TrackId, PDG, MomentumTMSStart, TruePathLengthInTMS). Both arrays share the
// same index/order for a given spill (verified empirically this session: the
// magic muon, vertex 1000000014002130, is index 714 in both Truth_Info[entry
// for slice 98] and Truth_Spill[spill 4]) -- so a slice's Truth_Info row is
// joined to its spill's Truth_Spill row via SpillNo, then walked by shared
// index.
//
// Convention (agreed with the user): plurality matching (a cluster matches
// whichever true track contributes the most of its space points), muons only
// (|PDG|==13), and a muon only counts as a "findable" target if it has at
// least min_points true hits in this slice -- reusing DBSCAN's own min_points
// floor rather than inventing a separate threshold, since nothing below that
// could ever form a cluster in the first place.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "TFile.h"
#include "TGeoManager.h"
#include "TTree.h"

#include "TMS_Geom.h"
#include "TMS_SpacePoint.h"
#include "TMS_SpacePointCluster.h"
#include "TMS_SpacePointDBScan.h"
#include "SpacePointLayerInput.h"

namespace {

const int kMaxSpacePoints = 10000;    // matches __TMS_MAX_SPACEPOINTS__
const int kMaxTrueParticles = 20000;  // matches __TMS_MAX_TRUE_PARTICLES__

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
  std::vector<int> parent_trackid;  // GetParent(): parent's TrackId within the same vertex, or -1
  std::vector<float> momentum;      // [n][4], flattened px,py,pz,E
  std::vector<float> path_length_tms;
  std::vector<bool> tms_fiducial_start;  // IsInsideTMS(birth position) -- true = born inside TMS
  std::vector<bool> lar_fiducial_start;  // IsInsideLarFiducial(birth position) -- true = born in ND-LAr fiducial
  std::unordered_map<TrueLabel, int, LabelHash> index_of;  // (vgid,trackid) -> particle index
  std::vector<int> collapsed_trackid;  // see CollapseTrackId()
};

// A raw G4 TrackId is per-hit-segment bookkeeping (edep-sim's own
// TG4HitSegment::PrimaryId is documented as "the track id of the most
// important particle for THIS hit segment") -- it is NOT collapsed to any
// ancestor. A muon's own delta-ray gets its own fresh TrackId the instant
// it's created, even though physically it's just a continuation of the
// muon's own ionization trail (the same fragmentation PR #317 fixed on the
// simulation side). Matching truth on the raw TrackId therefore treats a
// muon's own delta-ray hits as "a different particle" -- inflating both the
// X/Y-hit disagreement rate and depressing purity/completeness/efficiency
// for reasons that have nothing to do with clustering quality.
//
// Fix: walk the Parent chain (within the same vertex -- TrackIds are only
// unique per-vertex, not file-wide) and collapse to the first muon ancestor
// found (stopping there, so a decay-in-flight muon still claims its own
// descendants rather than being merged into its parent pion); if no muon
// ancestor exists, collapse to the topmost primary (Parent==-1) instead, so
// non-muon secondaries (e.g. a proton's own delta-rays) still group sensibly
// for the X/Y-agreement diagnostic.
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
  if (argc != 5 && argc != 6) {
    std::cerr << "Usage: " << argv[0]
              << " <edep_sim_geom_file> <input_reco_tree.root> <muons_output.csv> <clusters_output.csv>"
                 " [display_output_prefix]"
              << std::endl;
    std::cerr << "  display_output_prefix (optional): if given, also dumps <prefix>_points.csv and"
              << " <prefix>_pca.csv -- every space point + every cluster's PCA (not just track-like) for"
              << " the home slice of every ND-LAr-fiducial-origin muon candidate, for visualization."
              << std::endl;
    return -1;
  }

  const std::string geom_filename = argv[1];
  const std::string input_filename = argv[2];
  const std::string muons_csv_path = argv[3];
  const std::string clusters_csv_path = argv[4];
  const std::string display_prefix = argc == 6 ? argv[5] : "";

  // Clustering params: TMS_SpacePointDBScan::DefaultTolerance() (set below,
  // once the bar pitch is known) and these.
  const unsigned int min_points = 5;
  const double kLinearityThreshold = 0.8;
  const size_t kMinClusterSizeForTrack = 5;

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

  // --- Pass 1: load Truth_Spill entirely into memory, keyed by SpillNo. ---
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

  // --- Pass 2: loop Reco_Tree / Truth_Info entries in lockstep. ---
  int n_space_points = 0, spill_no = 0, slice_no = 0;
  static std::vector<float> sp_x(kMaxSpacePoints), sp_y(kMaxSpacePoints), sp_z(kMaxSpacePoints);
  static std::vector<float> sp_time(kMaxSpacePoints);
  static std::vector<long long> sp_x_vgid(kMaxSpacePoints), sp_y_vgid(kMaxSpacePoints);
  static std::vector<int> sp_x_trackid(kMaxSpacePoints), sp_y_trackid(kMaxSpacePoints);
  reco_tree->SetBranchAddress("nSpacePoints", &n_space_points);
  reco_tree->SetBranchAddress("SpacePointX", sp_x.data());
  reco_tree->SetBranchAddress("SpacePointY", sp_y.data());
  reco_tree->SetBranchAddress("SpacePointZ", sp_z.data());
  const SpacePointLayerInput sp_layer(reco_tree, kMaxSpacePoints);
  reco_tree->SetBranchAddress("SpacePointTime", sp_time.data());
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

  std::ofstream muons_csv(muons_csv_path);
  muons_csv << "entry,slice,spill,vertexglobalid,trackid,vertex_in_tms,vertex_in_lar_fiducial,momentum_mag,angle_deg,true_dir_x,"
               "true_dir_y,true_dir_z,"
               "path_length_mm,"
               "true_hits_in_slice,n_muons_in_slice,n_spacepoints_total,n_spacepoints_clustered,"
               "n_spacepoints_track_like,found,n_clusters_matched,best_cluster_size,"
               "best_cluster_matched_points,best_cluster_purity,best_cluster_linearity,"
               "best_cluster_dir_x,best_cluster_dir_y,best_cluster_dir_z\n";

  std::ofstream clusters_csv(clusters_csv_path);
  clusters_csv << "entry,slice,cluster_id,n_points,owner_vertexglobalid,owner_trackid,owner_pdg,"
                  "purity,is_muon_matched,linearity\n";

  // Optional (CTE_CLUSTER_DETAIL_CSV=<path>): one row per DBSCAN cluster,
  // track-like OR NOT (the main clusters CSV above only has track-like ones),
  // with what's needed to study (a) whether the PCA linearity cut could be
  // tightened to reject clusters that merge two real muons, and (b) whether
  // those muons are separable in time. Per cluster: all three PCA
  // eigenvalues, the per-layer transverse span (two parallel muons a few bar
  // pitches apart are nearly perfectly "linear" by (l1-l2)/l1, but show up
  // as a wide span within each layer), and the top-2 truth owners with the
  // mean time of each one's "pure" points (X-hit and Y-hit truth agree, so
  // the space-point time -- the X/Y hit-time average -- belongs to that one
  // particle and isn't a ghost mixing two particles' times).
  const char *detail_csv_env = std::getenv("CTE_CLUSTER_DETAIL_CSV");
  std::ofstream detail_csv;
  if (detail_csv_env) {
    detail_csv.open(detail_csv_env);
    detail_csv << "entry,slice,cluster_id,n_points,is_track_like,linearity,l1,l2,l3,n_layers,z_extent,"
                  "median_layer_span,frac_layers_span_gt150,"
                  "owner_vertexglobalid,owner_trackid,owner_pdg,owner_count,"
                  "second_vertexglobalid,second_trackid,second_pdg,second_count,"
                  "owner_pure_n,owner_pure_tmean,second_pure_n,second_pure_tmean,"
                  "n_muons_ge5pure,cluster_trms\n";
  }

  const bool dump_display = !display_prefix.empty();
  std::ofstream display_points_csv, display_pca_csv;
  if (dump_display) {
    display_points_csv.open(display_prefix + "_points.csv");
    display_points_csv << "entry,slice,muon_vertexglobalid,muon_trackid,x,y,z,cluster_id,linearity,"
                           "is_track_like,is_this_muon\n";
    display_pca_csv.open(display_prefix + "_pca.csv");
    display_pca_csv << "entry,slice,muon_vertexglobalid,muon_trackid,cluster_id,n,is_track_like,cx,cy,cz,"
                        "eval0,ex0,ey0,ez0,eval1,ex1,ey1,ez1,eval2,ex2,ey2,ez2\n";
  }

  Long64_t n_entries = reco_tree->GetEntries();
  long xy_mismatch_count_raw = 0, xy_agree_count_raw = 0;
  long xy_mismatch_count = 0, xy_agree_count = 0;  // after parent-chain collapse
  long xy_mismatch_diff_vertex = 0, xy_mismatch_same_vertex = 0;
  long n_slices_skipped_mismatch = 0;
  long n_muons_total = 0, n_muons_found = 0;
  long n_clusters_total = 0, n_clusters_muon_matched = 0;

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

    // Primary truth label per space point: X-hit link, falling back to Y.
    // Each raw label is collapsed via the particle's Parent chain (see
    // CollapseTrackId) before use, so a muon's own delta-ray hits are
    // attributed to the muon rather than counted as "a different particle".
    auto collapse = [&](const TrueLabel &raw) -> TrueLabel {
      if (!raw.Valid()) return raw;
      auto it = sp.index_of.find(raw);
      if (it == sp.index_of.end()) return raw;
      return TrueLabel{raw.vgid, sp.collapsed_trackid[it->second]};
    };
    std::vector<TrueLabel> point_label(n_space_points);
    // Set only when the X-hit and Y-hit truth labels agree (after collapse);
    // invalid otherwise. Used by the optional detail CSV's per-owner times.
    std::vector<TrueLabel> point_pure_label(n_space_points);
    for (int i = 0; i < n_space_points; ++i) {
      TrueLabel x_label_raw{sp_x_vgid[i], sp_x_trackid[i]};
      TrueLabel y_label_raw{sp_y_vgid[i], sp_y_trackid[i]};
      if (x_label_raw.Valid() && y_label_raw.Valid()) {
        if (y_label_raw == x_label_raw) ++xy_agree_count_raw; else ++xy_mismatch_count_raw;
      }
      const TrueLabel x_label = collapse(x_label_raw);
      const TrueLabel y_label = collapse(y_label_raw);
      if (x_label.Valid()) {
        point_label[i] = x_label;
        if (y_label.Valid()) {
          if (y_label == x_label) {
            ++xy_agree_count;
            point_pure_label[i] = x_label;
          } else {
            ++xy_mismatch_count;
            if (x_label.vgid != y_label.vgid) ++xy_mismatch_diff_vertex;
            else ++xy_mismatch_same_vertex;
          }
        }
      } else {
        point_label[i] = y_label;
      }
    }

    // Muon candidates in this slice: |PDG|==13, >= min_points true hits here.
    std::vector<int> muon_particle_idx;
    for (int i = 0; i < sp.n; ++i) {
      if (std::abs(sp.pdg[i]) == 13 && true_nhits_slice[i] >= (int)min_points) {
        muon_particle_idx.push_back(i);
      }
    }
    const int n_muons_in_slice = (int)muon_particle_idx.size();

    // Run the real clustering, exactly as SpillClusterDisplay does per slice.
    std::vector<TMS_SpacePoint> space_points;
    space_points.reserve(n_space_points);
    for (int i = 0; i < n_space_points; ++i) {
      space_points.emplace_back(sp_x[i], sp_y[i], sp_z[i], -1, -1, 0.0, sp_layer.Layer(i, sp_z[i]));
    }
    TMS_SpacePointDBScan dbscan(space_points, dbscan_min_points, dbscan_tolerance);
    std::vector<std::vector<int>> cluster_indices = dbscan.RunAndGetClusterIndices();
    std::vector<TMS_SpacePointCluster> clusters;
    clusters.reserve(cluster_indices.size());
    for (auto &indices : cluster_indices) clusters.emplace_back(space_points, indices);

    std::vector<int> point_cluster_id(n_space_points, 0);
    for (size_t c = 0; c < cluster_indices.size(); ++c)
      for (int idx : cluster_indices[c]) point_cluster_id[idx] = (int)c + 1;

    if (detail_csv_env) {
      for (size_t c = 0; c < clusters.size(); ++c) {
        const auto &cl = clusters[c];
        const std::vector<int> &idxs = cluster_indices[c];
        // Top-2 owners by plurality vote (same single-sided label as above).
        std::unordered_map<TrueLabel, int, LabelHash> votes;
        for (int idx : idxs)
          if (point_label[idx].Valid()) votes[point_label[idx]]++;
        TrueLabel first, second;
        int first_n = 0, second_n = 0;
        for (auto &kv : votes) {
          if (kv.second > first_n) {
            second = first; second_n = first_n;
            first = kv.first; first_n = kv.second;
          } else if (kv.second > second_n) {
            second = kv.first; second_n = kv.second;
          }
        }
        auto pdg_of = [&](const TrueLabel &l) {
          auto it = sp.index_of.find(l);
          return (l.Valid() && it != sp.index_of.end()) ? sp.pdg[it->second] : 0;
        };
        // Mean time of each owner's pure points, and how many distinct
        // muons have >= 5 pure points in the cluster.
        std::unordered_map<TrueLabel, std::pair<int, double>, LabelHash> pure;  // label -> (n, sum t)
        double tsum = 0, tsum2 = 0;
        for (int idx : idxs) {
          tsum += sp_time[idx];
          tsum2 += sp_time[idx] * sp_time[idx];
          if (point_pure_label[idx].Valid()) {
            auto &e = pure[point_pure_label[idx]];
            e.first++;
            e.second += sp_time[idx];
          }
        }
        int n_muons_ge5 = 0;
        for (auto &kv : pure)
          if (kv.second.first >= 5 && std::abs(pdg_of(kv.first)) == 13) ++n_muons_ge5;
        const double tmean = tsum / idxs.size();
        const double trms = std::sqrt(std::max(0.0, tsum2 / idxs.size() - tmean * tmean));
        auto pure_stats = [&](const TrueLabel &l, int &n, double &t) {
          n = 0; t = 0;
          auto it = pure.find(l);
          if (l.Valid() && it != pure.end()) { n = it->second.first; t = it->second.second / n; }
        };
        int first_pn, second_pn;
        double first_pt, second_pt;
        pure_stats(first, first_pn, first_pt);
        pure_stats(second, second_pn, second_pt);
        // Per-layer transverse span: max(x range, y range) of the cluster's
        // points in each point layer, then the median over layers with >= 2 points.
        std::map<int, std::array<double, 4>> layer_box;  // layer -> xmin,xmax,ymin,ymax
        std::map<int, int> layer_n;
        double zmin = 1e18, zmax = -1e18;
        for (int idx : idxs) {
          const int pl = space_points[idx].GetLayer();
          auto it = layer_box.find(pl);
          if (it == layer_box.end()) {
            layer_box[pl] = {sp_x[idx], sp_x[idx], sp_y[idx], sp_y[idx]};
          } else {
            auto &b = it->second;
            b[0] = std::min<double>(b[0], sp_x[idx]); b[1] = std::max<double>(b[1], sp_x[idx]);
            b[2] = std::min<double>(b[2], sp_y[idx]); b[3] = std::max<double>(b[3], sp_y[idx]);
          }
          layer_n[pl]++;
          zmin = std::min<double>(zmin, sp_z[idx]);
          zmax = std::max<double>(zmax, sp_z[idx]);
        }
        std::vector<double> spans;
        for (auto &kv : layer_box)
          if (layer_n[kv.first] >= 2)
            spans.push_back(std::max(kv.second[1] - kv.second[0], kv.second[3] - kv.second[2]));
        double median_span = 0, frac_wide = 0;
        if (!spans.empty()) {
          std::sort(spans.begin(), spans.end());
          median_span = spans[spans.size() / 2];
          for (double s : spans) if (s > 150.0) frac_wide += 1.0;
          frac_wide /= spans.size();
        }
        const auto &ev = cl.GetEigenvalues();
        detail_csv << entry << "," << slice_no << "," << (c + 1) << "," << cl.GetSize() << ","
                   << (cl.IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack) ? 1 : 0) << ","
                   << cl.GetLinearity() << "," << ev[0] << "," << ev[1] << "," << ev[2] << ","
                   << layer_box.size() << "," << (zmax - zmin) << "," << median_span << "," << frac_wide << ","
                   << first.vgid << "," << first.trackid << "," << pdg_of(first) << "," << first_n << ","
                   << second.vgid << "," << second.trackid << "," << pdg_of(second) << "," << second_n << ","
                   << first_pn << "," << first_pt << "," << second_pn << "," << second_pt << ","
                   << n_muons_ge5 << "," << trms << "\n";
      }
    }

    // Per track-like cluster: plurality vote -> owner label + owner_count.
    // muon_matches[owner_label] accumulates every track-like cluster owned by
    // that label, for the muon-side aggregation below.
    std::unordered_map<TrueLabel, std::vector<int>, LabelHash> muon_matches;  // label -> cluster indices (0-based)
    for (size_t c = 0; c < clusters.size(); ++c) {
      const auto &cl = clusters[c];
      const bool is_track_like = cl.IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack);
      if (!is_track_like) continue;
      ++n_clusters_total;

      std::unordered_map<TrueLabel, int, LabelHash> votes;
      for (int idx : cluster_indices[c]) {
        if (point_label[idx].Valid()) votes[point_label[idx]]++;
      }
      TrueLabel owner;
      int owner_count = 0;
      for (auto &kv : votes) {
        if (kv.second > owner_count) { owner = kv.first; owner_count = kv.second; }
      }
      const double purity = cl.GetSize() > 0 ? (double)owner_count / cl.GetSize() : 0.0;
      int owner_pdg = 0;
      auto owner_idx_it = sp.index_of.find(owner);
      if (owner.Valid() && owner_idx_it != sp.index_of.end()) owner_pdg = sp.pdg[owner_idx_it->second];
      const bool is_muon_matched = owner.Valid() && std::abs(owner_pdg) == 13;
      if (is_muon_matched) {
        ++n_clusters_muon_matched;
        muon_matches[owner].push_back((int)c);
      }

      clusters_csv << entry << "," << slice_no << "," << (c + 1) << "," << cl.GetSize() << "," << owner.vgid << ","
                   << owner.trackid << "," << owner_pdg << "," << purity << "," << (is_muon_matched ? 1 : 0) << ","
                   << cl.GetLinearity() << "\n";
    }

    // Per muon candidate: walk the funnel and aggregate over matched clusters.
    for (int pidx : muon_particle_idx) {
      ++n_muons_total;
      TrueLabel label{sp.vgid[pidx], sp.trackid[pidx]};

      int n_total = 0, n_clustered = 0, n_track_like = 0;
      for (int i = 0; i < n_space_points; ++i) {
        if (!(point_label[i] == label)) continue;
        ++n_total;
        if (point_cluster_id[i] > 0) {
          ++n_clustered;
          if (clusters[point_cluster_id[i] - 1].IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack)) ++n_track_like;
        }
      }

      auto match_it = muon_matches.find(label);
      const int n_clusters_matched = match_it != muon_matches.end() ? (int)match_it->second.size() : 0;
      const bool found = n_clusters_matched > 0;
      if (found) ++n_muons_found;

      int best_cluster_size = 0, best_cluster_matched_points = 0;
      double best_cluster_purity = 0, best_cluster_linearity = 0;
      double dir_x = 0, dir_y = 0, dir_z = 0;
      if (found) {
        int best_c = -1;
        size_t best_size = 0;
        for (int c : match_it->second) {
          if (best_c == -1 || clusters[c].GetSize() > best_size) { best_size = clusters[c].GetSize(); best_c = c; }
        }
        const auto &cl = clusters[best_c];
        best_cluster_size = cl.GetSize();
        best_cluster_linearity = cl.GetLinearity();
        int owner_count = 0;
        for (int idx : cluster_indices[best_c]) if (point_label[idx] == label) ++owner_count;
        best_cluster_matched_points = owner_count;
        best_cluster_purity = cl.GetSize() > 0 ? (double)owner_count / cl.GetSize() : 0.0;
        const auto &eigenvectors = cl.GetEigenvectors();
        dir_x = eigenvectors[0][0]; dir_y = eigenvectors[0][1]; dir_z = eigenvectors[0][2];
      }

      const float *p4 = &sp.momentum[pidx * 4];
      const double momentum_mag = std::sqrt(p4[0] * p4[0] + p4[1] * p4[1] + p4[2] * p4[2]);
      const double angle_deg = momentum_mag > 0 ? std::acos(p4[2] / momentum_mag) * 180.0 / M_PI : 0.0;
      const double true_dir_x = momentum_mag > 0 ? p4[0] / momentum_mag : 0.0;
      const double true_dir_y = momentum_mag > 0 ? p4[1] / momentum_mag : 0.0;
      const double true_dir_z = momentum_mag > 0 ? p4[2] / momentum_mag : 0.0;

      muons_csv << entry << "," << slice_no << "," << spill_no << "," << label.vgid << "," << label.trackid << ","
                << (sp.tms_fiducial_start[pidx] ? 1 : 0) << "," << (sp.lar_fiducial_start[pidx] ? 1 : 0) << ","
                << momentum_mag << "," << angle_deg << ","
                << true_dir_x << "," << true_dir_y << "," << true_dir_z
                << "," << sp.path_length_tms[pidx] << ","
                << true_nhits_slice[pidx] << "," << n_muons_in_slice << "," << n_total << "," << n_clustered << ","
                << n_track_like << "," << (found ? 1 : 0) << "," << n_clusters_matched << "," << best_cluster_size
                << "," << best_cluster_matched_points << "," << best_cluster_purity << "," << best_cluster_linearity
                << "," << dir_x << "," << dir_y << "," << dir_z << "\n";

      // Display dump: every space point + every cluster's PCA (not just
      // track-like -- a large non-track-like "shower" cluster's own shape is
      // exactly what we want to show) for this muon's home slice, if it's
      // ND-LAr-fiducial-origin (the actual population of interest).
      if (dump_display && sp.lar_fiducial_start[pidx]) {
        for (int i = 0; i < n_space_points; ++i) {
          const int cid = point_cluster_id[i];
          const bool is_track_like_pt = cid > 0 && clusters[cid - 1].IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack);
          const double lin = cid > 0 ? clusters[cid - 1].GetLinearity() : 0.0;
          display_points_csv << entry << "," << slice_no << "," << label.vgid << "," << label.trackid << ","
                              << sp_x[i] << "," << sp_y[i] << "," << sp_z[i] << "," << cid << "," << lin << ","
                              << (is_track_like_pt ? 1 : 0) << "," << (point_label[i] == label ? 1 : 0) << "\n";
        }
        for (size_t c = 0; c < clusters.size(); ++c) {
          const auto &cl = clusters[c];
          if (cl.GetSize() < 3) continue;  // matches TMS_SpacePointCluster's own PCA validity floor
          const bool is_track_like_c = cl.IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack);
          const auto &centroid = cl.GetCentroid();
          const auto &eigenvalues = cl.GetEigenvalues();
          const auto &eigenvectors = cl.GetEigenvectors();
          display_pca_csv << entry << "," << slice_no << "," << label.vgid << "," << label.trackid << "," << (c + 1)
                           << "," << cl.GetSize() << "," << (is_track_like_c ? 1 : 0) << "," << centroid[0] << ","
                           << centroid[1] << "," << centroid[2];
          for (int rank = 0; rank < 3; ++rank) {
            display_pca_csv << "," << eigenvalues[rank] << "," << eigenvectors[rank][0] << ","
                             << eigenvectors[rank][1] << "," << eigenvectors[rank][2];
          }
          display_pca_csv << "\n";
        }
      }
    }

    if (entry % 200 == 0) {
      std::cout << "  entry=" << entry << "/" << n_entries << " spill=" << spill_no << " slice=" << slice_no
                << " nSP=" << n_space_points << " muons_in_slice=" << n_muons_in_slice << std::endl;
    }
  }

  muons_csv.close();
  clusters_csv.close();
  if (dump_display) {
    display_points_csv.close();
    display_pca_csv.close();
    std::cout << "Wrote " << display_prefix << "_points.csv and " << display_prefix << "_pca.csv" << std::endl;
  }

  std::cout << "Done. " << n_muons_total << " muon candidates, " << n_muons_found << " found ("
            << (n_muons_total > 0 ? 100.0 * n_muons_found / n_muons_total : 0.0) << "%)" << std::endl;
  std::cout << n_clusters_total << " track-like clusters, " << n_clusters_muon_matched << " muon-matched ("
            << (n_clusters_total > 0 ? 100.0 * n_clusters_muon_matched / n_clusters_total : 0.0) << "%)"
            << std::endl;
  std::cout << "X/Y hit label agreement (raw TrackId): " << xy_agree_count_raw << " agree, "
            << xy_mismatch_count_raw << " mismatch" << std::endl;
  std::cout << "X/Y hit label agreement (after parent-chain collapse): " << xy_agree_count << " agree, "
            << xy_mismatch_count << " mismatch" << std::endl;
  std::cout << "  of the mismatches: " << xy_mismatch_diff_vertex << " are different vertices entirely (genuine "
            << "ghost pairing across unrelated interactions), " << xy_mismatch_same_vertex
            << " are the same vertex but a different collapsed track (correlated overlap within one interaction)"
            << std::endl;
  std::cout << "Slices skipped (nTrueParticles mismatch between Truth_Info/Truth_Spill): "
            << n_slices_skipped_mismatch << std::endl;
  std::cout << "Wrote " << muons_csv_path << " and " << clusters_csv_path << std::endl;

  return 0;
}
