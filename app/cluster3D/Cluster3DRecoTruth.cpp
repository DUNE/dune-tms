// Object-first, hit-level truth validation of the production Cluster3D
// reconstruction (TMS_Cluster3DReco::Run): every track it makes, and every
// true muon, scored by the energy shares of the hits the tracks used.
//
// The reconstruction under test is exactly the library stage conversion
// runs, given the slice's hits (hit-level fit, orphan pickup) -- no truth
// enters it. Truth only scores afterwards:
//   - a track's owner is the particle with the largest summed energy share
//     over the track's hits, if that share is more than half the hits;
//     otherwise the track is "mixed" (junk);
//   - a true muon (|PDG| 13 with >= 5 true hits in the slice, the muon-first
//     tools' population) is found if it owns at least one track; its best
//     track (most of its energy share) gives completeness and purity; more
//     than one owned track is a duplicate.
//
// Stage 2 (graph search in non-track-like clusters) is off by default, as in
// the library; CLUSTER3D_GRAPH=1 turns it on. Requires reco files with the per-hit
// table and per-hit energy shares (converted 2026-09-25 evening or later).
//
// Each muon row also says how its space points are spread over the DBSCAN
// clusters, and -- for a muon whose best track ends early -- where the rest of
// it went (its "tail": hits/points more than kTailMarginMM downstream of the
// best track's last hit): the track's own cluster, other clusters, DBSCAN
// noise, or another track. This is the ceiling study for linking clusters
// before the fit (2026-09-26).
//
// Usage: Cluster3DRecoTruth <edep_sim_geom_file> <reco.root> <tracks.csv> <muons.csv>

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

#include "TMS_Cluster3DReco.h"
#include "TMS_FieldModel.h"
#include "TMS_Geom.h"
#include "TMS_SpacePoint.h"
#include "SpacePointLayerInput.h"
#include "TruthLabels.h"

namespace {

const int kMaxSpacePoints = 10000;    // __TMS_MAX_SPACEPOINTS__
const int kMaxHits = 20000;           // __TMS_MAX_HITS__
const int kMaxTrueParticles = 20000;
// A muon hit/point counts as beyond its best track's end only when it is at
// least this far downstream of the track's last hit (so hits in the same
// plane or plane pair are not counted).
const double kTailMarginMM = 20.0;

struct SpillParticles {
  int n = 0;
  std::vector<long long> vgid;
  std::vector<int> trackid, pdg, parent_trackid;
  std::vector<bool> tms_fiducial_start, lar_fiducial_start, tms_fiducial_end;
  std::vector<float> momentum;  // MomentumTMSStart, 4 per particle
  // Starts inside the ND-LAr box the validation suite's Reco_Eff uses for its
  // ND-physics sample (is_lar_start_position) -- looser than LArFiducialStart.
  std::vector<bool> lar_box_start;
  std::unordered_map<TrueLabel, int, LabelHash> index_of;
  std::vector<int> collapsed_trackid;
};

// A particle's own secondaries are folded into the nearest muon ancestor, or
// the top primary -- same convention as the other Cluster3D truth tools.
int CollapseTrackId(const SpillParticles &sp, int start_idx) {
  int idx = start_idx;
  int fallback = sp.trackid[start_idx];
  int guard = 0;
  while (idx >= 0 && guard++ < 10000) {
    if (std::abs(sp.pdg[idx]) == 13) return sp.trackid[idx];
    fallback = sp.trackid[idx];
    const int parent = sp.parent_trackid[idx];
    if (parent < 0) break;
    auto it = sp.index_of.find(TrueLabel{sp.vgid[idx], parent});
    if (it == sp.index_of.end()) break;
    idx = it->second;
  }
  return fallback;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0] << " <edep_sim_geom_file> <reco.root> <tracks.csv> <muons.csv>" << std::endl;
    return -1;
  }
  const std::string geom_filename = argv[1], input_filename = argv[2];
  const unsigned int min_true_hits = 5;  // muon population, as the muon-first tools

  TFile geom_input(geom_filename.c_str());
  TGeoManager *geom = geom_input.IsZombie() ? nullptr : (TGeoManager *)geom_input.Get("EDepSimGeometry");
  if (!geom) {
    std::cerr << "No EDepSimGeometry in " << geom_filename << std::endl;
    return -1;
  }
  TMS_Geom::GetInstance().SetGeometry(geom);
  const double bar_pitch = TMS_Geom::GetInstance().GetMaxBarPitch();

  TMS_Cluster3DReco::Config config;
  if (const char *v = std::getenv("CLUSTER3D_GRAPH")) config.UseGraphSearch = std::atoi(v) != 0;
  if (const char *v = std::getenv("CLUSTER3D_GRAPH_MIN_LAYERS")) config.MinGraphPathLayers = std::atoi(v);
  if (const char *v = std::getenv("CLUSTER3D_GRAPH_MIN_HITS")) config.MinGraphTrackHits = std::atoi(v);
  if (const char *v = std::getenv("CLUSTER3D_MIN_HITS")) config.Split.MinHitsPerTrack = std::atoi(v);
  if (const char *v = std::getenv("CLUSTER3D_GRAPH_MIN_CLUSTER")) config.MinGraphClusterSize = std::atoi(v);
  const RegionFieldModel field;

  TFile input(input_filename.c_str());
  TTree *reco_tree = (TTree *)input.Get("Reco_Tree");
  TTree *truth_info = (TTree *)input.Get("Truth_Info");
  TTree *truth_spill = (TTree *)input.Get("Truth_Spill");
  if (!reco_tree || !truth_info || !truth_spill || !reco_tree->GetBranch("SpacePointHitTrueEnergyFrac")) {
    std::cerr << input_filename << ": needs Reco_Tree/Truth_Info/Truth_Spill with per-hit energy shares" << std::endl;
    return -1;
  }

  // --- Truth_Spill: every particle, per spill. ---
  int spill_no_ts = 0, n_tp_ts = 0;
  static std::vector<long long> vgid_ts(kMaxTrueParticles);
  static std::vector<int> trackid_ts(kMaxTrueParticles), pdg_ts(kMaxTrueParticles), parent_ts(kMaxTrueParticles);
  static std::vector<float> mom_ts(kMaxTrueParticles * 4), birth_ts(kMaxTrueParticles * 4);
  static bool fid_start_ts[kMaxTrueParticles], lar_start_ts[kMaxTrueParticles], fid_end_ts[kMaxTrueParticles];
  truth_spill->SetBranchAddress("SpillNo", &spill_no_ts);
  truth_spill->SetBranchAddress("nTrueParticles", &n_tp_ts);
  truth_spill->SetBranchAddress("VertexGlobalID", vgid_ts.data());
  truth_spill->SetBranchAddress("TrackId", trackid_ts.data());
  truth_spill->SetBranchAddress("PDG", pdg_ts.data());
  truth_spill->SetBranchAddress("Parent", parent_ts.data());
  truth_spill->SetBranchAddress("MomentumTMSStart", mom_ts.data());
  truth_spill->SetBranchAddress("BirthPosition", birth_ts.data());
  truth_spill->SetBranchAddress("TMSFiducialStart", fid_start_ts);
  truth_spill->SetBranchAddress("LArFiducialStart", lar_start_ts);
  truth_spill->SetBranchAddress("TMSFiducialEnd", fid_end_ts);
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
    sp.lar_box_start.resize(n_tp_ts);
    for (int k = 0; k < n_tp_ts; ++k) {
      const float *b = &birth_ts[k * 4];
      sp.lar_box_start[k] = b[0] >= -3478.48 && b[0] <= 3478.48 && b[1] >= -2166.71 && b[1] <= 829.282 &&
                            b[2] >= 4179.24 && b[2] <= 9135.88;
    }
    sp.tms_fiducial_start.assign(fid_start_ts, fid_start_ts + n_tp_ts);
    sp.lar_fiducial_start.assign(lar_start_ts, lar_start_ts + n_tp_ts);
    sp.tms_fiducial_end.assign(fid_end_ts, fid_end_ts + n_tp_ts);
    for (int i = 0; i < n_tp_ts; ++i) sp.index_of[{sp.vgid[i], sp.trackid[i]}] = i;
    sp.collapsed_trackid.resize(n_tp_ts);
    for (int i = 0; i < n_tp_ts; ++i) sp.collapsed_trackid[i] = CollapseTrackId(sp, i);
    spills[spill_no_ts] = std::move(sp);
  }

  // --- Reco_Tree: space points and the per-hit table. ---
  int n_sp = 0, n_hits = 0, spill_no = 0, slice_no = 0;
  static std::vector<float> sp_x(kMaxSpacePoints), sp_y(kMaxSpacePoints), sp_z(kMaxSpacePoints), sp_t(kMaxSpacePoints);
  static std::vector<int> sp_xi(kMaxSpacePoints), sp_yi(kMaxSpacePoints);
  static std::vector<float> h_t(kMaxHits), h_nz(kMaxHits), h_z(kMaxHits), h_f1(kMaxHits), h_f2(kMaxHits);
  static std::vector<int> h_view(kMaxHits), h_ped(kMaxHits), h_tk1(kMaxHits), h_tk2(kMaxHits);
  static std::vector<long long> h_vg1(kMaxHits), h_vg2(kMaxHits);
  reco_tree->SetBranchAddress("nSpacePoints", &n_sp);
  reco_tree->SetBranchAddress("SpacePointX", sp_x.data());
  reco_tree->SetBranchAddress("SpacePointY", sp_y.data());
  reco_tree->SetBranchAddress("SpacePointZ", sp_z.data());
  const SpacePointLayerInput sp_layer(reco_tree, kMaxSpacePoints);
  reco_tree->SetBranchAddress("SpacePointTime", sp_t.data());
  reco_tree->SetBranchAddress("SpacePointXHitIndex", sp_xi.data());
  reco_tree->SetBranchAddress("SpacePointYHitIndex", sp_yi.data());
  reco_tree->SetBranchAddress("nSpacePointHits", &n_hits);
  reco_tree->SetBranchAddress("SpacePointHitTime", h_t.data());
  reco_tree->SetBranchAddress("SpacePointHitNotZ", h_nz.data());
  reco_tree->SetBranchAddress("SpacePointHitZ", h_z.data());
  reco_tree->SetBranchAddress("SpacePointHitView", h_view.data());
  reco_tree->SetBranchAddress("SpacePointHitPedSup", h_ped.data());
  reco_tree->SetBranchAddress("SpacePointHitTrueVertexGlobalId", h_vg1.data());
  reco_tree->SetBranchAddress("SpacePointHitTrueTrackId", h_tk1.data());
  reco_tree->SetBranchAddress("SpacePointHitTrueEnergyFrac", h_f1.data());
  reco_tree->SetBranchAddress("SpacePointHitTrue2VertexGlobalId", h_vg2.data());
  reco_tree->SetBranchAddress("SpacePointHitTrue2TrackId", h_tk2.data());
  reco_tree->SetBranchAddress("SpacePointHitTrue2EnergyFrac", h_f2.data());
  reco_tree->SetBranchAddress("SpillNo", &spill_no);
  reco_tree->SetBranchAddress("SliceNo", &slice_no);
  int n_tp_ti = 0;
  static std::vector<int> true_nhits_slice(kMaxTrueParticles);
  truth_info->SetBranchAddress("nTrueParticles", &n_tp_ti);
  truth_info->SetBranchAddress("TrueNHitsInSlice", true_nhits_slice.data());

  std::ofstream tracks_csv(argv[3]), muons_csv(argv[4]);
  tracks_csv << "sourcefile,entry,slice,track,stage,cluster_size,iteration,hits_used,orphans,converged,"
                "owner_vgid,owner_trackid,owner_pdg,owner_is_muon,purity,duplicate,start_z,start_momentum_mev,"
                "start_charge\n";
  muons_csv << "sourcefile,entry,slice,vertexglobalid,trackid,vertex_in_tms,vertex_in_lar_fiducial,stops_in_tms,"
               "true_hits_in_slice,true_momentum_tms_mev,true_charge,tracks_owned,found,best_stage,"
               "hit_completeness_pct,hit_purity_pct,best_start_momentum_mev,best_start_charge,"
               "captor_share_pct,captor_owner_vgid,captor_owner_trackid,captor_stage,"
               "n_points,n_noise_points,n_clusters,n_clusters_2pt,n_tracklike_clusters,frac_in_largest_cluster,"
               "largest_cluster_muon_purity,true_last_hit_z,best_last_hit_z,best_cluster,tail_hits,"
               "tail_hits_other_tracks,tail_points,tail_points_best_cluster,tail_points_other_clusters,"
               "tail_other_clusters,tail_points_noise,vertex_in_lar_box\n";

  long n_tracks = 0, n_stage2 = 0;
  for (Long64_t entry = 0; entry < reco_tree->GetEntries(); ++entry) {
    reco_tree->GetEntry(entry);
    truth_info->GetEntry(entry);
    auto spill_it = spills.find(spill_no);
    if (spill_it == spills.end() || n_sp <= 0 || n_sp >= kMaxSpacePoints) continue;
    const SpillParticles &sp = spill_it->second;
    if (sp.n != n_tp_ti) continue;
    auto collapse = [&](const TrueLabel &raw) -> TrueLabel {
      if (!raw.Valid()) return raw;
      auto it = sp.index_of.find(raw);
      if (it == sp.index_of.end()) return raw;
      return TrueLabel{raw.vgid, sp.collapsed_trackid[it->second]};
    };

    // The slice's hits: fit measurements, truth, usability.
    std::vector<TMS_KalmanFollower::FitHit> hits(n_hits);
    std::vector<HitTruth> hit_truth(n_hits);
    std::vector<char> usable(n_hits, 0);
    for (int h = 0; h < n_hits; ++h) {
      hits[h].Z = h_z[h];
      hits[h].Coordinate = h_nz[h];
      hits[h].MeasuresX = h_view[h] == 1;
      hits[h].SigmaMM = bar_pitch / std::sqrt(12.0);
      hits[h].Time = h_t[h];
      usable[h] = !h_ped[h] && (h_view[h] == 0 || h_view[h] == 1);
      hits[h].Usable = usable[h] != 0;
      HitTruth &ht = hit_truth[h];
      ht.first = collapse(TrueLabel{h_vg1[h], h_tk1[h]});
      ht.first_frac = h_f1[h];
      ht.second = collapse(TrueLabel{h_vg2[h], h_tk2[h]});
      ht.second_frac = h_f2[h];
      if (ht.second.Valid() && ht.second == ht.first) {
        ht.first_frac += ht.second_frac;
        ht.second = TrueLabel();
        ht.second_frac = 0.0;
      }
    }
    std::vector<TMS_SpacePoint> points;
    points.reserve(n_sp);
    for (int i = 0; i < n_sp; ++i)
      points.emplace_back(sp_x[i], sp_y[i], sp_z[i], sp_xi[i], sp_yi[i], sp_t[i], sp_layer.Layer(i, sp_z[i]));

    // --- Reconstruction (no truth). ---
    TMS_Cluster3DReco::RunInfo run_info;
    const std::vector<TMS_Cluster3DReco::Track> tracks = TMS_Cluster3DReco::Run(points, hits, config, field, &run_info);

    // Which DBSCAN cluster each point is in (-1 = noise), each point's true
    // owner, and which track (if any) used each hit.
    std::vector<int> point_cluster(n_sp, -1);
    for (std::size_t c = 0; c < run_info.Clusters.size(); ++c)
      for (int idx : run_info.Clusters[c]) point_cluster[idx] = static_cast<int>(c);
    std::vector<TrueLabel> point_owner(n_sp);
    for (int i = 0; i < n_sp; ++i)
      if (sp_xi[i] >= 0 && sp_yi[i] >= 0 && sp_xi[i] < n_hits && sp_yi[i] < n_hits)
        point_owner[i] = PointOwner(hit_truth[sp_xi[i]], hit_truth[sp_yi[i]]);
    std::vector<int> hit_track(n_hits, -1);
    for (std::size_t t = 0; t < tracks.size(); ++t)
      for (int h : tracks[t].HitIndices) hit_track[h] = static_cast<int>(t);

    // --- Scoring: each track's owner by hit energy share. ---
    struct Owned {
      TrueLabel owner;
      double owner_share = 0.0;
      int used = 0;
      const TMS_Cluster3DReco::Track *track = nullptr;
    };
    std::vector<Owned> owned;
    std::set<TrueLabel, bool (*)(const TrueLabel &, const TrueLabel &)> owners_seen(
        [](const TrueLabel &a, const TrueLabel &b) { return a.vgid != b.vgid ? a.vgid < b.vgid : a.trackid < b.trackid; });
    for (std::size_t t = 0; t < tracks.size(); ++t) {
      const TMS_Cluster3DReco::Track &track = tracks[t];
      std::unordered_map<TrueLabel, double, LabelHash> share;
      for (int h : track.HitIndices) {
        const HitTruth &ht = hit_truth[h];
        if (ht.first.Valid()) share[ht.first] += ht.first_frac;
        if (ht.second.Valid()) share[ht.second] += ht.second_frac;
      }
      Owned o;
      o.used = static_cast<int>(track.HitIndices.size());
      o.track = &track;
      for (const auto &kv : share)
        if (kv.second > o.owner_share) {
          o.owner = kv.first;
          o.owner_share = kv.second;
        }
      if (!(o.owner_share > 0.5 * o.used)) o.owner = TrueLabel();  // mixed: nobody's
      int owner_pdg = 0;
      if (o.owner.Valid()) {
        auto it = sp.index_of.find(o.owner);
        if (it != sp.index_of.end()) owner_pdg = sp.pdg[it->second];
      }
      const bool duplicate = o.owner.Valid() && owners_seen.count(o.owner) > 0;
      if (o.owner.Valid()) owners_seen.insert(o.owner);
      const TMS_KalmanFollower::FitResult &fit = track.Fit;
      tracks_csv << input_filename << "," << entry << "," << slice_no << "," << t << "," << track.Stage << ","
                 << track.ClusterSize << "," << track.Iteration << "," << o.used << "," << fit.Orphans.size() << ","
                 << (fit.Converged ? 1 : 0) << "," << o.owner.vgid << "," << o.owner.trackid << "," << owner_pdg << ","
                 << (std::abs(owner_pdg) == 13 ? 1 : 0) << ","
                 << (o.used > 0 ? (o.owner.Valid() ? o.owner_share : 0.0) / o.used : 0.0) << ","
                 << (duplicate ? 1 : 0) << "," << (fit.HasStartState ? fit.StartZ : 0.0) << ","
                 << (fit.HasStartState ? fit.StartMomentumMeV : 0.0) << ","
                 << (fit.HasStartState ? fit.StartCharge : 0.0) << "\n";
      owned.push_back(o);
      ++n_tracks;
      if (track.Stage == 2) ++n_stage2;
    }

    // --- Muons: found if they own a track; the best one by energy share. ---
    for (int i = 0; i < sp.n; ++i) {
      if (std::abs(sp.pdg[i]) != 13 || true_nhits_slice[i] < (int)min_true_hits) continue;
      const TrueLabel label{sp.vgid[i], sp.trackid[i]};
      double total_share = 0.0;
      for (int h = 0; h < n_hits; ++h)
        if (usable[h]) total_share += hit_truth[h].Share(label);
      int n_owned = 0;
      const Owned *best = nullptr;
      for (const Owned &o : owned)
        if (o.owner.Valid() && o.owner == label) {
          ++n_owned;
          if (!best || o.owner_share > best->owner_share) best = &o;
        }
      // Diagnostic: the track holding the largest part of this muon's hit
      // energy, whoever owns it -- where a muon's hits went when it owns no track.
      double captor_share = 0.0;
      const Owned *captor = nullptr;
      for (const Owned &o : owned) {
        double s = 0.0;
        for (int h : o.track->HitIndices) s += hit_truth[h].Share(label);
        if (s > captor_share) {
          captor_share = s;
          captor = &o;
        }
      }
      // --- How this muon's points are spread over the DBSCAN clusters. ---
      std::map<int, int> muon_points_in_cluster;  // cluster -> this muon's points in it
      int n_points = 0, n_noise_points = 0;
      for (int k = 0; k < n_sp; ++k) {
        if (!(point_owner[k] == label)) continue;
        ++n_points;
        if (point_cluster[k] < 0) ++n_noise_points;
        else ++muon_points_in_cluster[point_cluster[k]];
      }
      int n_clusters_2pt = 0, n_tracklike_clusters = 0, largest_cluster = -1, largest_count = 0;
      for (const auto &kv : muon_points_in_cluster) {
        if (kv.second >= 2) ++n_clusters_2pt;
        if (run_info.ClusterTrackLike[kv.first]) ++n_tracklike_clusters;
        if (kv.second > largest_count) {
          largest_count = kv.second;
          largest_cluster = kv.first;
        }
      }
      const double largest_purity =
          largest_cluster >= 0 ? static_cast<double>(largest_count) / run_info.Clusters[largest_cluster].size() : 0.0;

      // --- Where the muon's tail went, beyond its best track's last hit. ---
      // (With no best track, the whole muon is tail.)
      double true_last_z = -1e9, best_last_z = -1e9;
      for (int h = 0; h < n_hits; ++h)
        if (usable[h] && hit_truth[h].Share(label) > 0.5) true_last_z = std::max(true_last_z, (double)h_z[h]);
      const int best_track = best ? static_cast<int>(best->track - tracks.data()) : -1;
      const int best_cluster = best ? best->track->ClusterIndex : -1;
      if (best)
        for (int h : best->track->HitIndices) best_last_z = std::max(best_last_z, (double)h_z[h]);
      const double tail_from = best ? best_last_z + kTailMarginMM : -1e9;
      int tail_hits = 0, tail_hits_other_tracks = 0;
      for (int h = 0; h < n_hits; ++h) {
        if (!usable[h] || !(hit_truth[h].Share(label) > 0.5) || h_z[h] <= tail_from) continue;
        ++tail_hits;
        if (hit_track[h] >= 0 && hit_track[h] != best_track) ++tail_hits_other_tracks;
      }
      int tail_points = 0, tail_points_best_cluster = 0, tail_points_other_clusters = 0, tail_points_noise = 0;
      std::set<int> tail_other_clusters;
      for (int k = 0; k < n_sp; ++k) {
        if (!(point_owner[k] == label) || sp_z[k] <= tail_from) continue;
        ++tail_points;
        const int c = point_cluster[k];
        if (c < 0) ++tail_points_noise;
        else if (c == best_cluster) ++tail_points_best_cluster;
        else {
          ++tail_points_other_clusters;
          tail_other_clusters.insert(c);
        }
      }

      const double p = std::sqrt(sp.momentum[i * 4] * sp.momentum[i * 4] + sp.momentum[i * 4 + 1] * sp.momentum[i * 4 + 1] +
                                 sp.momentum[i * 4 + 2] * sp.momentum[i * 4 + 2]);
      muons_csv << input_filename << "," << entry << "," << slice_no << "," << label.vgid << "," << label.trackid << ","
                << (sp.tms_fiducial_start[i] ? 1 : 0) << "," << (sp.lar_fiducial_start[i] ? 1 : 0) << ","
                << (sp.tms_fiducial_end[i] ? 1 : 0) << "," << true_nhits_slice[i] << "," << p << ","
                << (sp.pdg[i] > 0 ? -1 : 1) << "," << n_owned << "," << (best ? 1 : 0) << ","
                << (best ? best->track->Stage : 0) << ","
                << (best && total_share > 0 ? 100.0 * best->owner_share / total_share : 0.0) << ","
                << (best && best->used > 0 ? 100.0 * best->owner_share / best->used : 0.0) << ","
                << (best && best->track->Fit.HasStartState ? best->track->Fit.StartMomentumMeV : 0.0) << ","
                << (best && best->track->Fit.HasStartState ? best->track->Fit.StartCharge : 0.0) << ","
                << (captor && total_share > 0 ? 100.0 * captor_share / total_share : 0.0) << ","
                << (captor ? captor->owner.vgid : -1) << "," << (captor ? captor->owner.trackid : -999) << ","
                << (captor ? captor->track->Stage : 0) << "," << n_points << "," << n_noise_points << ","
                << muon_points_in_cluster.size() << "," << n_clusters_2pt << "," << n_tracklike_clusters << ","
                << (n_points > 0 ? static_cast<double>(largest_count) / n_points : 0.0) << "," << largest_purity << ","
                << true_last_z << "," << (best ? best_last_z : 0.0) << "," << best_cluster << "," << tail_hits << ","
                << tail_hits_other_tracks << "," << tail_points << "," << tail_points_best_cluster << ","
                << tail_points_other_clusters << "," << tail_other_clusters.size() << "," << tail_points_noise << ","
                << (sp.lar_box_start[i] ? 1 : 0) << "\n";
    }
  }
  std::cout << "Tracks: " << n_tracks << " (Stage 2: " << n_stage2 << ")" << std::endl;
  return 0;
}
