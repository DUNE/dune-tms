// Validates TMS_GraphTrackFinder against a single real, already-known-hard slice:
// by default, entry 101 of the 2026-09-07 clustering-benchmark reference
// file -- the shower-contaminated muon that originally motivated the
// Link-and-Tree proposal. DBSCAN puts all 3,820 of that slice's space points
// into one blob with PCA linearity 0.159 (well below the 0.8 track-like
// cutoff), so the whole cluster is discarded and the muon inside it is never
// recovered, even though a genuine 1.49 GeV muon runs cleanly through it.
//
// Reuses ClusterTruthEfficiency's truth machinery (Truth_Spill loading +
// Parent-chain collapse, so a muon's own delta-ray hits count as the muon's,
// not "a different particle") so "muon-owned" is defined exactly the same
// way here as in the DBSCAN benchmark -- results are directly comparable.
// Needs only Reco_Tree + Truth_Spill from a RecoCandidates file (no geometry
// file: unlike DBSCAN, Graph Track Finder doesn't need a clustering tolerance
// derived from bar/plane pitch).

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
#include "TTree.h"

#include "TMS_GraphTrackFinder.h"
#include "TMS_LayerGrouping.h"
#include "SpacePointLayerInput.h"
#include "TMS_SpacePoint.h"

namespace {

const int kMaxSpacePoints = 10000;    // matches __TMS_MAX_SPACEPOINTS__
const int kMaxTrueParticles = 20000;  // matches __TMS_MAX_TRUE_PARTICLES__

// A (vertex, track) pair identifying one true particle. Copied from
// ClusterTruthEfficiency.cpp so "muon-owned" means the same thing in both
// tools.
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

// One spill's worth of true particles, indexed for fast (vgid,trackid)
// lookup and pre-collapsed (see CollapseTrackId below).
struct SpillParticles {
  int n = 0;
  std::vector<long long> vgid;
  std::vector<int> trackid;
  std::vector<int> pdg;
  std::vector<int> parent_trackid;  // GetParent(): parent's TrackId within the same vertex, or -1
  std::unordered_map<TrueLabel, int, LabelHash> index_of;  // (vgid,trackid) -> particle index
  std::vector<int> collapsed_trackid;  // see CollapseTrackId()
};

// A raw G4 TrackId is per-hit-segment bookkeeping, not collapsed to any
// ancestor -- a muon's own delta-ray gets a fresh TrackId the instant it's
// created, even though physically it's just a continuation of the muon's
// own ionization trail. Walk the Parent chain and collapse to the first
// muon ancestor found (or the topmost primary if none), so a muon's own
// delta-ray hits count as "the muon" rather than "a different particle".
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

struct CandidateBreakdown {
  std::size_t index = 0;
  int layer = 0;
  bool target = false;
  double dxdz = 0.0;
  double dydz = 0.0;
  double cost = 0.0;
  int x_hit = -1;
  int y_hit = -1;
};

struct CandidateInspection {
  std::vector<CandidateBreakdown> candidates;
  std::size_t next_layer_total = 0;
  std::size_t next_layer_pass = 0;
};

CandidateInspection InspectNextLayerCandidates(
    const std::vector<TMS_SpacePoint> &points,
    const std::vector<TrueLabel> &labels,
    const std::vector<int> &z_layer,
    const TMS_GraphTrackFinder::Path &path,
    double max_abs_dxdz,
    double max_abs_dydz,
    double point_reward,
    double gap_penalty,
    double slope_penalty,
    double occupancy_penalty,
    double multiplicity_penalty) {
  CandidateInspection inspection;
  if (path.SpacePointIndices.empty()) return inspection;

  const std::size_t endpoint = path.SpacePointIndices.back();
  const int endpoint_layer = z_layer[endpoint];
  int next_layer = std::numeric_limits<int>::max();
  for (int layer : z_layer) {
    if (layer > endpoint_layer && layer < next_layer) next_layer = layer;
  }
  if (next_layer == std::numeric_limits<int>::max()) return inspection;

  std::unordered_map<int, std::size_t> x_multiplicity;
  std::unordered_map<int, std::size_t> y_multiplicity;
  for (const TMS_SpacePoint &point : points) {
    if (point.GetXHitIndex() >= 0) ++x_multiplicity[point.GetXHitIndex()];
    if (point.GetYHitIndex() >= 0) ++y_multiplicity[point.GetYHitIndex()];
  }

  const double endpoint_z = points[endpoint].GetZ();
  const double endpoint_time = points[endpoint].GetTime();
  std::size_t layer_size = 0;
  for (int layer : z_layer) if (layer == next_layer) ++layer_size;
  inspection.next_layer_total = layer_size;
  const double occupancy = layer_size > 0 ? static_cast<double>(layer_size - 1) : 0.0;

  for (std::size_t idx = 0; idx < points.size(); ++idx) {
    if (z_layer[idx] != next_layer) continue;
    const TMS_SpacePoint &candidate = points[idx];
    const double dz = candidate.GetZ() - endpoint_z;
    if (dz <= 0.0) continue;
    const double dxdz = (candidate.GetX() - points[endpoint].GetX()) / dz;
    const double dydz = (candidate.GetY() - points[endpoint].GetY()) / dz;
    if (std::abs(dxdz) > max_abs_dxdz || std::abs(dydz) > max_abs_dydz) continue;
    if (std::abs(candidate.GetTime() - endpoint_time) > 40.0) continue;
    ++inspection.next_layer_pass;
    const std::size_t mult_x = candidate.GetXHitIndex() >= 0 ? x_multiplicity[candidate.GetXHitIndex()] : 0;
    const std::size_t mult_y = candidate.GetYHitIndex() >= 0 ? y_multiplicity[candidate.GetYHitIndex()] : 0;
    const double multiplicity = static_cast<double>(std::max(mult_x, mult_y) > 0 ? std::max(mult_x, mult_y) - 1 : 0);
    CandidateBreakdown item;
    item.index = idx;
    item.layer = next_layer;
    item.target = labels[idx] == labels[endpoint];
    item.dxdz = dxdz;
    item.dydz = dydz;
    item.cost = -point_reward + gap_penalty * 0.0 + slope_penalty * (dxdz * dxdz + dydz * dydz) +
                occupancy_penalty * occupancy + multiplicity_penalty * multiplicity;
    item.x_hit = candidate.GetXHitIndex();
    item.y_hit = candidate.GetYHitIndex();
    inspection.candidates.push_back(item);
  }

  std::sort(inspection.candidates.begin(), inspection.candidates.end(), [](const CandidateBreakdown &a, const CandidateBreakdown &b) {
    if (a.cost != b.cost) return a.cost < b.cost;
    return a.index < b.index;
  });
  return inspection;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2 || argc > 14) {
    std::cerr << "Usage: " << argv[0]
              << " <input_reco_tree.root> [target_vertexglobalid=10000461] [target_trackid=0]"
                 " [occupancy_penalty] [hit_multiplicity_penalty] [max_seed_layer_occupancy]"
                 " [max_seed_hit_multiplicity] [beam_width] [max_links_per_target_layer]"
                 " [use_curvature_projection 0|1] [position_quantization_x] [position_quantization_y]"
                 " [output_json_path]"
              << std::endl;
    return 1;
  }
  const std::string input_filename = argv[1];
  const long long target_vgid = argc >= 3 ? std::stoll(argv[2]) : 10000461;
  const int target_trackid = argc >= 4 ? std::stoi(argv[3]) : 0;
  const TrueLabel target{target_vgid, target_trackid};
  // Optional overrides for the growth-phase costs and seed-gate thresholds
  // that scale with local occupancy -- left as command-line knobs (rather
  // than baked-in changes to TMS_GraphTrackFinder.h's Config defaults) since
  // retuning them for real TMS occupancy scales is exactly the open
  // question this tool exists to probe.
  const bool has_occ_override = argc >= 5;
  const bool has_mult_override = argc >= 6;
  const bool has_seed_occ_override = argc >= 7;
  const bool has_seed_mult_override = argc >= 8;
  const bool has_beam_width_override = argc >= 9;
  const bool has_max_links_override = argc >= 10;
  const double occ_override = has_occ_override ? std::stod(argv[4]) : 0.0;
  const double mult_override = has_mult_override ? std::stod(argv[5]) : 0.0;
  const std::size_t seed_occ_override = has_seed_occ_override
      ? static_cast<std::size_t>(std::stoul(argv[6])) : 0;
  const std::size_t seed_mult_override = has_seed_mult_override
      ? static_cast<std::size_t>(std::stoul(argv[7])) : 0;
  const std::size_t beam_width_override = has_beam_width_override
      ? static_cast<std::size_t>(std::stoul(argv[8])) : 0;
  const std::size_t max_links_override = has_max_links_override
      ? static_cast<std::size_t>(std::stoul(argv[9])) : 0;
  const bool has_curvature_override = argc >= 11;
  const bool curvature_override = has_curvature_override ? (std::stoi(argv[10]) != 0) : true;
  const bool has_quant_x_override = argc >= 12;
  const bool has_quant_y_override = argc >= 13;
  const double quant_x_override = has_quant_x_override ? std::stod(argv[11]) : 0.0;
  const double quant_y_override = has_quant_y_override ? std::stod(argv[12]) : 0.0;
  const std::string output_json_path = argc >= 14 ? argv[13] : "";

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
  Long64_t best_entry = -1;
  int best_count = 0;
  int best_total = 0;
  int best_spill = 0, best_slice = 0;
  // Cached copies of the winning entry's arrays -- reco_tree->GetEntry()
  // overwrites the buffers above on every call, so the entry we care about
  // has to be saved off as soon as we recognize it's the best one so far.
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
      // Primary truth label per space point: X-hit link, falling back to Y
      // -- same convention as ClusterTruthEfficiency, so results compare.
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
                                              sp_x_hitidx[i], sp_y_hitidx[i], sp_time[i],
                                              sp_layer.Layer(i, sp_z[i])));
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

  // Default config, same as the synthetic GraphTrackFinder_test -- real occupancy
  // is much higher (3800ish points vs. 330), so watch the diagnostics below
  // for ResourceLimitReached before trusting the result.
  TMS_GraphTrackFinder::Config config;
  if (has_occ_override) config.OccupancyPenalty = occ_override;
  if (has_mult_override) config.HitMultiplicityPenalty = mult_override;
  if (has_seed_occ_override) config.MaxSeedLayerOccupancy = seed_occ_override;
  if (has_seed_mult_override) config.MaxSeedHitMultiplicity = seed_mult_override;
  if (has_beam_width_override) config.BeamWidth = beam_width_override;
  if (has_max_links_override) config.MaxLinksPerTargetLayer = max_links_override;
  if (has_curvature_override) config.UseCurvatureProjection = curvature_override;
  if (has_quant_x_override) config.PositionQuantizationX = quant_x_override;
  if (has_quant_y_override) config.PositionQuantizationY = quant_y_override;
  const TMS_GraphTrackFinder::Result result = TMS_GraphTrackFinder::Finder(config).Find(best_points);

  // Distinct-plane bookkeeping (Section 15's real validation metric): a
  // combinatorial X/Y ghost can multiply how many *space points* carry the
  // target's label without the target ever visiting more than one real
  // plane per crossing, so raw point-count "completeness" overstates how
  // much of the muon's actual trajectory is missing.
  const std::vector<int> z_layer = TMS_LayerGrouping::GroupIndexOfEachPoint(best_points, config.LayerZTolerance);
  std::set<int> target_layers_in_slice;
  for (std::size_t i = 0; i < best_points.size(); ++i) {
    if (best_point_label[i] == target) target_layers_in_slice.insert(z_layer[i]);
  }
  std::cout << "Target particle touches " << target_layers_in_slice.size()
            << " distinct planes in this slice (out of " << result.Stats.Layers
            << " total planes with any activity).\n";

  std::cout << "\nGraph Track Finder on the real slice\n"
            << "  points: " << result.Stats.InputPoints << '\n'
            << "  z layers: " << result.Stats.Layers << '\n'
            << "  links tested/kept: " << result.Stats.LinksTested << "/"
            << result.Stats.LinksAccepted << '\n'
            << "  seeds generated/retained: " << result.Stats.SeedsGenerated
            << "/" << result.Stats.SeedsRetained << '\n'
            << "  hypotheses created/pruned: " << result.Stats.HypothesesCreated
            << "/" << result.Stats.HypothesesPruned << '\n'
            << "  paths before/after deduplication: "
            << result.Stats.PathsBeforeDeduplication << "/"
            << result.Stats.PathsAfterDeduplication << '\n'
            << "  resource limit reached: "
            << (result.Stats.ResourceLimitReached ? "YES (results may be truncated)" : "no") << '\n';

  std::size_t best_path_matched = 0;
  std::size_t best_path_planes = 0;
  for (std::size_t i = 0; i < result.Paths.size(); ++i) {
    const TMS_GraphTrackFinder::Path &path = result.Paths[i];
    std::size_t matched = 0;
    std::set<int> matched_layers;
    for (std::size_t idx : path.SpacePointIndices) {
      if (best_point_label[idx] == target) {
        ++matched;
        matched_layers.insert(z_layer[idx]);
      }
    }
    best_path_matched = std::max(best_path_matched, matched);
    best_path_planes = std::max(best_path_planes, matched_layers.size());
    const double purity = path.SpacePointIndices.empty() ? 0.0
        : 100.0 * matched / path.SpacePointIndices.size();
    std::cout << "  path " << i << ": points=" << path.SpacePointIndices.size()
              << ", layers=" << path.DistinctLayers << ", score=" << path.Score
              << ", target-matched=" << matched << " (" << purity << "% purity)"
              << ", target-planes=" << matched_layers.size() << "/"
              << target_layers_in_slice.size() << '\n';
  }

  std::cout << "\nFor comparison, DBSCAN+PCA on this same slice (2026-09-06 session):"
               " all " << best_total << " points land in one cluster with PCA linearity"
               " 0.159 (below the 0.8 track-like cutoff) -- the whole cluster is discarded,"
               " recovering 0/" << target_layers_in_slice.size() << " target planes.\n";
  std::cout << "Graph Track Finder recovered " << best_path_matched << "/" << best_count
            << " target-owned points (raw count, inflated by X/Y ghost pairing), covering "
            << best_path_planes << "/" << target_layers_in_slice.size()
            << " of the target's own distinct planes, in its best single path.\n";

  if (!result.Paths.empty()) {
    std::cout << "\nNext-layer candidate breakdown from the best path endpoint:" << std::endl;
    const CandidateInspection inspection = InspectNextLayerCandidates(
        best_points, best_point_label, z_layer, result.Paths.front(),
        1.25, 1.25, 3.0, 0.7, 0.05, 0.08, 0.35);
    std::cout << "  next-layer total=" << inspection.next_layer_total
              << ", pass-coarse-gates=" << inspection.next_layer_pass
              << ", kept-for-ranking=" << inspection.candidates.size() << std::endl;
    for (std::size_t i = 0; i < inspection.candidates.size() && i < 10; ++i) {
      const CandidateBreakdown &c = inspection.candidates[i];
      std::cout << "  cand " << i << ": idx=" << c.index
                << ", target=" << (c.target ? "yes" : "no")
                << ", xhit=" << c.x_hit
                << ", yhit=" << c.y_hit
                << ", dxdz=" << c.dxdz
                << ", dydz=" << c.dydz
                << ", cost=" << c.cost << std::endl;
    }
  }

  // Optional: dump the point cloud and every output path's geometry as JSON,
  // for an event display -- own/blob/other matches the same 3-way split
  // used by the earlier DBSCAN-truth-validation display (own = the target
  // particle; blob = a different track at the *same* vertex, e.g. a
  // co-produced shower; other = a different vertex entirely, or no truth).
  if (!output_json_path.empty()) {
    std::ofstream json(output_json_path);
    json << std::fixed;
    json << "{\"vgid\":" << target.vgid << ",\"trackid\":" << target.trackid
         << ",\"entry\":" << best_entry << ",\"spill\":" << best_spill
         << ",\"slice\":" << best_slice << ",\"n_total_slice\":" << best_total
         << ",\"n_own\":" << best_count
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
    writePoints("other_points", [&](std::size_t i) {
      return best_point_label[i].vgid != target.vgid;
    });

    json << ",\"paths\":[";
    for (std::size_t p = 0; p < result.Paths.size(); ++p) {
      if (p) json << ",";
      const TMS_GraphTrackFinder::Path &path = result.Paths[p];
      std::size_t matched = 0;
      for (std::size_t idx : path.SpacePointIndices) if (best_point_label[idx] == target) ++matched;
      const double purity = path.SpacePointIndices.empty() ? 0.0
          : 100.0 * matched / path.SpacePointIndices.size();
      json << "{\"score\":" << path.Score << ",\"n_matched\":" << matched
           << ",\"purity\":" << purity << ",\"points\":[";
      for (std::size_t k = 0; k < path.SpacePointIndices.size(); ++k) {
        if (k) json << ",";
        const std::size_t idx = path.SpacePointIndices[k];
        json << "{\"x\":" << best_points[idx].GetX() << ",\"y\":" << best_points[idx].GetY()
             << ",\"z\":" << best_points[idx].GetZ() << ",\"is_target\":"
             << (best_point_label[idx] == target ? "true" : "false") << "}";
      }
      json << "]}";
    }
    json << "]}";
    json.close();
    std::cout << "\nWrote event-display JSON to " << output_json_path << std::endl;
  }

  return 0;
}
