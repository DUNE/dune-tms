#include "TMS_SpacePointBuilder.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>

#include "TMS_Hit.h"
#include "TMS_SpacePointTiming.h"

namespace {

// (time, index into hits), sorted by time.
typedef std::vector<std::pair<double, int>> TimeSortedHits;

// Calls accept(x_idx, y_idx) for every X-bar hit / Y-bar hit pair whose
// times are within timing_window of each other, in X-bar-time order.
// Both inputs are time-sorted, so the scan over Y-bar hits can start from a
// forward-only pointer: as the X-bar time only increases, the "too early"
// cutoff (x_time - timing_window) only increases too, so a Y-bar hit that
// falls below it is permanently skipped instead of re-scanned each time.
template <typename Accept>
void PairWithinWindow(const TimeSortedHits &x_bar_hits, const TimeSortedHits &y_bar_hits,
                      double timing_window, Accept accept) {
  size_t y_start = 0;
  for (const auto &x_entry : x_bar_hits) {
    const double x_time = x_entry.first;
    while (y_start < y_bar_hits.size() && y_bar_hits[y_start].first - x_time < -timing_window) ++y_start;
    // The upper cutoff (x_time + timing_window) is NOT monotonic in the same
    // way, so this scan has to restart from y_start every time.
    for (size_t j = y_start; j < y_bar_hits.size(); ++j) {
      // No time-of-flight correction between the two planes: even the widest
      // pair (~130 mm apart) is ~0.4 ns at beta~1, negligible next to the window.
      if (y_bar_hits[j].first - x_time > timing_window) break;
      accept(x_entry.second, y_bar_hits[j].second);
    }
  }
}

}  // namespace

std::vector<TMS_SpacePoint> TMS_SpacePointBuilder::Build(
    const std::vector<TMS_Hit> &hits, double timing_window,
    const TMS_PlanePairing::Table &pairing, bool use_fallback, bool require_crossing,
    double crossing_slope, double corrected_timing_window) {
  std::vector<TMS_SpacePoint> space_points;

  // Bucket hits by plane, time-sorted. Only X-bar and Y-bar planes take part
  // (the pairing table leaves every other orientation out anyway).
  std::map<int, TimeSortedHits> hits_by_plane;
  for (size_t i = 0; i < hits.size(); ++i) {
    const TMS_Hit &hit = hits[i];
    // A pedestal-suppressed hit is noise-level and isn't treated as real
    // anywhere else in reconstruction (TMS_TrackFinder::FindTracks()
    // excludes them the same way, via GetHits()'s default
    // include_ped_sup=false). Pairing them here would inflate the ghost
    // space-point population with combinations that don't even correspond
    // to a real reconstructed hit.
    if (hit.GetPedSup()) continue;
    const TMS_Bar::BarType bar_type = hit.GetBar().GetBarType();
    if (bar_type != TMS_Bar::kXBar && bar_type != TMS_Bar::kYBar) continue;
    hits_by_plane[hit.GetBar().GetPlaneNumber()].push_back({hit.GetT(), static_cast<int>(i)});
  }
  for (auto &entry : hits_by_plane) std::sort(entry.second.begin(), entry.second.end());

  static const TimeSortedHits kNoHits;
  const auto plane_hits = [&](int plane) -> const TimeSortedHits & {
    auto it = hits_by_plane.find(plane);
    return it == hits_by_plane.end() ? kNoHits : it->second;
  };
  // With a corrected window the scan runs on raw times, so it has to be wide
  // enough to contain every pair the corrected test could accept: the raw
  // difference differs from the corrected one by at most the two hits' largest
  // possible transit delays. The corrected test in accept() then makes the cut.
  // The bound is exact even for pairs that bars_cross() admits slightly outside the
  // X-bar's ends (width/slope margin) or with crossing checks off: TransitDelayNs clamps
  // the along-bar position into the bar, so |delay| <= MaxTransitDelayNs(bar) always.
  const bool use_corrected = corrected_timing_window > 0.0;
  double scan_window = timing_window;
  if (use_corrected) {
    double max_x_delay = 0.0, max_y_delay = 0.0;
    for (const auto &entry : hits_by_plane) {
      for (const auto &timed_hit : entry.second) {
        const TMS_Bar &bar = hits[timed_hit.second].GetBar();
        double &max_delay = bar.GetBarType() == TMS_Bar::kXBar ? max_x_delay : max_y_delay;
        max_delay = std::max(max_delay, TMS_SpacePointTiming::MaxTransitDelayNs(bar));
      }
    }
    scan_window = corrected_timing_window + max_x_delay + max_y_delay;
  }
  // Does this pair pass the corrected time cut? (Always true with the raw
  // window, which the scan itself already applied.)
  const auto passes_time_cut = [&](int x_idx, int y_idx) {
    if (!use_corrected) return true;
    const TMS_Hit &x_bar_hit = hits[x_idx];
    const TMS_Hit &y_bar_hit = hits[y_idx];
    // The candidate point: x is where the Y-bar sits, y where the X-bar sits.
    const double dt = TMS_SpacePointTiming::PairedHitsTimeDifference(
        x_bar_hit.GetBar(), y_bar_hit.GetNotZ(), x_bar_hit.GetT(),
        y_bar_hit.GetBar(), x_bar_hit.GetNotZ(), y_bar_hit.GetT());
    return std::abs(dt) <= corrected_timing_window;
  };
  // Could one particle have crossed both bars? An X-bar spans x in
  // [center - length/2, center + length/2] (its own half of the detector when
  // it is one of the two halves split at x = 0); the Y-bar at x = GetNotZ()
  // is GetXw() wide. The two bars sit in different planes, so a track moves
  // sideways between them: allow crossing_slope times the planes' z gap on
  // top of the Y-bar's half width.
  const auto bars_cross = [&](int x_idx, int y_idx) {
    if (!require_crossing) return true;
    const TMS_Bar &x_bar = hits[x_idx].GetBar();
    const TMS_Bar &y_bar = hits[y_idx].GetBar();
    const double x_center = x_bar.GetAxisReadoutCenter();
    const double x_half_length = 0.5 * x_bar.GetBarLength();
    const double y_bar_x = y_bar.GetNotZ();
    const double margin = 0.5 * y_bar.GetXw() + crossing_slope * std::abs(x_bar.GetZ() - y_bar.GetZ());
    return y_bar_x + margin >= x_center - x_half_length &&
           y_bar_x - margin <= x_center + x_half_length;
  };
  const auto make_point = [&](int x_idx, int y_idx, const TMS_PlanePairing::PlanePair &pair) {
    // An X-bar measures y (GetNotZ() returns it); a Y-bar measures x.
    const TMS_Hit &x_bar_hit = hits[x_idx];
    const TMS_Hit &y_bar_hit = hits[y_idx];
    const double combined_time = (x_bar_hit.GetT() + y_bar_hit.GetT()) / 2.0;
    space_points.push_back(TMS_SpacePoint(y_bar_hit.GetNotZ(), x_bar_hit.GetNotZ(), pair.Z,
                                          x_idx, y_idx, combined_time, pair.Layer));
  };

  // Primary pairs, in the table's z order.
  std::vector<char> paired(hits.size(), 0);
  for (const TMS_PlanePairing::PlanePair &pair : pairing.Pairs) {
    if (pair.Fallback) continue;
    PairWithinWindow(plane_hits(pair.XBarPlane), plane_hits(pair.YBarPlane), scan_window,
                     [&](int x_idx, int y_idx) {
                       if (!bars_cross(x_idx, y_idx) || !passes_time_cut(x_idx, y_idx)) return;
                       make_point(x_idx, y_idx, pair);
                       paired[x_idx] = 1;
                       paired[y_idx] = 1;
                     });
  }
  if (!use_fallback) return space_points;

  // Fallback pairs: a combination is made if at least one of its two hits
  // was left unpaired by EVERY primary pair (its partner may be paired
  // already). The window scan visits each combination once, so a pair of two
  // unpaired hits is made once.
  for (const TMS_PlanePairing::PlanePair &pair : pairing.Pairs) {
    if (!pair.Fallback) continue;
    PairWithinWindow(plane_hits(pair.XBarPlane), plane_hits(pair.YBarPlane), scan_window,
                     [&](int x_idx, int y_idx) {
                       if (paired[x_idx] && paired[y_idx]) return;
                       if (!bars_cross(x_idx, y_idx) || !passes_time_cut(x_idx, y_idx)) return;
                       make_point(x_idx, y_idx, pair);
                     });
  }
  return space_points;
}
