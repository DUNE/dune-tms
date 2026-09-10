#include "TMS_SpacePointBuilder.h"

#include <algorithm>
#include <map>
#include <utility>

#include "TMS_Hit.h"

std::vector<TMS_SpacePoint> TMS_SpacePointBuilder::Build(
    const std::vector<TMS_Hit> &hits, double timing_window) {
  std::vector<TMS_SpacePoint> space_points;

  // Bucket hits by detector plane and by which axis the bar measures (X-bars
  // and Y-bars are laid out in alternating planes -- see TMS_Bar.h).
  std::map<int, std::vector<int>> x_hits_by_layer;  // plane -> indices into hits
  std::map<int, std::vector<int>> y_hits_by_layer;

  for (size_t i = 0; i < hits.size(); ++i) {
    const TMS_Hit &hit = hits[i];
    // A pedestal-suppressed hit is noise-level and isn't treated as real
    // anywhere else in reconstruction (TMS_TrackFinder::FindTracks()
    // excludes them the same way, via GetHits()'s default
    // include_ped_sup=false). Pairing them here would inflate the ghost
    // space-point population with combinations that don't even correspond
    // to a real reconstructed hit, on top of the expected real-hit-wrong-
    // particle ghosting.
    if (hit.GetPedSup()) continue;
    int layer = hit.GetBar().GetPlaneNumber();
    TMS_Bar::BarType bar_type = hit.GetBar().GetBarType();

    if (bar_type == TMS_Bar::kXBar) {
      x_hits_by_layer[layer].push_back(i);
    } else if (bar_type == TMS_Bar::kYBar) {
      y_hits_by_layer[layer].push_back(i);
    }
  }

  // Sort each plane's hits by time so the X/Y matching loop below can scan
  // forward instead of comparing every X hit against every Y hit.
  std::map<int, std::vector<std::pair<double, int>>> y_hits_by_layer_sorted;  // plane -> (time, index)
  for (const auto &y_layer_entry : y_hits_by_layer) {
    int y_layer = y_layer_entry.first;
    const std::vector<int> &y_indices = y_layer_entry.second;

    std::vector<std::pair<double, int>> time_index_pairs;
    for (int y_idx : y_indices) {
      time_index_pairs.push_back({hits[y_idx].GetT(), y_idx});
    }
    std::sort(time_index_pairs.begin(), time_index_pairs.end());
    y_hits_by_layer_sorted[y_layer] = time_index_pairs;
  }

  std::map<int, std::vector<std::pair<double, int>>> x_hits_by_layer_sorted;
  for (const auto &x_layer_entry : x_hits_by_layer) {
    int x_layer = x_layer_entry.first;
    const std::vector<int> &x_indices = x_layer_entry.second;

    std::vector<std::pair<double, int>> time_index_pairs;
    for (int x_idx : x_indices) {
      time_index_pairs.push_back({hits[x_idx].GetT(), x_idx});
    }
    std::sort(time_index_pairs.begin(), time_index_pairs.end());
    x_hits_by_layer_sorted[x_layer] = time_index_pairs;
  }

  // Create a space point for every X hit / Y hit pair in adjacent planes
  // (N-1, N+1) whose times land within timing_window of each other.
  for (const auto &x_layer_entry : x_hits_by_layer_sorted) {
    int x_layer = x_layer_entry.first;
    const auto &x_time_indices = x_layer_entry.second;

    for (int y_layer : {x_layer - 1, x_layer + 1}) {
      if (y_hits_by_layer_sorted.find(y_layer) == y_hits_by_layer_sorted.end()) {
        continue;  // No Y hits in this adjacent plane
      }

      const auto &y_time_indices = y_hits_by_layer_sorted[y_layer];

      // Forward-only pointer into y_time_indices. Both sequences are time-sorted,
      // so as x_time only increases across this loop, the "too early" cutoff
      // (x_time - timing_window) only increases too -- once a Y hit falls below
      // it, it falls below it for every later X hit as well, so it can be
      // permanently skipped instead of re-scanned from the start each time.
      size_t y_start = 0;

      for (const auto &x_time_idx : x_time_indices) {
        double x_time = x_time_idx.first;
        int x_idx = x_time_idx.second;
        const TMS_Hit &x_hit = hits[x_idx];
        // An X-bar is oriented along X, so the coordinate it actually measures is Y
        // (its Bar.x member is a sentinel -- see TMS_Bar.cpp's kXBar branch). GetNotZ()
        // already knows to return the real value for whichever axis the bar measures.
        double y_pos = x_hit.GetNotZ();
        double z_pos = x_hit.GetZ();

        // Advance past Y hits that are now permanently too early
        while (y_start < y_time_indices.size() &&
               y_time_indices[y_start].first - x_time < -timing_window) {
          ++y_start;
        }

        // Scan forward from y_start until Y hits become too late for this X hit.
        // The upper cutoff (x_time + timing_window) is NOT monotonic in the same
        // way, so this scan (unlike y_start) has to restart from y_start every time.
        for (size_t j = y_start; j < y_time_indices.size(); ++j) {
          double y_time = y_time_indices[j].first;
          int y_idx = y_time_indices[j].second;
          // No time-of-flight correction between the X and Y planes: even the
          // largest adjacent-plane gap (130mm, double-thick region) is only
          // ~0.43ns at beta~1, under 1.5% of the timing_window below -- negligible
          // next to the 30ns window and not worth the added complexity.
          double time_diff = y_time - x_time;

          if (time_diff > timing_window) {
            break;
          }

          const TMS_Hit &y_hit = hits[y_idx];
          // Symmetrically, a Y-bar measures X (its Bar.y member is the sentinel).
          double x_pos = y_hit.GetNotZ();

          double combined_time = (x_time + y_time) / 2.0;
          space_points.push_back(TMS_SpacePoint(x_pos, y_pos, z_pos, x_idx, y_idx, combined_time));
        }
      }
    }
  }

  return space_points;
}
