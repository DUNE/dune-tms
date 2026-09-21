#include "TMS_Passage.h"

#include <algorithm>
#include <cmath>
#include <map>

namespace TMS_Passage {

double Length(const Segment& s) {
  const double dx = s.stop[0] - s.start[0];
  const double dy = s.stop[1] - s.start[1];
  const double dz = s.stop[2] - s.start[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

namespace {

double Distance(const std::array<double, 3>& a, const std::array<double, 3>& b) {
  const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

}  // namespace

std::vector<std::vector<size_t>> BuildPassages(const std::vector<Segment>& segments,
                                               double max_gap_mm) {
  std::map<std::pair<int, int>, std::vector<size_t>> by_bar_and_trajectory;
  for (size_t i = 0; i < segments.size(); ++i) {
    by_bar_and_trajectory[{segments[i].bar_key, segments[i].trajectory_id}].push_back(i);
  }

  std::vector<std::vector<size_t>> passages;
  for (auto& entry : by_bar_and_trajectory) {
    std::vector<size_t>& idx = entry.second;
    std::stable_sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
      return segments[a].t_start < segments[b].t_start;
    });
    std::vector<size_t> current;
    for (size_t i : idx) {
      if (!current.empty() &&
          Distance(segments[current.back()].stop, segments[i].start) > max_gap_mm) {
        passages.push_back(current);
        current.clear();
      }
      current.push_back(i);
    }
    if (!current.empty()) passages.push_back(current);
  }

  std::stable_sort(passages.begin(), passages.end(),
                   [&](const std::vector<size_t>& a, const std::vector<size_t>& b) {
                     return segments[a.front()].t_start < segments[b.front()].t_start;
                   });
  return passages;
}

std::vector<OpticalDeposit> Resegment(const std::vector<Segment>& passage, double max_bin_length_mm) {
  std::vector<OpticalDeposit> bins;
  if (passage.empty()) return bins;

  // Arc-length coordinate of each segment along the passage.
  std::vector<double> s0(passage.size()), s1(passage.size());
  double total = 0;
  for (size_t j = 0; j < passage.size(); ++j) {
    s0[j] = total;
    total += Length(passage[j]);
    s1[j] = total;
  }

  const int nbins = (total > 0) ? std::max(1, static_cast<int>(std::ceil(total / max_bin_length_mm))) : 1;
  const double bin_len = (total > 0) ? total / nbins : 0.0;
  bins.resize(nbins);

  // Position/time at arc coordinate s, interpolated within whichever segment contains it.
  auto point_at = [&](double s, std::array<double, 3>& pos, double& t) {
    size_t j = 0;
    while (j + 1 < passage.size() && s > s1[j]) ++j;
    const double len = s1[j] - s0[j];
    const double f = (len > 0) ? std::min(1.0, std::max(0.0, (s - s0[j]) / len)) : 0.0;
    for (int k = 0; k < 3; ++k) pos[k] = passage[j].start[k] + f * (passage[j].stop[k] - passage[j].start[k]);
    t = passage[j].t_start + f * (passage[j].t_stop - passage[j].t_start);
  };

  for (int i = 0; i < nbins; ++i) {
    OpticalDeposit& b = bins[i];
    b.energy = 0;
    b.dx = bin_len;
    point_at((i + 0.5) * bin_len, b.position, b.t);
  }

  for (size_t j = 0; j < passage.size(); ++j) {
    const double seg_len = s1[j] - s0[j];
    if (seg_len <= 0) {
      // Zero-length step: no overlap to apportion, so the whole deposit goes to the bin
      // containing its position (clamped to the last bin) rather than being dropped.
      const int i = (bin_len > 0) ? std::min(nbins - 1, static_cast<int>(s0[j] / bin_len)) : 0;
      bins[i].energy += passage[j].energy;
      bins[i].shares.push_back({passage[j].trajectory_id, passage[j].energy});
      continue;
    }
    for (int i = 0; i < nbins; ++i) {
      const double lo = std::max(s0[j], i * bin_len);
      const double hi = std::min(s1[j], (i + 1) * bin_len);
      if (hi <= lo) continue;
      const double e = passage[j].energy * (hi - lo) / seg_len;
      bins[i].energy += e;
      bins[i].shares.push_back({passage[j].trajectory_id, e});
    }
  }
  return bins;
}

}  // namespace TMS_Passage
