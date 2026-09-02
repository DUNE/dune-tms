#ifndef __TMS_SPACEPOINTDBSCAN_H__
#define __TMS_SPACEPOINTDBSCAN_H__

#include <algorithm>
#include <array>
#include <vector>

#include "TMS_KDTree.h"
#include "TMS_SpacePoint.h"

// DBSCAN clustering over TMS_SpacePoints, backed by TMS_KDTree for
// KD-tree-accelerated neighbor lookups.
//
// This is a SEPARATE implementation from the existing src/TMS_DBScan.h, not a
// reuse of it: that class clusters on discretized (GetPlaneNumber(),
// GetBarNumber()) integers, which is fragile to sub-bar-width geometry shifts
// (see project history/memory "DBSCAN BarNumber sensitivity"). This class
// clusters on continuous physical (x,y,z) positions in mm instead, and never
// re-derives cluster membership from coordinates -- it carries space-point
// array indices through directly, so there's no analogous re-match fragility.
//
// The cluster-growth logic itself mirrors src/TMS_DBScan.h's GrowCluster()
// (a known-working classic DBSCAN expansion), adapted to operate on point
// indices via a parallel cluster-id array instead of embedding a ClusterID
// field in a point struct.
class TMS_SpacePointDBScan {
  public:
    TMS_SpacePointDBScan(const std::vector<TMS_SpacePoint> &sp, unsigned int min_points, double epsilon_mm)
      : _sp(sp),
        _positions(BuildPositions(sp)),
        _tree(_positions),
        _min_points(min_points),
        _epsilon(epsilon_mm) {}

    // Runs DBSCAN and returns clusters as vectors of SPACE-POINT INDICES
    // (into the `sp` vector passed to the constructor). Points classified as
    // noise are not included in any returned cluster.
    std::vector<std::vector<int>> RunAndGetClusterIndices() {
      const size_t n = _sp.size();
      std::vector<int> cluster_id_of_point(n, kUnclassified);
      int next_cluster_id = 1;

      for (size_t i = 0; i < n; ++i) {
        if (cluster_id_of_point[i] != kUnclassified) continue;
        if (GrowCluster(static_cast<int>(i), next_cluster_id, cluster_id_of_point)) {
          ++next_cluster_id;
        } else {
          cluster_id_of_point[i] = kNoise;
        }
      }

      const int n_clusters = next_cluster_id - 1;
      std::vector<std::vector<int>> clusters(n_clusters);
      for (size_t i = 0; i < n; ++i) {
        if (cluster_id_of_point[i] > 0) {
          clusters[cluster_id_of_point[i] - 1].push_back(static_cast<int>(i));
        }
      }
      return clusters;
    }

  private:
    static constexpr int kUnclassified = -1;
    static constexpr int kNoise = 0;

    static std::vector<std::array<double, 3>> BuildPositions(const std::vector<TMS_SpacePoint> &sp) {
      std::vector<std::array<double, 3>> positions;
      positions.reserve(sp.size());
      for (const auto &p : sp) {
        positions.push_back({p.GetX(), p.GetY(), p.GetZ()});
      }
      return positions;
    }

    std::vector<int> FindNeighbours(int point_index) const {
      return _tree.RadiusQuery(point_index, _epsilon);
    }

    // Mirrors src/TMS_DBScan.h's GrowCluster(): expand a cluster outward from
    // seed_index, absorbing unclassified/noise points reachable through
    // chains of dense (>= _min_points neighbours) points. Returns false (seed
    // should be marked noise by the caller) if the seed itself isn't a core
    // point.
    bool GrowCluster(int seed_index, int cluster_id, std::vector<int> &cluster_id_of_point) const {
      std::vector<int> neighbours = FindNeighbours(seed_index);
      if (neighbours.size() < _min_points) {
        return false;
      }

      // RadiusQuery includes the query point itself, so this also labels the seed.
      for (int idx : neighbours) {
        cluster_id_of_point[idx] = cluster_id;
      }
      neighbours.erase(std::remove(neighbours.begin(), neighbours.end(), seed_index), neighbours.end());

      size_t n = neighbours.size();
      for (size_t i = 0; i < n; ++i) {
        std::vector<int> nn = FindNeighbours(neighbours[i]);
        if (nn.size() < _min_points) continue;
        for (int j : nn) {
          if (cluster_id_of_point[j] == kUnclassified || cluster_id_of_point[j] == kNoise) {
            if (cluster_id_of_point[j] == kUnclassified) {
              neighbours.push_back(j);
              n = neighbours.size();
            }
            cluster_id_of_point[j] = cluster_id;
          }
        }
      }
      return true;
    }

    const std::vector<TMS_SpacePoint> &_sp;
    std::vector<std::array<double, 3>> _positions;
    TMS_KDTree _tree;
    unsigned int _min_points;
    double _epsilon;
};

#endif
