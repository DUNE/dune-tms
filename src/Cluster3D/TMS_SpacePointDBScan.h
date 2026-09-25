#ifndef __TMS_SPACEPOINTDBSCAN_H__
#define __TMS_SPACEPOINTDBSCAN_H__

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
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
// Neighbor test is ANISOTROPIC, not a plain isotropic radius: two points
// are neighbors if they are at most Tolerance::MaxDzMM apart in z, and their
// transverse (x/y) separation is within an allowance that grows with that z
// distance -- BaseTransverseMM + TransversePerDzMM * |dz|. The growth term is
// the transverse distance a straight track of slope up to TransversePerDzMM
// covers over |dz|: a flat transverse tolerance breaks links between genuine
// points on an inclined track exactly where they are farther apart in z
// (confirmed on event 487, where within-cluster transverse jumps tracked the
// z separation). The base term is about one bar pitch, the position
// quantization of a space point.
//
// History: until 2026-09-25 both limits were counted in PLANE INDEX gaps
// (max 3 planes; one extra bar pitch of transverse allowance per plane),
// because the old space points (BothNeighbors pairing) sat 130 mm apart in
// the front section but 390 mm apart in the back, and no single mm limit
// suited both. With NearestY pairing (TMS_PlanePairing) consecutive point
// layers are 130 mm apart in front and alternate 130/260 mm in the back, and
// a point's z is a midpoint between planes rather than a plane's z, so the
// limits are now in mm of z (decision of 2026-09-25).
//
// The caller supplies the tolerance (so this class stays independent of
// TMS_Geom/geometry-loading, and fully testable with synthetic data);
// DefaultTolerance() gives the standard values for a given bar pitch.
//
// The cluster-growth logic itself mirrors src/TMS_DBScan.h's GrowCluster()
// (a known-working classic DBSCAN expansion), adapted to operate on point
// indices via a parallel cluster-id array instead of embedding a ClusterID
// field in a point struct.
class TMS_SpacePointDBScan {
  public:
    struct Tolerance {
      // Neighbors are at most this far apart in z (mm).
      double MaxDzMM = 270.0;
      // Transverse allowance at dz = 0 (mm): about one bar pitch.
      double BaseTransverseMM = 36.0;
      // Extra transverse allowance per mm of |dz|: the largest track slope
      // (transverse/z) a link should survive.
      double TransversePerDzMM = 0.55;
      // Radius that contains every possible neighbor, for the KD-tree prefilter.
      double BroadPhaseRadiusMM() const {
        const double transverse = BaseTransverseMM + TransversePerDzMM * MaxDzMM;
        return std::sqrt(transverse * transverse + MaxDzMM * MaxDzMM);
      }
    };

    // Default core-point threshold (neighbors, including the point itself).
    // 5 until 2026-09-25, set on BothNeighbors points (~2.4 per muon crossing
    // in the front section); NearestY points carry ~1.3, and 3 did best in
    // the Phase 1 sweep (reports/2026-09-25_phase1_baselines/).
    static constexpr unsigned int kDefaultMinPoints = 3;

    // Standard tolerance for a detector with the given bar pitch (mm), e.g.
    // TMS_Geom::GetMaxBarPitch(). Starting values for NearestY points,
    // 2026-09-25: MaxDzMM 270 reaches the next point layer anywhere (130 mm
    // in front, up to 260 mm in the back); TransversePerDzMM 0.55 (~29 deg)
    // covers the observed muon angles (90th percentile ~20 deg, max ~33 deg)
    // and equals the old front-section allowance (one bar pitch per 65 mm plane).
    static Tolerance DefaultTolerance(double barPitchMM) {
      Tolerance tolerance;
      tolerance.BaseTransverseMM = barPitchMM;
      return tolerance;
    }

    TMS_SpacePointDBScan(const std::vector<TMS_SpacePoint> &sp, unsigned int min_points, const Tolerance &tolerance)
      : _sp(sp),
        _positions(BuildPositions(sp)),
        _tree(_positions),
        _min_points(min_points),
        _tolerance(tolerance),
        _broad_phase_radius(tolerance.BroadPhaseRadiusMM()) {}

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
      std::vector<int> candidates = _tree.RadiusQuery(point_index, _broad_phase_radius);
      std::vector<int> out;
      out.reserve(candidates.size());
      const TMS_SpacePoint &seed = _sp[point_index];
      for (int c : candidates) {
        const double dz = std::abs(_sp[c].GetZ() - seed.GetZ());
        if (dz > _tolerance.MaxDzMM) continue;
        const double allowed_transverse = _tolerance.BaseTransverseMM + _tolerance.TransversePerDzMM * dz;
        double dx = _sp[c].GetX() - seed.GetX();
        double dy = _sp[c].GetY() - seed.GetY();
        if (std::sqrt(dx * dx + dy * dy) > allowed_transverse) continue;
        out.push_back(c);
      }
      return out;
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

      // FindNeighbours includes the query point itself, so this also labels the seed.
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
    Tolerance _tolerance;
    double _broad_phase_radius;
};

#endif
