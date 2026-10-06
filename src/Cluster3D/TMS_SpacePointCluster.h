#ifndef __TMS_SPACEPOINTCLUSTER_H__
#define __TMS_SPACEPOINTCLUSTER_H__

#include <algorithm>
#include <array>
#include <vector>

#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"

#include "TMS_SpacePoint.h"

// A DBSCAN cluster of TMS_SpacePoints (see TMS_SpacePointDBScan), with a PCA
// shape analysis to classify it as "track-like" (elongated along one axis)
// vs. blob-like (diffuse / shower-like).
//
// Uses ROOT's TMatrixDSym + TMatrixDSymEigen for the 3x3 eigen-decomposition
// (already linked, already used identically throughout TMS_Kalman.cpp) rather
// than TPrincipal, which is a heavier API for a fixed 3x3 case.
class TMS_SpacePointCluster {
  public:
    // sp: the full space-point array this cluster's indices refer into.
    // indices: this cluster's space-point indices, e.g. one entry of
    // TMS_SpacePointDBScan::RunAndGetClusterIndices()'s result.
    TMS_SpacePointCluster(const std::vector<TMS_SpacePoint> &sp, std::vector<int> indices)
      : _indices(std::move(indices)) {
      ComputePCA(sp);
    }

    size_t GetSize() const { return _indices.size(); }
    const std::vector<int> &GetSpacePointIndices() const { return _indices; }

    const std::array<double, 3> &GetCentroid() const { return _centroid; }
    // Descending: eigenvalue[0] >= eigenvalue[1] >= eigenvalue[2] >= 0.
    const std::array<double, 3> &GetEigenvalues() const { return _eigenvalues; }
    // Unit vector along the largest-eigenvalue eigenvector (the cluster's estimated direction).
    const std::array<double, 3> &GetPrincipalDirection() const { return _principal_direction; }
    // All three unit eigenvectors, same descending eigenvalue order as GetEigenvalues();
    // _eigenvectors[0] == GetPrincipalDirection().
    const std::array<std::array<double, 3>, 3> &GetEigenvectors() const { return _eigenvectors; }

    // (lambda1 - lambda2) / lambda1 -- close to 1 for a line, close to 0 for a blob/disk.
    // 0 if the cluster is too small (<3 points) for a meaningful 3D PCA.
    double GetLinearity() const { return _linearity; }
    bool HasValidPCA() const { return _valid_pca; }

    // Default minimum size for IsTrackLike(). 5 until 2026-09-25, set on
    // BothNeighbors points; with NearestY points (about half the front-section
    // point density) short muons fell below it. Phase 1 benchmark: 4 matches
    // the old ND-LAr-fiducial finding efficiency (92.5% vs 92.7%) and gains
    // 4.6 pp overall, with 48% more non-muon track-like clusters (3: +6.9 pp,
    // but 2.4x). To revisit with hit-level purity after the hit-level fit.
    static constexpr size_t kDefaultMinTrackSize = 4;

    bool IsTrackLike(double linearity_threshold, size_t min_cluster_size) const {
      return _valid_pca && _linearity >= linearity_threshold && GetSize() >= min_cluster_size;
    }

  private:
    void ComputePCA(const std::vector<TMS_SpacePoint> &sp) {
      _centroid = {0.0, 0.0, 0.0};
      _eigenvalues = {0.0, 0.0, 0.0};
      _principal_direction = {0.0, 0.0, 0.0};
      _linearity = 0.0;
      _valid_pca = false;

      const size_t n = _indices.size();
      if (n < 3) return; // degenerate covariance -- not enough points for a meaningful 3D PCA

      for (int idx : _indices) {
        const TMS_SpacePoint &p = sp[idx];
        _centroid[0] += p.GetX();
        _centroid[1] += p.GetY();
        _centroid[2] += p.GetZ();
      }
      for (double &c : _centroid) c /= static_cast<double>(n);

      TMatrixDSym cov(3);
      cov.Zero();
      for (int idx : _indices) {
        const TMS_SpacePoint &p = sp[idx];
        double d[3] = {p.GetX() - _centroid[0], p.GetY() - _centroid[1], p.GetZ() - _centroid[2]};
        for (int j = 0; j < 3; ++j) {
          for (int k = 0; k < 3; ++k) {
            cov(j, k) += d[j] * d[k];
          }
        }
      }
      cov *= 1.0 / static_cast<double>(n);

      TMatrixDSymEigen eig(cov);
      const TVectorD &eigenvalues = eig.GetEigenValues();
      const TMatrixD &eigenvectors = eig.GetEigenVectors(); // column i = eigenvector for eigenvalues[i]

      // Sort descending explicitly -- don't rely on TMatrixDSymEigen's internal ordering convention.
      std::array<int, 3> order = {0, 1, 2};
      std::sort(order.begin(), order.end(), [&](int a, int b) { return eigenvalues[a] > eigenvalues[b]; });

      for (int i = 0; i < 3; ++i) _eigenvalues[i] = eigenvalues[order[i]];
      for (int rank = 0; rank < 3; ++rank) {
        for (int i = 0; i < 3; ++i) _eigenvectors[rank][i] = eigenvectors(i, order[rank]);
      }
      _principal_direction = _eigenvectors[0];

      _valid_pca = true;
      _linearity = (_eigenvalues[0] > 0.0) ? (_eigenvalues[0] - _eigenvalues[1]) / _eigenvalues[0] : 0.0;
    }

    std::vector<int> _indices;
    std::array<double, 3> _centroid{};
    std::array<double, 3> _eigenvalues{};
    std::array<double, 3> _principal_direction{};
    std::array<std::array<double, 3>, 3> _eigenvectors{};
    double _linearity = 0.0;
    bool _valid_pca = false;
};

#endif
