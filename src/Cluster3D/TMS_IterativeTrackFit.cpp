#include "TMS_IterativeTrackFit.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>

#include "TMS_SpacePointCluster.h"
#include "TMS_SpacePointDBScan.h"

namespace TMS_IterativeTrackFit {

bool IsFlaggedAsMerged(const std::vector<TMS_SpacePoint> &slicePoints, const std::vector<int> &planeIndex,
                       const std::vector<int> &clusterIndices, const Config &config) {
  if (clusterIndices.empty()) return false;

  // Per-plane bounding box of the cluster's points.
  std::map<int, std::array<double, 4> > box;  // plane -> xmin, xmax, ymin, ymax
  std::map<int, int> count;
  double tSum = 0.0, tSum2 = 0.0;
  for (int idx : clusterIndices) {
    const TMS_SpacePoint &p = slicePoints[idx];
    const int plane = planeIndex[idx];
    auto it = box.find(plane);
    if (it == box.end()) {
      box[plane] = {p.GetX(), p.GetX(), p.GetY(), p.GetY()};
    } else {
      std::array<double, 4> &b = it->second;
      b[0] = std::min(b[0], p.GetX());
      b[1] = std::max(b[1], p.GetX());
      b[2] = std::min(b[2], p.GetY());
      b[3] = std::max(b[3], p.GetY());
    }
    ++count[plane];
    tSum += p.GetTime();
    tSum2 += p.GetTime() * p.GetTime();
  }

  int nMulti = 0, nWide = 0;
  for (const auto &kv : box) {
    if (count[kv.first] < 2) continue;
    ++nMulti;
    const double span = std::max(kv.second[1] - kv.second[0], kv.second[3] - kv.second[2]);
    if (span > config.WideSpanMM) ++nWide;
  }
  const bool wide = nMulti > 0 && static_cast<double>(nWide) / nMulti > config.WideLayerFraction;

  const double n = static_cast<double>(clusterIndices.size());
  const double tMean = tSum / n;
  const double tRMS = std::sqrt(std::max(0.0, tSum2 / n - tMean * tMean));
  return wide || tRMS > config.TimeRMSNs;
}

namespace {

// Re-cluster a set of slice points and return the largest track-like
// sub-cluster (as slice indices), or empty if none qualifies.
std::vector<int> LargestTrackLikeSubCluster(const std::vector<TMS_SpacePoint> &slicePoints,
                                            const std::vector<int> &planeIndex, const std::vector<int> &indices,
                                            const Config &config) {
  if (indices.size() < config.MinClusterSizeForTrack) return {};
  std::vector<TMS_SpacePoint> subPoints;
  std::vector<int> subPlanes;
  subPoints.reserve(indices.size());
  subPlanes.reserve(indices.size());
  for (int idx : indices) {
    subPoints.push_back(slicePoints[idx]);
    subPlanes.push_back(planeIndex[idx]);
  }
  TMS_SpacePointDBScan dbscan(subPoints, subPlanes, config.DBScanMinPoints, config.BarPitchMM,
                              config.BaseTransverseBars, config.TransverseBarsPerPlaneGap, config.MaxPlaneGap,
                              config.BroadPhaseRadiusMM);
  const std::vector<std::vector<int> > subClusters = dbscan.RunAndGetClusterIndices();

  std::vector<int> best;
  for (const std::vector<int> &sub : subClusters) {
    TMS_SpacePointCluster cluster(subPoints, sub);
    if (!cluster.IsTrackLike(config.LinearityThreshold, config.MinClusterSizeForTrack)) continue;
    if (sub.size() <= best.size()) continue;
    best.clear();
    for (int local : sub) best.push_back(indices[local]);  // back to slice indices
  }
  return best;
}

// Fit one object over the slice minus every point using a claimed hit.
// Returns false if the object has no points left in the pool.
// allowed (optional): if non-null, only these slice points may enter the pool.
bool FitObject(const std::vector<TMS_SpacePoint> &slicePoints, const std::vector<int> &objectIndices,
               const TMS_KalmanFollower::Follower &follower, const ClaimedHits &claimed,
               const std::vector<char> *allowed, TMS_KalmanFollower::FitResult &fitOut) {
  std::vector<TMS_SpacePoint> pool;
  std::vector<std::size_t> poolToSlice;
  std::vector<long> sliceToPool(slicePoints.size(), -1);
  for (std::size_t i = 0; i < slicePoints.size(); ++i) {
    if (claimed.Uses(slicePoints[i])) continue;
    if (allowed && !(*allowed)[i]) continue;
    sliceToPool[i] = static_cast<long>(pool.size());
    pool.push_back(slicePoints[i]);
    poolToSlice.push_back(i);
  }
  std::vector<std::size_t> object;
  for (int idx : objectIndices)
    if (sliceToPool[idx] >= 0) object.push_back(static_cast<std::size_t>(sliceToPool[idx]));
  if (object.size() < 2) return false;

  fitOut = follower.RunBestSeed(pool, object);
  // Remap every index the fit reports back to the slice's own indexing.
  for (TMS_KalmanFollower::FollowedNode &node : fitOut.Nodes) {
    if (node.HasHit) node.ChosenSpacePointIndex = poolToSlice[node.ChosenSpacePointIndex];
    for (std::size_t &c : node.CandidateIndices) c = poolToSlice[c];
  }
  return true;
}

int CountHits(const TMS_KalmanFollower::FitResult &fit) {
  int n = 0;
  for (const TMS_KalmanFollower::FollowedNode &node : fit.Nodes)
    if (node.HasHit) ++n;
  return n;
}

}  // namespace

std::vector<Track> FitCluster(const std::vector<TMS_SpacePoint> &slicePoints, const std::vector<int> &planeIndex,
                              const std::vector<int> &clusterIndices, const TMS_KalmanFollower::Follower &follower,
                              const Config &config, ClaimedHits &claimed) {
  std::vector<Track> tracks;
  const bool flagged = IsFlaggedAsMerged(slicePoints, planeIndex, clusterIndices, config);
  const bool maySplit = config.Mode == Config::SplitMode::Always ||
                        (config.Mode == Config::SplitMode::Flagged && flagged);
  const int maxTracks = maySplit ? config.MaxTracksPerCluster : 1;
  std::vector<char> inCluster(slicePoints.size(), 0);
  for (int idx : clusterIndices) inCluster[idx] = 1;

  for (int iteration = 0; iteration < maxTracks; ++iteration) {
    // The cluster's points that no accepted track has claimed yet.
    std::vector<int> remaining;
    for (int idx : clusterIndices)
      if (!claimed.Uses(slicePoints[idx])) remaining.push_back(idx);

    // First pass: the cluster as the caller found it. Later passes: the
    // largest track-like piece of what is left, re-clustered from scratch,
    // since removing one track's points usually leaves the other particle
    // plus scattered leftovers.
    const std::vector<int> object =
        iteration == 0 ? remaining : LargestTrackLikeSubCluster(slicePoints, planeIndex, remaining, config);
    if (object.size() < config.MinClusterSizeForTrack) break;

    Track track;
    const std::vector<char> *allowed =
        (iteration > 0 && config.RestrictSplitFitToCluster) ? &inCluster : nullptr;
    if (!FitObject(slicePoints, object, follower, claimed, allowed, track.Fit)) break;
    if (CountHits(track.Fit) < (iteration == 0 ? config.MinHitsPerTrack : config.MinHitsPerSplitTrack)) break;

    for (const TMS_KalmanFollower::FollowedNode &node : track.Fit.Nodes) {
      if (!node.HasHit) continue;
      // Negative = no hit index recorded (see TMS_SpacePoint); never claim it,
      // or every other index-less point would look claimed too.
      const TMS_SpacePoint &chosen = slicePoints[node.ChosenSpacePointIndex];
      if (chosen.GetXHitIndex() >= 0) claimed.X.insert(chosen.GetXHitIndex());
      if (chosen.GetYHitIndex() >= 0) claimed.Y.insert(chosen.GetYHitIndex());
    }
    track.ObjectIndices.assign(object.begin(), object.end());
    track.Iteration = iteration;
    track.ClusterFlagged = flagged;
    tracks.push_back(std::move(track));
  }
  return tracks;
}

}  // namespace TMS_IterativeTrackFit
