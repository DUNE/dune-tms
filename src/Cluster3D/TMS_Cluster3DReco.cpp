#include "TMS_Cluster3DReco.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <utility>

#include "TMS_Geom.h"
#include "TMS_Hit.h"
#include "TMS_SpacePointCluster.h"
#include "TMS_SpacePointTiming.h"

namespace TMS_Cluster3DReco {

std::vector<TMS_KalmanFollower::FitHit> BuildFitHits(const std::vector<TMS_Hit> &hits, double barPitchMM) {
  std::vector<TMS_KalmanFollower::FitHit> fitHits(hits.size());
  for (std::size_t i = 0; i < hits.size(); ++i) {
    const TMS_Hit &hit = hits[i];
    const TMS_Bar::BarType type = hit.GetBar().GetBarType();
    TMS_KalmanFollower::FitHit &fitHit = fitHits[i];
    fitHit.Z = hit.GetZ();
    fitHit.Coordinate = hit.GetNotZ();
    fitHit.MeasuresX = type == TMS_Bar::kYBar;
    fitHit.SigmaMM = barPitchMM / std::sqrt(12.0);
    fitHit.Time = hit.GetT();
    fitHit.Usable = !hit.GetPedSup() && (type == TMS_Bar::kXBar || type == TMS_Bar::kYBar);
  }
  return fitHits;
}

namespace {

// Every hit a fit used: applied hits (hit-level fit) or both hits of each
// chosen point (space-point fit), plus orphans; sorted, unique.
std::vector<int> UsedHits(const TMS_KalmanFollower::FitResult &fit, const std::vector<TMS_SpacePoint> &points) {
  std::set<int> used;
  for (const TMS_KalmanFollower::FollowedNode &node : fit.Nodes) {
    if (!node.HasHit) continue;
    if (!node.Hits.empty()) {
      for (const auto &update : node.Hits)
        if (update.Applied) used.insert(update.HitIndex);
    } else {
      const TMS_SpacePoint &point = points[node.ChosenSpacePointIndex];
      if (point.GetXHitIndex() >= 0) used.insert(point.GetXHitIndex());
      if (point.GetYHitIndex() >= 0) used.insert(point.GetYHitIndex());
    }
  }
  for (const auto &orphan : fit.Orphans) used.insert(orphan.HitIndex);
  return std::vector<int>(used.begin(), used.end());
}

}  // namespace

std::vector<Track> Run(const std::vector<TMS_SpacePoint> &points,
                       const std::vector<TMS_KalmanFollower::FitHit> &hits, const Config &config,
                       const IFieldModel &field, RunInfo *info) {
  std::vector<Track> tracks;
  if (info) *info = RunInfo();
  if (points.empty()) return tracks;

  // DBSCAN tolerance: the geometry's bar pitch as the base, unless set.
  TMS_SpacePointDBScan::Tolerance tolerance = config.DBScanTolerance;
  if (tolerance.BaseTransverseMM <= 0.0) tolerance.BaseTransverseMM = TMS_Geom::GetInstance().GetMaxBarPitch();

  // The follower, with this slice's hits and X/Y time differences.
  TMS_KalmanFollower::Follower follower(config.Follower, field);
  if (!hits.empty()) follower.SetHits(&hits);
  std::map<std::pair<int, int>, double> xyTimeDifference;
  if (config.UseXYTime && !hits.empty()) {
    const int nHits = static_cast<int>(hits.size());
    for (const TMS_SpacePoint &point : points) {
      const int xi = point.GetXHitIndex(), yi = point.GetYHitIndex();
      if (xi < 0 || yi < 0 || xi >= nHits || yi >= nHits) continue;
      // X-bar hit: measures y (Coordinate), at its plane z; Y-bar hit likewise for x.
      double dt = 0.0;
      if (TMS_SpacePointTiming::CorrectedXYTimeDifference(point.GetX(), point.GetY(), hits[xi].Coordinate,
                                                           hits[xi].Z, hits[xi].Time, hits[yi].Coordinate,
                                                           hits[yi].Z, hits[yi].Time, dt))
        xyTimeDifference[{xi, yi}] = dt;
    }
    follower.SetXYTimeDifferenceSource([&xyTimeDifference](const TMS_SpacePoint &point, double &dt) {
      auto it = xyTimeDifference.find({point.GetXHitIndex(), point.GetYHitIndex()});
      if (it == xyTimeDifference.end()) return false;
      dt = it->second;
      return true;
    });
  }

  // DBSCAN; which objects are track-like.
  TMS_SpacePointDBScan dbscan(points, config.DBScanMinPoints, tolerance);
  const std::vector<std::vector<int>> dbscanClusters = dbscan.RunAndGetClusterIndices();

  // The objects to fit: DBSCAN's clusters, or with linking, each chain of
  // linked clusters merged into one object and every unlinked cluster as is.
  std::vector<std::vector<int>> clusters;
  std::vector<std::vector<int>> chains;
  if (config.UseClusterLinking) {
    chains = TMS_ClusterLinker::LinkClusters(points, dbscanClusters, config.Linker).Chains;
    std::vector<char> linked(dbscanClusters.size(), 0);
    for (const std::vector<int> &chain : chains) {
      std::vector<int> merged;
      for (int c : chain) {
        merged.insert(merged.end(), dbscanClusters[c].begin(), dbscanClusters[c].end());
        linked[c] = 1;
      }
      std::sort(merged.begin(), merged.end());
      clusters.push_back(std::move(merged));
    }
    for (std::size_t c = 0; c < dbscanClusters.size(); ++c)
      if (!linked[c]) clusters.push_back(dbscanClusters[c]);
  } else {
    clusters = dbscanClusters;
  }
  if (info) {
    info->DBScanClusters = dbscanClusters;
    info->Chains = chains;
  }
  std::vector<std::size_t> trackLike, other;
  for (std::size_t c = 0; c < clusters.size(); ++c) {
    TMS_SpacePointCluster cluster(points, clusters[c]);
    const bool isTrackLike = cluster.IsTrackLike(config.LinearityThreshold, config.MinClusterSizeForTrack);
    (isTrackLike ? trackLike : other).push_back(c);
    if (info) info->ClusterTrackLike.push_back(isTrackLike);
  }
  if (info) info->Clusters = clusters;
  // Largest first, so a big merged cluster claims its hits before any small
  // fragment beside it can.
  const auto largestFirst = [&](std::size_t a, std::size_t b) { return clusters[a].size() > clusters[b].size(); };
  std::sort(trackLike.begin(), trackLike.end(), largestFirst);
  std::sort(other.begin(), other.end(), largestFirst);

  TMS_IterativeTrackFit::Config split = config.Split;
  split.DBScanMinPoints = config.DBScanMinPoints;
  split.DBScanTolerance = tolerance;
  split.LinearityThreshold = config.LinearityThreshold;
  split.MinClusterSizeForTrack = config.MinClusterSizeForTrack;

  TMS_IterativeTrackFit::ClaimedHits claimed;

  // --- Stage 1: track-like clusters. ---
  for (std::size_t c : trackLike) {
    for (TMS_IterativeTrackFit::Track &fitted : TMS_IterativeTrackFit::FitCluster(points, clusters[c], follower, split, claimed)) {
      Track track;
      track.HitIndices = UsedHits(fitted.Fit, points);
      track.Fit = std::move(fitted.Fit);
      track.ObjectIndices = std::move(fitted.ObjectIndices);
      track.Stage = 1;
      track.ClusterIndex = static_cast<int>(c);
      track.ClusterSize = clusters[c].size();
      track.Iteration = fitted.Iteration;
      track.ClusterFlagged = fitted.ClusterFlagged;
      tracks.push_back(std::move(track));
    }
  }

  // --- Stage 2: graph search inside clusters that are not track-like. ---
  if (config.UseGraphSearch) {
    const TMS_GraphTrackFinder::Finder finder(config.Graph);
    for (std::size_t c : other) {
      // The cluster's unclaimed points.
      std::vector<int> free;
      for (int idx : clusters[c])
        if (!claimed.Uses(points[idx])) free.push_back(idx);
      if (free.size() < config.MinGraphClusterSize) continue;
      std::vector<TMS_SpacePoint> subset;
      subset.reserve(free.size());
      for (int idx : free) subset.push_back(points[idx]);

      const TMS_GraphTrackFinder::Result found = finder.Find(subset);
      for (const TMS_GraphTrackFinder::Path &path : found.Paths) {  // best first
        if (path.DistinctLayers < config.MinGraphPathLayers) continue;
        std::vector<int> object;
        for (std::size_t local : path.SpacePointIndices) {
          const int idx = free[local];
          if (!claimed.Uses(points[idx])) object.push_back(idx);  // an earlier path may have taken some
        }
        if (object.size() < config.MinGraphPathLayers) continue;
        Track track;
        if (!TMS_IterativeTrackFit::FitObject(points, object, follower, claimed, nullptr, track.Fit)) continue;
        if (TMS_IterativeTrackFit::CountHits(track.Fit) < config.MinGraphTrackHits) continue;
        TMS_IterativeTrackFit::ClaimHits(points, track.Fit, claimed);
        track.HitIndices = UsedHits(track.Fit, points);
        track.ObjectIndices.assign(object.begin(), object.end());
        track.Stage = 2;
        track.ClusterIndex = static_cast<int>(c);
        track.ClusterSize = clusters[c].size();
        tracks.push_back(std::move(track));
      }
    }
  }
  return tracks;
}

}  // namespace TMS_Cluster3DReco
