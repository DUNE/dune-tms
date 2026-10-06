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
  // --- Stitching of sequential pieces (Config::StitchSequentialTracks). ---
  if (config.StitchSequentialTracks && tracks.size() > 1) {
    // A track's end state (its last accepted node, or its single-hit
    // extension's end) and start state (the backward pass's, or its first
    // accepted node), as position + slopes at a z.
    struct Ends {
      bool ok = false;
      double sz = 0, sx = 0, sy = 0, sdx = 0, sdy = 0;  // start
      double ez = 0, ex = 0, ey = 0, edx = 0, edy = 0;  // end
    };
    auto endsOf = [](const TMS_KalmanFollower::FitResult &fit) {
      Ends e;
      const TMS_KalmanFollower::FollowedNode *first = nullptr, *last = nullptr;
      for (const TMS_KalmanFollower::FollowedNode &node : fit.Nodes)
        if (node.HasHit) {
          if (!first) first = &node;
          last = &node;
        }
      if (!first) return e;
      e.ok = true;
      e.sz = first->FilteredZ != 0.0 ? first->FilteredZ : first->Z;
      e.sx = first->FilteredX; e.sy = first->FilteredY; e.sdx = first->FilteredDXDZ; e.sdy = first->FilteredDYDZ;
      if (fit.HasStartState) {
        e.sz = fit.StartZ; e.sx = fit.StartX; e.sy = fit.StartY; e.sdx = fit.StartDXDZ; e.sdy = fit.StartDYDZ;
      }
      e.ez = last->FilteredZ != 0.0 ? last->FilteredZ : last->Z;
      e.ex = last->FilteredX; e.ey = last->FilteredY; e.edx = last->FilteredDXDZ; e.edy = last->FilteredDYDZ;
      if (fit.NExtensionHits > 0) { e.ez = fit.ExtensionEndZ; e.ex = fit.ExtensionEndX; e.ey = fit.ExtensionEndY; }
      return e;
    };
    // Score of B continuing A (lower = better); false if it does not.
    auto stitchScore = [&](const Ends &a, const Ends &b, double &score) {
      if (!a.ok || !b.ok) return false;
      if (b.sz < a.ez - config.StitchMaxOverlapMM || b.sz <= a.sz || b.ez <= a.ez) return false;
      const double gap = std::max(0.0, b.sz - a.ez);
      if (gap > config.StitchMaxGapMM) return false;
      const double tolY = config.StitchMissBaseMM + config.StitchMissPerMeterMM * gap / 1000.0, tolX = 2.0 * tolY;
      const double dz = b.sz - a.ez;
      const double mx1 = (a.ex + a.edx * dz - b.sx) / tolX, my1 = (a.ey + a.edy * dz - b.sy) / tolY;
      const double mx2 = (b.sx - b.sdx * dz - a.ex) / tolX, my2 = (b.sy - b.sdy * dz - a.ey) / tolY;
      if (std::abs(mx1) > 1 || std::abs(my1) > 1 || std::abs(mx2) > 1 || std::abs(my2) > 1) return false;
      const double c = (a.edx * b.sdx + a.edy * b.sdy + 1.0) /
                       std::sqrt((a.edx * a.edx + a.edy * a.edy + 1.0) * (b.sdx * b.sdx + b.sdy * b.sdy + 1.0));
      const double angle = std::acos(std::min(1.0, std::max(-1.0, c)));
      if (angle > config.StitchMaxAngleRad) return false;
      score = mx1 * mx1 + my1 * my1 + mx2 * mx2 + my2 * my2 + (angle / config.StitchMaxAngleRad) * (angle / config.StitchMaxAngleRad);
      return true;
    };
    auto release = [&](const TMS_KalmanFollower::FitResult &fit) {
      for (const TMS_KalmanFollower::FollowedNode &node : fit.Nodes) {
        if (!node.HasHit) continue;
        const TMS_SpacePoint &p = points[node.ChosenSpacePointIndex];
        claimed.X.erase(p.GetXHitIndex());
        claimed.Y.erase(p.GetYHitIndex());
      }
      for (const auto &orphan : fit.Orphans) {
        claimed.X.erase(orphan.HitIndex);
        claimed.Y.erase(orphan.HitIndex);
      }
    };
    bool merged = true;
    std::vector<char> failed(tracks.size() * tracks.size(), 0);  // pairs already tried without success
    while (merged) {
      merged = false;
      // The best-scoring (A, B) pair not yet tried.
      int bestA = -1, bestB = -1;
      double bestScore = 1e30;
      for (std::size_t a = 0; a < tracks.size(); ++a) {
        const Ends ea = endsOf(tracks[a].Fit);
        for (std::size_t b = 0; b < tracks.size(); ++b) {
          if (a == b || failed[a * tracks.size() + b]) continue;
          double score = 0.0;
          if (stitchScore(ea, endsOf(tracks[b].Fit), score) && score < bestScore) {
            bestScore = score;
            bestA = static_cast<int>(a);
            bestB = static_cast<int>(b);
          }
        }
      }
      if (bestA < 0) break;
      Track &A = tracks[bestA];
      Track &B = tracks[bestB];
      std::vector<int> object;
      for (const Track *t : {&A, &B})
        for (const TMS_KalmanFollower::FollowedNode &node : t->Fit.Nodes)
          if (node.HasHit) object.push_back(static_cast<int>(node.ChosenSpacePointIndex));
      std::sort(object.begin(), object.end());
      object.erase(std::unique(object.begin(), object.end()), object.end());
      release(A.Fit);
      release(B.Fit);
      TMS_KalmanFollower::FitResult fit;
      const Ends eb = endsOf(B.Fit);
      bool accept = TMS_IterativeTrackFit::FitObject(points, object, follower, claimed, nullptr, fit);
      if (accept) {
        const Ends em = endsOf(fit);
        accept = em.ok && em.ez >= eb.ez - config.StitchMaxOverlapMM &&
                 TMS_IterativeTrackFit::CountHits(fit) >=
                     std::max(TMS_IterativeTrackFit::CountHits(A.Fit), TMS_IterativeTrackFit::CountHits(B.Fit));
      }
      if (accept) {
        TMS_IterativeTrackFit::ClaimHits(points, fit, claimed);
        A.HitIndices = UsedHits(fit, points);
        A.Fit = std::move(fit);
        tracks.erase(tracks.begin() + bestB);
        failed.assign(tracks.size() * tracks.size(), 0);
        merged = true;
      } else {
        TMS_IterativeTrackFit::ClaimHits(points, A.Fit, claimed);
        TMS_IterativeTrackFit::ClaimHits(points, B.Fit, claimed);
        failed[bestA * tracks.size() + bestB] = 1;
        merged = true;  // try the next-best pair
      }
    }
  }

  // --- Shadow-track absorption (Config::AbsorbShadowTracks). ---
  if (config.AbsorbShadowTracks && !hits.empty() && tracks.size() > 1) {
    // A track's coordinate at z in one view: its fitted node nearest in z,
    // carried straight to z.
    auto coordinateAt = [](const TMS_KalmanFollower::FitResult &fit, double z, bool measuresX, double &out) {
      const TMS_KalmanFollower::FollowedNode *nearest = nullptr;
      double best = 0.0;
      for (const TMS_KalmanFollower::FollowedNode &node : fit.Nodes) {
        if (!node.HasHit) continue;
        const double nz = node.FilteredZ != 0.0 ? node.FilteredZ : node.Z;
        if (!nearest || std::abs(nz - z) < best) {
          nearest = &node;
          best = std::abs(nz - z);
        }
      }
      if (!nearest) return false;
      const double nz = nearest->FilteredZ != 0.0 ? nearest->FilteredZ : nearest->Z;
      out = measuresX ? nearest->FilteredX + nearest->FilteredDXDZ * (z - nz)
                      : nearest->FilteredY + nearest->FilteredDYDZ * (z - nz);
      return true;
    };
    auto zRange = [&](const Track &t, double &lo, double &hi) {
      lo = 1e30;
      hi = -1e30;
      for (int h : t.HitIndices) {
        lo = std::min(lo, hits[h].Z);
        hi = std::max(hi, hits[h].Z);
      }
    };
    std::vector<std::size_t> order(tracks.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    // Biggest first: a track can only be absorbed into a bigger one.
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return tracks[a].HitIndices.size() > tracks[b].HitIndices.size(); });
    std::vector<char> removed(tracks.size(), 0);
    for (std::size_t bi = order.size(); bi-- > 1;) {
      Track &b = tracks[order[bi]];
      if (b.HitIndices.empty()) continue;
      double b0, b1;
      zRange(b, b0, b1);
      for (std::size_t ai = 0; ai < bi; ++ai) {
        if (removed[order[ai]]) continue;
        Track &a = tracks[order[ai]];
        double a0, a1;
        zRange(a, a0, a1);
        if (std::min(a1, b1) - std::max(a0, b0) <= 0.0) continue;  // no z overlap
        std::vector<int> onA;
        for (int h : b.HitIndices) {
          const TMS_KalmanFollower::FitHit &hit = hits[h];
          double c = 0.0;
          if (!coordinateAt(a.Fit, hit.Z, hit.MeasuresX, c)) continue;
          const double pitch = hit.SigmaMM * std::sqrt(12.0);
          if (std::abs(hit.Coordinate - c) <= config.ShadowTolerancePitch * pitch) onA.push_back(h);
        }
        if (onA.size() < config.ShadowHitFraction * b.HitIndices.size()) continue;
        // B is A's shadow: A takes B's on-A hits; B goes.
        a.HitIndices.insert(a.HitIndices.end(), onA.begin(), onA.end());
        std::sort(a.HitIndices.begin(), a.HitIndices.end());
        a.HitIndices.erase(std::unique(a.HitIndices.begin(), a.HitIndices.end()), a.HitIndices.end());
        removed[order[bi]] = 1;
        break;
      }
    }
    std::vector<Track> kept;
    kept.reserve(tracks.size());
    for (std::size_t i = 0; i < tracks.size(); ++i)
      if (!removed[i]) kept.push_back(std::move(tracks[i]));
    tracks.swap(kept);
  }
  return tracks;
}

}  // namespace TMS_Cluster3DReco
