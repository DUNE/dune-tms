#ifndef _TMS_ITERATIVETRACKFIT_H_SEEN_
#define _TMS_ITERATIVETRACKFIT_H_SEEN_

#include <cstddef>
#include <unordered_set>
#include <vector>

#include "TMS_KalmanFollower.h"
#include "TMS_SpacePointDBScan.h"
#include "TMS_SpacePointCluster.h"
#include "TMS_SpacePoint.h"

// Iterative "fit, claim, refit the remainder" track extraction from one
// DBSCAN cluster -- the split step for clusters that merge more than one
// real particle.
//
// Motivation (2026-09-24, reports/2026-09-24_caseH_timing_pca/): DBSCAN can
// merge two real muons a few bar pitches apart into one cluster, and PCA
// still calls it track-like (linearity ~0.998 -- (l1-l2)/l1 is dominated by
// length, not width). The Kalman follower then fits ONE of the two, and the
// other is never found. Goal: find and fit every track in the cluster.
//
// Algorithm, per cluster (no truth used anywhere):
//   1. Fit the cluster's own points with Follower::RunBestSeed(), over a pool
//      of the whole slice minus any space point that uses an already-claimed
//      X or Y hit (so the fit can still extend past the cluster's edge, as
//      it does today, but never re-uses another track's hits).
//   2. Claim the X and Y hits of every point the fit chose. Removing every
//      space point that uses a claimed hit also removes the ghosts built from
//      them (a ghost pairs this track's hit with another particle's), which
//      leaves the other particle's own points behind.
//   3. If the cluster is allowed to split (see Config::Mode), re-run DBSCAN on
//      the cluster's remaining points, take the largest track-like
//      sub-cluster, and go back to 1 with it -- up to MaxTracksPerCluster.
//
// Works for any cause of the merge -- same-interaction pairs included, which
// timing alone cannot separate -- because it relies on the first fit
// following one particle, not on the two being separable a priori.
//
// Claimed hits are shared across every cluster of a slice (the caller owns
// the sets), so a later cluster cannot re-use hits a previous one claimed.
namespace TMS_IterativeTrackFit {

struct Config {
  // When to try extracting more than one track from a cluster.
  //   Off:     one fit per cluster (today's behavior).
  //   Flagged: only clusters that look like a merge (IsFlaggedAsMerged).
  //   Always:  every cluster; the remainder test alone decides.
  enum class SplitMode { Off, Flagged, Always };
  SplitMode Mode = SplitMode::Flagged;
  int MaxTracksPerCluster = 3;

  // A fit claims hits (and counts as a track) only with at least this many
  // accepted hits. Stops a remainder of stray ghosts or a delta ray from
  // turning into a "track".
  int MinHitsPerTrack = 4;
  // Stricter minimum for tracks fitted from a split REMAINDER (iteration >=
  // 1). File 13, 2026-09-24: with the plain 4-hit minimum, 30 of 40
  // remainder tracks were junk (15 ghost-only, 10 duplicates of an
  // already-found muon, 5 non-muon), almost all with <= 6 hits; 7 of the 10
  // genuinely new muons had >= 8.
  int MinHitsPerSplitTrack = 8;
  // If true, a split-remainder fit (iteration >= 1) walks only the
  // cluster's own unclaimed points, not the whole slice. The first fit
  // always walks the whole slice (as the follower does today). 15-file run,
  // 2026-09-24: with whole-slice remainder fits, 19 muons found with
  // splitting off were lost with it on -- every one in a slice where a
  // remainder track had been made, i.e. its hits were claimed by that track.
  bool RestrictSplitFitToCluster = true;

  // Merge flag: more than WideLayerFraction of the cluster's point layers
  // (those with >= 2 points) are wider than WideSpanMM in x or y, or the
  // cluster's space-point time RMS exceeds TimeRMSNs. 2026-09-24, track-like
  // clusters spanning >= 15 layers, 14 files: flags 90% of merged multi-muon
  // clusters and 18.6% of single-muon ones. Calibrated on BothNeighbors
  // points (~2.4 genuine points per muon per layer); NearestY layers hold
  // ~1.3, so these thresholds need re-checking on the new points.
  double WideSpanMM = 150.0;
  double WideLayerFraction = 0.3;
  double TimeRMSNs = 10.0;

  // DBSCAN + PCA settings for re-clustering a remainder. Must match the ones
  // the caller used for the original clustering, so "track-like" means the
  // same thing on both passes.
  unsigned int DBScanMinPoints = TMS_SpacePointDBScan::kDefaultMinPoints;
  // Required: set from TMS_SpacePointDBScan::DefaultTolerance(bar pitch) or
  // whatever the caller clustered with.
  TMS_SpacePointDBScan::Tolerance DBScanTolerance;
  double LinearityThreshold = 0.8;
  std::size_t MinClusterSizeForTrack = TMS_SpacePointCluster::kDefaultMinTrackSize;
};

// One extracted track.
struct Track {
  TMS_KalmanFollower::FitResult Fit;
  // Node's chosen points and candidates, remapped from the fit's filtered
  // pool back to indices into the slice's own point vector -- every index
  // in Fit.Nodes (ChosenSpacePointIndex, CandidateIndices) is already
  // remapped this way.
  std::vector<std::size_t> ObjectIndices;  // the (sub-)cluster this fit was seeded from
  int Iteration = 0;          // 0 = the cluster's first fit, >= 1 = from a split remainder
  bool ClusterFlagged = false;
};

// Claimed hit indices across one slice (X-view and Y-view hits are indexed
// separately by TMS_SpacePoint, so they are kept apart).
struct ClaimedHits {
  std::unordered_set<int> X;
  std::unordered_set<int> Y;
  bool Uses(const TMS_SpacePoint &p) const {
    return X.count(p.GetXHitIndex()) > 0 || Y.count(p.GetYHitIndex()) > 0;
  }
};

// Fit one object (indices into slicePoints) over the slice minus every point
// that uses a claimed hit, with Follower::RunBestSeed(); every index in the
// result is remapped to slicePoints. allowed (optional): only these slice
// points may enter the pool. False if fewer than two object points remain.
bool FitObject(const std::vector<TMS_SpacePoint> &slicePoints, const std::vector<int> &objectIndices,
               const TMS_KalmanFollower::Follower &follower, const ClaimedHits &claimed,
               const std::vector<char> *allowed, TMS_KalmanFollower::FitResult &fitOut);

// Nodes with an accepted point.
int CountHits(const TMS_KalmanFollower::FitResult &fit);

// Claim a fitted track's hits: both hits of every chosen point, and any
// orphan hits it picked up.
void ClaimHits(const std::vector<TMS_SpacePoint> &slicePoints, const TMS_KalmanFollower::FitResult &fit,
               ClaimedHits &claimed);

bool IsFlaggedAsMerged(const std::vector<TMS_SpacePoint> &slicePoints, const std::vector<int> &clusterIndices,
                       const Config &config);

// Extract up to Config::MaxTracksPerCluster tracks from one cluster.
// claimed is updated in place with every accepted track's hits.
std::vector<Track> FitCluster(const std::vector<TMS_SpacePoint> &slicePoints,
                              const std::vector<int> &clusterIndices, const TMS_KalmanFollower::Follower &follower,
                              const Config &config, ClaimedHits &claimed);

}  // namespace TMS_IterativeTrackFit

#endif
