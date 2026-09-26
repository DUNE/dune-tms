#ifndef _TMS_CLUSTERLINKER_H_SEEN_
#define _TMS_CLUSTERLINKER_H_SEEN_

#include <cstddef>
#include <vector>

#include "TMS_SpacePoint.h"

// Links the DBSCAN clusters of one slice that are pieces of the same track,
// before anything is fitted -- the truth-free counterpart of the truth tools'
// "merged_pca" step (which used truth to pick every cluster a muon touched).
//
// Motivation (2026-09-26, reports/2026-09-26_cluster_link_ceiling/): of the
// validation suite's ND-physics muons, 60 of 616 had their best track end more
// than 30 cm early with the rest of the muon in OTHER clusters (or DBSCAN
// noise), and 34 more ended correctly but as two or more tracks. 52 of those
// 60 stopped because the fit "ranged out": the follower seeds its momentum
// from the seed object's z-extent, a short first piece gives a seed ~0.6 of
// the true momentum, and energy loss runs it to the floor while the muon goes
// on. A cluster that spans the whole muon gives the right seed, fits once
// (one momentum, one charge), and claims its hits before any fragment can.
//
// The link graph: nodes are clusters; an edge A -> B (A upstream) is drawn
// when B carries on where A stops:
//   - ordered in z: B starts after A ends, up to MaxOverlapMM of overlap, and
//     B both starts and ends downstream of where A does -- so two clusters side
//     by side (overlapping in z) never link;
//   - a gap of at most MaxGapMM;
//   - colinear: each end's line, fitted to its last (first) EndLayers point
//     layers, extrapolated across the gap lands on the other end within
//     MissBaseMM + MissPerMeterMM per meter of gap (x tolerance scaled by
//     BendScaleX: the field bends tracks in x), and the two end directions
//     agree within MaxAngleRad;
//   - in time: the two ends' mean point times agree within MaxTimeDiffNs.
// A cluster too small for an end direction (< MinDirectionPoints points at
// that end) can link only to one that has a direction; an end whose points
// are not line-like (end-segment PCA linearity < MinEndLinearity: a shower)
// does not link at all.
//
// Each cluster keeps at most one downstream and one upstream link, and only
// mutual best ones (A's best downstream partner is B and B's best upstream
// partner is A), so a muon and a hadron fragment that both line up with the
// same piece cannot both take it. Linked clusters form chains, upstream first.
namespace TMS_ClusterLinker {

struct Config {
  double MaxGapMM = 1000.0;       // as TMS_KalmanFollower::Config::MaxGapMM
  double MaxOverlapMM = 150.0;    // ~one back-section plane pair
  int EndLayers = 4;              // point layers used for each end's line
  std::size_t MinDirectionPoints = 3;
  double MinEndLinearity = 0.7;
  double MissBaseMM = 100.0;      // ~3 bar pitches
  double MissPerMeterMM = 150.0;  // scattering and bending over the gap
  double BendScaleX = 2.0;
  double MaxAngleRad = 0.35;
  double MaxTimeDiffNs = 20.0;
};

struct Link {
  int From = -1;  // upstream cluster
  int To = -1;    // downstream cluster
  double Score = 0.0;  // lower = better; sum of squared normalized misses
};

struct Result {
  std::vector<Link> Links;                 // the accepted (mutual best) links
  std::vector<std::vector<int>> Chains;    // cluster indices, upstream first; only chains of >= 2
};

// clusters: point indices per cluster (TMS_SpacePointDBScan's output).
Result LinkClusters(const std::vector<TMS_SpacePoint> &points, const std::vector<std::vector<int>> &clusters,
                    const Config &config);

}  // namespace TMS_ClusterLinker

#endif
