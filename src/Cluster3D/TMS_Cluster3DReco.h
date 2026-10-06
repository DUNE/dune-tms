#ifndef _TMS_CLUSTER3DRECO_H_SEEN_
#define _TMS_CLUSTER3DRECO_H_SEEN_

#include <cstddef>
#include <vector>

#include "TMS_ClusterLinker.h"
#include "TMS_FieldModel.h"
#include "TMS_GraphTrackFinder.h"
#include "TMS_IterativeTrackFit.h"
#include "TMS_KalmanFollower.h"
#include "TMS_SpacePoint.h"
#include "TMS_SpacePointDBScan.h"

class TMS_Hit;

// The Cluster3D reconstruction of one slice, truth-free, as one call: the
// production counterpart of what the app/cluster3D truth tools orchestrate
// themselves.
//
//   Stage 1: DBSCAN over the space points; optionally (UseClusterLinking)
//     chains of clusters that are pieces of one track merged into one object
//     (TMS_ClusterLinker); every track-like object, largest first, through
//     TMS_IterativeTrackFit::FitCluster() (Kalman fit, and for objects that
//     look merged, further tracks from the unclaimed remainder).
//   Stage 2: for every cluster that is NOT track-like -- a muon inside hadron
//     activity -- the graph finder on the cluster's unclaimed points, and each
//     accepted path fitted over the slice's unclaimed points.
//
// Tracks claim their hits (TMS_IterativeTrackFit::ClaimHits) as they are
// accepted, so no two tracks share a hit. The Kalman follower is given the
// slice's hits (hit-level fit, orphan pickup) and the transit-corrected X/Y
// time differences, as the truth tools give it.
namespace TMS_Cluster3DReco {

struct Config {
  // DBSCAN and track-like selection (as every truth tool uses them).
  unsigned int DBScanMinPoints = TMS_SpacePointDBScan::kDefaultMinPoints;
  // BaseTransverseMM left at 0 = use the geometry's bar pitch (Run()).
  TMS_SpacePointDBScan::Tolerance DBScanTolerance = ZeroBase();
  double LinearityThreshold = 0.8;
  std::size_t MinClusterSizeForTrack = TMS_SpacePointCluster::kDefaultMinTrackSize;

  // Link DBSCAN clusters that are pieces of one track before fitting, and fit
  // each chain as one object. Off by default until validated (2026-09-26).
  bool UseClusterLinking = false;
  TMS_ClusterLinker::Config Linker;

  // Stage 1 per-cluster fitting and splitting. Its DBSCAN fields are
  // overwritten with the ones above.
  TMS_IterativeTrackFit::Config Split;

  // Stage 2. Off by default: as first designed (2026-09-25, file 7) it made
  // 20 tracks of which 3 were muons', and found 2 more muons.
  bool UseGraphSearch = false;
  // Only clusters with at least this many unclaimed points are searched.
  std::size_t MinGraphClusterSize = 8;
  // A graph path must span at least this many point layers to be fitted.
  std::size_t MinGraphPathLayers = 6;
  // And its fit must accept at least this many points to become a track.
  int MinGraphTrackHits = 6;
  // The validated real-data configuration the truth tools use.
  TMS_GraphTrackFinder::Config Graph = DefaultGraphConfig();

  // Stitching of sequential pieces, after all fits: a track B that starts
  // where a track A ends (B downstream, up to StitchMaxOverlapMM of z
  // overlap, a gap of at most StitchMaxGapMM), A's end extrapolating onto
  // B's start within StitchMissBaseMM + StitchMissPerMeterMM per meter of gap
  // (x tolerance x2: the bending plane) and their directions within
  // StitchMaxAngleRad, is refitted with A as one object; the merged fit
  // replaces both if it reaches B's end. Motivation (2026-09-27, 15 files):
  // ~80 of the extra tracks owned by an already-found muon are sequential
  // pieces (median gap 260 mm). The refit is the real test (it must reach
  // B's end with at least as many hits as either piece), so the geometric
  // pre-test is loose. Default on: with shadow absorption, 15 files, the
  // stitching adds 2 tracks ending correctly and removes 27 duplicates; most
  // sequential pieces fail the refit (their end directions disagree by ~0.3
  // rad) -- they are the muon itself, not descendants.
  bool StitchSequentialTracks = true;
  double StitchMaxGapMM = 1500.0;
  double StitchMaxOverlapMM = 150.0;
  double StitchMissBaseMM = 300.0;
  double StitchMissPerMeterMM = 300.0;
  double StitchMaxAngleRad = 0.6;

  // Shadow-track absorption, after all fits (needs the slice's hits). A
  // smaller track B that overlaps a bigger track A in z, with at least
  // ShadowHitFraction of its hits lying on A's fitted trajectory (within
  // ShadowTolerancePitch bar pitches of A's coordinate at the hit's plane, in
  // the view the hit measures), is A's shadow: its on-A hits join A and B is
  // dropped. Motivation (2026-09-27, 15 files): of 570 extra tracks owned by
  // a muon that already had a track, 484 overlap the main track in z and 302
  // of those hold the muon's hits in one view only -- a ghost built from the
  // muon's leftover hits in one view and foreign hits in the other.
  //
  // Default on at 0.5 (15 files, suite ND-physics muons 616): duplicate
  // tracks 621 -> 412, tracks ending correctly 492 -> 501, junk -19; all
  // muons found -9 of ~11600 (none from ND-LAr) -- near-parallel genuine
  // muons that match another track in one view. 0.4 / 0.3: duplicates 344 /
  // 268 but -28 / -54 muons found; a geometric test can't tell those apart.
  bool AbsorbShadowTracks = true;
  double ShadowHitFraction = 0.5;
  double ShadowTolerancePitch = 1.5;

  // Kalman follower settings (library defaults).
  TMS_KalmanFollower::Config Follower;
  // Give the follower the transit-corrected X/Y hit-time differences
  // (Config::UseXYTimeInSelection needs them).
  bool UseXYTime = true;

  static TMS_SpacePointDBScan::Tolerance ZeroBase() {
    TMS_SpacePointDBScan::Tolerance tolerance;
    tolerance.BaseTransverseMM = 0.0;
    return tolerance;
  }
  static TMS_GraphTrackFinder::Config DefaultGraphConfig() {
    TMS_GraphTrackFinder::Config graph;
    graph.MaxSeedLayerOccupancy = 150;
    graph.MaxSeedHitMultiplicity = 50;
    graph.OccupancyPenalty = 0.0;
    graph.HitMultiplicityPenalty = 0.0;
    graph.UseCurvatureProjection = false;
    return graph;
  }
};

struct Track {
  TMS_KalmanFollower::FitResult Fit;   // indices remapped to the slice's points
  std::vector<int> HitIndices;         // every hit the fit used (applied + orphans), into the slice's hit list
  std::vector<std::size_t> ObjectIndices;  // the object the fit was seeded from
  int Stage = 1;                       // 1 = track-like cluster, 2 = graph path in a non-track-like cluster
  int ClusterIndex = -1;               // into this slice's objects (RunInfo::Clusters)
  std::size_t ClusterSize = 0;
  int Iteration = 0;                   // Stage 1: 0 = first fit, >= 1 = split remainder
  bool ClusterFlagged = false;
};

// The slice's hits as Kalman-follower measurements, indexed like the space
// points' hit indices: an X-bar hit measures y, a Y-bar hit measures x, with
// a uniform-bar resolution of barPitchMM / sqrt(12); pedestal-suppressed
// hits and other bar orientations are marked unusable.
std::vector<TMS_KalmanFollower::FitHit> BuildFitHits(const std::vector<TMS_Hit> &hits, double barPitchMM);

// What Run() saw on the way, for diagnostics.
struct RunInfo {
  // The objects that were fitted (point indices; Track::ClusterIndex indexes
  // these) and which are track-like: the DBSCAN clusters, or with
  // UseClusterLinking, each chain of linked clusters as one object and every
  // unlinked cluster as itself.
  std::vector<std::vector<int>> Clusters;
  std::vector<bool> ClusterTrackLike;
  // DBSCAN's own clusters, and the chains the linker made (indices into them).
  std::vector<std::vector<int>> DBScanClusters;
  std::vector<std::vector<int>> Chains;
};

// Reconstruct one slice. hits: as BuildFitHits() makes them (or empty, for
// the space-point fit). Needs the TMS geometry loaded (TMS_Geom) for the
// material and the X/Y time transit correction. info, if given, is filled.
std::vector<Track> Run(const std::vector<TMS_SpacePoint> &points,
                       const std::vector<TMS_KalmanFollower::FitHit> &hits, const Config &config,
                       const IFieldModel &field, RunInfo *info = nullptr);

}  // namespace TMS_Cluster3DReco

#endif
