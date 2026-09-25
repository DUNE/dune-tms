#ifndef _TMS_CLUSTER3DRECO_H_SEEN_
#define _TMS_CLUSTER3DRECO_H_SEEN_

#include <cstddef>
#include <vector>

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
//   Stage 1: DBSCAN over the space points; every track-like cluster, largest
//     first, through TMS_IterativeTrackFit::FitCluster() (Kalman fit, and for
//     clusters that look merged, further tracks from the unclaimed remainder).
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

  // Stage 1 per-cluster fitting and splitting. Its DBSCAN fields are
  // overwritten with the ones above.
  TMS_IterativeTrackFit::Config Split;

  // Stage 2.
  bool UseGraphSearch = true;
  // Only clusters with at least this many unclaimed points are searched.
  std::size_t MinGraphClusterSize = 8;
  // A graph path must span at least this many point layers to be fitted.
  std::size_t MinGraphPathLayers = 6;
  // And its fit must accept at least this many points to become a track.
  int MinGraphTrackHits = 6;
  // The validated real-data configuration the truth tools use.
  TMS_GraphTrackFinder::Config Graph = DefaultGraphConfig();

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
  int ClusterIndex = -1;               // into the DBSCAN clusters of this slice
  std::size_t ClusterSize = 0;
  int Iteration = 0;                   // Stage 1: 0 = first fit, >= 1 = split remainder
  bool ClusterFlagged = false;
};

// The slice's hits as Kalman-follower measurements, indexed like the space
// points' hit indices: an X-bar hit measures y, a Y-bar hit measures x, with
// a uniform-bar resolution of barPitchMM / sqrt(12); pedestal-suppressed
// hits and other bar orientations are marked unusable.
std::vector<TMS_KalmanFollower::FitHit> BuildFitHits(const std::vector<TMS_Hit> &hits, double barPitchMM);

// Reconstruct one slice. hits: as BuildFitHits() makes them (or empty, for
// the space-point fit). Needs the TMS geometry loaded (TMS_Geom) for the
// material and the X/Y time transit correction.
std::vector<Track> Run(const std::vector<TMS_SpacePoint> &points,
                       const std::vector<TMS_KalmanFollower::FitHit> &hits, const Config &config,
                       const IFieldModel &field);

}  // namespace TMS_Cluster3DReco

#endif
