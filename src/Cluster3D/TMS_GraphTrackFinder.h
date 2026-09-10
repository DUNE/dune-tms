#ifndef _TMS_GRAPHTRACKFINDER_H_SEEN_
#define _TMS_GRAPHTRACKFINDER_H_SEEN_

#include "TMS_SpacePoint.h"

#include <cstddef>
#include <vector>

// A bounded, experimental 3D pattern-recognition stage.  It consumes a DBSCAN
// cluster of TMS_SpacePoints and returns candidate paths as indices into that
// input vector.  The returned indices preserve the X/Y hit provenance already
// stored by TMS_SpacePoint, so a later Kalman follower can recover native 2D
// measurements without rematching them.
namespace TMS_GraphTrackFinder {

struct Config {
  // Space points whose z positions differ by less than this are one graph layer.
  double LayerZTolerance = 1.0;       // mm
  int MaxLayerGap = 3;
  // 1.25 (~51 deg from the z-axis) was an uncalibrated placeholder. Checked
  // against the real ND-LAr-fiducial-origin muon population (n=31, the
  // actual target sample -- see muons.csv): median angle 12.5 deg, 90th
  // percentile 20.2 deg, max observed 33.3 deg. 0.8 (~39 deg) keeps a real
  // safety margin over that observed max (for scattering and the small
  // sample) while being far tighter than the old placeholder -- in a dense
  // shower layer this matters a lot: cone half-angle sets the *area* (not
  // linear extent) of admitted ghost combinations, so 51->39 deg cuts
  // admitted junk well more than the angle numbers alone suggest.
  double MaxAbsDXDZ = 0.8;
  double MaxAbsDYDZ = 0.8;
  double MaxTimeDifference = 40.0;    // ns; negative disables the gate

  // Hard combinatoric bounds.
  std::size_t MaxLinksPerTargetLayer = 16;
  std::size_t SeedLength = 4;
  std::size_t MaxSeedFrontier = 1024;
  std::size_t MaxSeeds = 96;
  double SeedOverlapFraction = 0.50;
  std::size_t BeamWidth = 48;
  std::size_t MaxPathsPerSeed = 8;
  std::size_t MaxOutputPaths = 32;
  std::size_t MaxHypotheses = 3000000;

  // Seeds are not allowed in shower-like layers or on highly ambiguous ghost
  // points. Growth (an already-established, >=2-point trajectory) handles
  // these layers differently instead of just being let through the same
  // gate: the static per-edge graph (MaxLinksPerTargetLayer, ranked by a
  // cost that knows nothing about this hypothesis's own curve) is skipped
  // for any layer whose occupancy exceeds this same threshold, and that
  // layer's points are searched directly, ranked purely by consistency
  // with the established trajectory (see DenseLayerCandidates below) --
  // this is what actually lets a real track punch through a shower core
  // instead of stopping at its edge.
  std::size_t MaxSeedLayerOccupancy = 24;
  std::size_t MaxSeedHitMultiplicity = 8;
  std::size_t MinPathPoints = 8;

  // Master switch for the dense-layer bypass described above. Default
  // false: validated on the full 26-case real-data benchmark (2026-09-09)
  // as a net regression relative to the occupancy-penalty fix alone --
  // 80.4% vs 86.4% mean plane-coverage, 2 new complete misses -- so it
  // should not be silently active just because MaxSeedLayerOccupancy was
  // raised for the (also-recommended) occupancy-penalty fix. Set true only
  // for deliberate, controlled A/B testing of the mechanism itself.
  bool UseDenseLayerSearch = false;

  // How many of a dense layer's points to keep, ranked by trajectory-
  // consistency (kink/curvature-projection cost), when growth bypasses the
  // static graph there. Small on purpose -- this is a final selection
  // after full-layer geometric search, not a coarse pre-filter.
  std::size_t DenseLayerCandidates = 6;
  // Absolute ceiling on that same kink/curvature-projection cost for a
  // dense-layer candidate -- "best of what's in this layer" isn't good
  // enough on its own: without a floor on actual quality, a rich enough
  // ghost population can present a chain of merely-plausible-looking
  // points all the way through a dense region, extending the path onto
  // noise instead of stopping at the trajectory's real edge. This value is
  // in the kink cost's own units (CurvatureProjectionPenaltyX/Y times a
  // deviation-from-projection in mm, squared) -- 80 corresponds to roughly
  // a 40-60mm deviation depending on axis, i.e. a bit more than one bar
  // pitch, enough margin for real scattering without accepting a genuine
  // mismatch.
  double DenseLayerMaxKinkCost = 80.0;

  // Lower score is better.  X curvature is expected in the TMS magnetic field,
  // so changes in dx/dz are penalized less than changes in dy/dz.
  double PointReward = 3.0;
  double GapPenalty = 0.7;
  double SlopePenalty = 0.05;
  double KinkXPenalty = 35.0;
  double KinkYPenalty = 100.0;
  double OccupancyPenalty = 0.08;
  double HitMultiplicityPenalty = 0.35;

  // Once a hypothesis has at least 3 already-accepted points, its transition
  // cost switches from penalizing raw kink (any bend at all, symmetrically)
  // to a curvature-projected model: fit the local bend rate from the last 3
  // accepted points, project forward assuming that rate continues, and
  // penalize a candidate's *deviation from the projection* instead. This is
  // the "handle on curvature" a longer accepted stretch gives you -- a track
  // curving smoothly and consistently (as expected in the TMS field) costs
  // almost nothing to continue, while a candidate that breaks the established
  // trend (a wrong/ghost point) costs a lot, even if its raw local kink looks
  // mild. Falls back to the raw-kink cost below when fewer than 3 points are
  // available yet (start of a seed). Set false to keep the old kink-only
  // behavior for comparison.
  bool UseCurvatureProjection = true;
  double CurvatureProjectionPenaltyX = 0.02;  // per mm^2 of (actual - projected) deviation
  double CurvatureProjectionPenaltyY = 0.05;  // per mm^2 -- Y still weighted stricter than X

  // A space point's coordinate in whichever axis a bar measures is quantized
  // to roughly one bar pitch, not continuous -- a bar's readout is centered
  // per-bar. A real, perfectly smooth trajectory still produces small
  // staircase jumps from this alone as it crosses bar boundaries, especially
  // over short plane-to-plane z gaps where that fixed position step
  // translates into a large apparent slope change. Confirmed on real data
  // (2026-09-08): even the objectively-correct hit sequence for a known
  // muon showed a Y kink cost 3-25x its X kink cost in the upstream half of
  // its own trajectory, from bar-pitch quantization alone, vanishing to
  // exactly zero downstream where the plane spacing is larger -- not real
  // curvature, since Y isn't expected to bend at all. These subtract the
  // expected worst-case quantization-driven noise from the raw kink/
  // projection-deviation before squaring it, so only genuine deviation
  // beyond what quantization can explain gets penalized. Set both to 0 to
  // recover the old, fully quadratic-from-zero behavior.
  double PositionQuantizationX = 36.0;  // mm, ~1 bar pitch
  double PositionQuantizationY = 36.0;  // mm, ~1 bar pitch

  // Paths sharing this fraction of the smaller path are treated as duplicates.
  double DuplicateOverlapFraction = 0.80;
};

struct Diagnostics {
  std::size_t InputPoints = 0;
  std::size_t Layers = 0;
  std::size_t LinksTested = 0;
  std::size_t LinksAccepted = 0;
  std::size_t SeedsGenerated = 0;
  std::size_t SeedsRetained = 0;
  std::size_t HypothesesCreated = 0;
  std::size_t HypothesesPruned = 0;
  std::size_t MaxLiveHypotheses = 0;
  std::size_t NativeHitConflicts = 0;
  std::size_t PathsBeforeDeduplication = 0;
  std::size_t PathsAfterDeduplication = 0;
  bool ResourceLimitReached = false;
};

struct Path {
  std::vector<std::size_t> SpacePointIndices;
  double Score = 0.0;
  std::size_t DistinctLayers = 0;
};

struct Result {
  std::vector<Path> Paths;
  Diagnostics Stats;
};

class Finder {
public:
  explicit Finder(const Config &config = Config());

  // Normally pass one DBSCAN cluster at a time.  No assumption is made about
  // which end is the interaction/entry end; clean stubs may seed anywhere and
  // each retained seed is grown toward both lower and higher z.
  Result Find(const std::vector<TMS_SpacePoint> &spacePoints) const;

private:
  Config fConfig;
};

} // namespace TMS_GraphTrackFinder

#endif
