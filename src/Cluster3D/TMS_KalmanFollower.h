#ifndef _TMS_KALMANFOLLOWER_H_SEEN_
#define _TMS_KALMANFOLLOWER_H_SEEN_

#include <cstddef>
#include <vector>

#include "TMatrixD.h"

#include "TMS_FieldModel.h"
#include "TMS_SpacePoint.h"

// A Kalman "follower": takes a topologically-found seed path (an ordered,
// z-increasing list of space-point indices, e.g. TMS_GraphTrackFinder::Path
// or a DBSCAN+PCA track-like cluster's own points) and re-walks it
// plane-by-plane with a physics-based filter that (a) fits position,
// direction, momentum and charge, and (b) resolves hit ambiguity at each
// plane -- when several candidate space points ("ghosts", see
// TMS_SpacePointBuilder.h) share a z-layer, picks the one consistent with
// the fit via a proper chi2 gate, instead of the legacy TMS_Kalman's
// behavior of silently keeping only the last hit per z and dropping the
// rest (TMS_Kalman.cpp:73-84).
//
// Unlike the legacy Kalman, this module actually applies magnetic-field
// bending during propagation (see TMS_FieldModel.h) -- legacy computes a
// deflection term and never adds it into the propagated state
// (TMS_Kalman.cpp:252-253).
//
// Forward-pass only for now (no RTS smoother yet -- see the project plan).
namespace TMS_KalmanFollower {

struct Config {
  // A candidate's chi2 (2 DoF: x, y position residual) must be below this
  // to be accepted at all. Tuned empirically 2026-09-15 against the 15-file
  // truth population (app/KalmanFollowerTruthEfficiency): the 2-DoF
  // statistical reference value (~9.21 @ 99% CL) was rejecting real
  // truth-matched hits outright -- 66.7% of all gap layers had the truth
  // point evaluated but chi2-gated out, not missing. Swept 9.21/15/25/40/60;
  // completeness rose monotonically the whole way (64.4%->71.2%) with
  // purity also improving slightly at every step (87.7%->88.5%, never
  // trading off) -- but gains halved at each step (+3.2/+1.8/+1.1/+0.7pp),
  // and by 60 the gate barely constrains anything statistically (ambiguous-
  // layer accuracy stayed flat at ~96.4% the entire sweep, since argmin-chi2
  // already picks correctly among candidates regardless of the gate value --
  // the gate only controls whether a hit is recorded at all). 25 keeps the
  // gate a real threshold while capturing 73% of the measured gain
  // (5.0/6.8pp). See kalman_follower memory for the full sweep table.
  double ChiSquareGateMax = 25.0;

  // Predict() sub-steps a layer-to-layer propagation into pieces no longer
  // than this, rather than one linearized jump. Needed for real numerical
  // stability, not just accuracy: the curvature Jacobian's position<->q/p
  // coupling term is quadratic in the step length, so one big ~100-200mm
  // plane-to-plane jump gives q/p an artificially outsized lever arm on
  // position -- discovered on the flagship real-data case, where a single
  // big jump let one ordinary ~45mm position residual (typical
  // quantization+scattering noise, not a real momentum signal) collapse
  // the fitted momentum to its floor in one update. ~1 bar pitch is a
  // reasonable scale (fine enough that the linearization stays valid,
  // coarse enough not to multiply the number of TMS_Geom::GetMaterials
  // calls unreasonably).
  double MaxSubstepLengthMM = 40.0;

  // How many consecutive layers with no accepted candidate the follower
  // tolerates before giving up (FitResult::Converged = false). Mirrors
  // TMS_GraphTrackFinder::Config::MaxLayerGap's role.
  int MaxConsecutiveGaps = 3;

  // How many layers past the seed path's OWN last point the follower will
  // keep walking. A generous margin (this lets the follower genuinely
  // recover real continuation the seed's finder missed, not just replay
  // it), but bounded: without this, a slice with a lot of unrelated
  // structure elsewhere (a big shower, other tracks) can keep presenting
  // some candidate that passes the chi2 gate at every subsequent layer
  // clear to the end of the slice's z-range, each one a real (slow)
  // TMS_Geom::GetMaterials navigation call -- discovered on the flagship
  // real-data case, which has activity across most of the detector's 82
  // planes even though the target muon only touches 14 of them.
  int MaxLayersBeyondSeed = 15;

  // MUST match whatever TMS_LayerGrouping::Build() tolerance the seed path
  // was grouped with (e.g. TMS_GraphTrackFinder::Config::LayerZTolerance),
  // or "gather every candidate at this layer" can silently mean a
  // different layer than what the seed path assumed at tolerance edges.
  double LayerZTolerance = 1.0;  // mm

  // Initial state covariance (diagonal), used only for the seed path's
  // first node. Position variance is deliberately wide -- standard Kalman
  // practice is to start uncertain and let the filter converge, rather
  // than bias early picks with an overconfident prior.
  double InitialCovXX = 200.0;       // mm^2
  double InitialCovYY = 1.0e3;       // mm^2
  // Slope variance is NOT deliberately wide -- unlike position, the
  // initial slope isn't a blind guess: SeedDirection() fits it from the
  // seed path's own first ~3 real points. TMS_Kalman.cpp:389-393's
  // hardcoded 1.5/2.5 (sigma~1.2-1.6, i.e. tens of degrees of "we have no
  // idea") were tried first here too, copied without checking they still
  // made sense for THIS module's different Update()/Jacobian structure --
  // combined with a comparatively tight InitialCovXX, that mismatch let
  // the Kalman gain badly over-correct direction on the very first real
  // position update (observed on the flagship real-data case: the fit
  // diverged after only 2 nodes, reproducing identically with the field
  // model swapped to zero, which is what isolated this from the separate
  // q/p-scale bug documented below). 0.1 rad-ish (sigma~0.1) is still
  // generous next to this project's real observed muon angular spread
  // (median ~12 deg, 90th percentile ~20 deg, max observed ~33 deg, i.e.
  // dxdz/dydz well under 0.6) without being so wide it swamps a real
  // 3-point seed fit's own information.
  double InitialCovDXDZDXDZ = 0.01;
  double InitialCovDYDZDYDZ = 0.01;
  // (e/MeV)^2. sigma_qp=0.01 matches this project's real observed muon
  // population's low-momentum end (q/p down to ~1/100 MeV for a ~100 MeV
  // muon) -- wide enough to comfortably cover the whole realistic range up
  // to a few GeV without biasing early picks, but NOT the literal value of
  // 1.0 first tried here: that's sigma_qp=1 e/MeV, i.e. "momentum could be
  // 1 MeV at 1-sigma" -- an unphysically wide prior that, combined with
  // the field-bending Jacobian's qp<->position coupling (Predict()'s
  // transfer(0,4)/transfer(2,4)), let a single real position measurement's
  // Kalman gain swing q/p by orders of magnitude on the very first real
  // update (observed: fitted momentum collapsed to ~0.0002 MeV on node 2
  // of the flagship real-data validation run before this fix).
  double InitialCovQPQP = 1.0e-4;

  // Momentum prior for the first node, used only to seed q/p before any
  // real fitting has happened. No range-based or other data-driven seed
  // exists anywhere in this repo (legacy's GetKEEstimateFromLength() is
  // dead in practice, see TMS_Kalman.cpp:15's hardcoded ForwardFitting =
  // false); a fixed, loosely-covaried prior is the standard fallback.
  double InitialMomentumSeedMeV = 1000.0;

  // Measurement uncertainty (mm) in a space point's own not-Z coordinate.
  // A per-hit lookup (TMS_Bar::GetNotZw()) would be more precise, but
  // needs a real, geometry-backed TMS_Hit for every candidate -- this
  // follower deliberately works from TMS_SpacePoint alone (see the
  // TMS_SpacePointBuilder ghosting comment above), so it uses the same
  // ~1-bar-pitch scale already validated empirically for this exact
  // purpose (TMS_GraphTrackFinder::Config::PositionQuantizationX/Y).
  double AssumedBarPitchMM = 36.0;
};

// One followed plane: which candidate (if any) was chosen, the filtered
// state there, and every candidate's chi2 (not just the chosen one) so
// validation tooling can measure ambiguity-resolution accuracy without
// re-running the fit.
struct FollowedNode {
  std::size_t Layer = 0;
  double Z = 0.0;
  bool HasHit = false;  // false = gap: no candidate passed the chi2 gate

  // Indices into the SAME allSpacePoints vector passed to Follower::Run(),
  // i.e. every space point seen at this layer, not just the seed's pick.
  std::vector<std::size_t> CandidateIndices;
  std::vector<double> CandidateChi2;  // parallel to CandidateIndices
  std::size_t ChosenSpacePointIndex = 0;  // valid only if HasHit
  double Chi2AtChosen = 0.0;

  double FilteredX = 0.0;
  double FilteredY = 0.0;
  double FilteredDXDZ = 0.0;
  double FilteredDYDZ = 0.0;
  double FilteredQP = 0.0;  // charge[e] / momentum[MeV/c]
  TMatrixD FilteredCovariance{5, 5};
};

struct FitResult {
  // Why the walk actually ended -- added to distinguish "ran out of search
  // budget" (GapLimit/RangeEnd, tunable via Config) from "the numerical
  // guards added during Phase 1's momentum-collapse debugging kicked in"
  // (Diverged, not a search-budget question at all). See kalman_follower
  // memory, "investigate the completeness ceiling" (2026-09-15).
  enum class StopReason { NotStarted, ReachedRangeEnd, GapLimitExceeded, Diverged };
  StopReason Stop = StopReason::NotStarted;

  bool Converged = false;
  std::vector<FollowedNode> Nodes;  // one per z-layer walked, low->high z

  double MomentumMeV = 0.0;  // from the final node's filtered q/p
  double Charge = 0.0;       // sign of the final node's filtered q/p
  double TotalChi2 = 0.0;
  int NDoF = 0;
  int NGapsFilled = 0;             // layers skipped for lack of a good candidate
  int NAmbiguousLayersResolved = 0;  // layers where >1 candidate existed
};

class Follower {
  public:
    Follower(const Config &config, const IFieldModel &field);

    // allSpacePoints: the FULL pool the seed was drawn from -- ambiguity
    //   resolution needs to see ghosts the seed path didn't pick.
    // seedPath: indices into allSpacePoints, low->high z (as produced by
    //   TMS_GraphTrackFinder::Path::SpacePointIndices, or a DBSCAN
    //   track-like cluster's own z-ordering). Only its first few points
    //   are used, to seed the initial position/direction -- everything
    //   after that is re-derived layer by layer from the full pool.
    FitResult Run(const std::vector<TMS_SpacePoint> &allSpacePoints,
                  const std::vector<std::size_t> &seedPath) const;

    // Multi-hypothesis seeding ("combinatorial Kalman filter" seeding, the
    // standard ATLAS/CMS/ACTS pattern for exactly this problem): for a
    // found track-like object whose points came from DBSCAN+PCA or a
    // merge-and-re-PCA (i.e. NOT already ordered by a directed search the
    // way TMS_GraphTrackFinder::Path is), z-sorting the object's points and
    // always starting from whichever one lands first can pick a bad anchor
    // when the object's own first z-layer has more than one point at
    // (near-)identical z -- discovered on real slices: one such object's
    // naive seed diverged after a single node, while a hand-picked
    // different first-layer point on the SAME object converged cleanly.
    // Run() itself can't distinguish these (it only ever sees one seedPath),
    // so this spawns one hypothesis per candidate at the object's own first
    // z-layer, fits each with Run(), and keeps the best -- reusing
    // Converged/NDoF/TotalChi2 as the ready-made selection signal rather
    // than inventing a new search. Prefer this over Run() for DBSCAN-direct
    // and merged-cluster seeds; TMS_GraphTrackFinder::Path seeds already
    // went through a directed graph search that resolved this same
    // ambiguity, so they should keep calling Run() directly.
    //
    // objectIndices: the found object's own point indices into
    // allSpacePoints, in ANY order (unlike seedPath above, this is not
    // expected to be pre-sorted).
    FitResult RunBestSeed(const std::vector<TMS_SpacePoint> &allSpacePoints,
                          const std::vector<std::size_t> &objectIndices) const;

  private:
    Config fConfig;
    const IFieldModel &fField;
};

}  // namespace TMS_KalmanFollower

#endif
