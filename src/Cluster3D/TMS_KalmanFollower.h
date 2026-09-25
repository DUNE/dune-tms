#ifndef _TMS_KALMANFOLLOWER_H_SEEN_
#define _TMS_KALMANFOLLOWER_H_SEEN_

#include <cstddef>
#include <functional>
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
  // truth population (app/cluster3D/KalmanFollowerTruthEfficiency): the 2-DoF
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

  // If > 0, sigma(q/p) of the first node is this fraction of the seed
  // |q/p| instead of the fixed sqrt(InitialCovQPQP). The fixed value is 10x
  // the seed q/p at 1000 MeV, i.e. an almost uninformative prior: ordinary
  // position noise then drags q/p to a few hundred MeV within one or two
  // layers (2026-09-21: a 2.4 GeV muon's fit went 2131 -> 231 MeV on a
  // chi2 of 3.7), after which real energy loss ranges the fit out while the
  // true track continues.
  //
  // Default 1.0 (with RangeSeedMargin below): 15-file sweep 2026-09-21
  // (reports/2026-09-21_kalman_prior_sweep/): completeness 75.2% -> 82.6%,
  // purity 86.4% -> 87.0%, ND-LAr-fiducial completeness 83.2% -> 92.4%.
  // A tight sigma on the plain 1000 MeV seed (0.3x) was WORSE (67.7% on file
  // 9) -- the tight prior is only safe once the seed itself is sensible.
  double InitialQPRelSigma = 1.0;

  // If > 0, seed the momentum with max(InitialMomentumSeedMeV, this *
  // p_range), where p_range is the smallest momentum that could carry a
  // muon across the seed object's own z-extent through the real material
  // budget (Bethe-Bloch energy loss walked backwards from the momentum
  // floor). A track that visibly spans N steel layers cannot have less
  // momentum than that; it is a lower bound, computed from reconstructed
  // quantities only.
  //
  // Default 1.5: the range is a lower bound and a candidate can be
  // truncated, so a margin above 1 compensates. File-9 sweep of the margin
  // (sigma(q/p)=1x): completeness 81.8 / 82.3 / 82.6 / 82.5% for margin
  // 1.0 / 1.5 / 2.0 / 3.0, but seeds more than 1.5x too high rise 9 / 17 /
  // 44 / 58%; 1.5 gives a median seed ~1.0x the true momentum.
  double RangeSeedMargin = 1.5;

  // End the walk when energy loss carries the fit to the momentum floor
  // (StopReason::RangedOut, counted as converged). The fitted momentum is
  // only as good as its seed, so this can fire while the true muon carries
  // on: with the defaults above 30% of ranged-out fits still have truth
  // planes beyond the stop. Set false to keep walking with the floor state.
  bool StopOnRangeOut = true;

  // RunBestSeed() also tries seeds that skip the object's first 1..MaxHeadSkip
  // layers, keeping whichever hypothesis IsBetterFit prefers. Guards against
  // an object whose head belongs to a different particle (see RunBestSeed).
  // 0 = only the original first-layer anchors.
  //
  // Default 2: near a vertex the first layers mix the muon with same-vertex
  // hadrons, so a seed built there can lock the fit onto the wrong particle
  // (case E: 0/16 -> 16/16 target planes). 15-file sweep 2026-09-21
  // (reports/2026-09-21_kalman_prior_sweep/), skip 0 / 1 / 2: completeness
  // 82.6 / 85.8 / 86.6%, purity 87.0 / 89.0 / 89.7%. Cost: up to 3x the fits
  // per muon (runtime not yet optimised), and a tail of 21 short tracks
  // (0.15%) lose >= 50 pp completeness because IsBetterFit counts hits
  // without checking they belong to one particle.
  int MaxHeadSkip = 2;

  // RunBestSeed() ranks its hypotheses (IsBetterFit) by: [converged, only if
  // this is true], then most hits, then lowest chi2/ndof. Default false.
  // With RangedOut and the gap-limit stop, "converged" stopped being a
  // quality signal: a prematurely ranged-out fit (counted converged) could
  // beat a longer fit that merely ended at the gap limit, and head-skip
  // hypotheses that filled gaps with high-chi2 wrong hits could win. Dropping
  // it (2026-09-21 hypothesis study, reports/2026-09-21_kalman_hypothesis_selection/,
  // 13,249 muons, files 1-7 vs 8-15 agree within 0.2 pp): completeness 88.6 ->
  // 89.2%, purity unchanged, muons losing >= 50 pp vs skip-0 21 -> 14, losing
  // >= 20 pp 115 -> 29. An oracle that sees the truth would reach 91.0% /
  // 95.7%, so the ranking still has headroom that reco-only hit counts and
  // chi2 do not reach.
  bool RankHypothesesByConvergence = false;

  // Measurement uncertainty (mm) in a space point's own not-Z coordinate.
  // A per-hit lookup (TMS_Bar::GetNotZw()) would be more precise, but
  // needs a real, geometry-backed TMS_Hit for every candidate -- this
  // follower deliberately works from TMS_SpacePoint alone (see the
  // TMS_SpacePointBuilder ghosting comment above), so it uses the same
  // ~1-bar-pitch scale already validated empirically for this exact
  // purpose (TMS_GraphTrackFinder::Config::PositionQuantizationX/Y).
  double AssumedBarPitchMM = 36.0;

  // RunBestSeed()'s head-skip hypotheses (above) only ever vary the ANCHOR
  // point (the object's own first-layer candidate); the next ~2 points that
  // SeedDirection() actually averages to get the initial slope are always
  // whichever sort first in global z order -- never explored as
  // alternatives, even when that layer has several candidates. Triplet
  // hypotheses fix that directly: at each skip level, enumerate every
  // (layer0, layer1, layer2) candidate combination, keep only the
  // MaxTripletHypotheses most nearly collinear (TripletCollinearityToleranceMM),
  // and run a full fit on each survivor. The collinearity prune is cheap
  // (arithmetic only); only the survivors pay for a real Kalman walk, so
  // this stays bounded even in a dense slice with hundreds of candidates
  // per layer (see the header comment above on TMS_SpacePointBuilder
  // ghosting for why a single layer can have that many). Purely additive to
  // the existing head-skip hypotheses -- worst case it finds nothing better
  // and IsBetterFit keeps the old winner.
  //
  // Default 5: 15-file sweep 2026-09-23 (reports/2026-09-21_kalman_prior_sweep/,
  // muons_triplets5.csv vs muons_rank_noconv.csv), on top of head-skip=2:
  // completeness 87.10 -> 88.95%, purity 89.68 -> 90.11%, ND-LAr-fiducial
  // completeness 94.53 -> 95.42% (purity also up), TMS-start completeness
  // 80.12 -> 83.20% (the largest single gain). Per muon: 996 better, 111
  // worse (27 lose >= 50pp completeness, almost all short dbscan_direct
  // tracks around 6 target planes -- the same known IsBetterFit short-track
  // weakness head-skip already has, not a new failure mode). Runtime cost
  // measured at ~1.4% (solo file-9 timing, 7m12s -> 7m18s) -- far cheaper
  // than head-skip's ~3x, for a larger completeness gain.
  int MaxTripletHypotheses = 5;

  // Max transverse deviation (mm) of the middle point from the straight
  // line through the first and third, for a (layer0,layer1,layer2)
  // candidate combination to be considered collinear enough to try. A few
  // bar pitches (AssumedBarPitchMM=36mm) -- wide enough to admit a real
  // muon's genuine multiple-scattering kink over 2 layers, tight enough to
  // reject combinations that are obviously not one particle.
  double TripletCollinearityToleranceMM = 100.0;

  // Time as a second discriminant in per-layer candidate selection. When on,
  // the follower keeps a running estimate of the track's time origin t0 --
  // the mean of (t - s/c) over accepted points, s = path length walked from
  // the seed, c = speed of light (muons treated as beta~1) -- and ranks each
  // candidate by position chi2 + time chi2, with time chi2 =
  // r^2 / (sigma_t^2 + sigma_t^2/n), r = (t - s/c) - t0, n = accepted points
  // so far. The acceptance gate itself stays the position-only
  // ChiSquareGateMax, unless TimeGateNSigma > 0 adds a separate |r| cut.
  //
  // Motivation (2026-09-24, reports/2026-09-24_caseH_timing_pca/): muons from
  // DIFFERENT interactions that DBSCAN merges into one cluster reach the TMS
  // a median 23 ns apart, and a ghost space point pairing one muon's X hit
  // with the other's Y hit carries the AVERAGE of the two hit times, so it
  // sits half that offset away. Position alone cannot separate them where the
  // two tracks come within a bar pitch of each other.
  //
  // Default on: 15-file truth run 2026-09-24 (reports/2026-09-24_kalman_timing/,
  // 16,369 muons), on vs off: +0.1 to +0.3 pp completeness and purity in every
  // population, strict (both-views) metrics included -- all: completeness
  // 88.95 -> 89.05%, purity 90.11 -> 90.16%; ND-LAr-fiducial 95.42 -> 95.60% /
  // 97.33 -> 97.59%. Per muon it is mixed (673 better, 544 worse, 36 lose
  // >= 50 pp strict completeness), and it cannot rescue a seed that started on
  // the wrong particle -- it keeps a fit consistent with its own start.
  bool UseTimeInSelection = true;
  // Per-space-point time resolution (ns) after the path-length TOF
  // correction. Measured 2026-09-24 on 120k X/Y-truth-agreeing muon space
  // points (files 1-4): pooled sd 5.85 ns, MAD-sigma 5.66 ns, only 0.18% of
  // points beyond 20 ns (i.e. close to Gaussian, no heavy tail to guard).
  double TimeSigmaNs = 5.8;
  // If > 0, also reject any candidate with |r| > TimeGateNSigma *
  // sqrt(sigma_t^2 + sigma_t^2/n). 0 = time only ranks, never gates.
  double TimeGateNSigma = 0.0;

  // X/Y hit-time agreement as a ghost discriminant in per-layer candidate
  // selection. A space point pairs an X-view and a Y-view hit; if both came
  // from one particle their times agree once each hit's light-transit delay
  // along its bar is removed (TMS_SpacePointTiming), while a ghost pairing
  // two particles' hits keeps their real time difference. When on, and a
  // source has been given with Follower::SetXYTimeDifferenceSource(), each
  // candidate's score gains dt^2 / XYTimeSigmaNs^2. Unlike the time term
  // above it needs no track context -- it is a property of the point itself.
  // Candidates whose dt is unavailable get no penalty.
  //
  // Motivation (2026-09-24, reports/2026-09-24_reco_hit_lookaside/):
  // transit-corrected |dt| <= 10 ns keeps 84% of genuine points but only 71%
  // of same-interaction ghosts and 37.5% of different-interaction ghosts.
  //
  // Default on: 15-file truth run 2026-09-24 (reports/2026-09-24_kalman_xytime/,
  // 16,544 muons), sigma 3.7 vs off: strict (both-views) completeness /
  // purity 74.94/75.69 -> 75.91/76.57% overall, +1.2 to +1.5 pp in
  // multi-muon slices and TMS-start tracks; loose metrics flat (+/-0.1).
  // Wider sigma (5.1, 8.9) gave monotonically less. Per muon: 904 better,
  // 432 worse, 19 lose >= 50 pp strict completeness. Only has an effect
  // when a source is set (SetXYTimeDifferenceSource) and a layer has more
  // than one passing candidate -- it ranks, it doesn't gate.
  bool UseXYTimeInSelection = true;
  // Width of the genuine-point transit-corrected dt distribution (ns): its
  // MAD-sigma with the physics transit correction (TMS_SpacePointTiming).
  // The tails are wider (sd 8.9 ns), so this may need relaxing.
  double XYTimeSigmaNs = 3.7;
  // If > 0, also reject any candidate with |dt| > XYTimeGateNSigma *
  // XYTimeSigmaNs. 0 = ranks only, never gates.
  double XYTimeGateNSigma = 0.0;
};

// Transit-corrected X-hit minus Y-hit time (ns) for a space point; returns
// false if it can't be computed for that point.
using XYTimeDifferenceFn = std::function<bool(const TMS_SpacePoint &, double &)>;

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
  std::vector<double> CandidateChi2;  // parallel to CandidateIndices (position-only chi2)
  // Parallel to CandidateIndices: each candidate's time chi2 against the
  // running track t0 (see Config::UseTimeInSelection). Filled only when
  // time is in use; empty otherwise.
  std::vector<double> CandidateTimeChi2;
  // Parallel to CandidateIndices: each candidate's X/Y time-agreement chi2
  // (see Config::UseXYTimeInSelection), 0 where unavailable. Filled only
  // when that term is in use; empty otherwise.
  std::vector<double> CandidateXYTimeChi2;
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
  // RangedOut: the muon's own energy loss through the material to the next
  // layer took it to the momentum floor, i.e. it physically stops before
  // reaching that layer. A normal termination (Converged stays true), not a
  // failure -- it is reported separately from ReachedRangeEnd only because
  // the walk ended before running out of layers. The result keeps every node
  // up to the last layer actually reached.
  enum class StopReason { NotStarted, ReachedRangeEnd, GapLimitExceeded, Diverged, RangedOut };
  StopReason Stop = StopReason::NotStarted;

  bool Converged = false;
  std::vector<FollowedNode> Nodes;  // one per z-layer walked, low->high z

  double MomentumMeV = 0.0;  // from the final node's filtered q/p
  double Charge = 0.0;       // sign of the final node's filtered q/p
  double TotalChi2 = 0.0;
  int NDoF = 0;
  int NGapsFilled = 0;             // layers skipped for lack of a good candidate
  int NAmbiguousLayersResolved = 0;  // layers where >1 candidate existed
  int HeadSkip = 0;  // RunBestSeed(): leading object layers this hypothesis's seed skipped
  // Running track time origin, mean of (t - s/c) over accepted points (ns).
  // Always computed, whether or not time is used in selection.
  double TrackT0Ns = 0.0;
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
    //
    // allHypotheses / bestIndex (both optional): every hypothesis' FitResult in
    // the order tried, and the index of the one returned -- for studying how
    // the hypotheses are ranked (see IsBetterFit in the .cpp).
    FitResult RunBestSeed(const std::vector<TMS_SpacePoint> &allSpacePoints,
                          const std::vector<std::size_t> &objectIndices,
                          std::vector<FitResult> *allHypotheses = nullptr,
                          std::size_t *bestIndex = nullptr) const;

    // Where Config::UseXYTimeInSelection gets each point's transit-corrected
    // X/Y time difference. The follower only sees TMS_SpacePoint, which keeps
    // the average of its two hit times; the caller, which has the hits,
    // supplies the difference (typically keyed on the point's hit indices).
    void SetXYTimeDifferenceSource(XYTimeDifferenceFn source) { fXYTimeDifference = std::move(source); }

  private:
    Config fConfig;
    const IFieldModel &fField;
    XYTimeDifferenceFn fXYTimeDifference;
};

}  // namespace TMS_KalmanFollower

#endif
