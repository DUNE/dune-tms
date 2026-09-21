#include "TMS_KalmanFollower.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>

#include "TVector3.h"

#include "BetheBloch.h"
#include "Material.h"
#include "MultipleScattering.h"
#include "TMS_Constants.h"
#include "TMS_Geom.h"
#include "TMS_LayerGrouping.h"

namespace TMS_KalmanFollower {
namespace {

// State vector order throughout this file: [x, y, dx/dz, dy/dz, q/p].
// q/p is charge in units of e over momentum magnitude in MeV/c, so its own
// sign carries the charge -- this is what lets curvature be linear in qp
// (see Kappa() below) instead of needing a separate sign variable.
struct StepState {
  double x = 0.0, y = 0.0, z = 0.0;
  double dxdz = 0.0, dydz = 0.0;
  double qp = 0.0;
  TMatrixD cov{5, 5};
  // Set by Predict() when a step's deterministic propagation lands outside
  // any physically sane bound (see the comment at its check site) --
  // signals the caller to stop rather than hand an absurd position/slope
  // to TMS_Geom::GetMaterials, which can spend a very long time (and a lot
  // of memory) trying to navigate between two wildly separated points.
  bool Diverged = false;
  // Set by PredictSubstep() when energy loss in the material stepped through
  // drives the momentum down to kMinMomentumMeV from above -- the particle
  // has stopped, which is a physical end of the track, not a fit failure.
  // Deliberately NOT set when q/p was already at the floor on entry (a
  // Kalman-update collapse rather than genuine range-out).
  bool RangedOut = false;
};

// The standard "0.3 rule" curvature constant: p[GeV/c] = 0.299792458 *
// B[T] * R[m] * q[e], rearranged to kappa[1/mm] = 1/R[mm] = 0.299792458 *
// B[T] * q[e]/p[MeV/c] = 0.299792458 * B[T] * qp. This is well-established
// physics, unlike the field MAGNITUDE itself (see TMS_FieldModel.h) --
// legacy TMS_Kalman.cpp's own "0.303*...*1.95" factor is undocumented and,
// tellingly, never actually applied (TMS_Kalman.cpp:252-253).
constexpr double kCurvatureConstant = 0.299792458;

double Kappa(double qp, double fieldTeslaY) { return kCurvatureConstant * fieldTeslaY * qp; }

// Wolin & Ho (Nucl Inst A329 1993 493-500) multiple-scattering covariance,
// same formula TMS_KalmanNode::FillUpdatedCovarianceMatrix uses
// (TMS_Kalman.h:220-250), reproduced here as a free function writing into a
// plain 5x5 rather than a class member. Always the forward-fitting sign
// convention (this follower only ever walks low->high z).
void FillMultipleScatteringCovariance(TMatrixD &cov, double pathLengthGcm2,
                                       double dxdz, double dydz, double ms) {
  const double norm = 1.0 + dxdz * dxdz + dydz * dydz;
  const double covAxAx = norm * ms * (1.0 + dxdz * dxdz);
  const double covAyAy = norm * ms * (1.0 + dydz * dydz);
  const double covAxAy = norm * ms * dxdz * dydz;
  const double pathLengthSq = pathLengthGcm2 * pathLengthGcm2;

  cov(0, 0) = covAxAx * pathLengthSq / 4.0;
  cov(1, 1) = covAyAy * pathLengthSq / 4.0;
  cov(2, 2) = covAxAx;
  cov(3, 3) = covAyAy;
  cov(1, 0) = cov(0, 1) = covAxAy * pathLengthSq / 4.0;
  cov(2, 0) = cov(0, 2) = covAxAx * pathLengthGcm2 / 2.0;
  cov(3, 1) = cov(1, 3) = covAyAy * pathLengthGcm2 / 2.0;
  cov(3, 0) = cov(0, 3) = covAxAy * pathLengthGcm2 / 2.0;
  cov(2, 1) = cov(1, 2) = covAxAy * pathLengthGcm2 / 2.0;
  cov(3, 2) = cov(2, 3) = covAxAy;
  // cov(4,4) (q/p variance) is filled separately by the caller from the
  // Bethe-Bloch straggling variance -- legacy sets this to a bare "qp"
  // value (TMS_Kalman.h:236), which is not a real variance in any unit
  // system; we derive one properly instead (see ApplyMaterialSteps).
}

// Physical floor on any state's momentum, applied after both the material-
// loss update (ApplyMaterialSteps) and the Kalman gain update (UpdateState)
// -- either one can in principle push q/p to an extreme value (a bad chi2-
// gated update on a real, messy real-data fit isn't impossible even with
// reasonable covariances), and an extreme q/p feeds straight into the next
// step's curvature kick (Kappa()), which is exactly what earlier caused a
// runaway TMS_Geom::GetMaterials navigation between wildly separated
// points. Preserves sign.
constexpr double kMinMomentumMeV = 20.0;
void ClampMomentum(double &qp) {
  if (qp == 0.0) return;
  const double momentum = 1.0 / std::abs(qp);
  if (momentum < kMinMomentumMeV) {
    const double chargeSign = (qp > 0.0) ? 1.0 : -1.0;
    qp = chargeSign / kMinMomentumMeV;
  }
}

// Walks the real material budget between two points (TMS_Geom::GetMaterials,
// already proven correct via app/ShootRay.cpp), applying mean Bethe-Bloch
// energy loss to qpInOut and accumulating Lynch-Dahl multiple-scattering
// variance. Returns the resulting process-noise covariance contribution;
// qpVarianceOut carries the straggling-derived q/p variance for the caller
// to add onto cov(4,4) separately (kept out of the 5x5 here since it needs
// the FINAL dxdz/dydz, which the Wolin-Ho formula above already accounts
// for through its own arguments).
TMatrixD ApplyMaterialSteps(const TVector3 &start, const TVector3 &end,
                             double dxdz, double dydz, double &qpInOut,
                             double &qpVarianceOut, bool &rangedOutOut) {
  TMatrixD scatterCov(5, 5);
  qpVarianceOut = 0.0;
  rangedOutOut = false;

  const double chargeSign = (qpInOut >= 0.0) ? 1.0 : -1.0;
  double momentum = (std::abs(qpInOut) > 1e-12) ? 1.0 / std::abs(qpInOut) : 1.0;
  double energy = std::sqrt(momentum * momentum + BetheBloch_Utils::Mm * BetheBloch_Utils::Mm);
  const double energyFloor = std::sqrt(kMinMomentumMeV * kMinMomentumMeV +
                                        BetheBloch_Utils::Mm * BetheBloch_Utils::Mm);

  const std::vector<std::pair<TGeoMaterial *, double> > materials =
      TMS_Geom::GetInstance().GetMaterials(start, end);

  // Placeholder material type, matching TMS_Kalman.cpp:6-7's convention --
  // .fMaterial is reassigned every step below before either calculator is
  // actually used.
  BetheBloch_Calculator bethe(Material::kPolyStyrene);
  MultipleScatter_Calculator msc(Material::kPolyStyrene);

  double totalPathLengthGcm2 = 0.0;
  double totalEnergyVarianceSq = 0.0;

  for (const auto &materialStep : materials) {
    // Same TGeoMaterial -> density -> Material(density) bridge legacy uses
    // (TMS_Kalman.cpp:304-328): fragile (Material(double) throws outside 3
    // tightly-toleranced hardcoded density windows) but matches this
    // detector's known-simple material budget (scintillator/steel/air).
    double density = materialStep.first->GetDensity() / (CLHEP::g / CLHEP::cm3);
    double thickness = materialStep.second / 10.0;  // mm -> cm
    const double scaleFactor = TMS_Geom::GetInstance().Scale(1.0);
    density /= std::pow(scaleFactor, 3);
    thickness = TMS_Geom::GetInstance().Scale(thickness);

    try {
      Material matter(density);
      bethe.fMaterial = matter;
      msc.fMaterial = matter;
    } catch (const std::invalid_argument &) {
      continue;  // unrecognised material at this step -- skip it, as legacy does
    }

    totalPathLengthGcm2 += density * thickness;

    // Always walking forward (low->high z): energy decreases.
    energy -= bethe.Calc_dEdx(energy) * density * thickness;
    // Floor at the energy of the kMinMomentumMeV momentum floor, NOT at the
    // bare rest mass. At E == Mm, beta is exactly 0 and Calc_dEdx /
    // Calc_dEdx_Straggling divide by beta^2 (and by MaximumEnergyTransfer,
    // also 0 there), so the NEXT material step returned inf/NaN. That NaN
    // then survived every later guard (`energy < Mm` and `momentum <
    // kMinMomentumMeV` are both false for NaN) and turned q/p, dx/dz and x
    // into NaN on the following gap nodes. Stopping the energy at the floor
    // keeps beta finite and agrees with the momentum floor applied below.
    if (energy < energyFloor) {
      energy = energyFloor;
      // Only a genuine range-out if we started this step above the floor;
      // a state already pinned there by a Kalman update is not.
      if (momentum > kMinMomentumMeV * 1.0001) rangedOutOut = true;
    }

    const double energyStragglingSigma = bethe.Calc_dEdx_Straggling(energy) * density * thickness;
    totalEnergyVarianceSq += energyStragglingSigma * energyStragglingSigma;

    msc.Calc_MS(energy, thickness * density);
  }

  momentum = BetheBloch_Utils::EnergyToMomentum(BetheBloch_Utils::Mm, energy);
  // A muon this deep into range-out is no longer something this v1 forward
  // fit should trust as a genuine mid-track state -- floor well above the
  // literal rest-mass edge case (which would let momentum -> 0) so q/p,
  // and therefore the curvature kick on the NEXT Predict() call, can't
  // blow up. (A 1 MeV floor here, tried first, let q/p reach ~1 e/MeV --
  // 100-1000x a real muon's q/p -- which sent the next step's position
  // to a wild extrapolation and TMS_Geom::GetMaterials into a very long,
  // very memory-hungry navigation between two absurdly distant points.)
  if (momentum < kMinMomentumMeV) momentum = kMinMomentumMeV;
  qpInOut = chargeSign / momentum;

  // qp = charge/p; d(qp)/dE = -charge/p^2 * dp/dE, and dp/dE = E/p for a
  // relativistic particle, so d(qp)/dE = -charge*E/p^3.
  const double dqpdE = -chargeSign * energy / (momentum * momentum * momentum);
  qpVarianceOut = dqpdE * dqpdE * totalEnergyVarianceSq;

  const double msSigma = msc.Calc_MS_Sigma();
  if (totalPathLengthGcm2 > 0.0) {
    FillMultipleScatteringCovariance(scatterCov, totalPathLengthGcm2, dxdz, dydz, msSigma);
  }
  return scatterCov;
}

// One sub-step of the swimmer: field-bent straight-line-plus-curvature step
// (the bend legacy computes and discards, see this file's header comment),
// then real energy loss + multiple scattering via ApplyMaterialSteps, over
// a single small dz. Kept separate from Predict() below because the
// quadratic-in-dz Jacobian term (transfer(0,4)) is only a safe
// linearization for a small dz -- see Predict()'s comment.
StepState PredictSubstep(const StepState &previous, double subDz, const IFieldModel &field,
                          bool stopOnRangeOut) {
  StepState predicted = previous;
  predicted.z = previous.z + subDz;

  const TVector3 midpoint(previous.x + 0.5 * previous.dxdz * subDz,
                           previous.y + 0.5 * previous.dydz * subDz,
                           previous.z + 0.5 * subDz);
  const double fieldY = field.GetField(midpoint).Y();
  const double kappa = Kappa(previous.qp, fieldY);

  predicted.dxdz = previous.dxdz + kappa * subDz;
  predicted.dydz = previous.dydz;  // B along y bends only dx/dz in this model
  predicted.x = previous.x + previous.dxdz * subDz + 0.5 * kappa * subDz * subDz;
  predicted.y = previous.y + previous.dydz * subDz;

  // Divergence guard: a slope or position this far outside anything
  // physical (TMS is a few meters across; even a hard-scattering muon
  // doesn't reach a multi-radian slope) means the fit has already gone
  // unphysical upstream -- stop here rather than hand these values to
  // TMS_Geom::GetMaterials, which can spend a very long time (and a lot
  // of memory) trying to navigate between two wildly separated points.
  // Written as !(|v| <= limit) rather than |v| > limit: every comparison
  // against NaN is false in C++, so the plain `>` form lets a NaN state
  // through as "not diverged".
  if (!(std::abs(predicted.dxdz) <= 10.0) || !(std::abs(predicted.dydz) <= 10.0) ||
      !(std::abs(predicted.x) <= 1.0e5) || !(std::abs(predicted.y) <= 1.0e5) ||
      !std::isfinite(predicted.qp)) {
    predicted.Diverged = true;
    return predicted;
  }

  // Transfer Jacobian d(new state)/d(old state); see this file's header
  // comment for the derivation. Off-diagonal q/p terms capture "how would
  // a different momentum have bent this step differently" -- physically
  // important, and exactly what a flat drift-only transfer matrix (like
  // legacy's, TMS_Kalman.cpp:112-119) cannot represent. The (0,4) term is
  // quadratic in subDz -- keeping subDz small (see Predict()) keeps this
  // linearization valid and keeps q/p from picking up an artificially
  // large lever arm on position over one step.
  TMatrixD transfer(5, 5);
  transfer.UnitMatrix();
  transfer(0, 2) = subDz;
  transfer(0, 4) = 0.5 * kCurvatureConstant * fieldY * subDz * subDz;
  transfer(1, 3) = subDz;
  transfer(2, 4) = kCurvatureConstant * fieldY * subDz;

  // NOTE: TMatrixD::T() transposes IN PLACE and returns *this -- it is NOT
  // a non-mutating "give me a transposed copy" like Eigen's/numpy's
  // .transpose(). Calling it inline inside this expression corrupted
  // `transfer` mid-evaluation (operand evaluation order across `*` isn't
  // sequenced, so `transfer` could already be transposed by the time
  // `transfer * previous.cov` itself runs) and silently broke the whole
  // covariance propagation -- discovered via KF_DEBUG instrumentation
  // showing cov(0,4)/cov(2,4) staying exactly zero after propagation
  // despite a manifestly nonzero curvature Jacobian, no matter how the
  // physics config (field magnitude, initial covariances, step size) was
  // tuned. Use the kTransposed constructor for an explicit, non-mutating
  // copy instead -- same reason legacy TMS_Kalman.h keeps a manually-
  // maintained separate TransferMatrixT member rather than ever calling
  // .T() inline.
  const TMatrixD transferT(TMatrixD::kTransposed, transfer);
  const TMatrixD propagatedCov = transfer * previous.cov * transferT;

  // Divergence guard: if repeated un-updated gaps have already inflated
  // the covariance past the detector's own transverse extent, the chi2
  // gate (which scales with this covariance) stops meaningfully rejecting
  // anything -- effectively any nearby point in a dense slice passes,
  // letting the walk wander indefinitely through unrelated material
  // instead of correctly running out of plausible candidates. Bound it
  // well above genuine values (this detector is a few meters across) so
  // real fits are never affected, but a runaway is caught here rather than
  // a hundred layers later.
  constexpr double kMaxPositionVarianceMM2 = 4.0e6;  // (2000mm)^2
  if (!(propagatedCov(0, 0) <= kMaxPositionVarianceMM2) || !(propagatedCov(1, 1) <= kMaxPositionVarianceMM2)) {
    predicted.Diverged = true;
    return predicted;
  }

  const TVector3 startPos(previous.x, previous.y, previous.z);
  const TVector3 endPos(predicted.x, predicted.y, predicted.z);
  double qpVariance = 0.0;
  bool rangedOut = false;
  const TMatrixD scatterCov = ApplyMaterialSteps(startPos, endPos, predicted.dxdz,
                                                  predicted.dydz, predicted.qp, qpVariance, rangedOut);
  if (rangedOut && stopOnRangeOut) {
    // The muon stops inside this sub-step: nothing beyond here is reachable,
    // so end the walk on the last layer actually reached instead of carrying
    // an ever-less-constrained state on to (typically unrelated) later layers.
    predicted.RangedOut = true;
    return predicted;
  }

  predicted.cov = propagatedCov + scatterCov;
  predicted.cov(4, 4) += qpVariance;

  // Last line of defence: any NaN/inf that still reaches the state (e.g. a
  // new degenerate material) stops the fit as Diverged instead of being
  // carried silently through gap nodes.
  if (!std::isfinite(predicted.qp) || !std::isfinite(predicted.cov(0, 0)) ||
      !std::isfinite(predicted.cov(4, 4))) {
    predicted.Diverged = true;
  }

  return predicted;
}

// Propagate one node forward to zTarget by sub-stepping through it (a real
// swimmer, not one big linearized jump). A single dz~100-200mm plane-to-
// plane gap makes PredictSubstep's transfer(0,4) term (quadratic in dz)
// large enough to give q/p an artificially outsized lever arm on position
// -- discovered on the flagship real-data case: a single ~130mm jump let
// one ordinary ~45mm position residual (typical quantization+scattering
// noise, not a real momentum signal) drive the Kalman gain to collapse
// momentum to its floor in one update. Sub-stepping keeps each
// linearization small and keeps the accumulated covariance honest.
StepState Predict(const StepState &previous, double zTarget, const IFieldModel &field,
                   double maxSubstepLengthMM, bool stopOnRangeOut) {
  const double totalDz = zTarget - previous.z;
  if (std::abs(totalDz) < 1e-9) return previous;

  const int nSubsteps = std::max(
      1, static_cast<int>(std::ceil(std::abs(totalDz) / maxSubstepLengthMM)));
  const double subDz = totalDz / nSubsteps;

  StepState current = previous;
  for (int i = 0; i < nSubsteps; ++i) {
    current = PredictSubstep(current, subDz, field, stopOnRangeOut);
    if (current.Diverged || current.RangedOut) return current;
  }
  return current;
}

// Measurement covariance for one candidate space point, from the project's
// own already-validated ~1-bar-pitch quantization scale (see this file's
// Config::AssumedBarPitchMM comment) rather than a per-hit bar-width
// lookup -- this follower works from TMS_SpacePoint alone, no native
// TMS_Hit round-trip. A uniform distribution of width w has variance
// w^2/12; both axes share the same pitch, matching
// TMS_GraphTrackFinder::Config::PositionQuantizationX/Y's precedent.
TMatrixD BuildMeasurementCovariance(double barPitchMM) {
  TMatrixD r(2, 2);
  const double sigma = barPitchMM / std::sqrt(12.0);
  r(0, 0) = sigma * sigma;
  r(1, 1) = sigma * sigma;
  return r;
}

double Chi2(const StepState &predicted, const TMS_SpacePoint &candidate, const TMatrixD &r) {
  const double rx = candidate.GetX() - predicted.x;
  const double ry = candidate.GetY() - predicted.y;
  TMatrixD innovationCov(2, 2);
  innovationCov(0, 0) = predicted.cov(0, 0) + r(0, 0) + 1e-6;  // tiny epsilon: numerical safety net only
  innovationCov(1, 1) = predicted.cov(1, 1) + r(1, 1) + 1e-6;
  innovationCov(0, 1) = predicted.cov(0, 1) + r(0, 1);
  innovationCov(1, 0) = predicted.cov(1, 0) + r(1, 0);
  innovationCov.Invert();
  return rx * rx * innovationCov(0, 0) + 2.0 * rx * ry * innovationCov(0, 1) +
         ry * ry * innovationCov(1, 1);
}

StepState UpdateState(const StepState &predicted, const TMS_SpacePoint &chosen, const TMatrixD &r) {
  StepState updated = predicted;
  const double residualX = chosen.GetX() - predicted.x;
  const double residualY = chosen.GetY() - predicted.y;

  TMatrixD innovationCov(2, 2);
  innovationCov(0, 0) = predicted.cov(0, 0) + r(0, 0) + 1e-6;
  innovationCov(1, 1) = predicted.cov(1, 1) + r(1, 1) + 1e-6;
  innovationCov(0, 1) = predicted.cov(0, 1) + r(0, 1);
  innovationCov(1, 0) = predicted.cov(1, 0) + r(1, 0);
  innovationCov.Invert();

  // Kalman gain K = Cov * H^T * S^-1 (5x2); H picks out (x,y), so H^T's
  // columns are just Cov's first two columns.
  double gain[5][2];
  for (int row = 0; row < 5; ++row) {
    const double covCol0 = predicted.cov(row, 0);
    const double covCol1 = predicted.cov(row, 1);
    gain[row][0] = covCol0 * innovationCov(0, 0) + covCol1 * innovationCov(1, 0);
    gain[row][1] = covCol0 * innovationCov(0, 1) + covCol1 * innovationCov(1, 1);
  }

  double stateVec[5] = {predicted.x, predicted.y, predicted.dxdz, predicted.dydz, predicted.qp};
  for (int row = 0; row < 5; ++row)
    stateVec[row] += gain[row][0] * residualX + gain[row][1] * residualY;

  if (std::getenv("KF_DEBUG")) {
    std::cerr << "[KF_DEBUG UpdateState] predicted.qp=" << predicted.qp
              << " cov(4,4)=" << predicted.cov(4, 4)
              << " cov(0,4)=" << predicted.cov(0, 4)
              << " cov(2,4)=" << predicted.cov(2, 4)
              << " cov(0,0)=" << predicted.cov(0, 0)
              << " residualX=" << residualX << " residualY=" << residualY
              << " gain[4][0]=" << gain[4][0] << " gain[4][1]=" << gain[4][1]
              << " dqp=" << (gain[4][0] * residualX + gain[4][1] * residualY)
              << " gain[2][0]=" << gain[2][0]
              << " ddxdz=" << (gain[2][0] * residualX + gain[2][1] * residualY)
              << std::endl;
  }

  updated.x = stateVec[0];
  updated.y = stateVec[1];
  updated.dxdz = stateVec[2];
  updated.dydz = stateVec[3];
  updated.qp = stateVec[4];
  ClampMomentum(updated.qp);

  // Cov_new = Cov_pred - K * H * Cov_pred; H*Cov_pred is just the top two
  // rows of Cov_pred.
  TMatrixD newCov = predicted.cov;
  for (int row = 0; row < 5; ++row) {
    for (int col = 0; col < 5; ++col) {
      newCov(row, col) -= gain[row][0] * predicted.cov(0, col) + gain[row][1] * predicted.cov(1, col);
    }
  }
  updated.cov = newCov;
  return updated;
}

struct GateResult {
  bool Accepted = false;
  std::size_t ChosenIndex = 0;
  std::vector<std::size_t> CandidateIndices;
  std::vector<double> CandidateChi2;
};

// The ambiguity-resolution core: score every candidate at this layer
// against the predicted state, accept the best one under the chi2 gate (or
// none, if nothing passes -- a gap, handled by the caller).
GateResult ResolveLayer(const StepState &predicted, const std::vector<std::size_t> &candidatesAtLayer,
                         const std::vector<TMS_SpacePoint> &allSpacePoints,
                         double barPitchMM, double chiSquareGateMax) {
  GateResult result;
  double bestChi2 = std::numeric_limits<double>::infinity();
  for (std::size_t index : candidatesAtLayer) {
    const TMS_SpacePoint &candidate = allSpacePoints[index];
    const TMatrixD measurementCov = BuildMeasurementCovariance(barPitchMM);
    const double chi2 = Chi2(predicted, candidate, measurementCov);
    result.CandidateIndices.push_back(index);
    result.CandidateChi2.push_back(chi2);
    if (chi2 <= chiSquareGateMax && chi2 < bestChi2) {
      bestChi2 = chi2;
      result.ChosenIndex = index;
      result.Accepted = true;
    }
  }
  return result;
}

// Local slope estimate from the seed's own first few points (up to 3),
// averaged over consecutive pairs -- just enough smoothing to not be
// thrown off by a single noisy first pair, without needing the whole path.
void SeedDirection(const std::vector<TMS_SpacePoint> &allSpacePoints,
                    const std::vector<std::size_t> &seedPath, double &dxdzOut, double &dydzOut) {
  dxdzOut = 0.0;
  dydzOut = 0.0;
  const std::size_t n = std::min<std::size_t>(seedPath.size(), 3);
  int nSlopes = 0;
  for (std::size_t i = 1; i < n; ++i) {
    const TMS_SpacePoint &a = allSpacePoints[seedPath[i - 1]];
    const TMS_SpacePoint &b = allSpacePoints[seedPath[i]];
    const double dz = b.GetZ() - a.GetZ();
    if (std::abs(dz) < 1e-6) continue;
    dxdzOut += (b.GetX() - a.GetX()) / dz;
    dydzOut += (b.GetY() - a.GetY()) / dz;
    ++nSlopes;
  }
  if (nSlopes > 0) {
    dxdzOut /= nSlopes;
    dydzOut /= nSlopes;
  }
}

// A self-contained re-implementation of TMS_ChargeID::ID_Track_Charge's
// bend-direction algorithm (src/TMS_ChargeID.cpp), operating directly on
// TMS_SpacePoint x/z instead of TMS_Hit::GetRecoX() -- avoids needing a
// native-hit round-trip (with its live-geometry dependency, see this
// file's header comment) just to seed a charge sign. Same physics: within
// each contiguous run of points inside one field region, draw a line from
// the run's first to last point and count how many interior points fall
// to each side -- more on one side means the track bent that way. Region
// 2's sign convention is flipped relative to region 1/3, matching the
// field pointing the opposite way in the central region (the same
// asymmetry RegionFieldModel encodes).
double SeedCharge(const std::vector<TMS_SpacePoint> &allSpacePoints,
                   const std::vector<std::size_t> &seedPath) {
  enum Region { kRegion1, kRegion2, kRegion3, kOutside };
  auto regionOf = [](double x) {
    if (x >= TMS_Const::TMS_Magnetic_region_1_outer_edge &&
        x <= TMS_Const::TMS_Magnetic_region_1_and_2_border)
      return kRegion1;
    if (x >= TMS_Const::TMS_Magnetic_region_1_and_2_border &&
        x <= TMS_Const::TMS_Magnetic_region_2_and_3_border)
      return kRegion2;
    if (x >= TMS_Const::TMS_Magnetic_region_2_and_3_border &&
        x <= TMS_Const::TMS_Magnetic_region_3_outer_edge)
      return kRegion3;
    return kOutside;
  };

  int nPlus = 0, nMinus = 0;
  std::size_t runStart = 0;
  while (runStart < seedPath.size()) {
    const Region region = regionOf(allSpacePoints[seedPath[runStart]].GetX());
    std::size_t runEnd = runStart + 1;
    while (runEnd < seedPath.size() && regionOf(allSpacePoints[seedPath[runEnd]].GetX()) == region) ++runEnd;

    if (region != kOutside && runEnd - runStart > 2) {
      const TMS_SpacePoint &front = allSpacePoints[seedPath[runStart]];
      const TMS_SpacePoint &back = allSpacePoints[seedPath[runEnd - 1]];
      const double dz = back.GetZ() - front.GetZ();
      if (std::abs(dz) > 1e-6) {
        const double slope = (back.GetX() - front.GetX()) / dz;
        const bool positiveMeansPlus = (region != kRegion2);
        for (std::size_t i = runStart + 1; i + 1 < runEnd; ++i) {
          const TMS_SpacePoint &p = allSpacePoints[seedPath[i]];
          const double interpolatedX = front.GetX() + slope * (p.GetZ() - front.GetZ());
          const double signedDist = p.GetX() - interpolatedX;
          if (signedDist == 0.0) continue;
          if ((signedDist > 0.0) == positiveMeansPlus) ++nPlus;
          else ++nMinus;
        }
      }
    }
    runStart = runEnd;
  }

  // Undecided (nPlus == nMinus, including the common 0-vs-0 case for a
  // short/straight seed) defaults to +1 -- the wide initial q/p covariance
  // (Config::InitialCovQPQP) already carries the "not really known"
  // uncertainty, per the project's open charge-seed-fallback decision.
  if (nPlus < nMinus) return -1.0;
  return 1.0;
}

// Ranks two hypotheses' fit results for RunBestSeed(): converged beats
// unconverged outright; among either group, more accepted (non-gap) nodes
// beats fewer (a hypothesis that only limped one step before diverging
// covers less of the object than one that walked its full length); ties
// broken by chi2/NDoF (lower is better), the standard goodness-of-fit
// comparison once coverage is equal.
bool IsBetterFit(const FitResult &a, const FitResult &b) {
  if (a.Converged != b.Converged) return a.Converged;
  int hitsA = 0, hitsB = 0;
  for (const FollowedNode &n : a.Nodes) if (n.HasHit) ++hitsA;
  for (const FollowedNode &n : b.Nodes) if (n.HasHit) ++hitsB;
  if (hitsA != hitsB) return hitsA > hitsB;
  const double chi2NDofA = a.NDoF > 0 ? a.TotalChi2 / a.NDoF : std::numeric_limits<double>::infinity();
  const double chi2NDofB = b.NDoF > 0 ? b.TotalChi2 / b.NDoF : std::numeric_limits<double>::infinity();
  return chi2NDofA < chi2NDofB;
}

// Smallest momentum (MeV/c) that lets a muon travel from start to end
// through the real material budget: energy loss walked BACKWARDS from the
// kMinMomentumMeV floor, so each step's dE/dx is evaluated at the energy the
// muon has after that step. A lower bound on the true momentum (a straight
// segment start->end, no scattering), used to seed q/p (Config::RangeSeedMargin).
double RangeMomentumMeV(const TVector3 &start, const TVector3 &end) {
  const std::vector<std::pair<TGeoMaterial *, double> > materials =
      TMS_Geom::GetInstance().GetMaterials(start, end);
  BetheBloch_Calculator bethe(Material::kPolyStyrene);
  double energy = std::sqrt(kMinMomentumMeV * kMinMomentumMeV +
                            BetheBloch_Utils::Mm * BetheBloch_Utils::Mm);
  for (auto it = materials.rbegin(); it != materials.rend(); ++it) {
    double density = it->first->GetDensity() / (CLHEP::g / CLHEP::cm3);
    double thickness = it->second / 10.0;  // mm -> cm
    const double scaleFactor = TMS_Geom::GetInstance().Scale(1.0);
    density /= std::pow(scaleFactor, 3);
    thickness = TMS_Geom::GetInstance().Scale(thickness);
    try {
      Material matter(density);
      bethe.fMaterial = matter;
    } catch (const std::invalid_argument &) {
      continue;
    }
    const double loss = bethe.Calc_dEdx(energy) * density * thickness;
    if (std::isfinite(loss)) energy += loss;
  }
  return BetheBloch_Utils::EnergyToMomentum(BetheBloch_Utils::Mm, energy);
}

}  // namespace

Follower::Follower(const Config &config, const IFieldModel &field) : fConfig(config), fField(field) {}

FitResult Follower::Run(const std::vector<TMS_SpacePoint> &allSpacePoints,
                         const std::vector<std::size_t> &seedPath) const {
  FitResult result;
  if (seedPath.size() < 2 || allSpacePoints.empty()) return result;

  const std::vector<std::vector<std::size_t> > zLayers =
      TMS_LayerGrouping::Build(allSpacePoints, fConfig.LayerZTolerance);

  std::size_t startLayer = zLayers.size();
  for (std::size_t layerIdx = 0; layerIdx < zLayers.size() && startLayer == zLayers.size(); ++layerIdx) {
    for (std::size_t index : zLayers[layerIdx]) {
      if (index == seedPath.front()) {
        startLayer = layerIdx;
        break;
      }
    }
  }
  if (startLayer == zLayers.size()) return result;  // seed's own point wasn't in the pool passed in

  std::size_t seedEndLayer = startLayer;
  for (std::size_t layerIdx = 0; layerIdx < zLayers.size(); ++layerIdx) {
    for (std::size_t index : zLayers[layerIdx]) {
      if (index == seedPath.back()) {
        seedEndLayer = layerIdx;
        break;
      }
    }
  }
  const std::size_t lastLayerToWalk =
      std::min(zLayers.size() - 1, seedEndLayer + static_cast<std::size_t>(fConfig.MaxLayersBeyondSeed));

  // --- Initial node: taken directly from the seed, not gated. ---
  double initialDxdz = 0.0, initialDydz = 0.0;
  SeedDirection(allSpacePoints, seedPath, initialDxdz, initialDydz);
  const double chargeSign = SeedCharge(allSpacePoints, seedPath);

  const TMS_SpacePoint &firstPoint = allSpacePoints[seedPath.front()];
  StepState current;
  current.x = firstPoint.GetX();
  current.y = firstPoint.GetY();
  current.z = firstPoint.GetZ();
  current.dxdz = initialDxdz;
  current.dydz = initialDydz;
  double seedMomentum = fConfig.InitialMomentumSeedMeV;
  if (fConfig.RangeSeedMargin > 0.0) {
    // z-extent of the seed object itself: its lowest-z and highest-z points.
    std::size_t lo = seedPath.front(), hi = seedPath.front();
    for (std::size_t idx : seedPath) {
      if (allSpacePoints[idx].GetZ() < allSpacePoints[lo].GetZ()) lo = idx;
      if (allSpacePoints[idx].GetZ() > allSpacePoints[hi].GetZ()) hi = idx;
    }
    const double pRange = RangeMomentumMeV(
        TVector3(allSpacePoints[lo].GetX(), allSpacePoints[lo].GetY(), allSpacePoints[lo].GetZ()),
        TVector3(allSpacePoints[hi].GetX(), allSpacePoints[hi].GetY(), allSpacePoints[hi].GetZ()));
    if (std::isfinite(pRange)) seedMomentum = std::max(seedMomentum, fConfig.RangeSeedMargin * pRange);
  }
  current.qp = chargeSign / seedMomentum;
  current.cov.Zero();
  current.cov(0, 0) = fConfig.InitialCovXX;
  current.cov(1, 1) = fConfig.InitialCovYY;
  current.cov(2, 2) = fConfig.InitialCovDXDZDXDZ;
  current.cov(3, 3) = fConfig.InitialCovDYDZDYDZ;
  current.cov(4, 4) = fConfig.InitialCovQPQP;
  if (fConfig.InitialQPRelSigma > 0.0) {
    const double sigmaQP = fConfig.InitialQPRelSigma / seedMomentum;
    current.cov(4, 4) = sigmaQP * sigmaQP;
  }

  FollowedNode firstNode;
  firstNode.Layer = startLayer;
  firstNode.Z = current.z;
  firstNode.HasHit = true;
  firstNode.ChosenSpacePointIndex = seedPath.front();
  firstNode.CandidateIndices = zLayers[startLayer];
  firstNode.CandidateChi2.assign(zLayers[startLayer].size(), 0.0);  // not evaluated for the seeded node
  firstNode.FilteredX = current.x;
  firstNode.FilteredY = current.y;
  firstNode.FilteredDXDZ = current.dxdz;
  firstNode.FilteredDYDZ = current.dydz;
  firstNode.FilteredQP = current.qp;
  firstNode.FilteredCovariance = current.cov;
  result.Nodes.push_back(firstNode);
  if (zLayers[startLayer].size() > 1) ++result.NAmbiguousLayersResolved;

  int consecutiveGaps = 0;
  result.Converged = true;
  result.Stop = FitResult::StopReason::ReachedRangeEnd;  // overridden below if the walk breaks early

  for (std::size_t layerIdx = startLayer + 1; layerIdx <= lastLayerToWalk; ++layerIdx) {
    const std::vector<std::size_t> &candidates = zLayers[layerIdx];
    if (candidates.empty()) continue;  // TMS_LayerGrouping never emits an empty layer; defensive only
    const double targetZ = allSpacePoints[candidates.front()].GetZ();

    const StepState predicted =
        Predict(current, targetZ, fField, fConfig.MaxSubstepLengthMM, fConfig.StopOnRangeOut);
    if (predicted.RangedOut) {
      // Physical end of the track (see StepState::RangedOut): keep every node
      // so far, and count it as a normal termination.
      result.Stop = FitResult::StopReason::RangedOut;
      break;  // result.Converged is still true from initialisation
    }
    if (predicted.Diverged) {
      result.Converged = false;
      result.Stop = FitResult::StopReason::Diverged;
      break;
    }
    const GateResult gate =
        ResolveLayer(predicted, candidates, allSpacePoints, fConfig.AssumedBarPitchMM, fConfig.ChiSquareGateMax);

    FollowedNode node;
    node.Layer = layerIdx;
    node.Z = targetZ;
    node.CandidateIndices = gate.CandidateIndices;
    node.CandidateChi2 = gate.CandidateChi2;

    if (gate.Accepted) {
      const TMS_SpacePoint &chosen = allSpacePoints[gate.ChosenIndex];
      const TMatrixD measurementCov = BuildMeasurementCovariance(fConfig.AssumedBarPitchMM);
      current = UpdateState(predicted, chosen, measurementCov);
      node.HasHit = true;
      node.ChosenSpacePointIndex = gate.ChosenIndex;
      for (std::size_t i = 0; i < gate.CandidateIndices.size(); ++i) {
        if (gate.CandidateIndices[i] == gate.ChosenIndex) {
          node.Chi2AtChosen = gate.CandidateChi2[i];
          break;
        }
      }
      result.TotalChi2 += node.Chi2AtChosen;
      result.NDoF += 2;
      consecutiveGaps = 0;
      if (candidates.size() > 1) ++result.NAmbiguousLayersResolved;
    } else {
      current = predicted;  // no update -- the inflated predicted covariance simply carries forward
      node.HasHit = false;
      ++result.NGapsFilled;
      ++consecutiveGaps;
      if (consecutiveGaps > fConfig.MaxConsecutiveGaps) {
        result.Converged = false;
        result.Stop = FitResult::StopReason::GapLimitExceeded;
      }
    }

    node.FilteredX = current.x;
    node.FilteredY = current.y;
    node.FilteredDXDZ = current.dxdz;
    node.FilteredDYDZ = current.dydz;
    node.FilteredQP = current.qp;
    node.FilteredCovariance = current.cov;
    result.Nodes.push_back(node);

    if (!result.Converged) break;
  }

  result.MomentumMeV = (std::abs(current.qp) > 1e-12) ? 1.0 / std::abs(current.qp) : 0.0;
  result.Charge = (current.qp >= 0.0) ? 1.0 : -1.0;
  result.NDoF -= 5;  // 5 fitted state parameters

  return result;
}

FitResult Follower::RunBestSeed(const std::vector<TMS_SpacePoint> &allSpacePoints,
                                 const std::vector<std::size_t> &objectIndices) const {
  FitResult best;
  if (objectIndices.empty()) return best;

  // Group just the object's own points to find ITS first z-layer -- not
  // allSpacePoints' first layer, which could belong to unrelated activity
  // elsewhere in a dense slice.
  std::vector<TMS_SpacePoint> objectPoints;
  objectPoints.reserve(objectIndices.size());
  for (std::size_t idx : objectIndices) objectPoints.push_back(allSpacePoints[idx]);
  const std::vector<std::vector<std::size_t> > objectLayers =
      TMS_LayerGrouping::Build(objectPoints, fConfig.LayerZTolerance);
  if (objectLayers.empty()) return best;

  // Hypotheses: for each number of leading layers to skip (0 .. MaxHeadSkip,
  // always leaving at least two layers to seed from), anchor the fit on every
  // point of the resulting first layer. Skipping matters when the object's
  // head belongs to something else -- e.g. a GraphTrackFinder path whose
  // first points sit on a companion particle from the same vertex (case E,
  // 2026-09-21): seeding there points the whole fit at the companion. The
  // hypotheses compete through IsBetterFit (converged, then most hits, then
  // chi2/ndof), all reco-only. A skip always gives up at least one layer's
  // hit, so it only wins when the full-head fit lost more than that.
  bool haveBest = false;
  const int maxSkip = std::max(0, fConfig.MaxHeadSkip);
  for (int skip = 0; skip <= maxSkip; ++skip) {
    if (objectLayers.size() < static_cast<std::size_t>(skip) + 2) break;

    const std::vector<std::size_t> &firstLayerLocal = objectLayers[skip];
    std::set<std::size_t> skippedOrFirstGlobal;
    for (int l = 0; l <= skip; ++l)
      for (std::size_t localIdx : objectLayers[l]) skippedOrFirstGlobal.insert(objectIndices[localIdx]);

    // The rest of the object, z-sorted (same (z,x,y) tie-break
    // TMS_LayerGrouping itself uses internally) -- shared across every anchor
    // of this skip level; only which point leads the seed path (and
    // therefore SeedDirection()'s first ~3-point average) changes.
    std::vector<std::size_t> restSorted;
    for (std::size_t idx : objectIndices)
      if (!skippedOrFirstGlobal.count(idx)) restSorted.push_back(idx);
    std::sort(restSorted.begin(), restSorted.end(), [&allSpacePoints](std::size_t a, std::size_t b) {
      const TMS_SpacePoint &pa = allSpacePoints[a];
      const TMS_SpacePoint &pb = allSpacePoints[b];
      if (pa.GetZ() != pb.GetZ()) return pa.GetZ() < pb.GetZ();
      if (pa.GetX() != pb.GetX()) return pa.GetX() < pb.GetX();
      return pa.GetY() < pb.GetY();
    });

    for (std::size_t localIdx : firstLayerLocal) {
      const std::size_t anchor = objectIndices[localIdx];
      std::vector<std::size_t> seedPath;
      seedPath.reserve(restSorted.size() + 1);
      seedPath.push_back(anchor);
      seedPath.insert(seedPath.end(), restSorted.begin(), restSorted.end());

      const FitResult candidate = Run(allSpacePoints, seedPath);
      if (!haveBest || IsBetterFit(candidate, best)) {
        best = candidate;
        haveBest = true;
      }
    }
  }
  return best;
}

}  // namespace TMS_KalmanFollower
