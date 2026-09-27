#include "TMS_KalmanFollower.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
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

// For the time-of-flight correction in the optional time term (muons treated
// as beta~1; see Config::UseTimeInSelection).
constexpr double kSpeedOfLightMMPerNs = 299.792458;

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

// Material steps along a straight line (TMS_Geom::GetMaterials, already
// proven correct via app/ShootRay.cpp).
typedef std::vector<std::pair<TGeoMaterial *, double> > MaterialSteps;

// Density (g/cm^3) of a material step, with the geometry's unit scaling --
// the same conversion ApplyMaterialSteps() uses.
double DensityGCm3(const TGeoMaterial *material) {
  const double scaleFactor = TMS_Geom::GetInstance().Scale(1.0);
  return material->GetDensity() / (CLHEP::g / CLHEP::cm3) / std::pow(scaleFactor, 3);
}

// Fraction of a step's path length in magnetized steel. The TMS field lives
// in the steel plates only (edep-sim's GDML field is attached to the steel
// volumes): measured 2026-09-25 from G4 truth, the field per unit steel is
// +-1.0-1.1 T at every |x|, and the effective field in each z section equals
// 1 T times its steel fraction (thin 15/65: 0.227 vs 0.231 T; thick 40/90:
// 0.424 vs 0.444 T; double 80/130: 0.597 vs 0.615 T). Steel is recognised by
// density (7.85 g/cm^3; nothing else in the TMS is above 5), not by name.
double SteelFraction(const MaterialSteps &materials) {
  double steel = 0.0, total = 0.0;
  for (const auto &step : materials) {
    total += step.second;
    if (DensityGCm3(step.first) > 5.0) steel += step.second;
  }
  return total > 0.0 ? steel / total : 0.0;
}

// Walks the real material budget of a step, applying mean Bethe-Bloch
// energy loss to qpInOut and accumulating Lynch-Dahl multiple-scattering
// variance. Returns the resulting process-noise covariance contribution;
// qpVarianceOut carries the straggling-derived q/p variance for the caller
// to add onto cov(4,4) separately (kept out of the 5x5 here since it needs
// the FINAL dxdz/dydz, which the Wolin-Ho formula above already accounts
// for through its own arguments).
// upstream: stepping toward lower z (the backward pass) -- the muon had MORE
// energy there, so the loss is added back instead of subtracted.
TMatrixD ApplyMaterialSteps(const MaterialSteps &materials,
                             double dxdz, double dydz, double &qpInOut,
                             double &qpVarianceOut, bool &rangedOutOut, bool upstream) {
  TMatrixD scatterCov(5, 5);
  qpVarianceOut = 0.0;
  rangedOutOut = false;

  const double chargeSign = (qpInOut >= 0.0) ? 1.0 : -1.0;
  double momentum = (std::abs(qpInOut) > 1e-12) ? 1.0 / std::abs(qpInOut) : 1.0;
  double energy = std::sqrt(momentum * momentum + BetheBloch_Utils::Mm * BetheBloch_Utils::Mm);
  const double energyFloor = std::sqrt(kMinMomentumMeV * kMinMomentumMeV +
                                        BetheBloch_Utils::Mm * BetheBloch_Utils::Mm);

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
    const double density = DensityGCm3(materialStep.first);
    double thickness = materialStep.second / 10.0;  // mm -> cm
    thickness = TMS_Geom::GetInstance().Scale(thickness);

    try {
      Material matter(density);
      bethe.fMaterial = matter;
      msc.fMaterial = matter;
    } catch (const std::invalid_argument &) {
      continue;  // unrecognised material at this step -- skip it, as legacy does
    }

    totalPathLengthGcm2 += density * thickness;

    // Walking forward (low->high z) energy decreases; walking upstream (the
    // backward pass) it is restored.
    if (upstream) {
      energy += bethe.Calc_dEdx(energy) * density * thickness;
    } else {
      energy -= bethe.Calc_dEdx(energy) * density * thickness;
    }
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
                          bool stopOnRangeOut, bool varianceGuard = true) {
  StepState predicted = previous;
  predicted.z = previous.z + subDz;

  const TVector3 midpoint(previous.x + 0.5 * previous.dxdz * subDz,
                           previous.y + 0.5 * previous.dydz * subDz,
                           previous.z + 0.5 * subDz);
  // Materials along the straight-line step (the curved path differs by well
  // under a mm over a substep), used both for the field -- which acts in the
  // steel only, so the step bends by the field times its steel fraction --
  // and below for energy loss and scattering.
  const TVector3 startPos(previous.x, previous.y, previous.z);
  const TVector3 straightEnd(previous.x + previous.dxdz * subDz, previous.y + previous.dydz * subDz,
                             previous.z + subDz);
  const MaterialSteps materials = TMS_Geom::GetInstance().GetMaterials(startPos, straightEnd);
  const double fieldY = field.GetField(midpoint).Y() * SteelFraction(materials);
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
  // (The backward pass turns this off: it refits a fixed list of
  // measurements, gating nothing, and near a stopping track's end -- tens of
  // MeV -- genuine multiple scattering exceeds this bound within one steel
  // layer; the next measurement sets the position anyway.)
  constexpr double kMaxPositionVarianceMM2 = 4.0e6;  // (2000mm)^2
  if (varianceGuard &&
      (!(propagatedCov(0, 0) <= kMaxPositionVarianceMM2) || !(propagatedCov(1, 1) <= kMaxPositionVarianceMM2))) {
    predicted.Diverged = true;
    return predicted;
  }

  double qpVariance = 0.0;
  bool rangedOut = false;
  const TMatrixD scatterCov = ApplyMaterialSteps(materials, predicted.dxdz,
                                                  predicted.dydz, predicted.qp, qpVariance, rangedOut,
                                                  /*upstream=*/subDz < 0.0);
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
                   double maxSubstepLengthMM, bool stopOnRangeOut, bool varianceGuard = true) {
  const double totalDz = zTarget - previous.z;
  if (std::abs(totalDz) < 1e-9) return previous;

  const int nSubsteps = std::max(
      1, static_cast<int>(std::ceil(std::abs(totalDz) / maxSubstepLengthMM)));
  const double subDz = totalDz / nSubsteps;

  StepState current = previous;
  for (int i = 0; i < nSubsteps; ++i) {
    current = PredictSubstep(current, subDz, field, stopOnRangeOut, varianceGuard);
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

// --- Hit-level measurement model (Config::MeasurementModel::Hits). ---

// Indices (into StepState) of the position and slope a hit measures.
int PositionIndex(const FitHit &hit) { return hit.MeasuresX ? 0 : 1; }
int SlopeIndex(const FitHit &hit) { return hit.MeasuresX ? 2 : 3; }

// Predicted coordinate at the hit's own plane, by straight-line transport
// from the state's z (field and material over the few cm between a point
// layer's z and its hits' planes are neglected here -- this is only used to
// SCORE candidates; the update itself steps properly), and the variance of
// the residual.
void ProjectToHit(const StepState &state, const FitHit &hit, double &predictedOut, double &residualVarOut) {
  const double dz = hit.Z - state.z;
  const int ip = PositionIndex(hit), is = SlopeIndex(hit);
  predictedOut = hit.MeasuresX ? state.x + state.dxdz * dz : state.y + state.dydz * dz;
  residualVarOut = state.cov(ip, ip) + 2.0 * dz * state.cov(ip, is) + dz * dz * state.cov(is, is) +
                   hit.SigmaMM * hit.SigmaMM + 1e-6;
}

double HitChi2(const StepState &state, const FitHit &hit) {
  double predicted = 0.0, var = 0.0;
  ProjectToHit(state, hit, predicted, var);
  const double r = hit.Coordinate - predicted;
  return r * r / var;
}

// Kalman update with one hit as a 1D measurement; predicted must already be
// at the hit's z. H picks out x or y, so the gain is that column of the
// covariance over the residual variance.
StepState UpdateWithHit(const StepState &predicted, const FitHit &hit, double &residualOut,
                        double &residualVarOut) {
  StepState updated = predicted;
  const int i = PositionIndex(hit);
  residualOut = hit.Coordinate - (hit.MeasuresX ? predicted.x : predicted.y);
  residualVarOut = predicted.cov(i, i) + hit.SigmaMM * hit.SigmaMM + 1e-6;
  double gain[5];
  for (int row = 0; row < 5; ++row) gain[row] = predicted.cov(row, i) / residualVarOut;
  double stateVec[5] = {predicted.x, predicted.y, predicted.dxdz, predicted.dydz, predicted.qp};
  for (int row = 0; row < 5; ++row) stateVec[row] += gain[row] * residualOut;
  updated.x = stateVec[0];
  updated.y = stateVec[1];
  updated.dxdz = stateVec[2];
  updated.dydz = stateVec[3];
  updated.qp = stateVec[4];
  ClampMomentum(updated.qp);
  for (int row = 0; row < 5; ++row)
    for (int col = 0; col < 5; ++col) updated.cov(row, col) -= gain[row] * predicted.cov(i, col);
  return updated;
}

// The two hits of a space point, if both indices are valid for this hit list.
bool PointHits(const TMS_SpacePoint &point, const std::vector<FitHit> &hits, int &xBarIndex, int &yBarIndex) {
  xBarIndex = point.GetXHitIndex();
  yBarIndex = point.GetYHitIndex();
  const int n = static_cast<int>(hits.size());
  return xBarIndex >= 0 && xBarIndex < n && yBarIndex >= 0 && yBarIndex < n;
}

struct GateResult {
  bool Accepted = false;
  std::size_t ChosenIndex = 0;
  std::vector<std::size_t> CandidateIndices;
  std::vector<double> CandidateChi2;
  std::vector<double> CandidateTimeChi2;  // filled only when time is in use
  std::vector<double> CandidateXYTimeChi2;  // filled only when X/Y time is in use
};

// The ambiguity-resolution core: score every candidate at this layer
// against the predicted state, accept the best one under the chi2 gate (or
// none, if nothing passes -- a gap, handled by the caller).
// Optional time information for ResolveLayer(): the running track t0 and
// its variance, plus where along the track this layer sits (path length).
struct TimeContext {
  bool Use = false;
  double T0 = 0.0;           // ns, mean of (t - s/c) over accepted points
  double ResidualVar = 0.0;  // ns^2, sigma_t^2 + var(T0)
  double PathLengthMM = 0.0; // s at this layer
  double GateNSigma = 0.0;   // 0 = no time gate
  // X/Y hit-time agreement of each candidate itself (Config::UseXYTimeInSelection);
  // null = not in use.
  const XYTimeDifferenceFn *XYTimeDifference = nullptr;
  double XYTimeVar = 0.0;         // ns^2
  double XYTimeGateNSigma = 0.0;  // 0 = no X/Y time gate
};

// hits: non-null for the hit-level model -- each candidate is then scored on
// its two hits at their own planes (2 DoF, as for the point) instead of on
// the point itself.
GateResult ResolveLayer(const StepState &predicted, const std::vector<std::size_t> &candidatesAtLayer,
                         const std::vector<TMS_SpacePoint> &allSpacePoints,
                         double barPitchMM, double chiSquareGateMax, const TimeContext &time,
                         const std::vector<FitHit> *hits) {
  GateResult result;
  double bestScore = std::numeric_limits<double>::infinity();
  for (std::size_t index : candidatesAtLayer) {
    const TMS_SpacePoint &candidate = allSpacePoints[index];
    const TMatrixD measurementCov = BuildMeasurementCovariance(barPitchMM);
    int xBarIndex = -1, yBarIndex = -1;
    const double chi2 = (hits != nullptr && PointHits(candidate, *hits, xBarIndex, yBarIndex))
        ? HitChi2(predicted, (*hits)[xBarIndex]) + HitChi2(predicted, (*hits)[yBarIndex])
        : Chi2(predicted, candidate, measurementCov);
    result.CandidateIndices.push_back(index);
    result.CandidateChi2.push_back(chi2);
    // Selection score: position chi2, plus the time chi2 when enabled. The
    // gate is still applied to position chi2 alone (and optionally a
    // separate time cut), so turning time on changes WHICH passing
    // candidate wins, not how permissive the gate is.
    double score = chi2;
    bool passes = chi2 <= chiSquareGateMax;
    if (time.Use) {
      const double r = (candidate.GetTime() - time.PathLengthMM / kSpeedOfLightMMPerNs) - time.T0;
      const double timeChi2 = r * r / time.ResidualVar;
      result.CandidateTimeChi2.push_back(timeChi2);
      score += timeChi2;
      if (time.GateNSigma > 0.0 && timeChi2 > time.GateNSigma * time.GateNSigma) passes = false;
    }
    if (time.XYTimeDifference != nullptr) {
      double dt = 0.0;
      const double xyChi2 = (*time.XYTimeDifference)(candidate, dt) ? dt * dt / time.XYTimeVar : 0.0;
      result.CandidateXYTimeChi2.push_back(xyChi2);
      score += xyChi2;
      if (time.XYTimeGateNSigma > 0.0 && xyChi2 > time.XYTimeGateNSigma * time.XYTimeGateNSigma) passes = false;
    }
    if (passes && score < bestScore) {
      bestScore = score;
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
bool IsBetterFit(const FitResult &a, const FitResult &b, bool rankByConvergence) {
  if (rankByConvergence && a.Converged != b.Converged) return a.Converged;
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
  return RunImpl(allSpacePoints, seedPath, 0.0, 0.0);
}

FitResult Follower::RunImpl(const std::vector<TMS_SpacePoint> &allSpacePoints, const std::vector<std::size_t> &seedPath,
                            double seedMomentumOverrideMeV, double qpRelSigmaOverride) const {
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
  // Last layer within MaxDistanceBeyondSeedMM of the seed's own last point.
  const double seedEndZ = allSpacePoints[zLayers[seedEndLayer].front()].GetZ();
  std::size_t lastLayerToWalk = seedEndLayer;
  while (lastLayerToWalk + 1 < zLayers.size() &&
         allSpacePoints[zLayers[lastLayerToWalk + 1].front()].GetZ() - seedEndZ <= fConfig.MaxDistanceBeyondSeedMM)
    ++lastLayerToWalk;

  // Hit-level model (Config::Measurement): the hits to update with, and the
  // ones already applied -- in the back section one y-measuring plane serves
  // the point layers on both sides of it, and its hit must enter the fit once.
  const std::vector<FitHit> *hits =
      (fConfig.Measurement == Config::MeasurementModel::Hits) ? fHits : nullptr;
  std::set<int> appliedHits;
  // Applies a point's two hits to state as 1D measurements in z order,
  // stepping field and material to each hit's own plane, and records each in
  // node.Hits. A hit already applied, or lying behind the state (possible
  // only if the walk skipped back), is recorded but not applied. Returns
  // false (stopOut set) if stepping to a hit ranged out or diverged.
  const auto applyPointHits = [&](int xBarIndex, int yBarIndex, StepState &state, FollowedNode &node,
                                  FitResult::StopReason &stopOut) -> bool {
    int order[2] = {xBarIndex, yBarIndex};
    if ((*hits)[yBarIndex].Z < (*hits)[xBarIndex].Z) std::swap(order[0], order[1]);
    for (int index : order) {
      const FitHit &hit = (*hits)[index];
      FollowedNode::HitUpdate record;
      record.HitIndex = index;
      record.Z = hit.Z;
      if (appliedHits.count(index) || hit.Z < state.z - 1e-3) {
        double predictedCoordinate = 0.0;
        ProjectToHit(state, hit, predictedCoordinate, record.ResidualVar);
        record.Residual = hit.Coordinate - predictedCoordinate;
        node.Hits.push_back(record);
        continue;
      }
      const StepState atHit = Predict(state, hit.Z, fField, fConfig.MaxSubstepLengthMM, fConfig.StopOnRangeOut);
      if (atHit.RangedOut) {
        stopOut = FitResult::StopReason::RangedOut;
        return false;
      }
      if (atHit.Diverged) {
        stopOut = FitResult::StopReason::Diverged;
        return false;
      }
      state = UpdateWithHit(atHit, hit, record.Residual, record.ResidualVar);
      record.Applied = true;
      appliedHits.insert(index);
      result.TotalChi2 += record.Residual * record.Residual / record.ResidualVar;
      result.NDoF += 1;
      node.Hits.push_back(record);
    }
    return true;
  };

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
  if (seedMomentumOverrideMeV > 0.0) seedMomentum = seedMomentumOverrideMeV;
  current.qp = chargeSign / seedMomentum;
  current.cov.Zero();
  current.cov(0, 0) = fConfig.InitialCovXX;
  current.cov(1, 1) = fConfig.InitialCovYY;
  current.cov(2, 2) = fConfig.InitialCovDXDZDXDZ;
  current.cov(3, 3) = fConfig.InitialCovDYDZDYDZ;
  current.cov(4, 4) = fConfig.InitialCovQPQP;
  const double qpRelSigma = qpRelSigmaOverride > 0.0 ? qpRelSigmaOverride : fConfig.InitialQPRelSigma;
  if (qpRelSigma > 0.0) {
    const double sigmaQP = qpRelSigma / seedMomentum;
    current.cov(4, 4) = sigmaQP * sigmaQP;
  }

  // Hit-level model: start at the seed point's first hit plane and apply both
  // its hits (the wide initial covariance lets them set the position).
  FollowedNode seedHits;
  int seedXBar = -1, seedYBar = -1;
  if (hits != nullptr && PointHits(firstPoint, *hits, seedXBar, seedYBar)) {
    current.z = std::min((*hits)[seedXBar].Z, (*hits)[seedYBar].Z);
    FitResult::StopReason unused = FitResult::StopReason::NotStarted;
    applyPointHits(seedXBar, seedYBar, current, seedHits, unused);
  }

  FollowedNode firstNode;
  firstNode.Layer = startLayer;
  firstNode.Z = firstPoint.GetZ();
  firstNode.Hits = seedHits.Hits;
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
  firstNode.FilteredZ = current.z;
  result.Nodes.push_back(firstNode);
  if (zLayers[startLayer].size() > 1) ++result.NAmbiguousLayersResolved;

  // Running time origin t0 = mean of (t - s/c) over accepted points, s =
  // path length from the seed point (see Config::UseTimeInSelection).
  // Tracked even when time is off, so FitResult::TrackT0Ns is always filled.
  double pathLengthMM = 0.0;
  double t0Sum = firstPoint.GetTime();
  int t0Count = 1;

  // z of the last accepted point, for the MaxGapMM limit.
  double lastAcceptedZ = firstPoint.GetZ();
  result.Converged = true;
  result.Stop = FitResult::StopReason::ReachedRangeEnd;  // overridden below if the walk breaks early

  for (std::size_t layerIdx = startLayer + 1; layerIdx <= lastLayerToWalk; ++layerIdx) {
    const std::vector<std::size_t> &candidates = zLayers[layerIdx];
    if (candidates.empty()) continue;  // TMS_LayerGrouping never emits an empty layer; defensive only
    const double targetZ = allSpacePoints[candidates.front()].GetZ();
    // Path length to this layer along the current direction estimate.
    pathLengthMM += std::abs(targetZ - current.z) *
                    std::sqrt(1.0 + current.dxdz * current.dxdz + current.dydz * current.dydz);

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
    TimeContext time;
    time.Use = fConfig.UseTimeInSelection;
    time.T0 = t0Sum / t0Count;
    time.ResidualVar = fConfig.TimeSigmaNs * fConfig.TimeSigmaNs * (1.0 + 1.0 / t0Count);
    time.PathLengthMM = pathLengthMM;
    time.GateNSigma = fConfig.TimeGateNSigma;
    if (fConfig.UseXYTimeInSelection && fXYTimeDifference) {
      time.XYTimeDifference = &fXYTimeDifference;
      time.XYTimeVar = fConfig.XYTimeSigmaNs * fConfig.XYTimeSigmaNs;
      time.XYTimeGateNSigma = fConfig.XYTimeGateNSigma;
    }
    const GateResult gate = ResolveLayer(predicted, candidates, allSpacePoints, fConfig.AssumedBarPitchMM,
                                         fConfig.ChiSquareGateMax, time, hits);

    FollowedNode node;
    node.Layer = layerIdx;
    node.Z = targetZ;
    node.CandidateIndices = gate.CandidateIndices;
    node.CandidateChi2 = gate.CandidateChi2;
    node.CandidateTimeChi2 = gate.CandidateTimeChi2;
    node.CandidateXYTimeChi2 = gate.CandidateXYTimeChi2;

    bool stopInsideLayer = false;
    if (gate.Accepted) {
      const TMS_SpacePoint &chosen = allSpacePoints[gate.ChosenIndex];
      node.HasHit = true;
      node.ChosenSpacePointIndex = gate.ChosenIndex;
      for (std::size_t i = 0; i < gate.CandidateIndices.size(); ++i) {
        if (gate.CandidateIndices[i] == gate.ChosenIndex) {
          node.Chi2AtChosen = gate.CandidateChi2[i];
          break;
        }
      }
      int xBarIndex = -1, yBarIndex = -1;
      if (hits != nullptr && PointHits(chosen, *hits, xBarIndex, yBarIndex)) {
        // Step from the last applied hit (not from the layer's z) to each hit.
        FitResult::StopReason stop = result.Stop;
        if (!applyPointHits(xBarIndex, yBarIndex, current, node, stop)) {
          stopInsideLayer = true;
          result.Stop = stop;
          if (stop == FitResult::StopReason::Diverged) result.Converged = false;
        }
      } else {
        const TMatrixD measurementCov = BuildMeasurementCovariance(fConfig.AssumedBarPitchMM);
        current = UpdateState(predicted, chosen, measurementCov);
        result.TotalChi2 += node.Chi2AtChosen;
        result.NDoF += 2;
      }
      t0Sum += chosen.GetTime() - pathLengthMM / kSpeedOfLightMMPerNs;
      ++t0Count;
      lastAcceptedZ = targetZ;
      if (candidates.size() > 1) ++result.NAmbiguousLayersResolved;
    } else {
      current = predicted;  // no update -- the inflated predicted covariance simply carries forward
      node.HasHit = false;
      ++result.NGapsFilled;
      if (targetZ - lastAcceptedZ > fConfig.MaxGapMM) {
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
    node.FilteredZ = current.z;
    result.Nodes.push_back(node);

    if (!result.Converged || stopInsideLayer) break;
  }

  // Extension past the walk's end on single hits (Config::ExtendOnHits).
  bool haveExtensionEnd = false;
  StepState extensionEndState;
  std::set<int> extensionHits;  // hit indices the extension took
  if (hits != nullptr && fConfig.ExtendOnHits && !appliedHits.empty()) {
    const FollowedNode *lastNode = nullptr;
    for (const FollowedNode &node : result.Nodes)
      if (node.HasHit) lastNode = &node;
    if (lastNode != nullptr) {
      StepState state;
      state.x = lastNode->FilteredX;
      state.y = lastNode->FilteredY;
      state.z = lastNode->FilteredZ != 0.0 ? lastNode->FilteredZ : lastNode->Z;
      state.dxdz = lastNode->FilteredDXDZ;
      state.dydz = lastNode->FilteredDYDZ;
      state.qp = lastNode->FilteredQP;
      state.cov = lastNode->FilteredCovariance;
      double lastHitZ = -std::numeric_limits<double>::infinity();
      for (int index : appliedHits) lastHitZ = std::max(lastHitZ, (*hits)[index].Z);
      // Usable, unused hits beyond the last applied hit, by plane (z, view).
      std::map<std::pair<long long, bool>, std::vector<int>> planes;
      for (int index = 0; index < static_cast<int>(hits->size()); ++index) {
        const FitHit &hit = (*hits)[index];
        if (hit.Usable && !appliedHits.count(index) && hit.Z > lastHitZ + 1.0)
          planes[{std::llround(hit.Z), hit.MeasuresX}].push_back(index);
      }
      std::vector<std::pair<double, std::vector<int>>> ordered;
      for (auto &kv : planes) ordered.push_back({(*hits)[kv.second.front()].Z, kv.second});
      std::stable_sort(ordered.begin(), ordered.end(),
                       [](const std::pair<double, std::vector<int>> &a, const std::pair<double, std::vector<int>> &b) {
                         return a.first < b.first;
                       });
      const double trackT0 = t0Sum / t0Count;
      const double seedZ = firstPoint.GetZ();
      double lastTakenZ = lastHitZ;
      StepState extensionEnd;
      for (const auto &plane : ordered) {
        if (plane.first - lastTakenZ > fConfig.ExtendMaxGapMM) break;
        const StepState predicted = Predict(state, plane.first, fField, fConfig.MaxSubstepLengthMM, false);
        if (predicted.Diverged) break;
        // Passing hits at this plane, best first.
        std::vector<std::pair<double, FitResult::OrphanHit>> passing;
        for (int index : plane.second) {
          const FitHit &hit = (*hits)[index];
          if (fConfig.ExtendTimeWindowNs > 0.0) {
            const double pathMM = (hit.Z - seedZ) *
                std::sqrt(1.0 + predicted.dxdz * predicted.dxdz + predicted.dydz * predicted.dydz);
            if (std::abs(hit.Time - (trackT0 + pathMM / kSpeedOfLightMMPerNs)) > fConfig.ExtendTimeWindowNs) continue;
          }
          double predictedCoordinate = 0.0, residualVar = 0.0;
          ProjectToHit(predicted, hit, predictedCoordinate, residualVar);
          const double residual = hit.Coordinate - predictedCoordinate;
          const double chi2 = residual * residual / residualVar;
          if (chi2 > fConfig.ExtendChi2Max) continue;
          FitResult::OrphanHit taken;
          taken.HitIndex = index;
          taken.Z = hit.Z;
          taken.Residual = residual;
          taken.ResidualVar = residualVar;
          passing.push_back({chi2, taken});
        }
        if (passing.empty()) {
          state = predicted;  // carry the prediction; the gap check above ends the walk
          continue;
        }
        std::sort(passing.begin(), passing.end(),
                  [](const std::pair<double, FitResult::OrphanHit> &a, const std::pair<double, FitResult::OrphanHit> &b) {
                    return a.first < b.first;
                  });
        const FitHit &best = (*hits)[passing.front().second.HitIndex];
        const double barPitch = best.SigmaMM * std::sqrt(12.0);
        bool ambiguous = false;
        for (const auto &candidate : passing)
          if (std::abs((*hits)[candidate.second.HitIndex].Coordinate - best.Coordinate) > 1.5 * barPitch) ambiguous = true;
        if (ambiguous) {
          state = predicted;
          continue;
        }
        double residual = 0.0, residualVar = 0.0;
        const StepState updated = UpdateWithHit(predicted, best, residual, residualVar);
        // A one-view update can drag the unmeasured direction; a muon here
        // never has |slope| above ~0.65 (33 degrees), so an update that
        // leaves the state steeper than 1.5 is following something else --
        // stop the extension rather than take it.
        if (!(std::abs(updated.dxdz) < 1.5 && std::abs(updated.dydz) < 1.5)) break;
        state = updated;
        result.Orphans.push_back(passing.front().second);
        appliedHits.insert(passing.front().second.HitIndex);  // keeps orphan pickup from re-taking it
        extensionHits.insert(passing.front().second.HitIndex);
        ++result.NExtensionHits;
        result.ExtensionEndX = state.x;
        result.ExtensionEndY = state.y;
        result.ExtensionEndZ = state.z;
        lastTakenZ = plane.first;
        extensionEnd = state;
      }
      if (result.NExtensionHits > 0) {
        haveExtensionEnd = true;
        extensionEndState = extensionEnd;
      }
    }
  }

  // Orphan-hit pickup (Config::PickUpOrphanHits): hits the track crosses
  // that are in no chosen space point. Each candidate is scored against the
  // filtered state nearest its plane, transported straight to the plane
  // (ProjectToHit), as candidates are scored during the walk.
  if (hits != nullptr && fConfig.PickUpOrphanHits && !appliedHits.empty()) {
    std::vector<StepState> states;
    for (const FollowedNode &node : result.Nodes) {
      if (!node.HasHit) continue;
      StepState state;
      state.x = node.FilteredX;
      state.y = node.FilteredY;
      state.z = node.FilteredZ;
      state.dxdz = node.FilteredDXDZ;
      state.dydz = node.FilteredDYDZ;
      state.qp = node.FilteredQP;
      state.cov = node.FilteredCovariance;
      states.push_back(state);
    }
    double zLow = std::numeric_limits<double>::infinity(), zHigh = -zLow;
    for (int index : appliedHits) {
      zLow = std::min(zLow, (*hits)[index].Z);
      zHigh = std::max(zHigh, (*hits)[index].Z);
    }
    zLow -= fConfig.OrphanZMarginMM;
    zHigh += fConfig.OrphanZMarginMM;
    // Passing candidates per plane (plane z, and which coordinate it measures).
    std::map<std::pair<long long, bool>, std::vector<std::pair<double, FitResult::OrphanHit>>> passing;
    // Expected track time at a plane: the running t0 (defined at the seed
    // point, path length 0) plus the path length from there over c.
    const double trackT0 = t0Sum / t0Count;
    const double seedZ = firstPoint.GetZ();
    for (int index = 0; index < static_cast<int>(hits->size()); ++index) {
      const FitHit &hit = (*hits)[index];
      if (!hit.Usable || appliedHits.count(index) || hit.Z < zLow || hit.Z > zHigh || states.empty()) continue;
      const StepState *nearest = &states.front();
      for (const StepState &state : states)
        if (std::abs(state.z - hit.Z) < std::abs(nearest->z - hit.Z)) nearest = &state;
      if (fConfig.OrphanTimeWindowNs > 0.0) {
        const double pathMM = (hit.Z - seedZ) *
            std::sqrt(1.0 + nearest->dxdz * nearest->dxdz + nearest->dydz * nearest->dydz);
        if (std::abs(hit.Time - (trackT0 + pathMM / kSpeedOfLightMMPerNs)) > fConfig.OrphanTimeWindowNs) continue;
      }
      double predictedCoordinate = 0.0, residualVar = 0.0;
      ProjectToHit(*nearest, hit, predictedCoordinate, residualVar);
      const double residual = hit.Coordinate - predictedCoordinate;
      const double chi2 = residual * residual / residualVar;
      if (chi2 > fConfig.OrphanChi2Max) continue;
      FitResult::OrphanHit orphan;
      orphan.HitIndex = index;
      orphan.Z = hit.Z;
      orphan.Residual = residual;
      orphan.ResidualVar = residualVar;
      passing[{std::llround(hit.Z), hit.MeasuresX}].push_back({chi2, orphan});
    }
    for (auto &plane : passing) {
      std::vector<std::pair<double, FitResult::OrphanHit>> &candidates = plane.second;
      std::sort(candidates.begin(), candidates.end(),
                [](const std::pair<double, FitResult::OrphanHit> &a, const std::pair<double, FitResult::OrphanHit> &b) {
                  return a.first < b.first;
                });
      // The best hit, plus any other passing hit in the next bar over from it.
      const FitHit &best = (*hits)[candidates.front().second.HitIndex];
      const double barPitch = best.SigmaMM * std::sqrt(12.0);
      std::vector<FitResult::OrphanHit> taken;
      bool ambiguous = false;
      for (const auto &candidate : candidates) {
        const FitHit &hit = (*hits)[candidate.second.HitIndex];
        if (&hit == &best || std::abs(hit.Coordinate - best.Coordinate) <= 1.5 * barPitch)
          taken.push_back(candidate.second);
        else
          ambiguous = true;
      }
      if (ambiguous && fConfig.OrphanSkipAmbiguousPlanes) continue;
      result.Orphans.insert(result.Orphans.end(), taken.begin(), taken.end());
    }
  }

  // Backward pass: the forward filter's first node only knows the seed, and
  // its last node -- the only one informed by every measurement -- sits at
  // the track's END. Refit the same measurements (plus any orphan hits)
  // from last to first, in z order, starting from the forward result with
  // its covariance inflated (so it acts as a weak prior), stepping upstream
  // through the material; the state at the first measurement is then the
  // track-start estimate (momentum and charge at, e.g., the TMS entrance).
  if (fConfig.BackwardPass && result.NDoF > 0) {
    // Measurements as (z, hit index), or (z, -1 - point index) for a node the
    // walk updated with its space point.
    std::vector<std::pair<double, int>> measurements;
    for (const FollowedNode &node : result.Nodes) {
      if (!node.HasHit) continue;
      if (hits != nullptr && !node.Hits.empty()) {
        for (const FollowedNode::HitUpdate &update : node.Hits)
          if (update.Applied) measurements.push_back({(*hits)[update.HitIndex].Z, update.HitIndex});
      } else {
        const int pointIndex = static_cast<int>(node.ChosenSpacePointIndex);
        measurements.push_back({allSpacePoints[pointIndex].GetZ(), -1 - pointIndex});
      }
    }
    for (const FitResult::OrphanHit &orphan : result.Orphans) measurements.push_back({orphan.Z, orphan.HitIndex});
    std::stable_sort(measurements.begin(), measurements.end(),
                     [](const std::pair<double, int> &a, const std::pair<double, int> &b) { return a.first > b.first; });

    // One backward pass over measurements[first..], from a given state. With
    // a range walker, it also steps the range momentum along the pass's
    // trajectory (rangeOk false if that walk fails -- which only loses the
    // range momentum, never the pass); afterFirst, if given, receives the
    // state updated at measurements[first] -- the track's last measurement
    // when first = 0.
    const ZeroFieldModel noField;
    auto runBackward = [&](StepState back, std::size_t first, StepState *range, bool &rangeOk,
                           StepState *afterFirst, StepState &out) {
      rangeOk = range != nullptr;
      for (std::size_t m = first; m < measurements.size(); ++m) {
        const auto &measurement = measurements[m];
        const StepState at = Predict(back, measurement.first, fField, fConfig.MaxSubstepLengthMM, /*stopOnRangeOut=*/false,
                                     /*varianceGuard=*/false);
        if (at.Diverged) return false;
        if (measurement.second >= 0) {
          double residual = 0.0, residualVar = 0.0;
          back = UpdateWithHit(at, (*hits)[measurement.second], residual, residualVar);
        } else {
          back = UpdateState(at, allSpacePoints[-1 - measurement.second],
                             BuildMeasurementCovariance(fConfig.AssumedBarPitchMM));
        }
        if (m == first && afterFirst) *afterFirst = back;
        if (rangeOk) {
          if (m == first) {
            // The range walk starts AT the last measurement, stopping there --
            // not at the forward walk's final state, which can lie up to
            // MaxGapMM past it after gap layers.
            *range = back;
            range->qp = (back.qp >= 0.0 ? 1.0 : -1.0) / fConfig.RangeStopMomentumMeV;
          } else {
            // Energy loss only (no field, no measurement updates on q/p), then
            // back onto the fitted trajectory -- position, direction and
            // covariance: the walk's own covariance means nothing (at tens of
            // MeV its scattering term alone would trip Predict()'s variance
            // guard within a few steps).
            StepState stepped = Predict(*range, measurement.first, noField, fConfig.MaxSubstepLengthMM, false,
                                        /*varianceGuard=*/false);
            if (stepped.Diverged) {
              rangeOk = false;
            } else {
              stepped.x = back.x;
              stepped.y = back.y;
              stepped.dxdz = back.dxdz;
              stepped.dydz = back.dydz;
              stepped.cov = back.cov;
              *range = stepped;
            }
          }
        }
      }
      out = back;
      return true;
    };

    // Start from the filtered state AT the track's last measurement -- the
    // last accepted node, or the single-hit extension's last hit -- not from
    // the walk's final state: that can lie up to MaxGapMM past the last
    // measurement after gap layers, its covariance grown by the predictions
    // carried through them. Scaled by BackwardCovScale on top, stepping back
    // over that stretch tripped Predict()'s position-variance guard and lost
    // the whole pass for 9% of stopping muons (2026-09-27).
    // The pass starts from the filtered state AT the track's last measurement
    // (withExtension: the single-hit extension's last hit, else the last
    // accepted node) -- not from the walk's final state, which can lie up to
    // MaxGapMM past it after gap layers, with the covariance grown by the
    // predictions carried through them (that lost the whole pass for 9% of
    // stopping muons, 2026-09-27). An orphan hit past that point (orphan
    // pickup looks OrphanZMarginMM beyond the applied hits) is reached by a
    // straight move, not by stepping downstream through steel at a
    // ranged-out momentum.
    auto backwardStart = [&](bool withExtension) {
      StepState start = current;
      for (const FollowedNode &node : result.Nodes) {
        if (!node.HasHit) continue;
        start.x = node.FilteredX;
        start.y = node.FilteredY;
        start.z = node.FilteredZ != 0.0 ? node.FilteredZ : node.Z;
        start.dxdz = node.FilteredDXDZ;
        start.dydz = node.FilteredDYDZ;
        start.qp = node.FilteredQP;
        start.cov = node.FilteredCovariance;
      }
      if (withExtension && haveExtensionEnd) start = extensionEndState;
      start.Diverged = start.RangedOut = false;
      if (!measurements.empty() && measurements.front().first > start.z) {
        const double dz = measurements.front().first - start.z;
        start.x += start.dxdz * dz;
        start.y += start.dydz * dz;
        start.z = measurements.front().first;
      }
      start.cov *= fConfig.BackwardCovScale;
      return start;
    };
    StepState back = backwardStart(true);
    StepState range, atLast, backOut;
    bool rangeOk = false;
    bool ok = runBackward(back, 0, &range, rangeOk, &atLast, backOut);
    if (!ok && haveExtensionEnd) {
      // Fallback: without the extension's hits, from the last walked node.
      // The track keeps them; only its momenta come from the rest.
      measurements.erase(std::remove_if(measurements.begin(), measurements.end(),
                                        [&](const std::pair<double, int> &m) {
                                          return extensionHits.count(m.second) > 0;
                                        }),
                         measurements.end());
      back = backwardStart(false);
      ok = runBackward(back, 0, &range, rangeOk, &atLast, backOut);
    }
    back = backOut;
    if (ok && rangeOk) result.RangeMomentumMeV = std::abs(range.qp) > 1e-12 ? 1.0 / std::abs(range.qp) : 0.0;

    if (ok && fConfig.RangeSeededBackwardPass) {
      // Same pass from the last measurement on, but starting from the range
      // hypothesis: |p| = RangeStopMomentumMeV there, with a tight q/p prior.
      StepState seeded = atLast;
      seeded.qp = (atLast.qp >= 0.0 ? 1.0 : -1.0) / fConfig.RangeStopMomentumMeV;
      for (int k = 0; k < 5; ++k) seeded.cov(4, k) = seeded.cov(k, 4) = 0.0;
      const double sigmaQP = fConfig.RangeSeedQPRelSigma * std::abs(seeded.qp);
      seeded.cov(4, 4) = sigmaQP * sigmaQP;
      StepState seededOut;
      bool unused = false;
      if (runBackward(seeded, 1, nullptr, unused, nullptr, seededOut))
        result.RangeSeededMomentumMeV = std::abs(seededOut.qp) > 1e-12 ? 1.0 / std::abs(seededOut.qp) : 0.0;
    }
    if (ok) {
      result.HasStartState = true;
      result.StartX = back.x;
      result.StartY = back.y;
      result.StartZ = back.z;
      result.StartDXDZ = back.dxdz;
      result.StartDYDZ = back.dydz;
      result.StartMomentumMeV = (std::abs(back.qp) > 1e-12) ? 1.0 / std::abs(back.qp) : 0.0;
      result.StartCharge = (back.qp >= 0.0) ? 1.0 : -1.0;
    }
  }

  result.MomentumMeV = (std::abs(current.qp) > 1e-12) ? 1.0 / std::abs(current.qp) : 0.0;
  result.Charge = (current.qp >= 0.0) ? 1.0 : -1.0;
  result.NDoF -= 5;  // 5 fitted state parameters
  result.TrackT0Ns = t0Sum / t0Count;

  return result;
}

FitResult Follower::RunBestSeed(const std::vector<TMS_SpacePoint> &allSpacePoints,
                                 const std::vector<std::size_t> &objectIndices,
                                 std::vector<FitResult> *allHypotheses,
                                 std::size_t *bestIndex) const {
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
  std::vector<std::size_t> bestSeedPath;
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

      FitResult candidate = Run(allSpacePoints, seedPath);
      candidate.HeadSkip = skip;
      if (allHypotheses) allHypotheses->push_back(candidate);
      if (!haveBest || IsBetterFit(candidate, best, fConfig.RankHypothesesByConvergence)) {
        best = candidate;
        bestSeedPath = seedPath;
        haveBest = true;
        if (bestIndex && allHypotheses) *bestIndex = allHypotheses->size() - 1;
      }
    }

    // Triplet hypotheses: explore genuinely different (layer0,layer1,layer2)
    // combinations, not just which point anchors the same fixed tail (see
    // the Config::MaxTripletHypotheses comment). Needs a third layer beyond
    // this skip level's first two.
    if (fConfig.MaxTripletHypotheses > 0 && objectLayers.size() >= static_cast<std::size_t>(skip) + 3) {
      const std::vector<std::size_t> &layer0 = firstLayerLocal;
      const std::vector<std::size_t> &layer1 = objectLayers[skip + 1];
      const std::vector<std::size_t> &layer2 = objectLayers[skip + 2];

      struct TripletCandidate {
        double residual;
        std::size_t a, b, c;  // global indices into allSpacePoints
      };
      std::vector<TripletCandidate> tripletCandidates;
      tripletCandidates.reserve(layer0.size() * layer1.size() * layer2.size());
      for (std::size_t la : layer0) {
        const std::size_t ga = objectIndices[la];
        const TMS_SpacePoint &pa = allSpacePoints[ga];
        for (std::size_t lc : layer2) {
          const std::size_t gc = objectIndices[lc];
          const TMS_SpacePoint &pc = allSpacePoints[gc];
          const double dz = pc.GetZ() - pa.GetZ();
          if (std::abs(dz) < 1e-6) continue;
          for (std::size_t lb : layer1) {
            const std::size_t gb = objectIndices[lb];
            const TMS_SpacePoint &pb = allSpacePoints[gb];
            // Predict the middle point by linear interpolation in z between
            // the outer two, and score by how far the actual candidate sits
            // from that prediction -- cheap, no fit involved.
            const double t = (pb.GetZ() - pa.GetZ()) / dz;
            const double predX = pa.GetX() + t * (pc.GetX() - pa.GetX());
            const double predY = pa.GetY() + t * (pc.GetY() - pa.GetY());
            const double residual = std::hypot(pb.GetX() - predX, pb.GetY() - predY);
            if (residual <= fConfig.TripletCollinearityToleranceMM) {
              tripletCandidates.push_back({residual, ga, gb, gc});
            }
          }
        }
      }
      std::sort(tripletCandidates.begin(), tripletCandidates.end(),
                [](const TripletCandidate &lhs, const TripletCandidate &rhs) {
                  return lhs.residual < rhs.residual;
                });
      const std::size_t nTriplets =
          std::min(tripletCandidates.size(), static_cast<std::size_t>(fConfig.MaxTripletHypotheses));
      for (std::size_t t = 0; t < nTriplets; ++t) {
        const TripletCandidate &tri = tripletCandidates[t];
        std::vector<std::size_t> seedPath;
        seedPath.reserve(restSorted.size() + 3);
        seedPath.push_back(tri.a);
        seedPath.push_back(tri.b);
        seedPath.push_back(tri.c);
        for (std::size_t idx : restSorted) {
          if (idx == tri.b || idx == tri.c) continue;  // already placed above
          seedPath.push_back(idx);
        }

        FitResult candidate = Run(allSpacePoints, seedPath);
        candidate.HeadSkip = skip;
        if (allHypotheses) allHypotheses->push_back(candidate);
        if (!haveBest || IsBetterFit(candidate, best, fConfig.RankHypothesesByConvergence)) {
          best = candidate;
          bestSeedPath = seedPath;
          haveBest = true;
          if (bestIndex && allHypotheses) *bestIndex = allHypotheses->size() - 1;
        }
      }
    }
  }

  // Range re-seed (Config::RangeReseedFactor): walk the best seed again with
  // its own range momentum (times a margin) as the seed, while that extends
  // the track.
  if (haveBest && fConfig.RangeReseedFactor > 0.0) {
    auto lastHitZ = [](const FitResult &fit) {
      double z = -1e30;
      for (const FollowedNode &node : fit.Nodes)
        if (node.HasHit) z = node.FilteredZ != 0.0 ? node.FilteredZ : node.Z;
      return z;
    };
    for (int pass = 0; pass < fConfig.RangeReseedMaxPasses; ++pass) {
      if (!(best.RangeMomentumMeV > 0.0)) break;
      FitResult again = RunImpl(allSpacePoints, bestSeedPath, fConfig.RangeReseedFactor * best.RangeMomentumMeV,
                                fConfig.RangeReseedQPRelSigma);
      again.HeadSkip = best.HeadSkip;
      // Keep it only if it reaches further downstream and does not rank worse.
      if (!(lastHitZ(again) > lastHitZ(best) + 1.0) || IsBetterFit(best, again, fConfig.RankHypothesesByConvergence)) break;
      best = again;
      if (allHypotheses) {
        allHypotheses->push_back(again);
        if (bestIndex) *bestIndex = allHypotheses->size() - 1;
      }
    }
  }
  return best;
}

}  // namespace TMS_KalmanFollower
