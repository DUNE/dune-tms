#include "TMS_Kalman.h"
#include "TMS_Geom.h"
#include "TMS_Manager.h"
#include <algorithm>
#include <limits>

TMS_Kalman::TMS_Kalman() : 
  Bethe(Material::kPolyStyrene),
  MSC(Material::kPolyStyrene),
  ForwardFitting(false),
  total_en(0.0),
  mass(0.0),
  momentum(0.0),
  charge_curvature(0.0),
  assumed_charge(0.0),
  AverageXSlope(0.0),
  AverageYSlope(0.0),
  Talk(false) {
}

// Helper to inflate the initial state covariance so the filter remains
// conservative about the starting state (large position/slope variances,
// modest curvature variance). Functional values unchanged.
namespace {
  inline void SetInflatedInitialCovariance(TMatrixD &P) {
    P.Zero();
    // Position (mm^2): large uncertainty to downweight initial estimate
    const double cov_x  = 1.0e4;  // std ~ 100 mm
    const double cov_y  = 1.0e4;  // std ~ 100 mm
    // Slopes ((dx/dz)^2, (dy/dz)^2): very uninformative initial slope
    const double cov_tx = 25.0;   // std ~ 5
    const double cov_ty = 25.0;   // std ~ 5
    // Curvature ((q/p)^2): large to avoid overconfidence on momentum
    const double cov_qp = 100.0;  // very conservative

    P(0,0) = cov_x;
    P(1,1) = cov_y;
    P(2,2) = cov_tx;
    P(3,3) = cov_ty;
    P(4,4) = cov_qp;
  }

  inline double AdaptiveChi2Threshold(const TMatrixD &S, int step_index, double base_thresh) {
    double factor = 1.0;
    const auto &mgr = TMS_Manager::GetInstance();
    double early_mult = mgr.Get_Reco_Kalman_Outlier_EarlyStepMultiplier();
    int    early_n    = mgr.Get_Reco_Kalman_Outlier_EarlyStepsCount();
    bool   ramp_en    = mgr.Get_Reco_Kalman_Outlier_EarlyRampEnable();
    int    ramp_n     = mgr.Get_Reco_Kalman_Outlier_EarlyRampSteps();
    double k_trace    = mgr.Get_Reco_Kalman_Outlier_InnovationTraceCoeff();
    double max_scale  = mgr.Get_Reco_Kalman_Outlier_MaxScale();

    if (step_index >= 0) {
      if (ramp_en && ramp_n > 0 && step_index < ramp_n) {
        double t = (ramp_n > 1) ? static_cast<double>(step_index) / static_cast<double>(ramp_n - 1) : 1.0;
        double ramp_factor = early_mult + (1.0 - early_mult) * std::clamp(t, 0.0, 1.0);
        factor = std::max(factor, ramp_factor);
      } else if (step_index < early_n) {
        factor = std::max(factor, early_mult);
      }
    }

    if (S.GetNrows() == 2 && S.GetNcols() == 2) {
      double traceS = S(0,0) + S(1,1);
      double conf_factor = 1.0 + k_trace * std::max(0.0, traceS);
      factor = std::max(factor, conf_factor);
    }
    if (std::isnan(factor) || std::isinf(factor)) factor = 1.0;
    factor = std::min(factor, max_scale);
    return base_thresh * factor;
  }
}

// Construct Kalman nodes from an initial set of hits.
TMS_Kalman::TMS_Kalman(std::vector<TMS_Hit> &Candidates, double charge) : 
  Bethe(Material::kPolyStyrene), 
  MSC(Material::kPolyStyrene),
  ForwardFitting(false),
  assumed_charge(charge),
  Talk(false)
{
//  TRandom3* RNG = new TRandom3(1337); // TODO: Seed properly sometime

  // Empty the KalmanStates
  KalmanNodes.clear();

  // Save the number of initial candidates
  int nCand = Candidates.size();
  
  // And muon mass
  mass = BetheBloch_Utils::Mm;

  std::sort(Candidates.begin(), Candidates.end(), TMS_Hit::SortByZ);

  // Make a new Kalman state for each hit
  // Each node i represents the plane at z_i, with dz_i = z_i - z_{i-1} (backward step).
  for (int i = 0; i < nCand; ++i) {

    TMS_Hit &hit = Candidates[i];
    double x_true = hit.GetTrueHit().GetX();
    double y_true = hit.GetTrueHit().GetY();
    double x = hit.GetRecoX();
    double y = hit.GetRecoY();
    double z = hit.GetZ();

    // Allow out-of-box hits silently; Kalman may still handle nudging
    TVector3 vecc(x,y,z);
    (void)vecc;

    // Determine dz as backward difference where possible; for first hit, use forward diff
    double dz = 0.0;
    if (i > 0) {
      dz = z - Candidates[i-1].GetZ();
    } else if (i+1 < nCand) {
      dz = Candidates[i+1].GetZ() - z;
    } else {
      dz = 0.0;
    }
    if (std::abs(dz) <= 1e-6) continue; // skip duplicate-z entries; TODO: merge hits per plane

    // Seed slopes from local difference (prefer backward)
    double dxdz = 0.0, dydz = 0.0;
    if (i > 0 && std::abs(dz) > 1e-6) {
      dxdz = (x - Candidates[i-1].GetRecoX()) / dz;
      dydz = (y - Candidates[i-1].GetRecoY()) / dz;
    } else if (i+1 < nCand) {
      double dz_f = Candidates[i+1].GetZ() - z;
      if (std::abs(dz_f) > 1e-6) {
        dxdz = (Candidates[i+1].GetRecoX() - x) / dz_f;
        dydz = (Candidates[i+1].GetRecoY() - y) / dz_f;
      }
    }

    TMS_KalmanNode Node(x, y, z, dz, dxdz, dydz);
    // Add truth for diagnostics and bar geometry for measurement/noise
    Node.SetTrueXY(x_true, y_true);
    Node.LayerOrientation = hit.GetBar().GetBarType();
    Node.LayerBarWidth    = hit.GetBar().GetBarWidth();
    Node.LayerBarLength   = hit.GetBar().GetBarLength();
    Node.LayerNotZ        = hit.GetBar().GetNotZ();
    Node.PlaneNumber      = hit.GetPlaneNumber();
    Node.Hit              = &hit;

    KalmanNodes.emplace_back(std::move(Node));
    KalmanNodes.back().StepIndex = KalmanNodes.size() - 1;
  }

  int N_LAYER_BACK = 10;
  // Can't look back further than the first element
  if (Candidates.size() < (unsigned)N_LAYER_BACK)
    N_LAYER_BACK = Candidates.size();

  double dz_x = (Candidates[Candidates.size() - N_LAYER_BACK].GetZ() - Candidates.back().GetZ());
  double dz_y = (Candidates.front().GetZ() - Candidates.back().GetZ());
  AverageXSlope = (std::abs(dz_x) > 1e-6)
                    ? (Candidates[Candidates.size() - N_LAYER_BACK].GetRecoX() - Candidates.back().GetRecoX()) / dz_x
                    : 0.0;
  AverageYSlope = (std::abs(dz_y) > 1e-6)
                    ? (Candidates.front().GetRecoY() - Candidates.back().GetRecoY()) / dz_y
                    : 0.0;


  // Set the momentum seed for the first hit from its length
  if (ForwardFitting) {
    double startx = Candidates.front().GetNotZ();
    double endx = Candidates.back().GetNotZ();
    double startz = Candidates.front().GetZ();
    double endz = Candidates.back().GetZ();
    double KEest = GetKEEstimateFromLength(startx, endx, startz, endz);
    double momest = sqrt((KEest+mass)*(KEest+mass)-mass*mass);
    if (Talk) std::cout << "momentum estimate from length: " << momest << std::endl;
    KalmanNodes.front().PreviousState.qp = 1./momest;
    KalmanNodes.front().CurrentState.qp = 1./momest;
    KalmanNodes.back().CurrentState.dxdz = 0.0;//AverageXSlope;
    KalmanNodes.back().CurrentState.dydz = 0.0;//AverageYSlope;
    KalmanNodes.back().PreviousState.dxdz = 0.0;//AverageXSlope;
    KalmanNodes.back().PreviousState.dydz = 0.0;//AverageYSlope;
  } else { // TODO check if 0 is sane
    KalmanNodes.back().CurrentState.dxdz = 0.0;//AverageXSlope;
    KalmanNodes.back().CurrentState.dydz = 0.0;//AverageYSlope;
    KalmanNodes.back().PreviousState.dxdz = 0.0;//AverageXSlope;
    KalmanNodes.back().PreviousState.dydz = 0.0;//AverageYSlope;
  }
  //InitializeMomentum();
  RunKalman();
  runRTSSmoother();
  BetheBloch();
  Runchi2();
  SignSelection();
  // Prune out nodes rejected as outliers and duplicates at the same z-layer
  PruneRejectedAndDuplicateZ();
}

void TMS_Kalman::InitializeMomentum(bool only_momentum) {
  // Sets the initial momentum of the first kalman node
  //SortNodesByRunOrder();
  int N_LAYER_BACK = 10;
  // Can't look back further than the first element
  if (KalmanNodes.size() < (unsigned)N_LAYER_BACK)
    N_LAYER_BACK = KalmanNodes.size();

  double dz_x = (KalmanNodes[KalmanNodes.size() - N_LAYER_BACK].z - KalmanNodes.back().z);
  double dz_y = (KalmanNodes.front().z - KalmanNodes.back().z);
  AverageXSlope = (std::abs(dz_x) > 1e-6)
                    ? (KalmanNodes[KalmanNodes.size() - N_LAYER_BACK].RecoX - KalmanNodes.back().RecoX) / dz_x
                    : 0.0;
  AverageYSlope = (std::abs(dz_y) > 1e-6)
                    ? (KalmanNodes.front().RecoY - KalmanNodes.back().RecoY) / dz_y
                    : 0.0;

  if (std::abs(mass) < 1e-3) { 
    std::cout<<"Fatal: particle mass is zero"<<std::endl;
    exit(-1);
  }
  // Set the momentum seed for the first hit from its length
  if (ForwardFitting) {
    // TODO GetNotZ might represent y position, which would give wrong numbers
    double startx = KalmanNodes.front().RecoX;
    double endx = KalmanNodes.back().RecoX;
    double startz = KalmanNodes.front().z;
    double endz = KalmanNodes.back().z;
    double KEest = GetKEEstimateFromLength(startx, endx, startz, endz);
    double momest = sqrt((KEest+mass)*(KEest+mass)-mass*mass);
    if (Talk) std::cout << "[InitializeMomentum] momentum estimate from length: " << momest << std::endl;
    KalmanNodes.front().PreviousState.qp = 1./momest;
    KalmanNodes.front().CurrentState.qp = 1./momest;
    if (Talk) std::cout << "[InitializeMomentum] momentum estimate from length as set: " << KalmanNodes[0].CurrentState.qp << std::endl;
    if (!only_momentum) {
      KalmanNodes.back().CurrentState.dxdz = 0.0;//AverageXSlope;
      KalmanNodes.back().CurrentState.dydz = 0.0;//AverageYSlope;
      KalmanNodes.back().PreviousState.dxdz = 0.0;//AverageXSlope;
      KalmanNodes.back().PreviousState.dydz = 0.0;//AverageYSlope;
    }
  } else {
    // todo: get better starting position in the future
    // But for now this is good enough since we're just calculating small correction
    // that's probably the same independent of x and y.
    double endz = KalmanNodes.back().z;
    TVector3 endpoint(0, 0, endz);
    double KEest = CalculateKEFromSteel(endpoint);
    double momest = std::max(sqrt((KEest+mass)*(KEest+mass)-mass*mass), 1e-3);
    if (Talk) std::cout << "[InitializeMomentum] KE est from steel " << KEest << std::endl;
    if (Talk) std::cout << "[InitializeMomentum] mass " << mass << std::endl;
    if (Talk) std::cout << "[InitializeMomentum] momentum estimate from half of steel: " << momest << std::endl;
    KalmanNodes.front().PreviousState.qp = 1./momest;
    KalmanNodes.front().CurrentState.qp = 1./momest;
    if (Talk) std::cout << "[InitializeMomentum] momentum estimate from length as set: " << KalmanNodes[0].CurrentState.qp << std::endl;
    // TODO check if 0 is sane
    if (!only_momentum) {
      KalmanNodes.back().CurrentState.dxdz = 0.0;//AverageXSlope;
      KalmanNodes.back().CurrentState.dydz = 0.0;//AverageYSlope;
      KalmanNodes.back().PreviousState.dxdz = 0.0;//AverageXSlope;
      KalmanNodes.back().PreviousState.dydz = 0.0;//AverageYSlope;
    }
  }
}

// Augment the track with additional candidate hits:
// - Seed from the smoothed state near the front (lowest z)
// - Propagate forward through a candidate pool
// - Accept candidates that pass gating/guards (one per plane)
// - Merge, repair/prune, and run the triple refit for consistency
void TMS_Kalman::AugmentWithCandidates(const std::vector<TMS_Hit> &candidate_pool, const size_t number_to_remove) {
  if (KalmanNodes.empty() || candidate_pool.empty()) return;

  // We will run a forward pass, so set this flag for physics (energy loss sign)
  bool prevForward = ForwardFitting;
  ForwardFitting = true;

  // Sort a local copy of candidates by z
  std::vector<TMS_Hit> pool(candidate_pool.begin(), candidate_pool.end());
  std::sort(pool.begin(), pool.end(), TMS_Hit::SortByZInc);

  // Ensure nodes are sorted by z ascending so back() is highest z
  SortNodesByZ();
  // Seed from the smoothed state at the highest-z node (back of vector)
  size_t new_length = KalmanNodes.size();
  size_t new_length_half = new_length / 2;
  for (size_t j = 0; j < number_to_remove && j < new_length_half; j++) {
    KalmanNodes.pop_back();
  }
  //KalmanNodes.resize(new_length);
  TMS_KalmanNode seed = KalmanNodes.back();
  TMS_KalmanNode prev_node(
      seed.SmoothState.x,
      seed.SmoothState.y,
      seed.SmoothState.z,
      seed.SmoothState.qp,
      seed.SmoothState.dxdz,
      seed.SmoothState.dydz);
  prev_node.SetTrueXY(seed.TrueX, seed.TrueY);
  prev_node.LayerOrientation = seed.LayerOrientation;
  prev_node.LayerBarWidth = seed.LayerBarWidth;
  prev_node.LayerBarLength = seed.LayerBarLength;
  prev_node.PlaneNumber = seed.PlaneNumber;
  prev_node.Hit = seed.Hit;
  prev_node.StepIndex = seed.StepIndex;
  // Carry covariance from the smoother as our starting uncertainty
  prev_node.CovarianceMatrix = KalmanNodes.back().SmoothCovarianceMatrix;
  prev_node.CurrentState = seed.SmoothState;
  prev_node.PreviousState = seed.SmoothState;
  // Boost initial momentum for forward pass using estimated energy gain
  // over a small forward z-window (range-based, not arbitrary scale).
  {
    double qp_seed = seed.SmoothState.qp;
    if (!std::isfinite(qp_seed) || std::abs(qp_seed) < 1e-12) qp_seed = 1.0/1000.0; // ~1 GeV
    double p_seed = 1.0 / qp_seed;
    double en_seed = std::sqrt(p_seed*p_seed + mass*mass);

    // Estimate additional KE the muon had earlier by integrating a short
    // straight segment ahead. This uses the same Bethe/geom
    // model as energy loss in Predict.
    // Skips any hits that are outside the window
    const double min_z = prev_node.SmoothState.z;
    const double max_window = TMS_Manager::GetInstance().Get_Reco_Kalman_Augment_MomentumWindowZ();
    double max_z = min_z;
    for (auto &hit : candidate_pool) {
      if (hit.GetZ() - min_z > max_window) continue; // outside window
      max_z = std::max(max_z, hit.GetZ());
    }
    const double window_z = max_z - min_z;
    TVector3 start(seed.SmoothState.x, seed.SmoothState.y, seed.SmoothState.z);
    TVector3 stop (seed.SmoothState.x, seed.SmoothState.y, seed.SmoothState.z + window_z);
    double dE_gain = 0.0;
    try {
      auto Materials = TMS_Geom::GetInstance().GetMaterials(start, stop);
      double en_tmp = en_seed;
      for (auto &mat : Materials) {
        double density = mat.first->GetDensity()/(CLHEP::g/CLHEP::cm3);
        double thickness = mat.second/10.;
        double scale_factor = TMS_Geom::GetInstance().Scale(1.0);
        density /= std::pow(scale_factor, 3);
        thickness = TMS_Geom::GetInstance().Scale(thickness);
        try {
          Material matter(density);
          Bethe.fMaterial = matter;
        } catch (...) {}
        double dE = Bethe.Calc_dEdx(en_tmp)*density*thickness;
        if (std::isfinite(dE) && dE > 0) dE_gain += dE;
      }
    } catch (...) {
      if (Talk) std::cout << "[Augment] Warning: Caught exception. Setting dE_gain to zero" << std::endl;
      dE_gain = 0.0;
    }
    // Convert back to momentum
    double en_aug = en_seed + dE_gain;
    if (!std::isfinite(en_aug) || en_aug <= mass) en_aug = en_seed;
    if (Talk) std::cout << "[Augment] Set en_aug to " << en_aug << ", up from en_seed " << en_seed << std::endl;
    double p2_aug = en_aug*en_aug - mass*mass;
    double p_aug = (p2_aug > 0.0) ? std::sqrt(p2_aug) : p_seed;
    // Clamp to sane muon range [0.2, 20] GeV
    double p_min = 200.0;   // MeV
    double p_max = 20000.0; // MeV
    p_aug = std::clamp(p_aug, p_min, p_max);
    double qp_aug = 1.0 / p_aug;
    prev_node.PreviousState.qp = qp_aug;
    prev_node.CurrentState.qp  = qp_aug;
  }
  prev_node.StepIndex = 0;
  // Strongly inflate the entire covariance so the augmentation pass
  // is humble about the smoothed seed. This reduces gain-induced
  // zig-zagging when adding planes with alternating orientations.
  SetInflatedInitialCovariance(prev_node.CovarianceMatrix);

  // Reset outlier streak and per-z acceptance map for this pass
  ConsecutiveOutliers = 0;
  ZLayerAccepted.clear();
  AugmentationMode = true;
  // Configure freeze steps from manager
  AugmentFreezeSteps = TMS_Manager::GetInstance().Get_Reco_Kalman_Augment_FreezeSteps();

  // Track current z to enforce monotonic forward propagation beyond the current end
  double current_z = prev_node.CurrentState.z;
  int step_index = 0;

  std::vector<TMS_KalmanNode> accepted_new;
  accepted_new.reserve(pool.size());

  int nCandidates = 0;
  bool EnableOutlierRejectionOld = EnableOutlierRejection;
  // Keep outlier rejection enabled with adaptive thresholding
  // to avoid mis-associations that can kick the filter.
  EnableOutlierRejection = true;

  // Debug accumulators
  int dbg_considered = 0, dbg_accepted = 0, dbg_rejected = 0;
  double dbg_sum_chi2_acc = 0.0, dbg_sum_chi2_rej = 0.0;
  double dbg_sum_rnorm_acc = 0.0, dbg_sum_rnorm_rej = 0.0;
  double dbg_sum_dp_acc = 0.0, dbg_sum_dp_rej = 0.0;

  // Highest-z candidate in the pool (pool is sorted ascending by z)
  double z_max_pool = pool.empty() ? -std::numeric_limits<double>::infinity() : pool.back().GetZ();

  for (auto &hit : pool) {
    double z = hit.GetZ();
    double dz = z - prev_node.CurrentState.z;
    if (z <= current_z + 1e-6) continue; // only forward propagation

    // Form a node at the candidate's z with a measurement (x,y) derived
    // from bar geometry and the current track projection.
    double x = 0.0, y = 0.0;
    BuildBarMeasurement(hit, prev_node.CurrentState, z, x, y);

    // Construct the node with plane z = hit z and dz = z - prev.z.
    // The node constructor defines PreviousState at z-dz and CurrentState at z,
    // matching the Kalman step semantics and avoiding z overshoot.
    TMS_KalmanNode node(x, y, z, dz, prev_node.CurrentState.dxdz, prev_node.CurrentState.dydz);
    double x_true = hit.GetTrueHit().GetX();
    double y_true = hit.GetTrueHit().GetY();
    node.SetTrueXY(x_true, y_true);
    node.LayerOrientation = hit.GetBar().GetBarType();
    node.LayerBarWidth    = hit.GetBar().GetBarWidth();
    node.LayerBarLength   = hit.GetBar().GetBarLength();
    node.LayerNotZ        = hit.GetBar().GetNotZ();
    node.PlaneNumber      = hit.GetPlaneNumber();
    node.Hit              = &hit;
    node.StepIndex        = ++step_index;
    
    // Skip if over a certain number of planes, but allow the very last pool plane
    bool is_last_plane = (std::abs(z - z_max_pool) < 1e-6);
    int maximum_n_planes = TMS_Manager::GetInstance().Get_Reco_Kalman_Augment_MaximumPlaneGap();
    bool over_gap = (maximum_n_planes > 0 && node.PlaneNumber > prev_node.PlaneNumber + maximum_n_planes);
    if (over_gap && !is_last_plane) continue;

    nCandidates += 1;
    dbg_considered += 1;

    // Propagate from prev_node state to this node and perform update/gating
    Update(prev_node, node);
    Predict(node);

    // If accepted by gating, keep it and advance the prev_node state
    if (node.Accepted) {
      accepted_new.emplace_back(node);
      prev_node = node; // advance
      current_z = prev_node.CurrentState.z;
      dbg_accepted += 1;
    } else {
      // Loosen constraints: if this is the last pool plane, allow an
      // acceptance override based on relaxed absolute residuals.
      if (is_last_plane) {
        const auto &mgr = TMS_Manager::GetInstance();
        double max_x = mgr.Get_Reco_Kalman_Outlier_AbsoluteResidualMaxX();
        double max_y = mgr.Get_Reco_Kalman_Outlier_AbsoluteResidualMaxY();
        double relax = mgr.Get_Reco_Kalman_Augment_FinalPlaneRelaxFactor();
        // Relaxed absolute residual ceilings
        double base_x = mgr.Get_Reco_Kalman_Augment_FinalPlaneRelaxAbsX();
        double base_y = mgr.Get_Reco_Kalman_Augment_FinalPlaneRelaxAbsY();
        double lax_x = (max_x > 0 ? relax*max_x : base_x);
        double lax_y = (max_y > 0 ? relax*max_y : base_y);
        double r0 = node.rVec.GetNrows() >= 1 ? std::abs(node.rVec[0]) : 0.0;
        double r1 = node.rVec.GetNrows() >= 2 ? std::abs(node.rVec[1]) : 0.0;
        if (r0 <= lax_x && r1 <= lax_y) {
          node.Accepted = true;
          accepted_new.emplace_back(node);
          prev_node = node; // advance
          current_z = prev_node.CurrentState.z;
          dbg_accepted += 1;
          if (Talk) std::cout << "[Augment] Override-accepted final plane at z=" << node.z << std::endl;
        } else {
          dbg_rejected += 1;
        }
      } else {
        dbg_rejected += 1;
      }
    }

    // Debug logging per candidate
    if (Talk) {
      // Reconstruct B sign used consistently with Predict()
      double Bsign = 0.0;
      {
        const double RegionBoundaryX = TMS_Const::TMS_Magnetic_region_2_and_3_border;
        double MagneticField = 0.0;
        if (std::abs(node.PreviousState.x) <= RegionBoundaryX) MagneticField = -1.0;
        else MagneticField = 1.0;
        Bsign = MagneticField;
        if (AugmentationMode && node.StepIndex <= AugmentFreezeSteps) Bsign = 0.0;
      }
      // Momentum before/after
      double p_prev = 1.0 / std::max(std::abs(node.PreviousState.qp), 1e-12);
      double p_curr = 1.0 / std::max(std::abs(node.CurrentState.qp), 1e-12);
      double dp = std::abs(p_curr - p_prev);
      // Residual info
      double r0 = node.rVec.GetNrows() >= 1 ? node.rVec[0] : 0.0;
      double r1 = node.rVec.GetNrows() >= 2 ? node.rVec[1] : 0.0;
      double rnorm = std::sqrt(r0*r0 + r1*r1);
      double chi2 = node.chi2;
      // Accumulate stats
      if (node.Accepted) {
        dbg_sum_chi2_acc += std::isfinite(chi2) ? chi2 : 0.0;
        dbg_sum_rnorm_acc += rnorm;
        dbg_sum_dp_acc += dp;
      } else {
        dbg_sum_chi2_rej += std::isfinite(chi2) ? chi2 : 0.0;
        dbg_sum_rnorm_rej += rnorm;
        dbg_sum_dp_rej += dp;
      }
      std::cout << "[Augment] z=" << node.z
                << ", plane=" << node.PlaneNumber
                << ", Bsign=" << Bsign
                << ", p_prev(MeV)=" << p_prev
                << ", p_curr(MeV)=" << p_curr
                << ", chi2=" << chi2
                << ", |r|=" << rnorm
                << ", accepted=" << (node.Accepted ? 1 : 0)
                << std::endl;
    }
  }
  // Change back to old value
  EnableOutlierRejection = EnableOutlierRejectionOld;
  if (Talk) {
    auto avg = [](double sum, int n){ return (n>0)? sum/static_cast<double>(n) : 0.0; };
    std::cout<<"Checked "<<pool.size()<<" hits. Found "<<nCandidates<<" candidates. Augmented "<<accepted_new.size()<<" nodes."<<std::endl;
    std::cout<<"[Augment Summary] considered="<<dbg_considered
             <<", accepted="<<dbg_accepted
             <<", rejected="<<dbg_rejected
             <<", avg_chi2_acc="<<avg(dbg_sum_chi2_acc, dbg_accepted)
             <<", avg_|r|_acc="<<avg(dbg_sum_rnorm_acc, dbg_accepted)
             <<", avg_dp_acc(MeV)="<<avg(dbg_sum_dp_acc, dbg_accepted)
             <<", avg_chi2_rej="<<avg(dbg_sum_chi2_rej, dbg_rejected)
             <<", avg_|r|_rej="<<avg(dbg_sum_rnorm_rej, dbg_rejected)
             <<", avg_dp_rej(MeV)="<<avg(dbg_sum_dp_rej, dbg_rejected)
             << std::endl;
  }

  // Merge: append accepted new nodes, then re-sort by z
  if (!accepted_new.empty()) {
    wasAugmented = true;
    nAugmentedNodes = accepted_new.size();
    // Set the smooth states. runRTSSmoother is unhappy with new nodes
    for (size_t i = 0; i < accepted_new.size(); ++i) {
      accepted_new[i].SmoothState = accepted_new[i].CurrentState;
      accepted_new[i].SmoothCovarianceMatrix = accepted_new[i].CovarianceMatrix;
    }
    KalmanNodes.insert(KalmanNodes.end(), accepted_new.begin(), accepted_new.end());
    SortNodesByZ();

    // Reindex step indices
    for (size_t i = 0; i < KalmanNodes.size(); ++i) {
      KalmanNodes[i].StepIndex = static_cast<int>(i);
    }
    InitializeMomentum(false);
    // Repair suspicious (x,y) measurements using bar geometry, then prune any
    // nodes that still look dynamically inconsistent, then refit.
    RepairBadNodeMeasurements();
    PruneInconsistentNodes();
    // Perform backward → forward → backward triple refit and smoothing
    TripleRefitSmooth();
  }

  // Restore the previous direction flag
  ForwardFitting = prevForward;
  AugmentationMode = false;
}

  // Static helper: compute a 2D measurement (x,y) for a hit's bar at z using
  // the previous state projected to that z. Mapping mirrors FillNoiseMatrix():
  //  - X bars: precise y (use bar y), x from projection
  //  - Y bars: precise x (use bar x), y from projection
  //  - U/V bars: ±tilt around Y; adjust x by tan(tilt) * (y_proj - y_mid)
void TMS_Kalman::BuildBarMeasurement(const TMS_Hit &hit,
                                     const TMS_KalmanState &state_at_prev,
                                     double meas_z,
                                     double &meas_x,
                                     double &meas_y) {
  const auto &bar = hit.GetBar();
  // Project previous state to measurement z
  double dz = meas_z - state_at_prev.z;
  double x_proj = state_at_prev.x + state_at_prev.dxdz * dz;
  double y_proj = state_at_prev.y + state_at_prev.dydz * dz;

  // Configured stereo slope (tan(tilt)) and midplane reference (mm)
  const double t = TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TiltAngle();
  const double y_mid_mm = TMS_Manager::GetInstance().Get_Geometry_YMIDDLE() * 1000.0;
  const double y_rel = y_proj - y_mid_mm;

  switch (bar.GetBarType()) {
    case TMS_Bar::kXBar:
      // X bar: precise y; use bar's y center, x from projection
      meas_y = bar.GetNotZ();
      meas_x = x_proj;
      break;
    case TMS_Bar::kYBar:
      // Y bar: precise x; use bar's x center, y from projection
      meas_x = bar.GetNotZ();
      meas_y = y_proj;
      break;
    case TMS_Bar::kUBar:
      // U bar: +tilt; use y relative to detector midplane
      meas_x = bar.GetNotZ() - t * y_rel;
      meas_y = y_proj;
      break;
    case TMS_Bar::kVBar:
      // V bar: -tilt; use y relative to detector midplane
      meas_x = bar.GetNotZ() + t * y_rel;
      meas_y = y_proj;
      break;
    default:
      // Fallback: use projected values
      meas_x = x_proj;
      meas_y = y_proj;
      break;
  }
  // Basic sanity guards
  if (!std::isfinite(meas_x)) meas_x = x_proj;
  if (!std::isfinite(meas_y)) meas_y = y_proj;
}

void TMS_Kalman::BuildMeasurementsFromBar() {
  for (size_t i = 1; i < KalmanNodes.size(); i++) {
    auto& node = KalmanNodes[i];
    if (node.Hit == NULL) continue;
    auto prev = KalmanNodes[i-1].SmoothState;
    double z = node.z;
    double x = node.x;
    double y = node.y;
    TMS_Hit& hit = (*node.Hit);
    BuildBarMeasurement(hit, prev, z, x, y);
    node.x = x;
    node.y = y;
  }
}

// Seed the starting KE for backtracking from steel thickness at a position.
// Assumes ~half steel penetration on average.
double TMS_Kalman::CalculateKEFromSteel(TVector3 position) {
  if (Talk) std::cout << "[CalculateKEFromSteel] Checking for position(" << position.X()<<","<<position.Y()<<","<<position.Z() << ")" << std::endl;
  // Sanity checks
  bool bad = false;
  if (std::abs(position.X()) > TMS_Const::TMS_Start[0]) { 
    position.SetX(0);
    bad = true;
  }
  if (position.Y() > TMS_Const::TMS_End[1] || position.Y() < TMS_Const::TMS_Start[1]) {
    position.SetY(0);
    bad = true;
  }
  if (Talk && bad) std::cout << "[CalculateKEFromSteel] Position was out of bounds. Changed to (" << position.X()<<","<<position.Y()<<","<<position.Z() << ")" << std::endl;
  // Needs to be at least the particle mass, but it'll fail if it's exactly the mass
  double en = mass + 1e-2;
  // Assume half steel thickness penetration depth on average
  double steel_depth = TMS_Geom::GetInstance().GetSteelThickness(position) * 0.5;
  TVector3 endpoint = position + TVector3(0, 0, steel_depth);
  // Get the materials between the two points
  std::vector<std::pair<TGeoMaterial*, double> > Materials = TMS_Geom::GetInstance().GetMaterials(position, endpoint);

  if (Talk) std::cout << "[CalculateKEFromSteel] Looping over " << Materials.size() << " materials" << std::endl;
  // Loop over the materials between the two projection points
  for (auto material : Materials) {
    // Read these directly from a TGeoManager
    // If the geometry is in mm (CLHEP units), then want to scale density to g/cm3 and thickness to cm
    // Otherwise, assume it's like that and then fix it with geometry scaling functions
    double density = material.first->GetDensity()/(CLHEP::g/CLHEP::cm3); 
    double thickness = material.second/10.; 
    // Potentially need to scale from g/cm3 to g/mm3, so find the scale factor and scale by 1/that^3.
    double scale_factor = TMS_Geom::GetInstance().Scale(1.0);
    density /= std::pow(scale_factor, 3);
    thickness = TMS_Geom::GetInstance().Scale(thickness);

    // Update the Bethe Bloch calculator to use this material
    try {
      Material matter(density);
      Bethe.fMaterial = matter;
      // Set the material for the multiple scattering
      MSC.fMaterial = matter;
    }
    catch (std::invalid_argument const& ex) {
      std::cout<<"[CalculateKEFromSteel] Could not make a material using density "<<density<<", is position within tms bounds?"<<std::endl;
    }

    double old_en = en;
    // Bethe.Calc_dEdx(mass) is very high, so move alittle away from mass
    // so we're not oversampling the highest point
    double b = Bethe.Calc_dEdx(std::max(mass + 5, en));
    if (!std::isfinite(b) || b < 0) {
      std::cout << "[CalculateKEFromSteel] Warning: Found Bethe(" << en << " MeV) = " << b << ", setting it to 1" << std::endl;
      b = 1;
    }
    en += b*density*thickness;
    

    // Check that new energy is valid
    if (en < 0 || std::isnan(en) || std::isinf(en)) {
      if (Talk) { 
        std::cout<<"[CalculateKEFromSteel] Got invalid energy "<<en<<" from thickness ";
        std::cout<<thickness<<" and density "<<density;
        std::cout<<" and Bethe.Calc_dEdx(en) "<<Bethe.Calc_dEdx(old_en);
        std::cout<<". Switching to old energy: "<<old_en<<std::endl;
      }
      en = old_en;
    }
  }
  double ke = en - mass;
  if (Talk) std::cout << "[CalculateKEFromSteel] Found KE " << ke << " MeV for a steel thickness of " << steel_depth << "mm" << std::endl;
  return ke;
}

// Estimate starting KE from straight-line length (start→end), with separate
// fits for thin/thick regions.
double TMS_Kalman::GetKEEstimateFromLength(double startx, double endx, double startz, double endz) {
  // if in thin and thick target there's a different relationship
  double distx = endx-startx;
  double distz = endz-startz;
  double dist = sqrt(distx*distx+distz*distz);

  /* pol0 fit to xz distance of start and death point using truth, for track ending in THICK (death > 739 cm in z)
     p0                        =      101.506   +/-   4.05823     
     p1                        =     0.132685   +/-   0.00354222  
     */

  /* pol0 fit to xz distance of start and death point using truth, for track ending in THIN (death < 739 cm in z)
     p0                        =     -1.13305   +/-   0.129217    
     p1                        =     0.234404   +/-   0.00218364  
     */

  double KEest = 0;
  // If in thick region
  if (endz > TMS_Const::TMS_Thick_Start) KEest = 101.5+0.133*dist;
  else KEest = -1.13+0.234*dist;

  return KEest;
}

// Run a single-direction Kalman filter across nodes in the current run order
void TMS_Kalman::RunKalman() {

  // Ensure all measurements are inside bounds
  //RepairBadNodeMeasurements(false);
 
  int nCand = KalmanNodes.size();
  if (nCand == 0) return;
  // Assume nodes already sorted for the chosen run direction
  ConsecutiveOutliers = 0;
  ZLayerAccepted.clear();
  KalmanNodes[0].SetRecoXY(KalmanNodes[0].PreviousState);
  // Ensure the very first node has an inflated initial covariance.
  // This keeps the filter humble about the starting state.
  if (KalmanNodes[0].CovarianceMatrix.GetNrows() == 5 && KalmanNodes[0].CovarianceMatrix.GetNcols() == 5) {
    SetInflatedInitialCovariance(KalmanNodes[0].CovarianceMatrix);
  }
  KalmanNodes[0].StepIndex = 0;
  for (int i = 1; i < nCand; ++i) {
    // Carry forward prior state/covariance and predict
    Update(KalmanNodes[i-1], KalmanNodes[i]);
    KalmanNodes[i].StepIndex = i; // keep index synced for adaptive gating
    Predict(KalmanNodes[i]);
  }

  // Update the momentum, and start/end position
  UpdateParameters();
}

void TMS_Kalman::UpdateParameters() {
  if (KalmanNodes.empty()) return;

  // Do NOT reorder KalmanNodes here; later components assume run-order
  // adjacency and transfer matrices aligned to the current order. Instead,
  // scan to find the min/max z nodes and the second-closest to the max z
  // for a stable end direction estimate.
  auto pick_state = [](const TMS_KalmanNode &n) -> const TMS_KalmanState& {
    // Prefer smoothed, then current, then previous if needed
    const TMS_KalmanState *s = &n.SmoothState;
    if (!std::isfinite(s->x) || !std::isfinite(s->y) || !std::isfinite(s->z)) s = &n.CurrentState;
    if (!std::isfinite(s->x) || !std::isfinite(s->y) || !std::isfinite(s->z)) s = &n.PreviousState;
    return *s;
  };

  size_t idx_min = static_cast<size_t>(-1);
  size_t idx_max = static_cast<size_t>(-1);
  double z_min = std::numeric_limits<double>::infinity();
  double z_max = -std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < KalmanNodes.size(); ++i) {
    const auto &st = pick_state(KalmanNodes[i]);
    if (!std::isfinite(st.z)) continue;
    if (st.z < z_min) { z_min = st.z; idx_min = i; }
    if (st.z > z_max) { z_max = st.z; idx_max = i; }
  }
  if (idx_min == static_cast<size_t>(-1) || idx_max == static_cast<size_t>(-1)) return; // nothing sane

  const auto &FrontState = pick_state(KalmanNodes[idx_min]);
  const auto &EndState   = pick_state(KalmanNodes[idx_max]);

  // Momentum from front-most qp if finite; else fallback to end qp
  double qp_front = FrontState.qp;
  double qp_end   = EndState.qp;
  double qp_use = (std::isfinite(qp_front) && std::abs(qp_front) > 1e-12) ? qp_front : qp_end;
  if (!std::isfinite(qp_use) || std::abs(qp_use) < 1e-12) qp_use = 1.0/1000.0; // ~1 GeV
  SetMomentum(1.0 / qp_use);

  // Set start/end pos/dir using finite states
  SetStartPosition(FrontState.x, FrontState.y, FrontState.z);
  SetStartDirection(FrontState.dxdz, FrontState.dydz);

  SetEndPosition(EndState.x, EndState.y, EndState.z);
  SetEndDirection(EndState.dxdz, EndState.dydz);

  if (Talk) {
    std::cout << "[UpdateParameters] nNodes=" << KalmanNodes.size()
              << ", idx_min=" << idx_min << ", z_min=" << FrontState.z
              << ", idx_max=" << idx_max << ", z_max=" << EndState.z
              << ", start=(" << FrontState.x << "," << FrontState.y << "," << FrontState.z << ")"
              << ", end=(" << EndState.x << "," << EndState.y << "," << EndState.z << ")"
              << ", p(MeV)=" << (1.0 / qp_use)
              << std::endl;
  }

  if (std::isnan(momentum) || std::isinf(momentum)){
    //std::cerr << "[TMS_Kalmann.cpp] Weirdness -- Momentum from fitter isn't a sane number: " << momentum << std::endl;
  }
}

void TMS_Kalman::Update(TMS_KalmanNode &PreviousNode, TMS_KalmanNode &CurrentNode) {
  CurrentNode.PreviousState = PreviousNode.CurrentState;
  CurrentNode.CovarianceMatrix = PreviousNode.CovarianceMatrix;
}

// Recompute dz and transfer matrices after (re)sorting nodes according to
// current run order. Ensures consistent stepping and z bookkeeping for
// subsequent filtering/smoothing passes.
void TMS_Kalman::RebuildRunOrderSteps() {
  const size_t n = KalmanNodes.size();
  if (n == 0) return;
  for (size_t i = 0; i < n; ++i) {
    double zi = KalmanNodes[i].z;
    double dz = 0.0;
    if (i > 0) dz = zi - KalmanNodes[i-1].z; // signed step consistent with order
    KalmanNodes[i].dz = dz;

    // Ensure state z bookkeeping is consistent with the plane
    KalmanNodes[i].CurrentState.z  = zi;
    KalmanNodes[i].SmoothState.z   = zi;
    KalmanNodes[i].PreviousState.z = zi - dz;

    // Refresh transfer matrices for this step
    TMatrixD &F  = KalmanNodes[i].TransferMatrix;
    TMatrixD &FT = KalmanNodes[i].TransferMatrixT;
    F.Zero(); FT.Zero();
    for (int j = 0; j < KALMAN_DIM; ++j) { F(j,j) = 1.0; FT(j,j) = 1.0; }
    F(0,2) = dz; F(1,3) = dz;
    FT(2,0) = dz; FT(3,1) = dz;
  }
}

// Backward → Forward → Backward refit with RTS smoother in each pass,
// then update physics summaries.
void TMS_Kalman::TripleRefitSmooth() {
  if (Talk) std::cout<<"[TripleRefitSmooth] Running triple refit"<<std::endl;
  if (KalmanNodes.empty()) return;
  bool saved_dir = ForwardFitting;

  // 1) Backward pass
  ForwardFitting = false;
  SortNodesByRunOrder();
  RebuildRunOrderSteps();
  InitializeMomentum(true);
  if (Talk) std::cout<<"[TripleRefitSmooth] 1st backfit p / MeV = "<<(1/KalmanNodes.at(0).PreviousState.qp)<<std::endl;
  RunKalman();
  runRTSSmoother();
  PruneRejectedAndDuplicateZ();
  //BuildMeasurementsFromBar();

  // 2) Forward pass, seeded with momentum from the backward-smoothed front
  ForwardFitting = true;
  SortNodesByRunOrder();
  RebuildRunOrderSteps();
  if (!KalmanNodes.empty()) {
    // Seed qp at the first node using its smoothed qp from the previous pass
    double qp_seed = KalmanNodes[0].SmoothState.qp;
    if (!std::isfinite(qp_seed) || std::abs(qp_seed) < 1e-12) qp_seed = 1.0/1000.0; // ~1 GeV fallback
    KalmanNodes[0].PreviousState.qp = qp_seed;
    KalmanNodes[0].CurrentState.qp  = qp_seed;
    if (Talk) std::cout<<"[TripleRefitSmooth] Forward fit p / MeV = "<<(1/qp_seed)<<std::endl;
  }
  RunKalman();
  runRTSSmoother();

  // 3) Final backward pass, seeded with momentum corrected by final angle
  //    to account for increased path through steel at non-normal incidence.
  //    Approximation: scale the steel KE by 1/cos(theta_z).
  double qp_back_seed = 1.0/1000.0; // ~1 GeV fallback
  if (!KalmanNodes.empty()) {
    // After forward smoothing, the highest-z node is at back()
    const TMS_KalmanState &end_s = KalmanNodes.back().SmoothState;
    double ax = end_s.dxdz;
    double ay = end_s.dydz;
    double mag = std::sqrt(1.0 + ax*ax + ay*ay);
    double cosz = (mag > 0.0) ? (1.0/mag) : 1.0;
    double scale = 1.0 / std::max(cosz, 1e-3); // path length scale
    TVector3 endpoint(end_s.x, end_s.y, end_s.z);
    double KE = CalculateKEFromSteel(endpoint);
    double KE_corr = std::isfinite(KE) ? KE * scale : 0.0;
    double en = KE_corr + mass;
    double p2 = en*en - mass*mass;
    double p = (p2 > 0.0 && std::isfinite(p2)) ? std::sqrt(p2) : 1.0; // MeV
    qp_back_seed = 1.0 / std::max(p, 1e-3);
  }
  ForwardFitting = false;
  SortNodesByRunOrder();
  RebuildRunOrderSteps();
  if (!KalmanNodes.empty()) {
    KalmanNodes[0].PreviousState.qp = qp_back_seed;
    KalmanNodes[0].CurrentState.qp  = qp_back_seed;

    if (Talk) std::cout<<"[TripleRefitSmooth] 2nd backfit p / MeV = "<<(1/KalmanNodes.at(0).PreviousState.qp)<<std::endl;
    if (Talk) std::cout<<"[TripleRefitSmooth] 2nd backfit p / MeV back = "<<(1/KalmanNodes.back().PreviousState.qp)<<std::endl;
  }
  RunKalman();
  runRTSSmoother();

  ForwardFitting = saved_dir;

  // Refresh downstream physics
  BetheBloch();
  Runchi2();
  SignSelection();
  UpdateParameters();
  if (Talk) std::cout<<"[TripleRefitSmooth] Done running triple refit"<<std::endl;
}

// Recompute measurements for nodes whose input (x,y) is suspicious. Uses the
// previous state's projection to the node z and bar geometry to form a
// consistent (x,y) per BuildBarMeasurement logic. Intended for cases where
// 3D-tracker-provided x/y are unreliable.
void TMS_Kalman::RepairBadNodeMeasurements(bool extrapolate) {
  if (KalmanNodes.size() < 2) return;

  const auto &mgr = TMS_Manager::GetInstance();
  const double abs_max_x = mgr.Get_Reco_Kalman_Outlier_AbsoluteResidualMaxX();
  const double abs_max_y = mgr.Get_Reco_Kalman_Outlier_AbsoluteResidualMaxY();
  const double y_mid_mm  = mgr.Get_Geometry_YMIDDLE() * 1000.0;
  const double t = mgr.Get_Reco_TRACKMATCH_TiltAngle(); // tan(tilt)

  // Assume current order is fine (we already sorted earlier in the flow).

  auto in_bounds = [](double v, double lo, double hi){ return (v >= lo && v <= hi); };
  
  bool check_initial = true;
  if (check_initial) {
    bool bad = false;
    auto &node = KalmanNodes[0];
    if (!in_bounds(node.x, TMS_Const::TMS_Start[0], TMS_Const::TMS_End[0])) bad = true;
    if (!in_bounds(node.y, TMS_Const::TMS_Start[1], TMS_Const::TMS_End[1])) bad = true;
    if (bad) {
      double x_proj = node.x;
      double y_proj = node.y;
      // First try bringing into bounds
      if (node.x < TMS_Const::TMS_Start[0]) x_proj = TMS_Const::TMS_Start[0];
      if (node.x > TMS_Const::TMS_End[0]) x_proj = TMS_Const::TMS_End[0];
      if (node.y < TMS_Const::TMS_Start[1]) y_proj = TMS_Const::TMS_Start[1];
      if (node.y > TMS_Const::TMS_End[1]) y_proj = TMS_Const::TMS_End[1];
      
      // Now try to fix with bar information
      double mx = x_proj, my = y_proj;
      if (extrapolate) {
        const double y_rel = y_proj - y_mid_mm;
        switch (node.LayerOrientation) {
          case TMS_Bar::kXBar:
            my = std::isfinite(node.LayerNotZ) ? node.LayerNotZ : y_proj; // y center of bar
            mx = x_proj;         // weakly measured axis from projection
            break;
          case TMS_Bar::kYBar:
            mx = std::isfinite(node.LayerNotZ) ? node.LayerNotZ : x_proj; // x center of bar
            my = y_proj;
            break;
          case TMS_Bar::kUBar:
            mx = (std::isfinite(node.LayerNotZ) ? node.LayerNotZ : x_proj) - t * y_rel;
            my = y_proj;
            break;
          case TMS_Bar::kVBar:
            mx = (std::isfinite(node.LayerNotZ) ? node.LayerNotZ : x_proj) + t * y_rel;
            my = y_proj;
            break;
          default:
            // Fallback to projection
            mx = x_proj;
            my = y_proj;
            break;
        }
      }
      if (Talk) {
        std::cout << "[RepairMeas] Warning: plane=" << node.PlaneNumber
                  << " z=" << node.z
                  << " out of bounds(x,y)=(" << node.x << "," << node.y << ")"
                  << " updated to(x,y)=(" << mx << "," << my << ")"
                  << std::endl;
      }
      
      node.x = mx;
      node.y = my;
    }
  }

  int repaired = 0;
  for (size_t i = 1; i < KalmanNodes.size(); ++i) {
    auto &node = KalmanNodes[i];
    const auto &prev = KalmanNodes[i-1];

    // Use smoothed previous state if available, else current
    const TMS_KalmanState &ps = (std::isfinite(prev.SmoothState.x) && std::isfinite(prev.SmoothState.y)) ? prev.SmoothState : prev.CurrentState;

    double dz = node.z - ps.z;
    double x_proj = ps.x + ps.dxdz * dz;
    double y_proj = ps.y + ps.dydz * dz;
    
    if (x_proj < TMS_Const::TMS_Start[0]) x_proj = TMS_Const::TMS_Start[0];
    if (x_proj > TMS_Const::TMS_End[0]) x_proj = TMS_Const::TMS_End[0];
    if (y_proj < TMS_Const::TMS_Start[1]) y_proj = TMS_Const::TMS_Start[1];
    if (y_proj > TMS_Const::TMS_End[1]) y_proj = TMS_Const::TMS_End[1];

    // Residuals compared to current node measurement
    double dx = node.x - x_proj;
    double dy = node.y - y_proj;

    // Heuristics: flag as bad if outside detector bounds or excessive residual
    bool bad = false;
    if (!in_bounds(node.x, TMS_Const::TMS_Start[0], TMS_Const::TMS_End[0])) bad = true;
    if (!in_bounds(node.y, TMS_Const::TMS_Start[1], TMS_Const::TMS_End[1])) bad = true;

    // Orientation-specific residual checks (robust to missing metadata)
    if (node.LayerOrientation == TMS_Bar::kXBar) {
      // X bar measures y precisely; dx can be large, dy should be close to bar center
      double w = std::isfinite(node.LayerBarWidth) ? node.LayerBarWidth : 10.0;
      if (!std::isfinite(abs_max_y) || abs_max_y <= 0) {
        if (std::abs(dy) > 3.0 * w) bad = true;
      } else if (std::abs(dy) > std::max(abs_max_y, 3.0 * w)) bad = true;
    } else if (node.LayerOrientation == TMS_Bar::kYBar) {
      // Y bar measures x precisely
      double w = std::isfinite(node.LayerBarWidth) ? node.LayerBarWidth : 10.0;
      if (!std::isfinite(abs_max_x) || abs_max_x <= 0) {
        if (std::abs(dx) > 3.0 * w) bad = true;
      } else if (std::abs(dx) > std::max(abs_max_x, 3.0 * w)) bad = true;
    } else {
      // U/V bars: x is correlated with y via stereo tilt
      double w = std::isfinite(node.LayerBarWidth) ? node.LayerBarWidth : 10.0;
      if (!std::isfinite(abs_max_x) || abs_max_x <= 0) {
        if (std::abs(dx) > 3.0 * w) bad = true;
      } else if (std::abs(dx) > std::max(abs_max_x, 3.0 * w)) bad = true;
    }

    if (!bad) continue;
    
    double mx = x_proj, my = y_proj;
    if (extrapolate) {

      // Recompute measurement using bar geometry and projected seed
      const double y_rel = y_proj - y_mid_mm;
      switch (node.LayerOrientation) {
        case TMS_Bar::kXBar:
          my = std::isfinite(node.LayerNotZ) ? node.LayerNotZ : y_proj; // y center of bar
          mx = x_proj;         // weakly measured axis from projection
          break;
        case TMS_Bar::kYBar:
          mx = std::isfinite(node.LayerNotZ) ? node.LayerNotZ : x_proj; // x center of bar
          my = y_proj;
          break;
        case TMS_Bar::kUBar:
          mx = (std::isfinite(node.LayerNotZ) ? node.LayerNotZ : x_proj) - t * y_rel;
          my = y_proj;
          break;
        case TMS_Bar::kVBar:
          mx = (std::isfinite(node.LayerNotZ) ? node.LayerNotZ : x_proj) + t * y_rel;
          my = y_proj;
          break;
        default:
          // Fallback to projection
          mx = x_proj;
          my = y_proj;
          break;
      }

      if (!std::isfinite(mx)) mx = x_proj;
      if (!std::isfinite(my)) my = y_proj;
    
    }

    if (Talk) {
      std::cout << "[RepairMeas] plane=" << node.PlaneNumber
                << " z=" << node.z
                << " was(x,y)=(" << node.x << "," << node.y << ")"
                << " proj(x,y)=(" << x_proj << "," << y_proj << ")"
                << " -> new(x,y)=(" << mx << "," << my << ")"
                << std::endl;
    }

    node.x = mx;
    node.y = my;
    repaired++;
  }
  
  // Double check that everything is within bounds
  if (extrapolate) RepairBadNodeMeasurements(false);

  if (Talk) std::cout << "[RepairMeas] Repaired nodes: " << repaired << std::endl;
}

// Drop nodes that look dynamically inconsistent, to prevent residual
// oscillations from a few bad planes.
void TMS_Kalman::PruneInconsistentNodes() {
  const size_t n = KalmanNodes.size();
  if (n < 3) return;

  // Work in increasing z to make projection logic straightforward
  SortNodesByZ();

  const auto &mgr = TMS_Manager::GetInstance();
  const double abs_max_x = mgr.Get_Reco_Kalman_Outlier_AbsoluteResidualMaxX();
  const double abs_max_y = mgr.Get_Reco_Kalman_Outlier_AbsoluteResidualMaxY();

  std::vector<char> keep(n, 1);

  // First pass: flag midpoints in oscillatory triplets
  for (size_t i = 2; i < n; ++i) {
    const auto &n0 = KalmanNodes[i-2].SmoothState;
    const auto &n1 = KalmanNodes[i-1].SmoothState;
    const auto &n2 = KalmanNodes[i].SmoothState;
    if (!(std::isfinite(n0.x) && std::isfinite(n1.x) && std::isfinite(n2.x))) continue;
    double dx0 = n1.x - n0.x;
    double dx1 = n2.x - n1.x;
    double dy0 = n1.y - n0.y;
    double dy1 = n2.y - n1.y;
    // Look for alternating sign with sizeable magnitude
    if (dx0*dx1 < 0.0 && std::min(std::abs(dx0), std::abs(dx1)) > 50.0) keep[i-1] = 0; // remove middle
    if (dy0*dy1 < 0.0 && std::min(std::abs(dy0), std::abs(dy1)) > 50.0) keep[i-1] = 0;
  }

  // Second pass: large residual vs projection from previous smoothed state
  for (size_t i = 1; i < n; ++i) {
    if (!keep[i]) continue; // already flagged
    const auto &prev = KalmanNodes[i-1].SmoothState;
    const auto &cur  = KalmanNodes[i].SmoothState;
    if (!(std::isfinite(prev.x) && std::isfinite(prev.y) && std::isfinite(prev.dxdz) && std::isfinite(prev.dydz))) continue;
    double dz = KalmanNodes[i].z - prev.z;
    double x_proj = prev.x + prev.dxdz * dz;
    double y_proj = prev.y + prev.dydz * dz;
    double rx = std::abs(cur.x - x_proj);
    double ry = std::abs(cur.y - y_proj);
    double thrx = (abs_max_x > 0) ? abs_max_x : 50.0; // mm defaults
    double thry = (abs_max_y > 0) ? abs_max_y : 50.0;
    if (rx > thrx || ry > thry) keep[i] = 0;
  }

  // Build pruned list
  std::vector<TMS_KalmanNode> pruned;
  pruned.reserve(n);
  int removed = 0;
  for (size_t i = 0; i < n; ++i) {
    if (keep[i]) pruned.emplace_back(KalmanNodes[i]); else removed++;
  }
  if (removed == 0) return;

  KalmanNodes.swap(pruned);
  for (size_t i = 0; i < KalmanNodes.size(); ++i) KalmanNodes[i].StepIndex = static_cast<int>(i);
  if (Talk) std::cout << "[Prune] Removed inconsistent nodes: " << removed << std::endl;
}

// Predict the next step: propagate state/covariance from PreviousState to
// CurrentState, accounting for energy loss, multiple scattering, and magnetic
// bending. Functional logic unchanged; comments tightened for clarity.
void TMS_Kalman::Predict(TMS_KalmanNode &Node) {

  TMS_KalmanState &PreviousState = Node.PreviousState;
  TMS_KalmanState &CurrentState = Node.CurrentState;

  // Matrices to propagate the current state
  TMatrixD &Transfer = Node.TransferMatrix;
  TMatrixD &TransferT = Node.TransferMatrixT;
  // Initialise to something sane(-ish)
  if (PreviousState.dxdz ==  -999.9) // Only on initialisation?
  {
    //PreviousState.dxdz = AverageXSlope;
    if ( std::abs(AverageXSlope) > 2.0 )
    {
      //std::cerr << "[TMS_Kalman.cpp] Excessive average X slope = " << TMS_Kalman::AverageXSlope << " of track (first to last hit), setting to 0" << std::endl;
      PreviousState.dxdz = 0.0;
    } else {
      PreviousState.dxdz = AverageXSlope;
    }
  }
  if (PreviousState.dydz ==  -999.9) // Only on initialisation?
  {
    if ( std::abs(AverageYSlope) > 1.5 )
    {
      //std::cerr << "[TMS_Kalman.cpp] Excessive average Y slope = " << TMS_Kalman::AverageYSlope << " of track (first to last hit), setting to 0" << std::endl;
      PreviousState.dydz = 0.0;
    } else {
      PreviousState.dydz = AverageYSlope;
    }
  }
  

   // Modification begins here: introduce magnetic field and deflection based on regions
  // Determine magnetic field based on x-position
  double MagneticField = 0.0;
  const double RegionBoundaryX = TMS_Const::TMS_Magnetic_region_2_and_3_border; // mm
  if (std::abs(PreviousState.x) <= RegionBoundaryX) {
      MagneticField = -1.0; // Central region: field +y or -y encoded by sign
  } else if (PreviousState.x > RegionBoundaryX) {
      MagneticField = 1.0; // Right side region
  } else {
      MagneticField = 1.0; // Left side region
  }

  // Calculate Lorentz deflection (linearized change in slope)
  double p = 1.0 / PreviousState.qp;  // momentum magnitude (MeV)
  double dz = Node.dz;                // step in mm
  // During augmentation, freeze B for the first few nodes and clamp magnitude,
  // to avoid oscillations from overly strong impulses.
  double magnetic_deflection_tx = 0.0;
  {
    double Bsign = MagneticField;
    if (AugmentationMode && Node.StepIndex <= AugmentFreezeSteps) {
      Bsign = 0.0; // freeze field effect for initial augmented steps
    }
    double raw_defl = 0.303 * assumed_charge * Bsign * 1.95 * dz * 0.5 / std::max(p, 1e-9);
    // Clamp to a small change in slope to prevent ringing
    const double max_dtx = TMS_Manager::GetInstance().Get_Reco_Kalman_Augment_MaxDeflectionPerStep();
    magnetic_deflection_tx = std::clamp(raw_defl, -max_dtx, max_dtx);
  }
  // Modification ends here
  

  TVectorD PreviousVec(5);
  PreviousVec[0] = PreviousState.x;
  PreviousVec[1] = PreviousState.y;
  PreviousVec[2] = PreviousState.dxdz;
  PreviousVec[3] = PreviousState.dydz;
  PreviousVec[4] = PreviousState.qp;

  TVectorD UpdateVec = Transfer*(PreviousVec);
  // Add magnetic deflection to slope and position (x only for By field)
  UpdateVec[2] += magnetic_deflection_tx;
  UpdateVec[0] += 0.5 * magnetic_deflection_tx * dz; // small-angle approx

  if (Talk) {
    std::cout << "\nPrevious vector: " << std::endl;
    PreviousState.Print();
  }

  // ENERGY LOSS PART
  // Update the energy
  double mom = 1./PreviousState.qp;
  //if (std::isinf(mom)) mom = 4800; // set to 1 GeV
  // Initial energy before energy loss
  double en_initial = sqrt(mom*mom+mass*mass);
  // The energy we'll be changing
  double en = en_initial;
  if (en < 0 || std::isnan(en) || std::isinf(en)) {
    if (Talk) { 
      std::cout<<"Got invalid energy "<<en<<" from momentum ";
      std::cout<<mom<<" and mass "<<mass;
      std::cout<<". Switching to energy: 0.25"<<std::endl;
    }
    en = 0.25;
  }

  // Read the position between current point and extrapolated into next bar
  double xval = PreviousState.x;
  double yval = PreviousState.y;
  double zval = PreviousState.z;

  double xval2 = CurrentState.x;
  double yval2 = CurrentState.y;
  double zval2 = PreviousState.z + Transfer(0,2); // Probably a nicer way to do this (:

  //TODO: When it is blowed up, temporary
  if (std::abs(xval)>TMS_Const::TMS_End_Exact[0]) xval=0;
  if (std::abs(xval2)>TMS_Const::TMS_End_Exact[0]) xval2=0;
  if (yval>TMS_Const::TMS_End_Exact[1]||yval<TMS_Const::TMS_Start_Exact[1]) yval=0;
  if (yval2>TMS_Const::TMS_End_Exact[1] || yval2 < TMS_Const::TMS_Start_Exact[1]) yval2=0;

  // Make TVector3s of the two points
  TVector3 start(xval,yval,zval); // Start
  TVector3 stop(xval2,yval2,zval2); // Stop

  if (Talk) {
    std::cout << "Going from " << start.X() << " " << start.Y() << " " << start.Z() << std::endl;
    std::cout << "To " << stop.X() << " " << stop.Y() << " " << stop.Z() << std::endl;
  }

  // Get the materials between the two points
  std::vector<std::pair<TGeoMaterial*, double> > Materials = TMS_Geom::GetInstance().GetMaterials(start, stop);

  if (Talk) std::cout << "Looping over " << Materials.size() << " materials" << std::endl;
  double TotalPathLength = 0;
  double TotalLength = 0;

  // Loop over the materials between the two projection points
  int counter = 0;
  double total_en_var = 0;
  if (Talk) std::cout << "Energy before " << Materials.size() << " materials: " << en_initial << std::endl;
  for (auto material : Materials) {

    // Read these directly from a TGeoManager
    // If the geometry is in mm (CLHEP units), then want to scale density to g/cm3 and thickness to cm
    // Otherwise, assume it's like that and then fix it with geometry scaling functions
    double density = material.first->GetDensity()/(CLHEP::g/CLHEP::cm3); 
    double thickness = material.second/10.; 
    // Potentially need to scale from g/cm3 to g/mm3, so find the scale factor and scale by 1/that^3.
    double scale_factor = TMS_Geom::GetInstance().Scale(1.0);
    density /= std::pow(scale_factor, 3);
    thickness = TMS_Geom::GetInstance().Scale(thickness);

    TotalPathLength += density*thickness;
    TotalLength += thickness;

    // Update the Bethe Bloch calculator to use this material
    try {
      Material matter(density);
      Bethe.fMaterial = matter;
      // Set the material for the multiple scattering
      MSC.fMaterial = matter;
    }
    catch (std::invalid_argument const& ex) {
      std::cout<<"Could not make a material using density "<<density<<", is position within tms bounds?"<<std::endl;
    }

    double old_en = en;
    // Subtract off the energy loss for this material
    if (ForwardFitting) en -= Bethe.Calc_dEdx(en)*density*thickness;
    else                en += Bethe.Calc_dEdx(en)*density*thickness;


    // Variance assuming Gaussian straggling
    double en_var = Bethe.Calc_dEdx_Straggling(en)*density*thickness;
    total_en_var += en_var*en_var;

    // Calculate this before or after the energy subtraction/addition?
    MSC.Calc_MS(en, thickness*density);
    
    // Check that new energy is valid
    if (en < 0 || std::isnan(en) || std::isinf(en)) {
      if (Talk) { 
        std::cout<<"Got invalid energy "<<en<<" from thickness ";
        std::cout<<thickness<<" and density "<<density;
        std::cout<<". Switching to old energy: "<<old_en<<std::endl;
      }
      en = old_en;
    }

    counter++;
  }
  if (Talk) {
    std::cout << " - energy after " << counter << ": " << en << std::endl;
    std::cout << "Total path length (g/cm2): " << TotalPathLength << std::endl;
    std::cout << "Total length (cm): " << TotalLength << std::endl;
  }

  // Updated momentum^2
  double p_2_up = en*en-BetheBloch_Utils::Mm*BetheBloch_Utils::Mm;
  double p_up;
  if (p_2_up > 0) p_up = sqrt(p_2_up);
  else {
    //std::cerr << "[TMS_Kalman.cpp] negative momentum squared, setting momentum to 1 MeV" << std::endl;
    //p_up = 1;
    p_up = sqrt(en*en);
  }

  // Update the state's q/p with sane clamping and optional augmentation smoothing
  if (Talk) std::cout << "setting current p = " << p_up << std::endl;
  // Guard against zero/NaN momentum; clamp to [0.2, 20] GeV
  double p_min = 200.0;   // MeV
  double p_max = 20000.0; // MeV
  if (!std::isfinite(p_up) || p_up <= 0.0) p_up = 1000.0; // 1 GeV fallback
  // In augmentation, optionally freeze momentum updates for first few nodes
  if (AugmentationMode && Node.StepIndex <= AugmentFreezeSteps) {
    p_up = 1.0 / std::max(std::abs(PreviousState.qp), 1e-9);
  }
  p_up = std::clamp(p_up, p_min, p_max);
  double qp_up = 1.0 / p_up;
  if (!std::isfinite(qp_up)) qp_up = 1.0/1000.0; // fallback ~1 GeV
  // During augmentation, blend qp to avoid oscillations from step-to-step
  if (AugmentationMode) {
    double beta = 0.2; // small update fraction
    qp_up = (1.0 - beta) * PreviousState.qp + beta * qp_up;
  }
  CurrentState.qp = qp_up;


  double p_var = (2*en/p_up)*(2*en/p_up) * total_en_var;
  // Guard qp variance against zero/NaN
  double p4 = p_up*p_up*p_up*p_up;
  double qp_var = (p4 > 0.0 && std::isfinite(p4)) ? (p_var / p4) : 0.0;
  if (!std::isfinite(qp_var)) qp_var = 0.0;


  // Covariance matrices
  TMatrixD &NoiseMatrix = Node.NoiseMatrix;
  TMatrixD &CovarianceMatrix = Node.CovarianceMatrix;
  TMatrixD &UpdatedCovarianceMatrix = Node.UpdatedCovarianceMatrix;
  TMatrixD &EstimatedCovarianceMatrix = Node.EstimatedCovarianceMatrix;//For smooth
  
  if (Talk) std::cout << "Initial cov\n" << std::flush;
  if (Talk) CovarianceMatrix.Print();

  double sigma = MSC.Calc_MS_Sigma();
  
  if (sigma < 0 || std::isnan(sigma) || std::isinf(sigma)) {
    if (Talk) { 
      std::cout<<"Got invalid sigma "<<sigma<<", setting it to 0.0001"<<std::endl;
    }
    sigma = 0.0001;
  }

  if (TotalPathLength >= 1.0) { // If path is 0 we're in the first Node, set initial cov
    Node.FillUpdatedCovarianceMatrix(TotalPathLength, UpdateVec[2], UpdateVec[3], qp_var, sigma, ForwardFitting); // Use direction for proper sign
  } else { // Initialise cov to something 'sane'-ish
    UpdatedCovarianceMatrix.Zero(); // zero this out to be sure
    // Inflate the initial covariance strongly to avoid overconfidence
    SetInflatedInitialCovariance(CovarianceMatrix);
    if (Talk) std::cout << "Initialising covariance!" << std::endl;
  }

  Node.FillNoiseMatrix(); // Build measurement noise matrix

  CovarianceMatrix = Transfer*CovarianceMatrix*TransferT;

  CovarianceMatrix += UpdatedCovarianceMatrix;

  EstimatedCovarianceMatrix = CovarianceMatrix;

  // Proper Kalman update using H and R
  TMatrixD H(2,5); H.Zero(); H(0,0) = 1.0; H(1,1) = 1.0;
  TMatrixD HT = TMatrixD(TMatrixD::kTransposed, H);
  TMatrixD R(2,2); R.Zero();
  // Use the top-left 2x2 of the measurement noise matrix
  R(0,0) = NoiseMatrix(0,0);
  R(0,1) = NoiseMatrix(0,1);
  R(1,0) = NoiseMatrix(1,0);
  R(1,1) = NoiseMatrix(1,1);

  TMatrixD S = H * CovarianceMatrix * HT; S += R; // innovation covariance (2x2)
  // Symmetrize S to reduce numerical asymmetry
  {
    TMatrixD ST(TMatrixD::kTransposed, S);
    S += ST; S *= 0.5;
  }
  // Guard inversion
  bool invertible = true;
  double detS = S.Determinant();
  if (!std::isfinite(detS) || std::abs(detS) < 1e-18) invertible = false;
  TMatrixD S_inv = S;
  if (invertible) {
    S_inv.Invert();
    // Basic sanity: check result isn't exploding
    if (!std::isfinite(S_inv(0,0)) || !std::isfinite(S_inv(1,1))) invertible = false;
  }
  // Only compute gain with a valid inverse; otherwise keep zero
  TMatrixD GainMatrix(5,2); GainMatrix.Zero();
  if (invertible) {
    GainMatrix = CovarianceMatrix * HT * S_inv;
  }

  if (Talk) std::cout << "Final cov\n" << std::flush;
  if (Talk) CovarianceMatrix.Print();

  TVectorD FilteredVec = TVectorD(5);
  TVectorD Measurement = TVectorD(2);

  // Measurement is the observed hit position of this node
  Measurement[0] = Node.x;
  Measurement[1] = Node.y;


  if (Talk) std::cout << "Gain" << std::flush;
  if (Talk) GainMatrix.Print();

  // Innovation and chi2
  TVectorD innovation = Measurement - H * UpdateVec; // size 2
  // Compute per-node chi2 ahead of update (only if inversion succeeded)
  double node_chi2 = std::numeric_limits<double>::infinity();
  if (invertible) {
    TVectorD tmp_for_chi2 = S_inv * innovation;
    node_chi2 = innovation * tmp_for_chi2; // r^T S^{-1} r (2 dof)
  }

  // Outlier rejection: skip measurement update if chi2 too large
  bool use_measurement = true;
  double effective_chi2_thresh = Chi2RejectThreshold;
  if (EnableOutlierRejection) {
    effective_chi2_thresh = AdaptiveChi2Threshold(S, Node.StepIndex, Chi2RejectThreshold);
    if (!invertible || node_chi2 > effective_chi2_thresh || std::isnan(node_chi2) || std::isinf(node_chi2)) {
      use_measurement = false;
    }
    // Absolute residual guardrail (per-axis, in mm)
    if (use_measurement && TMS_Manager::GetInstance().Get_Reco_Kalman_Outlier_UseAbsoluteResidual() \
        && Node.StepIndex > TMS_Manager::GetInstance().Get_Reco_Kalman_Outlier_EarlyStepsCount()) {
      const auto &mgr = TMS_Manager::GetInstance();
      double max_x = mgr.Get_Reco_Kalman_Outlier_AbsoluteResidualMaxX();
      double max_y = mgr.Get_Reco_Kalman_Outlier_AbsoluteResidualMaxY();
      if (std::abs(innovation[0]) > max_x || std::abs(innovation[1]) > max_y) use_measurement = false;
    }
  }
  // Duplicate-at-same-z guard: if a measurement at this z was already
  // accepted, reject others (always enforced, independent of outlier switch)
  {
    long long zkey = QuantizeZKey(Node.z);
    if (use_measurement && ZLayerAccepted.count(zkey) && ZLayerAccepted[zkey]) {
      use_measurement = false;
    }
  }

  if (use_measurement) {
    FilteredVec = UpdateVec + GainMatrix * innovation;
    ConsecutiveOutliers = 0; // reset streak when we accept a hit
    // Mark this z-layer as having an accepted measurement
    long long zkey = QuantizeZKey(Node.z);
    ZLayerAccepted[zkey] = true;
  } else {
    // prediction only
    FilteredVec = UpdateVec;
    // Nudge towards raw measurement in position space (configurable alpha)
    const auto &mgr = TMS_Manager::GetInstance();
    double alpha = mgr.Get_Reco_Kalman_Outlier_NudgeAlpha();
    alpha = std::clamp(alpha, 0.0, 1.0);
    FilteredVec[0] += alpha * (Measurement[0] - UpdateVec[0]);
    FilteredVec[1] += alpha * (Measurement[1] - UpdateVec[1]);
    if (mgr.Get_Reco_Kalman_Outlier_NudgeSlopes()) {
      double as = mgr.Get_Reco_Kalman_Outlier_NudgeSlopesAlpha();
      as = std::clamp(as, 0.0, 1.0);
      double dz_step = std::max(1e-6, Node.dz);
      // Approximate slope residuals from position innovation over dz
      FilteredVec[2] += as * (innovation[0] / dz_step);
      FilteredVec[3] += as * (innovation[1] / dz_step);
    }
    ConsecutiveOutliers++;
    if (ConsecutiveOutliers > OutlierResetThreshold) {
      // Too many outliers in a row: reset covariance to a humble prior
      SetInflatedInitialCovariance(CovarianceMatrix);
      EstimatedCovarianceMatrix = CovarianceMatrix;
      if (Talk) std::cout << "[Kalman] Resetting covariance after " << ConsecutiveOutliers << " consecutive outliers" << std::endl;
      ConsecutiveOutliers = 0; // avoid repeated resets
    }
  }

  TMatrixD IdentityMatrix(5,5); IdentityMatrix.UnitMatrix();
  TMatrixD KH = GainMatrix * H;
  // Joseph form to preserve symmetry/positivity: P = (I-KH) P (I-KH)^T + K R K^T
  TMatrixD IminusKH = IdentityMatrix; IminusKH -= KH;
  if (use_measurement) {
    CovarianceMatrix = IminusKH * CovarianceMatrix * TMatrixD(TMatrixD::kTransposed, IminusKH);
    TMatrixD KRKt = GainMatrix * R * TMatrixD(TMatrixD::kTransposed, GainMatrix);
    CovarianceMatrix += KRKt;
  }
  // Symmetrize P
  {
    TMatrixD PT(TMatrixD::kTransposed, CovarianceMatrix);
    CovarianceMatrix += PT; CovarianceMatrix *= 0.5;
  }

  // Store innovation-based chi2 per node (2D residual normalized by innovation covariance)
  Node.rVec[0] = innovation[0];
  Node.rVec[1] = innovation[1];
  // Store a symmetrized inverse for stability when available
  if (invertible) {
    TMatrixD S_inv_sym = S_inv;
    TMatrixD S_inv_T(TMatrixD::kTransposed, S_inv);
    S_inv_sym += S_inv_T; S_inv_sym *= 0.5;
    Node.RMatrix = S_inv_sym;
  } else {
    Node.RMatrix.Zero();
  }
  Node.chi2 = node_chi2;
  Node.Accepted = use_measurement;

  CurrentState.x    = FilteredVec[0];
  CurrentState.y    = FilteredVec[1];
  CurrentState.dxdz = FilteredVec[2];
  CurrentState.dydz = FilteredVec[3];

  if ( (CurrentState.x < TMS_Const::TMS_Start[0]) || (CurrentState.x > TMS_Const::TMS_End[0]) ) // point outside x region
  {
    //std::cerr << "[TMS_Kalman.cpp] x value outside TMS: " << CurrentState.x << "\tTMS: [" << TMS_Const::TMS_Start[0] << ", "<< TMS_Const::TMS_End[0] << "]" << std::endl;
    if (Talk)
    {
      Node.PrintTrueReco();
      PreviousState.Print();
      CurrentState.Print();
    }
  }
  if ( (CurrentState.y < TMS_Const::TMS_Start[1]) || (CurrentState.y > TMS_Const::TMS_End[1]) ) // point outside y region
  {
    //std::cerr << "[TMS_Kalman.cpp] y value outside TMS: " << CurrentState.y << "\tTMS: [" << TMS_Const::TMS_Start[1] << ", "<< TMS_Const::TMS_End[1] << "]" << std::endl;
    if (Talk)
    {
      Node.PrintTrueReco();
      PreviousState.Print();
      CurrentState.Print();
    }
  }

  Node.SetRecoXY(CurrentState);
  if (Talk) Node.PrintTrueReco();
  if (Talk) CurrentState.Print();
}

// Generate a toy noise vector from the node's noise covariance (used for tests)
TVectorD TMS_Kalman::GetNoiseVector(TMS_KalmanNode Node) {
  TVectorD rand_vec;
  rand_vec.ResizeTo(KALMAN_DIM);
  for(int i = 0; i < KALMAN_DIM; i++)
    rand_vec[i] = RNG.Gaus();


  TMatrixD &cov = Node.NoiseMatrix;
  TVectorD toy;
  toy.ResizeTo(5);
  toy.Zero();

  toy = cov*rand_vec;
  if (Talk) rand_vec.Print();

  return toy;
}

// Remove nodes rejected as outliers and any duplicates at the same z-layer
void TMS_Kalman::PruneRejectedAndDuplicateZ() {
  if (KalmanNodes.empty()) return;

  std::vector<TMS_KalmanNode> pruned;
  pruned.reserve(KalmanNodes.size());
  std::unordered_map<long long, bool> keptZ;
  for (auto &node : KalmanNodes) {
    if (!node.Accepted) continue; // drop outlier-rejected nodes
    long long key = QuantizeZKey(node.z);
    if (keptZ.count(key)) continue; // drop duplicates at same z
    pruned.emplace_back(node);
    keptZ[key] = true;
  }

  KalmanNodes.swap(pruned);

  // Reindex and refresh start/end for downstream users
  for (size_t i = 0; i < KalmanNodes.size(); ++i) {
    KalmanNodes[i].StepIndex = static_cast<int>(i);
  }

  if (!KalmanNodes.empty()) {
    // Momentum already set; keep it. Refresh geometry endpoints to match pruned nodes
    SetStartPosition(KalmanNodes.back().CurrentState.x, KalmanNodes.back().CurrentState.y, KalmanNodes.back().CurrentState.z);
    SetStartDirection(KalmanNodes.back().CurrentState.dxdz, KalmanNodes.back().CurrentState.dydz);
    if (KalmanNodes.size() > 1) {
      SetEndPosition(KalmanNodes.front().PreviousState.x, KalmanNodes.front().PreviousState.y, KalmanNodes.front().PreviousState.z);
      if (KalmanNodes.size() > 1) SetEndDirection(KalmanNodes[1].CurrentState.dxdz, KalmanNodes[1].CurrentState.dydz);
    } else {
      SetEndPosition(KalmanNodes.front().CurrentState.x, KalmanNodes.front().CurrentState.y, KalmanNodes.front().CurrentState.z);
      SetEndDirection(KalmanNodes.front().CurrentState.dxdz, KalmanNodes.front().CurrentState.dydz);
    }
  }
}


void TMS_Kalman::SetStartDirection(double ax, double ay)
{
  // Defining a vector by two slopes is a bit odd;
  // it's not possible to define vectors in the x/y plane for example.
  // Consider the dector (dx/dz, dy/dz, dz/dz), the z component is 1 by construction,
  // and thus we normalise this
  double mag  = sqrt(ax*ax + ay*ay + 1);

  StartDirection[0]=ax/mag;
  StartDirection[1]=ay/mag;
  StartDirection[2]= 1/mag;
}

void TMS_Kalman::SetEndDirection(double ax, double ay)
{
  // Defining a vector by two slopes is a bit odd;
  // it's not possible to define vectors in the x/y plane for example.
  // Consider the dector (dx/dz, dy/dz, dz/dz), the z component is 1 by construction,
  // and thus we normalise this
  double mag  = sqrt(ax*ax + ay*ay + 1);

  EndDirection[0]=ax/mag;
  EndDirection[1]=ay/mag;
  EndDirection[2]= 1/mag;
}

// RTS smoother (Rauch–Tung–Striebel). Reference:
// https://jwmi.github.io/ASM/6-KalmanFilter.pdf
void TMS_Kalman::runRTSSmoother() {
  int nCand = KalmanNodes.size();
  if (nCand == 0) return;
  // Initialize last node
  KalmanNodes[nCand-1].SmoothState = KalmanNodes[nCand-1].CurrentState;
  KalmanNodes[nCand-1].SmoothCovarianceMatrix = KalmanNodes[nCand-1].CovarianceMatrix;
  // Backward pass
  for (int t = nCand-2; t >= 0; --t) {
      // Filtered covariance (at t) and predicted covariance (t+1|t)
      TMatrixD V_t = KalmanNodes[t].CovarianceMatrix;
      TMatrixD P_t1_pred = KalmanNodes[t+1].EstimatedCovarianceMatrix;
      TMatrixD SmoothCov_t1 = KalmanNodes[t+1].SmoothCovarianceMatrix;
      TMatrixD F = KalmanNodes[t+1].TransferMatrix; // transition from t -> t+1
      TMatrixD FT(F); FT.T();

      if (P_t1_pred.Determinant() == 0) {
          continue;
      }

      // Smoothing gain
      TMatrixD P_t1_pred_inv = P_t1_pred; P_t1_pred_inv.Invert();
      TMatrixD C_t = V_t * FT * P_t1_pred_inv;

      // States
      TVectorD u_t(5);
      u_t[0] = KalmanNodes[t].CurrentState.x;
      u_t[1] = KalmanNodes[t].CurrentState.y;
      u_t[2] = KalmanNodes[t].CurrentState.dxdz;
      u_t[3] = KalmanNodes[t].CurrentState.dydz;
      u_t[4] = KalmanNodes[t].CurrentState.qp;

      TVectorD u_t1_smooth(5);
      u_t1_smooth[0] = KalmanNodes[t + 1].SmoothState.x;
      u_t1_smooth[1] = KalmanNodes[t + 1].SmoothState.y;
      u_t1_smooth[2] = KalmanNodes[t + 1].SmoothState.dxdz;
      u_t1_smooth[3] = KalmanNodes[t + 1].SmoothState.dydz;
      u_t1_smooth[4] = KalmanNodes[t + 1].SmoothState.qp;

      // Predicted state at t+1|t
      TVectorD u_t1_pred = F * u_t;

      TVectorD smoothed_state = u_t + C_t * (u_t1_smooth - u_t1_pred);
      TMatrixD Cov_smooth = V_t + C_t * (SmoothCov_t1 - P_t1_pred) * TMatrixD(TMatrixD::kTransposed, C_t);

      KalmanNodes[t].SmoothCovarianceMatrix = Cov_smooth;
      KalmanNodes[t].SmoothState.x    = smoothed_state[0];
      KalmanNodes[t].SmoothState.y    = smoothed_state[1];
      KalmanNodes[t].SmoothState.dxdz = smoothed_state[2];
      KalmanNodes[t].SmoothState.dydz = smoothed_state[3];
      KalmanNodes[t].SmoothState.qp   = KalmanNodes[t].CurrentState.qp;
      KalmanNodes[t].SetRecoXY(KalmanNodes[t].SmoothState);
  }
}

void TMS_Kalman::BetheBloch() {
    int nCand = KalmanNodes.size();
    for (int i = 1; i < nCand; ++i) {
        TMS_KalmanState &PreviousState = KalmanNodes[i-1].SmoothState;
        TMS_KalmanState &CurrentState = KalmanNodes[i].SmoothState;


        TVectorD PreviousVec(5);
        PreviousVec[0] = PreviousState.x;
        PreviousVec[1] = PreviousState.y;
        PreviousVec[2] = PreviousState.dxdz;
        PreviousVec[3] = PreviousState.dydz;
        PreviousVec[4] = PreviousState.qp;

        TVectorD CurrentVec(5);
        CurrentVec[0] = CurrentState.x;
        CurrentVec[1] = CurrentState.y;
        CurrentVec[2] = CurrentState.dxdz;
        CurrentVec[3] = CurrentState.dydz;
        CurrentVec[4] = CurrentState.qp;


        double mom = 1./PreviousState.qp;
        double en_initial = sqrt(mom*mom+mass*mass);
        double en = en_initial;

        //
        // Read the position between current point and extrapolated into next bar
        double xval = PreviousState.x;
        double yval = PreviousState.y;
        double zval = PreviousState.z;

        double xval2 = CurrentState.x;
        double yval2 = CurrentState.y;
        double zval2 = CurrentState.z; // Probably a nicer way to do this (:

        //TODO: When it is blowed up, temporary
        if (std::abs(xval)>TMS_Const::TMS_End_Exact[0]) xval=0;
        if (std::abs(xval2)>TMS_Const::TMS_End_Exact[0]) xval2=0;
        if (yval>TMS_Const::TMS_End_Exact[1]||yval<TMS_Const::TMS_Start_Exact[1]) yval=0;
        if (yval2>TMS_Const::TMS_End_Exact[1] || yval2 < TMS_Const::TMS_Start_Exact[1]) yval2=0;

        TVector3 start(xval,yval,zval); // Start
        TVector3 stop(xval2,yval2,zval2); // Stop

        // Get the materials between the two points
        std::vector<std::pair<TGeoMaterial*, double> > Materials = TMS_Geom::GetInstance().GetMaterials(start, stop);

        double TotalPathLength = 0;
        double TotalLength = 0;

        // Loop over the materials between the two projection points
        int counter = 0;
        double total_en_var = 0;
        for (auto material : Materials) {

            // Read these directly from a TGeoManager
            // If the geometry is in mm (CLHEP units), then want to scale density to g/cm3 and thickness to cm
            // Otherwise, assume it's like that and then fix it with geometry scaling functions
            double density = material.first->GetDensity()/(CLHEP::g/CLHEP::cm3); 
            double thickness = material.second/10.; 
            // Potentially need to scale from g/cm3 to g/mm3, so find the scale factor and scale by 1/that^3.
            double scale_factor = TMS_Geom::GetInstance().Scale(1.0);
            density /= std::pow(scale_factor, 3);
            thickness = TMS_Geom::GetInstance().Scale(thickness);

            TotalPathLength += density*thickness;
            TotalLength += thickness;

            // Update the Bethe Bloch calculator to use this material
            try {
                Material matter(density);
                Bethe.fMaterial = matter;
                // Set the material for the multiple scattering
                MSC.fMaterial = matter;
            }
            catch (std::invalid_argument const& ex) {
                std::cout<<"Could not make a material using density "<<density<<", is position within tms bounds?"<<std::endl;
            }

            if (ForwardFitting) en -= Bethe.Calc_dEdx(en)*density*thickness;
            else                en += Bethe.Calc_dEdx(en)*density*thickness;

            // Variance assuming Gaussian straggling
            double en_var = Bethe.Calc_dEdx_Straggling(en)*density*thickness;
            total_en_var += en_var*en_var;

            // Calculate this before or after the energy subtraction/addition?
            MSC.Calc_MS(en, thickness*density);

            counter++;
        }

        // Updated momentum^2
        double p_2_up = en*en-BetheBloch_Utils::Mm*BetheBloch_Utils::Mm;
        double p_up;
        if (p_2_up > 0) p_up = sqrt(p_2_up);
        else {
            //std::cerr << "[TMS_Kalman.cpp] negative momentum squared, setting momentum to 1 MeV" << std::endl;
            p_up = sqrt(en*en);
        }


        CurrentState.qp = 1./p_up;

    }
    SetMomentum(1./KalmanNodes.back().SmoothState.qp);
}
void TMS_Kalman::Runchi2() {
    int nCand = KalmanNodes.size();
    for (int i = 0; i < nCand; ++i) {
        // If innovation and S^{-1} were stored during filtering, compute normalized chi2
        if (KalmanNodes[i].RMatrix.GetNrows() == 2 && KalmanNodes[i].RMatrix.GetNcols() == 2) {
            TVectorD tmp = KalmanNodes[i].RMatrix * KalmanNodes[i].rVec;
            KalmanNodes[i].chi2 = KalmanNodes[i].rVec * tmp; // r^T S^{-1} r
        }
    }
}
void TMS_Kalman::SignSelection() {
    int nCand = KalmanNodes.size();
    double sum_signals=0.0;
    for (int i = 1; i < nCand; ++i) {
        TMS_KalmanState &PreviousState = KalmanNodes[i-1].SmoothState;
        TMS_KalmanState &CurrentState = KalmanNodes[i].SmoothState;


        TVector3 B;
        const double RegionBoundaryX = TMS_Const::TMS_Magnetic_region_2_and_3_border; //1860; // mm
        if (std::abs(0.5 * (PreviousState.x + CurrentState.x)) <= RegionBoundaryX)
            B = TVector3(0.0, 1.0, 0.0); // Central region
        else
            B = TVector3(0.0, -1.0, 0.0); // Central region

        TVector3 p1 = TVector3(PreviousState.dxdz,PreviousState.dydz,1);
        TVector3 p2 = TVector3(CurrentState.dxdz,CurrentState.dydz,1);
        p1=p1.Unit();
        p2=p2.Unit();
        p1=p1*(1/PreviousState.qp);
        p2=p2*(1/CurrentState.qp);
        TVector3 dp = p2 - p1;//It supposed to be p1-p2 but backword and keep the B-field, So p2-p1 is right, 
        TVector3 p_avg = 0.5 * (p1 + p2);
        TVector3 v_cross_B = p_avg.Cross(B);

        double signal = v_cross_B.Dot(dp);
        //std::cout << "p_avg: (" << p_avg.X() << ", " << p_avg.Y() << ", " << p_avg.Z() << ")"
        //  << ", v_cross_B: (" << v_cross_B.X() << ", " << v_cross_B.Y() << ", " << v_cross_B.Z() << ")"
        //  << ", dp: (" << dp.X() << ", " << dp.Y() << ", " << dp.Z() << ")"
        //  << std::endl;
        sum_signals += signal;
    }
    SetCharge_curvature(sum_signals);
}
