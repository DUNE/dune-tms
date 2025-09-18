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

// Helper to inflate the initial state covariance so the filter
// is humble about the initial state estimate.
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

// Take a collection of hits, return collection of Kalman nodes
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
  //KalmanNodes.reserve(nCand);
  for (int i = 0; i < nCand; ++i) {

    TMS_Hit hit = Candidates[i];
    double x_true = hit.GetTrueHit().GetX();
    double y_true = hit.GetTrueHit().GetY();
    double x = hit.GetRecoX();
    double y = hit.GetRecoY();
    double z = hit.GetZ();

    TVector3 vecc = TVector3(x,y,z);

    if (! (TMS_Geom::GetInstance().IsInsideBox(vecc, TMS_Const::TMS_Start_Exact, TMS_Const::TMS_End_Exact)))
    {
      //std::cerr << "[TMS_Kalman.cpp] Hit " << i << "/" << nCand << " position not within TMS before Kalman filter, (x,y,z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
      //std::cerr << "[TMS_Kalman.cpp] Were reco x and y values set before running Kalman?" << std::endl;
      //throw; // yeet it
    }

    int j;
    for (j=i; j<nCand; j++)
      if (std::abs(Candidates[j].GetZ() - z) > 1E-3)
      {
        break;
      }

    // Use the next different-z hit as the future reference, if it exists
    int next_idx = (j < nCand) ? j : i;
    double future_z = (next_idx == i) ? z : Candidates[next_idx].GetZ();
    double DeltaZ = future_z - z;

    double future_x = (next_idx == i) ? x : Candidates[next_idx].GetRecoX();
    double DeltaX = future_x - x;

    double future_y = (next_idx == i) ? y : Candidates[next_idx].GetRecoY();
    double DeltaY = future_y - y;


    // This also initialises the state vectors in each of the nodes
    if (std::abs(DeltaZ) > 1E-3) // TODO: Only add one hit per z for now, noise breaks
    {
      // TODO: Combine multiple hits into a single 'node' <-> 'measurement'
      TMS_KalmanNode Node(x, y, z, DeltaZ, DeltaX/DeltaZ, DeltaY/DeltaZ);
      Node.SetTrueXY(x_true, y_true); // Add truth to enable reco to truth comparison
      Node.LayerOrientation = hit.GetBar().GetBarType();
      Node.LayerBarWidth    = hit.GetBar().GetBarWidth();
      Node.LayerBarLength   = hit.GetBar().GetBarLength();
      Node.PlaneNumber      = hit.GetPlaneNumber();

      KalmanNodes.emplace_back(std::move(Node));
      KalmanNodes.back().StepIndex = KalmanNodes.size() - 1;
//    } else { // TODO: Handle layers with more than one hit, waiting on Asa to confirm potential structures
//      //std::cout << "more than one hit per layer? Kalman unhappy " << i << "\t " << j-i << std::endl;
    }
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
  RunKalman();
  runRTSSmoother();
  BetheBloch();
  Runchi2();
  SignSelection();
  // Prune out nodes rejected as outliers and duplicates at the same z-layer
  PruneRejectedAndDuplicateZ();
}

// Minimal second pass: start from the smoothed state near the front (lowest z),
// propagate forward in z through an external pool of candidates, and accept any
// non-outlier hits. Merge accepted nodes, re-sort by z, and re-smooth.
void TMS_Kalman::AugmentWithCandidates(const std::vector<TMS_Hit> &candidate_pool) {
  if (KalmanNodes.empty() || candidate_pool.empty()) return;

  // We will run a forward pass, so set this flag for physics (energy loss sign)
  bool prevForward = ForwardFitting;
  ForwardFitting = true;

  // Sort a local copy of candidates by z
  std::vector<TMS_Hit> pool(candidate_pool.begin(), candidate_pool.end());
  std::sort(pool.begin(), pool.end(), TMS_Hit::SortByZ);

  // Ensure nodes are sorted by z ascending so back() is highest z
  SortNodesByZ();
  // Seed from the smoothed state at the highest-z node (back of vector)
  TMS_KalmanNode seed = KalmanNodes.back();
  TMS_KalmanNode prev_node(
      seed.SmoothState.x,
      seed.SmoothState.y,
      seed.SmoothState.z,
      0.0,
      seed.SmoothState.dxdz,
      seed.SmoothState.dydz);
  prev_node.SetTrueXY(seed.TrueX, seed.TrueY);
  prev_node.LayerOrientation = seed.LayerOrientation;
  prev_node.LayerBarWidth = seed.LayerBarWidth;
  prev_node.LayerBarLength = seed.LayerBarLength;
  prev_node.PlaneNumber = seed.PlaneNumber;
  prev_node.StepIndex = seed.StepIndex;
  // Carry covariance from the smoother as our starting uncertainty
  prev_node.CovarianceMatrix = KalmanNodes.back().SmoothCovarianceMatrix;
  prev_node.CurrentState = seed.SmoothState;
  prev_node.PreviousState = seed.SmoothState;
  // Boost initial momentum for forward pass (reduce qp) on BOTH states
  {
    double qp_seed = seed.SmoothState.qp;
    if (!std::isfinite(qp_seed) || std::abs(qp_seed) < 1e-12) {
      // Fallback to a conservative finite qp (~1/1 GeV)
      qp_seed = 1.0/1000.0;
    }
    double scale = 0.25; // 4x momentum to avoid starvation
    prev_node.PreviousState.qp = qp_seed * scale;
    prev_node.CurrentState.qp  = qp_seed * scale;
  }
  prev_node.StepIndex = 0;
  // Be conservative about momentum during augmentation
  if (prev_node.CovarianceMatrix.GetNrows() == 5 && prev_node.CovarianceMatrix.GetNcols() == 5) {
    prev_node.CovarianceMatrix(4,4) *= 25.0; // inflate qp variance modestly
  }

  // Reset outlier streak and per-z acceptance map for this pass
  ConsecutiveOutliers = 0;
  ZLayerAccepted.clear();

  // Track current z to enforce monotonic forward propagation beyond the current end
  double current_z = prev_node.CurrentState.z;
  int step_index = 0;

  std::vector<TMS_KalmanNode> accepted_new;
  accepted_new.reserve(pool.size());

  for (const auto &hit : pool) {
    double z = hit.GetZ();
    double dz = z - prev_node.CurrentState.z;
    if (z <= current_z + 1e-6) continue; // only forward propagation

    // Form a node at the candidate's z with a measurement (x,y) derived
    // from bar geometry and the current track projection.
    double x = 0.0, y = 0.0;
    BuildBarMeasurement(hit, prev_node.CurrentState, z, x, y);

    TMS_KalmanNode node(x, y, z, dz, prev_node.CurrentState.dxdz, prev_node.CurrentState.dydz);
    double x_true = hit.GetTrueHit().GetX();
    double y_true = hit.GetTrueHit().GetY();
    node.SetTrueXY(x_true, y_true);
    node.LayerOrientation = hit.GetBar().GetBarType();
    node.LayerBarWidth    = hit.GetBar().GetBarWidth();
    node.LayerBarLength   = hit.GetBar().GetBarLength();
    node.PlaneNumber      = hit.GetPlaneNumber();
    node.StepIndex        = ++step_index;
    
    // Skip if over a certain number of planes
    int maximum_n_planes = 5;
    if (maximum_n_planes > 0 && node.PlaneNumber > prev_node.PlaneNumber + maximum_n_planes) continue;

    // Propagate from prev_node state to this node and perform update/gating
    Update(prev_node, node);
    Predict(node);

    // If accepted by gating, keep it and advance the prev_node state
    if (node.Accepted) {
      accepted_new.emplace_back(node);
      prev_node = node; // advance
      current_z = prev_node.CurrentState.z;
    }
  }

  // Merge: append accepted new nodes, then re-sort by z
  if (!accepted_new.empty()) {
    KalmanNodes.insert(KalmanNodes.end(), accepted_new.begin(), accepted_new.end());
    SortNodesByZ();

    // Reindex step indices
    for (size_t i = 0; i < KalmanNodes.size(); ++i) {
      KalmanNodes[i].StepIndex = static_cast<int>(i);
    }
    // Skip global smoothing to keep interior stable; optionally we could
    // perform a trailing-window smoother/update only on the last N nodes.
    //runRTSSmoother();
    BetheBloch();
    Runchi2();
    SignSelection();
    UpdateParameters();
  }

  // Restore the previous direction flag
  ForwardFitting = prevForward;
}

// Static helper: compute a 2D measurement for a hit's bar at z using the
// previous state projected to that z. Coordinate mapping follows the existing
// measurement/noise model used in FillNoiseMatrix():
//  - X bars constrain y (x is weak), so use bar y and projected x
//  - Y bars constrain x (y is weak), so use bar x and projected y
//  - U/V bars are ±3° stereo around Y; use projected y to adjust x by tan(±3°)
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

// Used for seeding the starting KE for the Kalman filter from the start and end point of a track
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

// Update to the next step
void TMS_Kalman::RunKalman() {
 
  int nCand = KalmanNodes.size();
  if (nCand == 0) return;
  ConsecutiveOutliers = 0; // reset outlier streak at start
  ZLayerAccepted.clear();  // clear per-z acceptance map
  KalmanNodes[0].SetRecoXY(KalmanNodes[0].PreviousState);
  // Ensure the very first node has an inflated initial covariance.
  // This keeps the filter humble about the starting state.
  if (KalmanNodes[0].CovarianceMatrix.GetNrows() == 5 && KalmanNodes[0].CovarianceMatrix.GetNcols() == 5) {
    SetInflatedInitialCovariance(KalmanNodes[0].CovarianceMatrix);
  }
  KalmanNodes[0].StepIndex = 0;
  for (int i = 1; i < nCand; ++i) {
    // Perform the update from the (i-1)th node's predicted to the ith node's previous
    Update(KalmanNodes[i-1], KalmanNodes[i]);
    KalmanNodes[i].StepIndex = i; // keep index synced for adaptive gating
    Predict(KalmanNodes[i]);
  }

  // Update the momentum, and start/end position
  UpdateParameters();
}

void TMS_Kalman::UpdateParameters() {
  // Make sure of a consistent sorting
  SortNodesByZ();
  auto FrontNode = KalmanNodes.front();
  auto BackNode = KalmanNodes.back();
  auto DirectionNode = BackNode;
  if (KalmanNodes.size() > 1) DirectionNode = KalmanNodes.at(KalmanNodes.size() - 2);

  SetMomentum(1./FrontNode.CurrentState.qp);

  // Set start pos/dir
  SetStartPosition(FrontNode.CurrentState.x, FrontNode.CurrentState.y, FrontNode.CurrentState.z);
  SetStartDirection(FrontNode.CurrentState.dxdz, FrontNode.CurrentState.dydz);
  // Set end pos/dir
  SetEndPosition(BackNode.CurrentState.x, BackNode.CurrentState.y, BackNode.CurrentState.z);
  SetEndDirection(DirectionNode.CurrentState.dxdz, DirectionNode.CurrentState.dydz);

  if (std::isnan(momentum) || std::isinf(momentum)){
    //std::cerr << "[TMS_Kalmann.cpp] Weirdness -- Momentum from fitter isn't a sane number: " << momentum << std::endl;
  }
}

void TMS_Kalman::Update(TMS_KalmanNode &PreviousNode, TMS_KalmanNode &CurrentNode) {
  CurrentNode.PreviousState = PreviousNode.CurrentState;
  CurrentNode.CovarianceMatrix = PreviousNode.CovarianceMatrix;
}

// Predict the next step
// Here we use the previous state to predict the current state
// We also calculate the updated noise matrix from multiple scattering and energy loss
//jAccount for energy loss, multiple scattering, and bending due to the magnetic field
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

  // Calculate Lorentz force (deflection in x only)
  double p = 1.0 / PreviousState.qp;  // Total momentum
  //double px = p * PreviousState.dxdz; // Momentum in the x direction (proportional to slope dxdz)
                      

  // a crude calculation. delta px(momentum increase in the x direction)= f*delta_t = q*v*B*delta_z/ v, roughly 30 MeV per layer 
  // in natural unit, q= 0.303, 1T = 1.95*10^-10 MeV^2, 1mm = 5*10^9MeV^-1
  // Linearized deflection of slope in x due to B ~ By; units chosen to be consistent with internal code
  // Use the node's configured step size to avoid accidental
  // inconsistencies in CurrentState.z during augmentation.
  double dz = Node.dz;
  double magnetic_deflection_tx = 0.303 * assumed_charge * MagneticField * 1.95 * dz * 0.5 / std::max(p, 1e-9);
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

  if (Talk) std::cout << "\nPrevious vector: " << std::endl;
  if (Talk) PreviousState.Print();

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

  //std::cout << "Going from " << start.X() << " " << start.Y() << " " << start.Z() << std::endl;
  //std::cout << "To " << stop.X() << " " << stop.Y() << " " << stop.Z() << std::endl;
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

  // Update the state's q/p
  if (Talk) std::cout << "setting current p = " << p_up << std::endl;
  // Guard against zero/NaN momentum
  if (!std::isfinite(p_up) || p_up <= 0.0) {
    p_up = 1e-3; // clamp to small but finite momentum to avoid INF
  }
  CurrentState.qp = 1.0 / p_up;
  if (!std::isfinite(CurrentState.qp)) {
    CurrentState.qp = 1.0; // fallback to modest curvature if something went wrong
  }


  double p_var = (2*en/p_up)*(2*en/p_up) * total_en_var;
  // Guard qp variance against zero/NaN
  double p4 = p_up*p_up*p_up*p_up;
  double qp_var = (p4 > 0.0 && std::isfinite(p4)) ? (p_var / p4) : 0.0;
  if (!std::isfinite(qp_var)) qp_var = 0.0;


  // Set pointers to the noise matrix
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

  Node.FillNoiseMatrix(); // Full the matrix for multiple scattering

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

  TMatrixD S = H * CovarianceMatrix * HT; S += R; // innovation covariance
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

// Use the NoiseMatrix from the Kalman state to throw a vector of random noise // Liam: ??????????????
// to add to the state measurement
TVectorD TMS_Kalman::GetNoiseVector(TMS_KalmanNode Node) {
  // Clarence my man... what is this function even for?
  TVectorD rand_vec;
  rand_vec.ResizeTo(KALMAN_DIM);
  for(int i = 0; i < KALMAN_DIM; i++)
    rand_vec[i] = RNG.Gaus();


  TMatrixD &cov = Node.NoiseMatrix;
  TVectorD toy;
  toy.ResizeTo(5);
  toy.Zero();

  toy = cov*rand_vec;
  rand_vec.Print();

  return toy;
}

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

//reference:https://jwmi.github.io/ASM/6-KalmanFilter.pdf
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
