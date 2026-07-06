#ifndef _TMS_KALMAN_H_SEEN_
#define _TMS_KALMAN_H_SEEN_

#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_map>

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TRandom3.h"

// Include the physics for propgataion
#include "BetheBloch.h"
#include "MultipleScattering.h"
#include "Material.h"

// And the geometry for the materials
#include "TMS_Geom.h"
// Need to understand what a hit is
#include "TMS_Hit.h"
#include "TMS_Bar.h"

// Define the number of dimensions for a Kalman Node
#ifndef KALMAN_DIM
#define KALMAN_DIM 5
#endif


// State vector for the Kalman filter (not directly observed)
class TMS_KalmanState {
  public:
    TMS_KalmanState() = delete;
    TMS_KalmanState(double xvar, double yvar, double zvar, double dxdzvar, double dydzvar, double qpvar)
      : x(xvar), y(yvar), dxdz(dxdzvar), dydz(dydzvar), qp(qpvar), z(zvar) {
    };

    double x;
    double y;
    double dxdz;
    double dydz;
    double qp;
    double z; // the dependent variable of the state vector

    //TMatrixD &cov;//[KALMAN_DIM*KALMAN_DIM];

    void Print() {
      std::cout << "  {x, y, dx/dz, dy/dz, q/p, z} = {" << x << ", " << y << ", " << dxdz << ", " << dydz << ", " << qp << ", " << z << "}" << std::endl;
    }
};

// One Kalman node containing measurement, state, and matrices for a step
class TMS_KalmanNode {
  public:
  TMS_KalmanNode() = delete;

  // x,y at plane zvar; dzvar = zvar - z_prev (step size). PreviousState is
  // defined at z_prev = zvar - dzvar; Current/Smooth at zvar.
  TMS_KalmanNode(double xvar, double yvar, double zvar, double dzvar,double dxdzvar, double dydzvar) :
    x(xvar), y(yvar), z(zvar), dz(dzvar), dxdz(dxdzvar), dydz(dydzvar),
    StepIndex(-1),
    RecoX(xvar), RecoY(yvar),
    TrueX(-999999999.0), TrueY(-999999999.0),
    LayerOrientation(TMS_Bar::kError),
    LayerBarWidth(0.0),
    LayerBarLength(0.0),
    PlaneNumber(-1),
    LayerNotZ(-999999999.0),
    Hit(NULL),
    PreviousState(x, y, zvar - dzvar, dxdzvar, dydzvar, 1./20.),
    CurrentState(x, y, zvar,          dxdzvar, dydzvar, 1./20.), // Initialise at plane z
    SmoothState(x, y, zvar,           -999.9, -999.9, -1./20.),  // Initialise at plane z
    TransferMatrix(KALMAN_DIM,KALMAN_DIM),
    TransferMatrixT(KALMAN_DIM,KALMAN_DIM),
    NoiseMatrix(KALMAN_DIM,KALMAN_DIM),
    CovarianceMatrix(KALMAN_DIM,KALMAN_DIM),
    UpdatedCovarianceMatrix(KALMAN_DIM,KALMAN_DIM),
    EstimatedCovarianceMatrix(KALMAN_DIM, KALMAN_DIM),
    SmoothCovarianceMatrix(KALMAN_DIM, KALMAN_DIM),
    MeasurementMatrix(KALMAN_DIM,KALMAN_DIM),
    MeasurementVec(5),
    chi2(0.0),
    Accepted(true)
//    DeflectedVec(5)
  {
    TransferMatrix.ResizeTo(KALMAN_DIM, KALMAN_DIM);
    TransferMatrixT.ResizeTo(KALMAN_DIM, KALMAN_DIM);
    NoiseMatrix.ResizeTo(KALMAN_DIM, KALMAN_DIM);
    CovarianceMatrix.ResizeTo(KALMAN_DIM, KALMAN_DIM);
    UpdatedCovarianceMatrix.ResizeTo(KALMAN_DIM, KALMAN_DIM);
    EstimatedCovarianceMatrix.ResizeTo(KALMAN_DIM, KALMAN_DIM);
    SmoothCovarianceMatrix.ResizeTo(KALMAN_DIM, KALMAN_DIM);
    MeasurementMatrix.ResizeTo(KALMAN_DIM, KALMAN_DIM);
    rVec.ResizeTo(2);
    rVecT.ResizeTo(2);
    RMatrix.ResizeTo(2, 2);
    MeasurementVec.ResizeTo(5);
//    DeflectedVec.ResizeTo(5);

    // Build identity transfer and set dz terms for x/y from slopes
    TransferMatrix.Zero();
    TransferMatrixT.Zero(); // Transposed
    rVec.Zero();
    rVecT.Zero();
    RMatrix.Zero();
    MeasurementVec.Zero();
//    DeflectedVec.Zero();


    // Diagonal element
    for (int j = 0; j < KALMAN_DIM; ++j)
    {
      TransferMatrix(j,j)  = 1.;
      TransferMatrixT(j,j) = 1.;
    }

    // Could put in energy loss into transfer matrix?
    // dz for the slope
    TransferMatrix(0,2)  = TransferMatrix(1,3)  = dzvar; // Clarence matrix was this and 1 diagonal
    TransferMatrixT(2,0) = TransferMatrixT(3,1) = dzvar; // Transpose of previous
  }


  double x;
  double y;
  double z;
  double dz;
  double dxdz;
  double dydz;

  // Index of this node within the track (0-based). Used for adaptive gating.
  int StepIndex;

  double RecoX; // Reco X and Reco Y get updated with Kalman prediction info
  double RecoY;

  double TrueX; // True X and True Y to compare reco to truth
  double TrueY;

  TMS_Bar::BarType LayerOrientation;
  double LayerBarWidth;
  double LayerBarLength;
  int PlaneNumber;
  double LayerNotZ; // Bar center along the measured axis (mm)
  TMS_Hit *Hit;

  // The state vectors carry information about the covariance matrices etc
  TMS_KalmanState PreviousState;
  TMS_KalmanState CurrentState;
  TMS_KalmanState SmoothState;

  // Propagator matrix (from k-1 to k)
  TMatrixD TransferMatrix;
  TMatrixD TransferMatrixT; // Transposed

  // Process/measurement covariance storage
  TMatrixD NoiseMatrix;
  TMatrixD CovarianceMatrix;
  TMatrixD UpdatedCovarianceMatrix;
  TMatrixD EstimatedCovarianceMatrix;
  TMatrixD SmoothCovarianceMatrix;


  // Measurement matrix
  TMatrixD MeasurementMatrix;
  // Innovation and chi2 bookkeeping
  TVectorD rVec;
  TVectorD rVecT;
  TMatrixD RMatrix;
  TVectorD MeasurementVec;
//  TVectorD DeflectedVec;
  double chi2;
  bool Accepted; // true if the measurement was used (passed outlier cuts)


  void SetRecoXY(TMS_KalmanState& State)
  {
    RecoX = State.x;
    RecoY = State.y;
  }

  void SetTrueXY(double xarg, double yarg)
  {
    TrueX = xarg;
    TrueY = yarg;
  }

  void PrintTrueReco()
  {
    std::cout << "True x: " << TrueX << ",\t" << "Reco x: " << RecoX << ",\t" << "True - Reco: " << TrueX - RecoX << std::endl;
    std::cout << "True y: " << TrueY << ",\t" << "Reco y: " << RecoY << ",\t" << "True - Reco: " << TrueY - RecoY << std::endl;
  }

  void FillNoiseMatrix()
  {
    double H = 0.00274576; // (tan(3 deg))^2
    double A = LayerBarWidth;
    double B = 2000;       // bar-length uncertainty (mm)

    int sign;
    if (       LayerOrientation == TMS_Bar::kUBar) {
      sign = -1;
    } else if (LayerOrientation == TMS_Bar::kVBar) {
      sign =  1;
    } else if (LayerOrientation == TMS_Bar::kXBar) {
      NoiseMatrix(0,0) = B*B;
      NoiseMatrix(1,1) = A*A;
      NoiseMatrix(1,0) = NoiseMatrix(0,1) = 0.0;
      return;
    } else if (LayerOrientation == TMS_Bar::kYBar) {
      NoiseMatrix(0,0) = A*A;
      NoiseMatrix(1,1) = B*B;
      NoiseMatrix(1,0) = NoiseMatrix(0,1) = 0.0;
      return;
    } else {
      // Unknown orientation: fall back to diagonal measurement noise
      NoiseMatrix(0,0) = A*A;
      NoiseMatrix(1,1) = B*B;
      NoiseMatrix(1,0) = NoiseMatrix(0,1) = 0.0;
      return;
    }
    H *= sign;

    NoiseMatrix(0,0) = A*A;
    NoiseMatrix(1,1) = B*B;
    NoiseMatrix(1,0) = NoiseMatrix(0,1) = H*A*B;
  }

  void FillUpdatedCovarianceMatrix(double pathLength, double dxdz, double dydz, double qp, double ms, bool ForwardFitting=false)
  {
    // Now proceed with Wolin and Ho (Nucl Inst A329 1993 493-500)
    // covariance for multiple scattering
    // Also see MINOS note on Kalman filter (John Marshall, Nov 15 2005)
    double norm = 1 + dxdz*dxdz + dydz*dydz; // 1+P3^2+P4^2 in eq 16, 17, 18 in Wolin and Ho
    if (std::isnan(norm) || std::isinf(norm)) {
      std::cout<<"Norm is incorrect: "<<norm<<" with dxdz: "<<dxdz<<" and dydz: "<<dydz<<std::endl;
      norm = 1;
    }
    double covAxAx = norm*ms*(1 + dxdz*dxdz);// eq 16 Wolin and Ho
    if (std::isnan(covAxAx) || std::isinf(covAxAx)) {
      std::cout<<"covAxAx is incorrect: "<<norm<<" with dxdz: "<<dxdz<<" and dydz: "<<dydz<<" and ms: "<<ms<<std::endl;
      norm = 1;
    }
    double covAyAy = norm*ms*(1 + dydz*dydz);// eq 17 Wolin and Ho
    double covAxAy = norm*ms*dxdz*dydz;// eq 18 Wolin and Ho

    double TotalPathLengthSq = pathLength*pathLength;
    UpdatedCovarianceMatrix(0,0) = covAxAx * TotalPathLengthSq / 4.;
    UpdatedCovarianceMatrix(1,1) = covAyAy * TotalPathLengthSq / 4.;
    UpdatedCovarianceMatrix(2,2) = covAxAx;
    UpdatedCovarianceMatrix(3,3) = covAyAy;
    UpdatedCovarianceMatrix(4,4) = qp;

    // Negative signs depend on if we're doing backward or forward fitting (- sign for decreasing z, + sign for increasing z)
    int Sign = +1;
    if (!ForwardFitting) Sign *= -1;

    UpdatedCovarianceMatrix(1,0) = UpdatedCovarianceMatrix(0,1) = covAxAy * TotalPathLengthSq/4.;

    UpdatedCovarianceMatrix(2,0) = UpdatedCovarianceMatrix(0,2) = (Sign)*covAxAx * pathLength/2.;
    UpdatedCovarianceMatrix(3,1) = UpdatedCovarianceMatrix(1,3) = (Sign)*covAyAy * pathLength/2.;

    UpdatedCovarianceMatrix(3,0) = UpdatedCovarianceMatrix(0,3) = (Sign)*covAxAy * pathLength/2.;
    UpdatedCovarianceMatrix(2,1) = UpdatedCovarianceMatrix(1,2) = (Sign)*covAxAy * pathLength/2.;
    UpdatedCovarianceMatrix(3,2) = UpdatedCovarianceMatrix(2,3) = (Sign)*covAxAy;
  }

  bool operator<(const TMS_KalmanNode &other) const {
    return z < other.z;
  }
  bool operator>(const TMS_KalmanNode &other) const {
    return z > other.z;
  }
};

class TMS_Kalman {
  public:
    TRandom3 RNG;
    TMS_Kalman();
    TMS_Kalman(std::vector<TMS_Hit> &Candidates, double charge);
    // Minimal API to add a forward augment pass using a pool of extra hits.
    // Call this after constructing with the seed hits above. It will start
    // from the smoothed state at the front (lowest z) and propagate forward,
    // considering each candidate in z order and accepting those that pass
    // outlier gating. Accepted hits are merged into the node list and the
    // track is re-smoothed.
    void AugmentWithCandidates(const std::vector<TMS_Hit> &candidate_pool, const size_t number_to_remove = 5);

    // Snap downstream hits within a large cylinder aligned to +z whose axis
    // passes through the last accepted node's (x0,y0) position. Creates up to
    // max_planes new nodes, forces their measurements onto the cylinder axis
    // projection for that plane, merges, and refits. Intended to recover the
    // final 2–3 planes that are commonly missing after augment.
    // Defaults: radius_mm = 500 mm (50 cm), max_planes = 5.
    // Snap any downstream hits within a forward z-window [z_last, z_last + dz_mm]
    // aligned to +z. Default dz is 500 mm (50 cm). Outlier rejection is
    // disabled during the refit, and duplicates at identical z are avoided
    // by nudging z slightly.
    void SnapDownstreamHitsAndRefit(const std::vector<TMS_Hit> &pool,
                                    double dz_mm = 500.0,
                                    int max_nodes = 0,
                                    double max_axis_distance = -1.0);

    // Print a concise per-node diagnostic summary. Includes z, view, acceptance,
    // chi2, residuals, and measured vs filtered positions. Only prints if Talk is true.
    void PrintNodesDiagnostics(const char* tag = "KalmanNodes");

    // Scan the node sequence for expected X-U-V view coverage. Reports windows
    // of 3 consecutive nodes that are missing one or more views. Only prints if Talk is true.
    void PrintMissingXUVTriplets(const char* tag = "XUVPattern");

    double Start[3];
    double End[3];
    double StartDirection[3];
    double EndDirection[3];

    void InitializeMomentum(bool only_momentum = true);
    double GetKEEstimateFromLength(double startx, double endx, double startz, double endz);
    double CalculateKEFromSteel(TVector3 position);

    void SetMomentum(double mom) {momentum = mom;}
    void SetCharge_curvature(double charge) {
        if (charge > 0) charge_curvature= +1;
        else charge_curvature= -1;
    }

    // Set direction unit vectors from only x and y slope
    void SetStartDirection(double ax, double ay);// {StartDirection[0]=ax; StartDirection[1]=ay; StartDirection[2]=sqrt(1 - ax*ax - ay*ay);};
    void SetEndDirection  (double ax, double ay);// {EndDirection[0]=ax;   EndDirection[1]=ay;   EndDirection[2]=sqrt(1 - ax*ax - ay*ay);};

    // Set position unit vectors
    void SetStartPosition(double ax, double ay, double az) {Start[0]=ax; Start[1]=ay; Start[2]=az;};
    void SetEndPosition  (double ax, double ay, double az) {End[0]=ax;   End[1]=ay;   End[2]=az;};

    double GetMomentum() {return momentum;}
    double GetCharge_curvature() {return charge_curvature;}

    double GetTrackChi2()
    {
      double tmp_chi2 = 0.0;
      for (auto node : KalmanNodes)
        tmp_chi2 += node.chi2;

      return tmp_chi2;
    }


    std::vector<TMS_KalmanNode> GetKalmanNodes() {return KalmanNodes;}
    bool GetWasAugmented() { return wasAugmented; }
    int GetNAugmentedNodes() { return nAugmentedNodes; }

    // Baseline (pre-snap) diagnostics
    bool   HasPreSnapBaseline() const { return hasPreSnapBaseline; }
    double GetPreSnapChi2() const { return preSnapChi2; }
    void   GetPreSnapEnd(double &x, double &y, double &z) const { x = preSnapEnd[0]; y = preSnapEnd[1]; z = preSnapEnd[2]; }

    TVectorD GetNoiseVector(TMS_KalmanNode Node);

    void SetTalk(bool value) { Talk = value; };

  private:
    // Energy-loss calculator
    BetheBloch_Calculator Bethe;
    MultipleScatter_Calculator MSC;

    void Predict(TMS_KalmanNode &Node);
    void Update(TMS_KalmanNode &PreviousNode, TMS_KalmanNode &CurrentNode);
    void RunKalman();
    void PruneRejectedAndDuplicateZ();
    void runRTSSmoother();
    void BetheBloch();
    void SignSelection();
    void Runchi2();
    void UpdateParameters();
    // Refit helper: backward → forward → backward passes with smoothing
    void TripleRefitSmooth();

    // Repair obviously bad node measurements (x,y) by recomputing them
    // from bar geometry using the previous state's projection as seed.
    // Only adjusts nodes that exceed residual/bounds heuristics.
    void RepairBadNodeMeasurements(bool extrapolate = true);

    // Remove nodes that are inconsistent after rebuilding/refit.
    // Heuristics:
    //  - Large residual vs projection from previous smoothed state
    //  - High-frequency oscillation pattern across consecutive nodes
    void PruneInconsistentNodes();
    int PruneNodesOutsideTMS();
    void RepairStatesOutsideTMS();

    void SortNodesByZ() { std::sort(KalmanNodes.begin(), KalmanNodes.end(), [](const TMS_KalmanNode &a, const TMS_KalmanNode &b){return a.z < b.z;}); };
    // Sorts so that first node is first node processed by RunKalman. ie. based on ForwardFitting flag.
    void SortNodesByRunOrder() {
        const bool fit_direction = ForwardFitting;
        std::sort(KalmanNodes.begin(), KalmanNodes.end(), [fit_direction](const TMS_KalmanNode &a, const TMS_KalmanNode &b)
        {return fit_direction ? a.z < b.z : a.z > b.z;});
      };

    // After sorting in the current run order, rebuild each node's dz and
    // transfer matrices to reflect adjacency in that order. Also ensure
    // PreviousState.z = z - dz and CurrentState.z = z at each node.
    void RebuildRunOrderSteps();

    // Build a 2D measurement (x,y) at hit z using bar geometry and
    // the current track state projected to that z. For bars that do not
    // directly measure one axis, use the projected value for that axis.
    // For U/V stereo bars, apply a ±3° rotation correction using the
    // projected y to adjust the effective x measurement.
    static void BuildBarMeasurement(const TMS_Hit &hit,
                                    const TMS_KalmanState &state_at_prev,
                                    double meas_z,
                                    double &meas_x,
                                    double &meas_y);
    void BuildMeasurementsFromBar();



    // State vector
    // x, y, dx/dz, dy/dz, q/p
    std::vector<TMS_KalmanNode> KalmanNodes;

    // Remember if we're forward fitting
    bool ForwardFitting;

    double total_en;
    double mass;
    double momentum;
    double charge_curvature;
    double assumed_charge;
    double AverageXSlope; // Seeding initial X slope in Kalman
    double AverageYSlope; // Seeding initial Y slope in Kalman

    bool wasAugmented;
    int nAugmentedNodes;

    bool Talk;

    // Diagnostics captured before snap augmentation
    bool   hasPreSnapBaseline = false;
    double preSnapEnd[3] = {0.0, 0.0, 0.0};
    double preSnapChi2 = 0.0;

    // True during augmentation pass so we can tame qp updates
    bool AugmentationMode = false;
    int  AugmentFreezeSteps = 2; // freeze physics for first few augmented nodes

    // Outlier handling state
    int ConsecutiveOutliers = 0;           // running count during filtering
    int OutlierResetThreshold = TMS_Manager::GetInstance().Get_Reco_Kalman_Outlier_ResetThreshold();
    // Track if a z-layer already has an accepted measurement, to reject duplicates at same z
    std::unordered_map<long long, bool> ZLayerAccepted;
    static long long QuantizeZKey(double z) { return llround(z * 1000.0); } // quantize to micron to group same-z hits

    // Outlier rejection control
    bool EnableOutlierRejection = TMS_Manager::GetInstance().Get_Reco_Kalman_Use_Outlier_Rejection();
    // 2-DOF chi2 threshold (tunable)
    double Chi2RejectThreshold = TMS_Manager::GetInstance().Get_Reco_Kalman_Outlier_Rejection_Chi2_Threshold();
};

#endif
