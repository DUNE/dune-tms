#ifndef _TMS_KALMAN_H_SEEN_
#define _TMS_KALMAN_H_SEEN_

#include <iostream>

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


// State vector for Kalman filter
// Not actually ever truly observed (not a measurement!)
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

// One node in the Kalman code
// Has information about the two states i and i+1
class TMS_KalmanNode {
  public:
  TMS_KalmanNode() = delete;

  // x,y,z, delta_z = distance from previous hit to current in z
  TMS_KalmanNode(double xvar, double yvar, double zvar, double dzvar) :
    x(xvar), y(yvar), z(zvar), dz(dzvar),
    RecoX(xvar), RecoY(yvar),
    CurrentState(x, y, z+dz, -999.9, -999.9, -1./20.), // Initialise the state vectors
    PreviousState(x, y, z, -999.9, -999.9, -1./20.),
    TransferMatrix(KALMAN_DIM,KALMAN_DIM),
    TransferMatrixT(KALMAN_DIM,KALMAN_DIM),
    NoiseMatrix(KALMAN_DIM,KALMAN_DIM),
    CovarianceMatrix(KALMAN_DIM,KALMAN_DIM),
    UpdatedCovarianceMatrix(KALMAN_DIM,KALMAN_DIM),
    MeasurementMatrix(KALMAN_DIM,KALMAN_DIM)
  {
    TransferMatrix.ResizeTo(KALMAN_DIM, KALMAN_DIM);
    TransferMatrixT.ResizeTo(KALMAN_DIM, KALMAN_DIM);
    NoiseMatrix.ResizeTo(KALMAN_DIM, KALMAN_DIM);
    CovarianceMatrix.ResizeTo(KALMAN_DIM, KALMAN_DIM);
    UpdatedCovarianceMatrix.ResizeTo(KALMAN_DIM, KALMAN_DIM);
    MeasurementMatrix.ResizeTo(KALMAN_DIM, KALMAN_DIM);
    rVec.ResizeTo(2);
    rVecT.ResizeTo(2);
    RMatrix.ResizeTo(2, 2);

    // Make the transfer matrix for each of the states
    // Initialise to zero
    TransferMatrix.Zero();
    TransferMatrixT.Zero(); // Transposed
    rVec.Zero();
    rVecT.Zero();
    RMatrix.Zero();


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

  double RecoX; // Reco X and Reco Y get updated with Kalman prediction info
  double RecoY;

  double TrueX; // True X and True Y to compare reco to truth
  double TrueY;

  TMS_Bar::BarType LayerOrientation;
  double LayerBarWidth;
  double LayerBarLength;

  // The state vectors carry information about the covariance matrices etc
  TMS_KalmanState CurrentState;
  TMS_KalmanState PreviousState;

  // Propagator matrix
  // Takes us from detector k-1 to detector k
  TMatrixD TransferMatrix;
  TMatrixD TransferMatrixT; // Transposed

  // Random variable w(k-1) includes random disturbances of track between z(k-1) and z(k) from multiple scattering
  // Noise matrix
  TMatrixD NoiseMatrix;
  TMatrixD CovarianceMatrix;
  TMatrixD UpdatedCovarianceMatrix;
  // Measurement matrix
  TMatrixD MeasurementMatrix;
  // For chi2 stuff
  TVectorD rVec;
  TVectorD rVecT;
  TMatrixD RMatrix;
  double chi2;


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
    double H = 0.00274576; // ( tan(3 deg) )**2
    double A = LayerBarWidth; //10.0; //10.0 mm bar width based uncert
    //double B = LayerBarLength;//4000.0; //4000.0 mm bar length based uncert
    double B = 2000;//4000.0; //4000.0 mm bar length based uncert

    int sign;
    if (       LayerOrientation == TMS_Bar::kUBar) {
      sign = -1;
    } else if (LayerOrientation == TMS_Bar::kVBar) {
      sign =  1;
    } else if (LayerOrientation == TMS_Bar::kXBar) { // this should just work right?
      NoiseMatrix(0,0) = B*B;
      NoiseMatrix(1,1) = A*A;
      NoiseMatrix(1,0) = NoiseMatrix(0,1) = 0.0;
      return;
    } else if (LayerOrientation == TMS_Bar::kYBar) { // this should just work right?
      NoiseMatrix(0,0) = A*A;
      NoiseMatrix(1,1) = B*B;
      NoiseMatrix(1,0) = NoiseMatrix(0,1) = 0.0;
      return;
    } else {
      throw; // xd haha TODO tho
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
    //double norm = 1+dxdz+dydz; // 1+P3^2+P4^2 in eq 16, 17, 18 in Wolin and Ho
    double norm = 1 + dxdz*dxdz + dydz*dydz; // 1+P3^2+P4^2 in eq 16, 17, 18 in Wolin and Ho
    double covAxAx = norm*ms*(1 + dxdz*dxdz);// eq 16 Wolin and Ho
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
    
    double Start[3];
    double End[3];
    double StartDirection[3];
    double EndDirection[3];
   
    double GetKEEstimateFromLength(double startx, double endx, double startz, double endz);

    void SetMomentum(double mom) {momentum = mom;}
    // Set direction unit vectors from only x and y slope
    void SetStartDirection(double ax, double ay);// {StartDirection[0]=ax; StartDirection[1]=ay; StartDirection[2]=sqrt(1 - ax*ax - ay*ay);};
    void SetEndDirection  (double ax, double ay);// {EndDirection[0]=ax;   EndDirection[1]=ay;   EndDirection[2]=sqrt(1 - ax*ax - ay*ay);};

    // Set position unit vectors
    void SetStartPosition(double ax, double ay, double az) {Start[0]=ax; Start[1]=ay; Start[2]=az;};
    void SetEndPosition  (double ax, double ay, double az) {End[0]=ax;   End[1]=ay;   End[2]=az;};

    double GetMomentum() {return momentum;}

    double GetTrackChi2()
    {
      double tmp_chi2 = 0.0;
      for (auto node : KalmanNodes)
        tmp_chi2 += node.chi2;

      return tmp_chi2;
    }


    std::vector<TMS_KalmanNode> GetKalmanNodes() {return KalmanNodes;}

    TVectorD GetNoiseVector(TMS_KalmanNode Node);

  private:
    // Energy-loss calculator
    BetheBloch_Calculator Bethe;
    MultipleScatter_Calculator MSC;

    void Predict(TMS_KalmanNode &Node);
    void Update(TMS_KalmanNode &PreviousNode, TMS_KalmanNode &CurrentNode);
    void RunKalman();


    // State vector
    // x, y, dx/dz, dy/dz, q/p
    std::vector<TMS_KalmanNode> KalmanNodes;

    // Remember if we're forward fitting
    bool ForwardFitting;

    double total_en;
    double mass;
    double momentum;
    double assumed_charge;
    double AverageXSlope; // Seeding initial X slope in Kalman
    double AverageYSlope; // Seeding initial Y slope in Kalman
    
    bool Talk;
};

#endif
