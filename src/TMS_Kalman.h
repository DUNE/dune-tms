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

    void Print() {
      //std::cout << "Printing Kalman node: " << std::endl;
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
    CurrentState(x, y, z+dz, -999.9, -999.9, -1./20.), // Initialise the state vectors
    PreviousState(x, y, z, -999.9, -999.9, -1./20.),
    TransferMatrix(KALMAN_DIM,KALMAN_DIM),
    NoiseMatrix(KALMAN_DIM,KALMAN_DIM),
    MeasurementMatrix(KALMAN_DIM,KALMAN_DIM) {

    TransferMatrix.ResizeTo(KALMAN_DIM, KALMAN_DIM);
    NoiseMatrix.ResizeTo(KALMAN_DIM, KALMAN_DIM);
    MeasurementMatrix.ResizeTo(KALMAN_DIM, KALMAN_DIM);

    // Make the transfer matrix for each of the states
    // Initialise to zero
    TransferMatrix.Zero();
    // Diagonal element
    for (int j = 0; j < KALMAN_DIM; ++j) TransferMatrix(j,j) = 1.;
    // Could put in energy loss into transfer matrix?
    // dz for the slope
    TransferMatrix(0,2) = TransferMatrix(1,3) = dzvar; // Clarence matrix was this and 1 diagonal
  }


  double x;
  double y;
  double z;
  double dz;

  double RecoX; // Reco X and Reco Y get updated with Kalman prediction info
  double RecoY;

  TMS_Bar::BarType LayerOrientation;

  // The state vectors carry information about the covariance matrices etc
  TMS_KalmanState CurrentState;
  TMS_KalmanState PreviousState;

  // Propagator matrix
  // Takes us from detector k-1 to detector k
  TMatrixD TransferMatrix;

  // Random variable w(k-1) includes random disturbances of track between z(k-1) and z(k) from multiple scattering
  // Noise matrix
  TMatrixD NoiseMatrix;

  // Measurement matrix
  TMatrixD MeasurementMatrix;

  void SetRecoXY(TMS_KalmanState& State)
  {
    RecoX = State.x;
    RecoY = State.y;
  }

  TMatrixD GetRecoNoiseMatrix()
  {
    // vvv   args: n_row, n_col, elements
//    double mat[25] = {9.0, 0, 0, 0, 0,
//                      0, 900.0, 0, 0, 0,
//                      0, 0, 0.05, 0, 0,
//                      0, 0, 0, 0.05, 0,
//                      0, 0, 0, 0, 0.01};
    double H = 0.00274576; // ( tan(3 deg) )**2
    int sign;
    if (LayerOrientation == TMS_Bar::kUBar) {
      sign = -1;
    } else if (LayerOrientation == TMS_Bar::kVBar) {
      sign =  1;
    } else {
      throw; // xd
    }
    H *= sign;

    double A = 0.0001;
    double B = 3.5;
    double mat[25] = {A*A, H*A*B, 0, 0, 0,
                      H*A*B, B*B, 0, 0, 0,
                      0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0};
    return TMatrixD(5,5, mat);
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
    TMS_Kalman(std::vector<TMS_Hit> &Candidates);

    double GetKEEstimateFromLength(double startx, double endx, double startz, double endz);

    void SetMomentum(double mom) {momentum = mom;}

    double GetMomentum() {return momentum;}

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

    double AverageXSlope; // Seeding initial X slope in Kalman
    double AverageYSlope; // Seeding initial Y slope in Kalman

    bool Talk;
};

#endif
