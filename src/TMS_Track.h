#ifndef _TMS_TRACK_H_SEEN_

#include "TMS_Hit.h"
#include "TMS_Kalman.h"
#include "TMS_TrueParticle.h"

#define _TMS_TRACK_H_SEEN_

// General 3D-Track class
class TMS_Track {

  public:
    TMS_Track();
    // Phase III: the hand-rolled copy constructor/assignment operator that used to live here
    // existed solely to shepherd the now-removed fTrueParticle member around (confirmed dead:
    // never set anywhere in the codebase). Every other member is a plain value type or STL
    // container, so the compiler-generated copy constructor/assignment operator (member-wise
    // copy) behaves identically -- no need to hand-write them anymore.


    void Print();

    int    Charge;
    int    Charge_Kalman;
    int    Charge_Kalman_curvature;
    double Start[4];     // Start point in x,y,z,t
    double End[4];       // End point in x,y,z,t
    double StartDirection[3]; // Unit vector in track direction at start
    double EndDirection[3]; // Unit vector in track direction at end
    double Length;
    double Occupancy;
    double EnergyDeposit;
    double EnergyRange;
    double Momentum;
    double Time;         // TODO: Fill this in a sensible way
    double Chi2;
    double Chi2_minus;
    double Chi2_plus;
    // Fit-quality counts, filled by the Cluster3D reconstruction (-1 = not
    // filled, e.g. legacy tracks): the fit's degrees of freedom (Chi2 / NDoF
    // is its chi2 per degree of freedom), the point layers the Kalman walk
    // stepped through, how many of those had no accepted hit, and the hits
    // taken outside any space point (orphan pickup and the single-hit
    // extension). 2026-09-27, 15 files: tracks whose hits no single particle
    // owns have a median chi2/ndf of 6.1 (muon tracks 1.4) and gap fraction
    // NGapLayers / NLayersWalked of 0.33 (muons 0); chi2/ndf > 6 flags half of
    // them at a cost of 5% of muon tracks (and 36% of real pion/proton tracks,
    // which look much like them).
    int NDoF = -1;
    int NLayersWalked = -1;
    int NGapLayers = -1;
    int NOrphanHits = -1;
    


    double GetEnergyDeposit() {return EnergyDeposit;};
    double GetEnergyRange()   {return EnergyRange;};
    double GetMomentum()      {return Momentum;};
    double GetChi2_minus()          {return Chi2_minus;};
    double GetChi2_plus()          {return Chi2_plus;};

    // Manually set variables
    void SetEnergyDeposit (double val) {EnergyDeposit = val;};
    void SetEnergyRange   (double val) {EnergyRange   = val;};
    void SetMomentum      (double val) {Momentum      = val;};
    void SetChi2_minus          (double val) {Chi2_minus          = val;};
    void SetChi2_plus          (double val) {Chi2_plus          = val;};
    void SetChi2               (double val) {Chi2               = val;};
   

    // Set direction unit vectors from only x and y slope
    void SetStartDirection(double ax, double ay, double az);// {StartDirection[0]=ax; StartDirection[1]=ay; StartDirection[2]=az;};
    void SetEndDirection  (double ax, double ay, double az);// {EndDirection[0]=ax;   EndDirection[1]=ay;   EndDirection[2]=az;};

    // Set position unit vectors
    void SetStartPosition(double ax, double ay, double az) {Start[0]=ax; Start[1]=ay; Start[2]=az;};
    void SetEndPosition  (double ax, double ay, double az) {End[0]=ax;   End[1]=ay;   End[2]=az;};

    int nHits;
    std::vector<TMS_Hit> Hits;

    // Kalman filter track info
    std::vector<TMS_KalmanNode> KalmanNodes;
    std::vector<TMS_KalmanNode> KalmanNodes_minus;
    std::vector<TMS_KalmanNode> KalmanNodes_plus;

    void Compare()
    {
      std::cout << "Lengths: Hits " << Hits.size() << " KalmanNodes " << KalmanNodes.size() << std::endl;
      double PrevPlane = -1.0;
      for (long unsigned int i=0; i<KalmanNodes.size(); i++)
      {
        if (PrevPlane == Hits[i].GetZ())
          continue; // TODO: Skip duplicate planes, Kalman currently ignores them

        PrevPlane = Hits[i].GetZ();

        std::cout << "Layer " << Hits[i].GetZ()
                  << "\t" << KalmanNodes[i].z << std::endl;
      }
    }

    void ApplyTrackSmoothing();
    double CalculateTrackSmoothnessY();
    void LookForHitsOutsideTMS();


  // a lot of the vars from above can be moved into this in future
  private:
    void setDefaultUncertainty();
    std::vector<size_t> findYTransitionPoints();
    double getAvgYSlopeBetween(size_t ia, size_t ib) const;
    double getMaxAllowedSlope(size_t ia, size_t ib) const;
    void simpleTrackSmoothing();

};


#endif
