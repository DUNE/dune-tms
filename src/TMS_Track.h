#ifndef _TMS_TRACK_H_SEEN_

#include "TMS_Hit.h"
#include "TMS_Kalman.h"
#include "TMS_TrueParticle.h"

#define _TMS_TRACK_H_SEEN_

// General 3D-Track class
class TMS_Track {

  public:
    TMS_Track() //std::vector<TMS_Hit>& OneTrack, std::vector<TMS_Hit>& OtherTrack)
    {
      // TODO: Take first and last hits and do the maffs
    };
    void Print();

    int    Charge;
    double Start[3];     // Start point in x,y,z
    double End[3];       // End point in x,y,z
    double Direction[3]; // Unit vector in track direction
    double Length;
    double Occupancy;
    double EnergyDeposit;
    double EnergyRange;
    double Momentum;
    double Time;         // TODO: Fill this in a sensible way

    double GetEnergyDeposit() {return EnergyDeposit;};
    double GetEnergyRange()   {return EnergyRange;};
    double GetMomentum()      {return Momentum;};

    TMS_TrueParticle GetTrueParticle() {return fTrueParticle;};

    // Manually set variables
    void SetEnergyDeposit(double val) {EnergyDeposit = val;};
    void SetEnergyRange  (double val) {EnergyRange   = val;};
    void SetMomentum     (double val) {Momentum      = val;};

    int nHits;
    std::vector<TMS_Hit> Hits;
    int nKalmanNodes;
    std::vector<TMS_KalmanNode> KalmanNodes;

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


  // a lot of the vars from above can be moved into this in future
  private:
    TMS_TrueParticle fTrueParticle;

};


#endif
