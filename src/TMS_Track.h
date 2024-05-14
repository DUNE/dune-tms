#include "TMS_Hit.h"

#ifndef _TMS_TRACK_H_SEEN_
#define _TMS_TRACK_H_SEEN_

// General 3D-Track class
class TMS_Track {

  public:
    TMS_Track() //std::vector<TMS_Hit>& OneTrack, std::vector<TMS_Hit>& OtherTrack)
    {
      // TODO: Take first and last hits and do the maffs
      //0;
    };
    void Print();

    int    Charge;
    double Start[3];     // Start point in x,y,z, in mm
    double End[3];       // End point in x,y,z, in mm
    double Direction[3]; // Unit vector in track direction
    double Length;       // in mm
    double Occupancy;
    double EnergyDeposit;// MeV
    double EnergyRange;  // MeV
    double Time;         // TODO: Fill this in a sensible way

    double GetEnergyDeposit(){return EnergyDeposit;}; // MeV
    double GetEnergyRange(){return EnergyRange;};     // MeV

    int nHits;
    std::vector<TMS_Hit> Hits;
};


#endif
