#include "TMS_Hit.h"

#ifndef _TMS_TRACK_H_SEEN_
#define _TMS_TRACK_H_SEEN_

// General 3D-Track class
class TMS_Track {

  public:
    TMS_Track() //std::vector<TMS_Hit>& OneTrack, std::vector<TMS_Hit>& OtherTrack)
    {
      // TODO: Take first and last hits and do the maffs
      0;
    };
    void Print();

    double Start[3];      // Start point in x,y,z
    double End[3];        // End point in x,y,z
    double Direction[3];  // Unit vector in track direction
    double Length;
    double EnergyDeposit;
    double EnergyRange;
    double ActivityStart; // Measure of activity near start of track
    double ActivityEnd;   // Measure of activity near end of track
    double Time;          // TODO: Fill this in a sensible way

    double GetEnergyDeposit(){return EnergyDeposit;};
    double GetEnergyRange(){return EnergyRange;};

    int nHits;
    std::vector<TMS_Hit*> Hits;
};


#endif
