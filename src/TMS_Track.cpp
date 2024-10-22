#include "TMS_Track.h"

void TMS_Track::Print()
{
  //0x90; // TODO: add a function here
};

// Set the start direction of the track object, normalised so magnitude == 1
void TMS_Track::SetStartDirection(double ax, double ay, double az) {
  double mag = ax*ax + ay*ay + az*az;
 
  Start[0] = ax/mag;
  Start[1] = ay/mag;
  Start[2] = az/mag;
};

// Set the end direction of the track object, normalised so magnitude == 1
void TMS_Track::SetEndDirection(double ax, double ay, double az) {
  double mag = ax*ax + ay*ay + az*az;

  End[0] = ax/mag;
  End[1] = ay/mag;
  End[2] = az/mag;
};
