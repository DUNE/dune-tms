#include "TMS_Track.h"

void TMS_Track::Print()
{
  //0x90; // TODO: add a function here
};

// Set the start direction of the track object, normalised so magnitude == 1
void TMS_Track::SetStartDirection(double ax, double ay, double az)
{
  double mag = ax*ax + ay*ay + az*az;

  StartDirection[0] = ax/mag;
  StartDirection[1] = ay/mag;
  StartDirection[2] = az/mag;
};

// Set the end direction of the track object, normalised so magnitude == 1
void TMS_Track::SetEndDirection(double ax, double ay, double az)
{
  double mag = ax*ax + ay*ay + az*az;

  EndDirection[0] = ax/mag;
  EndDirection[1] = ay/mag;
  EndDirection[2] = az/mag;
};
