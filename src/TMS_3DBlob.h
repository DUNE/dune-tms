#ifndef _TMS_3DBLOB_H_SEEN_
#define _TMS_3DBLOB_H_SEEN_

#include <string>

// Include the constants
#include "TMS_Constants.h"

#include "TMS_Bar.h"
#include "TMS_Hit.h"

// Class describing a '3D Blob'; A Blob is an imaginary 3-dimensional object representing
// the region in which the hits from a track are expected to belong to.
// A Blob is constructed by fitting a track to the hits, and then using this track to
// determine a 'best-fit' value of sorts for the coordinates.
// The extent of the Blob in a dimension is the positional uncertainty.

class TMS_3DBlob {

  public:
    // The constructor for the TMS hit
    TMS_3DBlob() { // Initialise the object
       X = 0;  Y = 0;  Z = 0; // Central blob position
      dX = 0; dY = 0; dZ = 0; // Blob extent. ~= uncertainty
      Energy     = -999;
      TrueEnergy = -999;
    };
    virtual ~TMS_3DBlob() {};

    void Print() const;

    //std::vector<TMS_Hit> Hits; // Vector of actual detector hits in one layer that contribute
    //std::vector<TMS_Bar> Bars; // Vector of actual bars hit and associated with the track

    double X, Y, Z;  // Central blob position
    double dX,dY,dZ; // Blob extent. ~= uncertainty

    // Energy deposited that is associated to the track in this layer
    double Energy;
    double TrueEnergy; // Can we even fill this with something useful? (:
};

#endif
