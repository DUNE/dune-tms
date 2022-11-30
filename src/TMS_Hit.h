#ifndef _TMS_HIT_H_SEEN_
#define _TMS_HIT_H_SEEN_

#include <string>

// Include the constants
#include "TMS_Constants.h"
#include "TMS_Bar.h"
#include "TMS_TrueHit.h"
// To get the geometry
#include "TMS_Geom.h"

// Include the edep-sim
#include "EDepSim/TG4HitSegment.h"

// Maybe should be in singleton?
#include "TGeoManager.h"
#include "TGeoNavigator.h"

// A low-level hit
class TMS_Hit {

  public:
    void Print() const;
    // The constructor for the TMS hit
    TMS_Hit(TG4HitSegment &edep_seg);

    const TMS_Bar &GetBar() const { return Bar; };
    void SetBar(TMS_Bar bar) { Bar = bar; };

    // Sort by decreasing Z
    static bool SortByZ(TMS_Hit &a, TMS_Hit &b) {
      return ( a.GetBar().GetPlaneNumber() > b.GetBar().GetPlaneNumber() );
    }

    // Sort by increasing Z
    static bool SortByZInc(TMS_Hit &a, TMS_Hit &b) {
      return ( a.GetBar().GetPlaneNumber() < b.GetBar().GetPlaneNumber() );
    }

    // A helper function to determine if a hit is close to a gap
    bool NextToGap();

    // Sort by increasing T
    static bool SortByT(TMS_Hit &a, TMS_Hit &b) {
      return ( a.GetT() < b.GetT() );
    }
    // Sort by increasing Z
    static bool SortByZThenT(TMS_Hit &a, TMS_Hit &b) {
      if ( a.GetBar().GetPlaneNumber() == b.GetBar().GetPlaneNumber() ) return a.GetT() < b.GetT();
      return ( a.GetBar().GetPlaneNumber() < b.GetBar().GetPlaneNumber() );
    }

    // The true hit
    const TMS_TrueHit &GetTrueHit() const { return TrueHit; };

    // Over-riders (maybe delete in future)
    void SetTrueHit(TMS_TrueHit hit) {TrueHit = hit;};

    void SetE(double E) {EnergyDeposit = E;};
    void SetT(double t) {Time = t;};
    
    void SetPedSup(bool isPedSup) { PedSuppressed = isPedSup;};
    bool GetPedSup() { return PedSuppressed; };
    
    void SetPE(double pe) { PE = pe; };
    double GetPE() { return PE; };

    double GetE() const {return EnergyDeposit;};
    double GetT() const {return Time;};
    
    void SetSlice(int slice) { Slice = slice; };
    int GetSlice() { return Slice; };

    double GetX() const { return Bar.GetX(); };
    double GetY() const { return Bar.GetY(); };
    double GetZ() const { return Bar.GetZ(); };
    double GetNotZ() const { return Bar.GetNotZ(); };

    double GetXw() const { return Bar.GetXw(); };
    double GetYw() const { return Bar.GetYw(); };
    double GetZw() const { return Bar.GetZw(); };
    double GetNotZw() const { return Bar.GetNotZw(); };

    int GetPlaneNumber() const {return Bar.GetPlaneNumber(); };
    int GetBarNumber() const {return Bar.GetBarNumber(); };

  private:
    // The true hit (x,y,z,t) --- does not quantise hit into bars
    TMS_TrueHit TrueHit;
    // The true particle that created this hit
    // The bar that registered the hit
    TMS_Bar Bar;
    // The energy deposited
    double EnergyDeposit;
    // The timing of the hit
    double Time;
    
    int Slice;
    
    
    bool PedSuppressed;
    double PE;
};

inline bool operator==(const TMS_Hit &a, const TMS_Hit &b) {
  if (a.GetZ()    == b.GetZ() && 
      a.GetNotZ() == b.GetNotZ() && 
      a.GetE()    == b.GetE() &&
      a.GetT()    == b.GetT() ) return true;
  return false;
}

#endif
