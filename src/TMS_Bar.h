#ifndef _TMS_BAR_H_SEEN_
#define _TMS_BAR_H_SEEN_

// Include the constants
#include "TMS_Constants.h"
// To get the geometry
#include "TMS_Geom.h"

#include "EDepSim/TG4HitSegment.h"

#include "TGeoBBox.h"

class TMS_Bar {
  public:

    TMS_Bar(TG4HitSegment &edep_seg);
    TMS_Bar(double x, double y, double z);

    // Enum for the x, y, U, V bar orientation
    enum BarType { kXBar, kYBar, kUBar, kVBar, kError };
    std::string BarType_ToString(BarType bar) const;

    // Getter functions
    int GetBarNumber() const { return BarNumber; }; // Contiguous transverse bar index within its plane, from TMS_Geom's node-path lookup (not a coordinate offset)
    int GetBarWidth() const { return BarWidth; }; // Get Bar width in mm
    int GetBarLength() const { return BarLength; }; // Get Bar length in mm
    int GetPlaneNumber() const { return PlaneNumber; }; // Plane or layer number through the detector starting at smallest z
    int GetGlobalBarNumber() const { return GlobalBarNumber; }; // Number of hit Scintillator Module (4 modules)

    BarType GetBarType() const { return BarOrient; };
    int GetBarTypeNumber() const;

    double GetX() const { return x; };
    double GetY() const { return y; };
    double GetZ() const { return z; };

    // The bar's own global-frame center along its readout axis (x for X-bars, y for
    // U/V/Y-bars), captured in FindModules() before GetX()/GetY() gets sentinelled out for
    // that axis. Needed by the optical model to place each bar's readout end in the global
    // frame -- bars are not centered on the detector's coordinate origin, so barLength alone
    // isn't enough to locate a readout end.
    double GetAxisReadoutCenter() const { return AxisReadoutCenter; };

    // Get the dimension which is not z (i.e. x or y dependent on the bar
    double GetNotZ() const { 
      if (BarOrient == kXBar) return y;
      else if (BarOrient == kYBar) return x;
      else if (BarOrient == kVBar) return x;
      else if (BarOrient == kUBar) return x;
      return -9999999999;
    }

    double GetXw() const { return xw; };
    double GetYw() const { return yw; };
    double GetZw() const { return zw; };

    double GetNotZw() const {
      // FindModules swaps xw/yw for X bars so xw remains the transverse
      // (Y) width used by their bar-numbering convention.  The X-bar
      // not-Z coordinate is Y, so use that transverse width here too.
      if (BarOrient == kXBar) return xw;
      else if (BarOrient == kYBar) return xw;
      else if (BarOrient == kVBar) return xw;
      else if (BarOrient == kUBar) return xw;
      return -9999999999;
    }

    void Print() const;

    int FindBar(double x_pos, double y_pos, double z_pos);
    bool FindModules(double x_pos, double y_pos, double z_pos);

    double FindYbar(double yval);

    bool CheckBar();

    // Find if a 2D point is inside the bar
    // x_pos here denotes the other view than z
    // can be both x and y views (depending on bar type)!
    bool Contains(double x_pos, double z_pos) {

      // Get the maxium and minimum of the bar
      double zmin = GetZ()-GetZw()/2;
      double zmax = GetZ()+GetZw()/2;

      // Check the 2D point is inside in z
      if (z_pos > zmax || z_pos < zmin) return false;

      // Now check 2D point is inside in not z
      double xmin = -9999999999999;
      double xmax = 9999999999999;
      //if (BarOrient == kXBar) {
      //  xmin = GetY() - GetYw() / 2;
      //  xmax = GetY() + GetYw() / 2;
      //} else if (BarOrient == kYBar || BarOrient == kVBar || BarOrient == kUBar) {
      //  xmin = GetX() - GetXw() / 2;
      //  xmax = GetX() + GetXw() / 2;
      //}
      xmin = GetNotZ() - GetNotZw() / 2;
      xmax = GetNotZ() + GetNotZw() / 2;

      if (x_pos > xmax || x_pos < xmin) return false;

      return true;
    }

  private:
    // Plane that the bar belongs in
    int PlaneNumber;
    // The bar number in this plane
    int BarNumber;
    // Bar Width and Length in mm, Length is always the long dimension
    int BarWidth;
    int BarLength;
    // The global bar number (0-100) 
    int GlobalBarNumber;
    // All in mm units!
    // The bar start positions
    double x;
    double y;
    double z;
    // The bar widths
    double xw;
    double yw;
    double zw;
    // Which type of bar
    BarType BarOrient;
    // See GetAxisReadoutCenter()
    double AxisReadoutCenter;
};

inline bool operator==(const TMS_Bar &a, const TMS_Bar &b) {
  if (a.GetZ() == b.GetZ() && 
      a.GetNotZ() == b.GetNotZ() )
    return true;

  return false;
}


#endif
