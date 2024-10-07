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

    // Enum for the x, y, U, V bar orientation
    enum BarType { kXBar, kYBar, kUBar, kVBar, kError };
    std::string BarType_ToString(BarType bar) const;

    // Getter functions
    int GetBarNumber() const { return BarNumber; }; // Number of bar (start counting from -3520mm onwards)
    int GetBarWidth() const { return BarWidth; }; // Get Bar width in mm
    int GetBarLength() const { return BarLength; }; // Get Bar length in mm
    int GetPlaneNumber() const { return PlaneNumber; }; // Plane or layer number through the detector starting at smallest z
    int GetGlobalBarNumber() const { return GlobalBarNumber; }; // Number of hit Scintillator Module (4 modules)

    BarType GetBarType() const { return BarOrient; };

    double GetX() const { return x; };
    double GetY() const { return y; };
    double GetZ() const { return z; };

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
      if (BarOrient == kXBar) return yw;
      else if (BarOrient == kYBar) return xw;
      else if (BarOrient == kVBar) return xw;
      else if (BarOrient == kUBar) return xw;
      return -9999999999;
    }

    void Print() const;

    int FindBar(double x, double y, double z);
    bool FindModules(double x, double y, double z);

    double FindYbar(double yval);

    bool CheckBar();

    // Find if a 2D point is inside the bar
    // x here denotes the other view than z
    // can be both x and y views (depending on bar type)!
    bool Contains(double x, double z) {

      // Get the maxium and minimum of the bar
      double zmin = GetZ()-GetZw()/2;
      double zmax = GetZ()+GetZw()/2;

      // Check the 2D point is inside in z
      if (z > zmax || z < zmin) return false;

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

      if (x > xmax || x < xmin) return false;

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
};

inline bool operator==(const TMS_Bar &a, const TMS_Bar &b) {
  if (a.GetZ() == b.GetZ() && 
      a.GetNotZ() == b.GetNotZ() )
    return true;

  return false;
}


#endif
