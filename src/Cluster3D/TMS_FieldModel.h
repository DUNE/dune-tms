#ifndef _TMS_FIELDMODEL_H_SEEN_
#define _TMS_FIELDMODEL_H_SEEN_

#include <cmath>

#include "TVector3.h"

#include "TMS_Constants.h"

// A swappable magnetic-field lookup for the Kalman follower's swimmer. The
// follower only ever calls GetField() -- swapping RegionFieldModel for a
// real interpolated field map later is a one-line change at construction,
// not a redesign.
class IFieldModel {
  public:
    virtual ~IFieldModel() = default;

    // Field INSIDE THE MAGNETIZED STEEL at a lab-frame position's x/y, in
    // Tesla. The TMS field exists in the steel plates only; the Kalman
    // follower scales this by the fraction of each propagation step spent in
    // steel (see SteelFraction() in TMS_KalmanFollower.cpp).
    virtual TVector3 GetField(const TVector3 &position) const = 0;
};

// No field at all -- lets Predict()/ResolveLayer() be exercised and
// debugged (synthetic tests, straight-track sanity checks) independently
// of whether the field-bending math is right.
class ZeroFieldModel : public IFieldModel {
  public:
    TVector3 GetField(const TVector3 & /*position*/) const override {
      return TVector3(0.0, 0.0, 0.0);
    }
};

// v1 field model: the same 3-zone piecewise-constant region split TMS_Kalman
// already uses (TMS_Const::TMS_Magnetic_region_1_and_2_border /
// _2_and_3_border), field along y only (matching the follower's assumption
// that only dx/dz bends -- see TMS_KalmanFollower.cpp), sign flipping
// between the central and outer regions same as legacy's SignSelection().
//
// The per-region magnitude was a total placeholder (0.15T) until confirmed
// two independent ways on 2026-09-10: (1) the actual production GDML used
// for this MC (nd_hall_with_lar_tms_sand_drift1_v2026.03.06.txt) defines a
// real field via the standard edep-sim <auxiliary auxtype="BField"> GDML
// convention on TMS's own steel volumes -- thinvolTMS/thickvolTMS/
// doublevolTMS at (0,+1.0T,0), thinvol2TMS/thickvol2TMS/doublevol2TMS at
// (0,-1.0T,0) -- i.e. exactly this two-sided-sign-flip structure, at 1.0
// Tesla, field along y (bends only dx/dz, matching this file's model).
// (2) Verified empirically too: the flagship muon's own true (both-sides-
// verified) trajectory shows dx/dz genuinely, monotonically steepening
// along the track (~-0.14 near the entrance to ~-0.56 near the exit) while
// dy/dz stays flat at zero -- the real signature of a systematic Y-field
// bend, not multiple-scattering noise (which would random-walk, not trend
// monotonically in one axis only). Back-solving this file's own curvature
// formula with that muon's known true momentum (~1.49 GeV) gives an
// implied field of ~0.89T, matching the GDML's 1.0T within the accuracy of
// a crude segment-averaged slope estimate over real bar-pitch-quantized
// data. Still worth confirming against a second production geometry (the
// user's caution: newer geometries may use a different, progressively-
// varying field structure across thickness regions -- this hasn't been
// checked yet) before treating 1.0T as universally correct.
//
// Correction 2026-09-25: the field is in the STEEL only (as the GDML says),
// not along the whole path -- until then the follower applied 1.0 T over
// every mm, bending ~2-4x too much and fitting momentum 2.2-2.7x too high.
// Measured from G4 truth (1,507 muons): 1.0-1.1 T per unit steel at every
// |x|, sign flip near |x| = 1860 mm (see SignBorderMM). The "~0.89T" one-muon
// estimate above averaged over steel AND gaps with a crude slope estimate.
struct RegionFieldModelConfig {
  // Tesla. Central region (|x| < region_2_and_3_border) and outer regions
  // share a magnitude but opposite sign, matching legacy's
  // TMS_Kalman.cpp:225-230 region split AND the GDML's thinvolTMS (+) vs.
  // thinvol2TMS (-) volume pairing.
  double CentralRegionFieldY = 1.0;
  double OuterRegionFieldY = -1.0;
  // |x| (mm) where the field flips sign. Measured 2026-09-25 from G4 truth
  // (reports/2026-09-25_phase2_hitfit/true_field_vs_x.out, 1,507 muons): the
  // field per unit steel is +1.0-1.1 T for |x| < 1750 mm and -1.1 T beyond
  // 2000 mm, i.e. the 1860 mm border of TMS_Constants.h's IS_PDR (7 m wide)
  // branch. The build uses the non-PDR branch (1467.5 mm), which does NOT
  // match this production's geometry -- hence a setting here rather than
  // TMS_Const::TMS_Magnetic_region_2_and_3_border.
  double SignBorderMM = 1860.0;
};

class RegionFieldModel : public IFieldModel {
  public:
    using Config = RegionFieldModelConfig;

    explicit RegionFieldModel(const Config &config = Config()) : fConfig(config) {}

    TVector3 GetField(const TVector3 &position) const override {
      const double fieldY =
          (std::abs(position.X()) <= fConfig.SignBorderMM)
              ? fConfig.CentralRegionFieldY
              : fConfig.OuterRegionFieldY;
      return TVector3(0.0, fieldY, 0.0);
    }

  private:
    Config fConfig;
};

#endif
