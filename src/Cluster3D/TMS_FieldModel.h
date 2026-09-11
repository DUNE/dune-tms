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

    // Field at a lab-frame position, in Tesla.
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
// The per-region magnitude has NO source-of-truth value anywhere in this
// repo (checked: no field map, and legacy's own deflection term is computed
// but never applied -- see TMS_Kalman.cpp:252-253). RegionFieldModelConfig
// below defaults to a clearly-flagged placeholder; treat any momentum this
// produces as "plausibly shaped," not calibrated, until a real value is
// confirmed against TMS design documentation.
struct RegionFieldModelConfig {
  // Tesla. Central region (|x| < region_2_and_3_border) and outer regions
  // share a magnitude but opposite sign, matching legacy's
  // TMS_Kalman.cpp:225-230 region split.
  double CentralRegionFieldY = 0.15;  // PLACEHOLDER -- not a confirmed TMS design value
  double OuterRegionFieldY = -0.15;   // PLACEHOLDER -- not a confirmed TMS design value
};

class RegionFieldModel : public IFieldModel {
  public:
    using Config = RegionFieldModelConfig;

    explicit RegionFieldModel(const Config &config = Config()) : fConfig(config) {}

    TVector3 GetField(const TVector3 &position) const override {
      const double fieldY =
          (std::abs(position.X()) <= TMS_Const::TMS_Magnetic_region_2_and_3_border)
              ? fConfig.CentralRegionFieldY
              : fConfig.OuterRegionFieldY;
      return TVector3(0.0, fieldY, 0.0);
    }

  private:
    Config fConfig;
};

#endif
