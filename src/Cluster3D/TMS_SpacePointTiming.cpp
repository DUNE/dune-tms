#include "TMS_SpacePointTiming.h"

#include "TMS_Geom.h"
#include "TMS_Readout_Manager.h"

namespace TMS_SpacePointTiming {

double TransitDelayNs(const TMS_Bar &bar, double alongBarMM) {
  // Constants as in TMS_DetectorSimulation::SimulateTimingModel().
  const double kSpeedOfLightMPerNs = 0.2998;
  const double kFiberIndex = 1.5;
  const double speedInFiberMPerNs = kSpeedOfLightMPerNs / kFiberIndex;
  const double lengthMultiplier = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_WSFLengthMultiplier();

  // Distance to the readout end, as GetTrueDistanceFromReadout() in
  // TMS_DetectorSimulation.cpp: X-bars are read out at whichever side of
  // the detector the bar half sits on, Y/U/V-bars from the top.
  const double barLength = bar.GetBarLength();
  const double barCenter = bar.GetAxisReadoutCenter();
  const TMS_Geom &geom = TMS_Geom::GetInstance();
  double distanceFromReadout;
  if (bar.GetBarType() == TMS_Bar::kXBar) {
    if (alongBarMM < 0) distanceFromReadout = alongBarMM - geom.XBarNegReadoutLocation(barCenter, barLength);
    else distanceFromReadout = geom.XBarPosReadoutLocation(barCenter, barLength) - alongBarMM;
  } else {
    distanceFromReadout = geom.YBarReadoutLocation(barCenter, barLength) - alongBarMM;
  }
  const double distanceFromMiddleM = (distanceFromReadout - 0.5 * barLength) * 1e-3;
  return distanceFromMiddleM * lengthMultiplier / speedInFiberMPerNs;
}

bool CorrectedXYTimeDifference(double x, double y,
                               double xHitNotZ, double xHitZ, double tXHit,
                               double yHitNotZ, double yHitZ, double tYHit,
                               double &dtOut) {
  // The X-view hit's bar runs along x at transverse y = xHitNotZ; the
  // Y-view hit's bar runs along y at transverse x = yHitNotZ. Look each up
  // at the space point's position along it (on the bar's own center line,
  // so the geometry lookup lands inside the bar).
  TMS_Bar xBar(x, xHitNotZ, xHitZ);
  TMS_Bar yBar(yHitNotZ, y, yHitZ);
  // A failed lookup leaves the bar/plane/module indices at -1 (and the
  // orientation unset), so check those before trusting GetBarType(); then
  // require each hit to sit on the orientation its view implies.
  const auto found = [](const TMS_Bar &bar) {
    return bar.GetPlaneNumber() >= 0 && bar.GetBarNumber() >= 0 && bar.GetGlobalBarNumber() >= 0;
  };
  if (!found(xBar) || !found(yBar)) return false;
  if (xBar.GetBarType() != TMS_Bar::kXBar) return false;
  if (yBar.GetBarType() == TMS_Bar::kXBar || yBar.GetBarType() == TMS_Bar::kError) return false;
  const double tX = tXHit - TransitDelayNs(xBar, x);
  const double tY = tYHit - TransitDelayNs(yBar, y);
  dtOut = tX - tY;
  return true;
}

}  // namespace TMS_SpacePointTiming
