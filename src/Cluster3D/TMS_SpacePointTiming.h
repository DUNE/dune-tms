#ifndef _TMS_SPACEPOINTTIMING_H_SEEN_
#define _TMS_SPACEPOINTTIMING_H_SEEN_

// Light-transit correction for a space point's two component hit times.
//
// Each hit's time includes the time its light spent travelling along the
// wavelength-shifting fiber to the bar's readout end (see
// TMS_DetectorSimulation::SimulateTimingModel()). That delay depends on
// where along the bar the particle crossed -- which a single hit does not
// know, but a space point does: an X-bar hit's along-bar position is the
// space point's x, and a Y-bar hit's is its y. Removing the delay from both
// hits leaves tX - tY centred on zero for a genuine (one-particle) space
// point, while a ghost pairing two different particles' hits keeps the real
// time difference between them.
//
// Measured 2026-09-24 (reports/2026-09-24_reco_hit_lookaside/, 15 files):
// uncorrected, genuine and ghost points have the same tX - tY spread (MAD
// ~14.6 ns) -- the transit delay dominates. Corrected (empirical fit),
// genuine MAD 5.1 ns vs 20.5 ns for ghosts pairing different interactions.

#include "TMS_Bar.h"

namespace TMS_SpacePointTiming {

// Transit delay (ns) that SimulateTimingModel() adds for light from
// along-bar coordinate alongBarMM (global x for X-bars, global y for
// Y/U/V-bars) travelling the short way to the bar's readout end, relative
// to light from the bar's centre. Same formula and constants as the
// simulation: (distance to readout - half bar length) * WSF length
// multiplier / speed of light in fiber.
double TransitDelayNs(const TMS_Bar &bar, double alongBarMM);

// Transit-corrected (tX - tY) for a space point at (x, y) built from an
// X-view hit (bar centre transverse coordinate xHitNotZ, z xHitZ, time
// tXHit) and a Y-view hit (yHitNotZ, yHitZ, tYHit). Each hit's bar is
// looked up in the loaded geometry at the space point's position along it.
// Returns false (dtOut untouched) if either bar can't be found, e.g. an x
// that falls in the X-bars' central gap.
bool CorrectedXYTimeDifference(double x, double y,
                               double xHitNotZ, double xHitZ, double tXHit,
                               double yHitNotZ, double yHitZ, double tYHit,
                               double &dtOut);

}  // namespace TMS_SpacePointTiming

#endif
