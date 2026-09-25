#ifndef _TMS_SPACEPOINTBUILDER_H_SEEN_
#define _TMS_SPACEPOINTBUILDER_H_SEEN_

#include <vector>

#include "TMS_PlanePairing.h"
#include "TMS_SpacePoint.h"

class TMS_Hit;

// Builds 3D space points by pairing X-bar hits (which measure y) with Y-bar
// hits (which measure x) from the planes a TMS_PlanePairing::Table pairs,
// when they land close together in time. Pulled out of TMS_Event (which
// still owns and calls it, see TMS_Event::BuildSpacePoints()) so the
// space-point-building logic can be read, tested, and changed on its own.
class TMS_SpacePointBuilder {
  public:
    // hits: an event's (or slice's) reconstructed hits. Pedestal-suppressed
    //   hits are skipped automatically -- they're noise-level and are not
    //   treated as real hits anywhere else in reconstruction.
    // timing_window: how close in time (ns) the two hits have to land to be
    //   paired into one space point.
    // pairing: which planes pair with which, each pair's z and point layer
    //   (TMS_PlanePairing::BuildFromGeometry() for the loaded geometry).
    // use_fallback: after the primary pairs, let a hit that found no partner
    //   in any of its primary pairs pair across its fallback pairs instead
    //   (see TMS_PlanePairing.h). Ignored if the table has no fallback pairs.
    // Returns one TMS_SpacePoint per accepted pair of hits, with z and layer
    // from its plane pair. A single hit can end up in more than one space
    // point if several partners fall inside the window -- that combinatorial
    // ambiguity is resolved later by clustering, not here.
    static std::vector<TMS_SpacePoint> Build(const std::vector<TMS_Hit> &hits,
                                              double timing_window,
                                              const TMS_PlanePairing::Table &pairing,
                                              bool use_fallback);
};

#endif
