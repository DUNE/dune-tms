#ifndef _TMS_SPACEPOINTBUILDER_H_SEEN_
#define _TMS_SPACEPOINTBUILDER_H_SEEN_

#include <vector>

#include "TMS_SpacePoint.h"

class TMS_Hit;

// Builds 3D space points by pairing X-bar and Y-bar hits from adjacent
// detector planes that land close together in time. Pulled out of TMS_Event
// (which still owns and calls it, see TMS_Event::BuildSpacePoints()) so the
// space-point-building logic can be read, tested, and changed on its own.
class TMS_SpacePointBuilder {
  public:
    // hits: an event's (or slice's) reconstructed hits. Pedestal-suppressed
    //   hits are skipped automatically -- they're noise-level and are not
    //   treated as real hits anywhere else in reconstruction.
    // timing_window: how close in time (ns) an X hit and a Y hit in an
    //   adjacent plane have to land to be paired into one space point.
    // Returns one TMS_SpacePoint per accepted X/Y pair. A single hit can end
    // up in more than one space point if several partners fall inside the
    // window -- that combinatorial ambiguity is resolved later by clustering,
    // not here.
    static std::vector<TMS_SpacePoint> Build(const std::vector<TMS_Hit> &hits,
                                              double timing_window);
};

#endif
