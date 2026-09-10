#ifndef _TMS_LAYERGROUPING_H_SEEN_
#define _TMS_LAYERGROUPING_H_SEEN_

#include <vector>

#include "TMS_SpacePoint.h"

// Groups space points into z-layers -- detector planes, in practice -- so
// graph search (TMS_GraphTrackFinder) and, later, a Kalman follower agree on
// exactly where one plane ends and the next begins. Pulled out of
// TMS_GraphTrackFinder (which still owns and calls it) so a follower can
// reuse the identical layer boundaries a seed path was built against,
// instead of silently re-deriving a slightly different grouping.
class TMS_LayerGrouping {
  public:
    // spacePoints: the pool to group (order in the input doesn't matter).
    // zTolerance: two points belong to the same layer if their z values are
    //   within this distance (mm) of each other.
    // Returns groups of ORIGINAL indices into spacePoints, one group per
    // layer, in increasing z; within a layer, indices are ordered by
    // (z, then x, then y) of the point they refer to. Grouping uses
    // "first-anchor" tolerance: a layer's boundary is the z of the FIRST
    // point placed in that layer, not the immediately preceding point --
    // so a layer's total z-span can exceed zTolerance if points drift
    // gradually, by design (matches what TMS_GraphTrackFinder has always
    // done; changing this would change its output).
    static std::vector<std::vector<std::size_t>> Build(
        const std::vector<TMS_SpacePoint> &spacePoints, double zTolerance);
};

#endif
