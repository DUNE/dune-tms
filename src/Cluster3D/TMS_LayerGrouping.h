#ifndef _TMS_LAYERGROUPING_H_SEEN_
#define _TMS_LAYERGROUPING_H_SEEN_

#include <vector>

#include "TMS_SpacePoint.h"

// Groups space points into point layers so graph search
// (TMS_GraphTrackFinder) and the Kalman follower agree on exactly where one
// layer ends and the next begins. Pulled out of TMS_GraphTrackFinder (which
// still owns and calls it) so a follower can reuse the identical layer
// boundaries a seed path was built against, instead of silently re-deriving
// a slightly different grouping.
//
// Only layers that contain points appear: the i-th group is the i-th
// POPULATED layer, so group-index differences are not layer counts. Use the
// groups' z (every point in a group shares one z) for distances.
class TMS_LayerGrouping {
  public:
    // spacePoints: the pool to group (order in the input doesn't matter).
    // If every point carries a point layer (TMS_SpacePoint::GetLayer() >= 0,
    // set by TMS_SpacePointBuilder), groups are exactly those layers. With
    // NearestY pairing a point's z is a midpoint between two planes, so
    // layers must come from the builder, not from z.
    // Otherwise (points from older files or synthetic tests) points are
    // grouped by z: two points belong to the same layer if their z values
    // are within zTolerance (mm), using "first-anchor" tolerance -- a layer's
    // boundary is the z of the FIRST point placed in it, not the immediately
    // preceding point (matches what TMS_GraphTrackFinder has always done).
    // Returns groups of ORIGINAL indices into spacePoints, one group per
    // populated layer, in increasing z; within a layer, indices are ordered
    // by (z, then x, then y) of the point they refer to.
    static std::vector<std::vector<std::size_t>> Build(
        const std::vector<TMS_SpacePoint> &spacePoints, double zTolerance);

    // The same grouping as Build(), as a per-point group index (parallel to
    // spacePoints). For tools that need "which layer is point i in".
    static std::vector<int> GroupIndexOfEachPoint(
        const std::vector<TMS_SpacePoint> &spacePoints, double zTolerance);
};

#endif
