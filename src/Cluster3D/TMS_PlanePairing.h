#ifndef _TMS_PLANEPAIRING_H_SEEN_
#define _TMS_PLANEPAIRING_H_SEEN_

#include <string>
#include <vector>

// Which detector planes' hits TMS_SpacePointBuilder pairs into space points,
// and which "point layer" each pair belongs to.
//
// Naming follows the code's bar types: an X-bar plane (TMS_Bar::kXBar) has
// bars running along x, so it MEASURES y; a Y-bar plane (kYBar) measures x.
// The TMS alternates the two in its front section (65 mm pitch) and uses
// Y-Y-X triplets -- two x-measuring planes per y-measuring plane -- in the
// back (130 mm pitch).
//
// Schemes:
//  NearestY: every Y-bar (x-measuring) plane pairs with its nearest X-bar
//    (y-measuring) plane; a tie (front section) goes to the downstream one,
//    giving disjoint pairs there. One point layer per Y-bar plane, at the
//    midpoint z of the pair. Each Y-bar hit is used once; in the back, an
//    X-bar plane serves the Y-bar plane on each side. Chosen 2026-09-25
//    (reports/2026-09-25_spacepoint_pairing/): vs BothNeighbors, genuine
//    muon points per layer 2.38 -> 1.29 and nearby competing points 1.61 ->
//    0.73, at the same crossing coverage (with the fallback below).
//  BothNeighbors: the original scheme -- every X-bar plane pairs with the
//    Y-bar planes directly before and after it, and both pairs share one
//    point layer at the X-bar plane's z. Kept only for comparison until
//    NearestY is validated.
//
// Fallback pairs (NearestY only): each adjacent X-bar/Y-bar plane pair that
// isn't already a primary pair. A hit with no partner in its own primary
// pair(s) may pair across one of these instead (see TMS_SpacePointBuilder).
// They get their own point layers; their midpoints are distinct too.
namespace TMS_PlanePairing {

enum class Scheme { NearestY, BothNeighbors };

// Parses the [Recon.SpacePoints] Pairing config string; throws on anything else.
Scheme SchemeFromString(const std::string &name);

struct PlanePair {
  int XBarPlane = -1;    // plane index (TMS_Bar::GetPlaneNumber()) of the y-measuring plane
  int YBarPlane = -1;    // plane index of the x-measuring plane
  double Z = 0.0;        // z (mm) given to space points from this pair
  int Layer = -1;        // point layer index: z-ordered, shared by pairs with the same Z
  bool Fallback = false;
};

struct Table {
  std::vector<PlanePair> Pairs;  // sorted by Z
  // Per plane index: indices into Pairs of that plane's primary pairs, and
  // of the fallback pairs it may use when a hit finds no primary partner.
  std::vector<std::vector<int>> PrimaryPairsOfPlane;
  std::vector<std::vector<int>> FallbackPairsOfPlane;
  int NLayers = 0;
};

// planeZ / planeOrientation: indexed by plane, as TMS_Geom::GetPlaneZs() and
// GetPlaneOrientations() give them. Planes of any other orientation (U, V,
// unknown) are left out of every pair, as the builder has always done.
Table Build(const std::vector<double> &planeZ, const std::vector<int> &planeOrientation, Scheme scheme);

// The table for the currently loaded geometry and the configured scheme.
Table BuildFromGeometry();

}  // namespace TMS_PlanePairing

#endif
