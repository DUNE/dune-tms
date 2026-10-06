#include "TMS_PlanePairing.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>

#include "TMS_Geom.h"
#include "TMS_Manager.h"

namespace TMS_PlanePairing {

namespace {
// Orientation numbers as TMS_Geom::GetPlaneOrientations() / TMS_Bar::GetBarTypeNumber() give them.
const int kXBar = 0;  // measures y
const int kYBar = 1;  // measures x

// Two distances within this (mm) count as a tie when choosing the nearest plane.
const double kTieToleranceMM = 1.0;
}  // namespace

Scheme SchemeFromString(const std::string &name) {
  if (name == "NearestY") return Scheme::NearestY;
  if (name == "BothNeighbors") return Scheme::BothNeighbors;
  throw std::invalid_argument("TMS_PlanePairing: unknown [Recon.SpacePoints] Pairing '" + name +
                              "' (expected NearestY or BothNeighbors)");
}

Table Build(const std::vector<double> &planeZ, const std::vector<int> &planeOrientation, Scheme scheme) {
  Table table;
  const int nPlanes = static_cast<int>(std::min(planeZ.size(), planeOrientation.size()));
  const auto isXBar = [&](int p) { return p >= 0 && p < nPlanes && planeOrientation[p] == kXBar; };
  const auto isYBar = [&](int p) { return p >= 0 && p < nPlanes && planeOrientation[p] == kYBar; };

  std::set<std::pair<int, int>> primary;  // (XBarPlane, YBarPlane)
  if (scheme == Scheme::NearestY) {
    for (int y = 0; y < nPlanes; ++y) {
      if (!isYBar(y)) continue;
      // Nearest X-bar plane on each side.
      int below = -1, above = -1;
      for (int p = y - 1; p >= 0; --p)
        if (isXBar(p)) { below = p; break; }
      for (int p = y + 1; p < nPlanes; ++p)
        if (isXBar(p)) { above = p; break; }
      int x = -1;
      if (below < 0) x = above;
      else if (above < 0) x = below;
      else {
        const double dBelow = planeZ[y] - planeZ[below];
        const double dAbove = planeZ[above] - planeZ[y];
        // Tie (the alternating front section): downstream, so pairs there are disjoint.
        x = (std::abs(dBelow - dAbove) <= kTieToleranceMM || dAbove < dBelow) ? above : below;
      }
      if (x < 0) continue;
      primary.insert({x, y});
      PlanePair pair;
      pair.XBarPlane = x;
      pair.YBarPlane = y;
      pair.Z = 0.5 * (planeZ[x] + planeZ[y]);
      table.Pairs.push_back(pair);
    }
    // Fallback: every directly adjacent X-bar/Y-bar plane pair not already primary.
    for (int y = 0; y < nPlanes; ++y) {
      if (!isYBar(y)) continue;
      for (int x : {y - 1, y + 1}) {
        if (!isXBar(x) || primary.count({x, y})) continue;
        PlanePair pair;
        pair.XBarPlane = x;
        pair.YBarPlane = y;
        pair.Z = 0.5 * (planeZ[x] + planeZ[y]);
        pair.Fallback = true;
        table.Pairs.push_back(pair);
      }
    }
  } else {
    // BothNeighbors: each X-bar plane with the Y-bar planes directly before
    // and after it, both at the X-bar plane's own z.
    for (int x = 0; x < nPlanes; ++x) {
      if (!isXBar(x)) continue;
      for (int y : {x - 1, x + 1}) {
        if (!isYBar(y)) continue;
        PlanePair pair;
        pair.XBarPlane = x;
        pair.YBarPlane = y;
        pair.Z = planeZ[x];
        table.Pairs.push_back(pair);
      }
    }
  }

  std::sort(table.Pairs.begin(), table.Pairs.end(), [](const PlanePair &a, const PlanePair &b) {
    if (a.Z != b.Z) return a.Z < b.Z;
    return a.YBarPlane < b.YBarPlane;
  });
  // Point layers: one per distinct z.
  int layer = -1;
  double previousZ = 0.0;
  for (std::size_t i = 0; i < table.Pairs.size(); ++i) {
    if (i == 0 || std::abs(table.Pairs[i].Z - previousZ) > 1e-3) {
      ++layer;
      previousZ = table.Pairs[i].Z;
    }
    table.Pairs[i].Layer = layer;
  }
  table.NLayers = layer + 1;

  table.PrimaryPairsOfPlane.assign(nPlanes, {});
  table.FallbackPairsOfPlane.assign(nPlanes, {});
  for (std::size_t i = 0; i < table.Pairs.size(); ++i) {
    const PlanePair &pair = table.Pairs[i];
    auto &target = pair.Fallback ? table.FallbackPairsOfPlane : table.PrimaryPairsOfPlane;
    target[pair.XBarPlane].push_back(static_cast<int>(i));
    target[pair.YBarPlane].push_back(static_cast<int>(i));
  }
  return table;
}

Table BuildFromGeometry() {
  TMS_Geom &geom = TMS_Geom::GetInstance();
  const Scheme scheme = SchemeFromString(TMS_Manager::GetInstance().Get_RECO_SPACEPOINTS_Pairing());
  return Build(geom.GetPlaneZs(), geom.GetPlaneOrientations(), scheme);
}

}  // namespace TMS_PlanePairing
