#include "TMS_LayerGrouping.h"

#include <algorithm>
#include <cmath>

std::vector<std::vector<std::size_t>> TMS_LayerGrouping::Build(
    const std::vector<TMS_SpacePoint> &spacePoints, double zTolerance) {
  std::vector<std::size_t> order(spacePoints.size());
  for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&spacePoints](std::size_t a,
                                                        std::size_t b) {
    if (spacePoints[a].GetZ() != spacePoints[b].GetZ())
      return spacePoints[a].GetZ() < spacePoints[b].GetZ();
    if (spacePoints[a].GetX() != spacePoints[b].GetX())
      return spacePoints[a].GetX() < spacePoints[b].GetX();
    return spacePoints[a].GetY() < spacePoints[b].GetY();
  });

  const bool useLayers = std::all_of(spacePoints.begin(), spacePoints.end(),
                                     [](const TMS_SpacePoint &p) { return p.GetLayer() >= 0; });

  std::vector<std::vector<std::size_t>> layers;
  for (std::size_t inputIndex : order) {
    const TMS_SpacePoint &point = spacePoints[inputIndex];
    bool newLayer = layers.empty();
    if (!newLayer) {
      const TMS_SpacePoint &anchor = spacePoints[layers.back().front()];
      // Layer indices are z-ordered (one z per layer), so after sorting by z
      // a change of layer index always starts a new group.
      newLayer = useLayers ? point.GetLayer() != anchor.GetLayer()
                           : std::abs(point.GetZ() - anchor.GetZ()) > zTolerance;
    }
    if (newLayer) layers.push_back(std::vector<std::size_t>());
    layers.back().push_back(inputIndex);
  }
  return layers;
}

std::vector<int> TMS_LayerGrouping::GroupIndexOfEachPoint(
    const std::vector<TMS_SpacePoint> &spacePoints, double zTolerance) {
  std::vector<int> group(spacePoints.size(), -1);
  const std::vector<std::vector<std::size_t>> layers = Build(spacePoints, zTolerance);
  for (std::size_t g = 0; g < layers.size(); ++g)
    for (std::size_t index : layers[g]) group[index] = static_cast<int>(g);
  return group;
}
