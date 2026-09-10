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

  std::vector<std::vector<std::size_t>> layers;
  for (std::size_t inputIndex : order) {
    const TMS_SpacePoint &point = spacePoints[inputIndex];
    if (layers.empty() ||
        std::abs(point.GetZ() -
                 spacePoints[layers.back().front()].GetZ()) > zTolerance) {
      layers.push_back(std::vector<std::size_t>());
    }
    layers.back().push_back(inputIndex);
  }
  return layers;
}
