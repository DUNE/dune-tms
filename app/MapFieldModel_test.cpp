// Closure test for MapFieldModel: a synthetic map holding the same +-1 T y-field as RegionFieldModel must give
// the same answer as RegionFieldModel away from the sign border, and interpolate linearly across it.
#include <cmath>
#include <iostream>
#include <vector>

#include "TMS_MapFieldModel.h"

int main() {
  RegionFieldModel region;

  // Grid: x every 10 mm over +-4000 (includes the 1860 border only between nodes), y and z coarse.
  std::vector<TVector3> positions, fields;
  for (double x = -4000; x <= 4000; x += 10)
    for (double y : {-3000.0, 0.0, 1000.0})
      for (double z : {11000.0, 14000.0, 18000.0}) {
        positions.emplace_back(x, y, z);
        fields.push_back(region.GetField(TVector3(x, y, z)));
      }
  const MapFieldModel map(positions, fields);

  int failures = 0;
  auto check = [&](bool ok, const char *what) {
    if (!ok) {
      std::cerr << "FAIL: " << what << std::endl;
      ++failures;
    }
  };

  // Away from the border (and at a mix of off-node positions) the map reproduces the region model.
  for (double x : {-3500.3, -1000.7, 0.0, 123.4, 1500.0, 2500.9, 3999.0})
    for (double y : {-2000.0, 500.5})
      for (double z : {11500.0, 15000.2}) {
        const TVector3 p(x, y, z);
        check((map.GetField(p) - region.GetField(p)).Mag() < 1e-9, "map != region away from sign border");
      }

  // Across the sign flip the map interpolates linearly between its two nodes (1860 is on a node: x=1860 -> +1, 1870 -> -1).
  check(std::abs(map.GetField(TVector3(1865, 0, 14000)).Y()) < 1e-9, "midpoint of sign flip should be 0");
  // Outside the grid the field is zero.
  check(map.GetField(TVector3(5000, 0, 14000)).Mag() == 0.0, "outside x range should be zero");
  check(map.GetField(TVector3(0, 0, 9000)).Mag() == 0.0, "outside z range should be zero");

  std::cout << "MapFieldModel_test: " << map.NPoints() << " grid points, " << failures << " failure(s)" << std::endl;
  return failures == 0 ? 0 : 1;
}
