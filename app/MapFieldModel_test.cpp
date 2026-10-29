// Tests for MapFieldModel: it must read the edep-sim ArbBField map format, reproduce RegionFieldModel when
// the map holds the same field, and average a plate-by-plate map along a step.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

#include "TMS_MapFieldModel.h"

static int failures = 0;
static void check(bool ok, const char *what) {
  if (!ok) {
    std::cerr << "FAIL: " << what << std::endl;
    ++failures;
  }
}

// Write an edep-sim map (header, then "x y z Bx By Bz Bmag", z fastest then y then x) of an x/z-dependent By.
template <typename ByFunction>
static void WriteMap(const std::string &path, double x0, double x1, double dx, double zLo, double zHi, double dz,
                     ByFunction by) {
  std::ofstream out(path);
  out << "# test map\n# units position=mm field=tesla\n";
  out << x0 << " " << -500 << " " << zLo << " " << dx << " " << 500 << " " << dz << "\n";
  for (double x = x0; x <= x1 + 1e-9; x += dx)
    for (double y : {-500.0, 0.0})  // spacing 500 matches the header's dy
      for (double z = zLo; z <= zHi + 1e-9; z += dz) {
        const double b = by(x, z);
        out << x << " " << y << " " << z << " 0 " << b << " 0 " << std::abs(b) << "\n";
      }
}

int main() {
  const std::string path = "MapFieldModel_test_map.txt";
  const RegionFieldModel region;

  // 1. A map holding RegionFieldModel's field reproduces it away from the sign border.
  WriteMap(path, -4000, 4000, 10, 11000, 18000, 100, [&](double x, double) { return region.GetField(TVector3(x, 0, 0)).Y(); });
  {
    const MapFieldModel map(path);
    check(map.NPoints() == 801 * 2 * 71, "edep-format file read with the right number of points");
    for (double x : {-3500.3, -1000.7, 0.0, 123.4, 1500.0, 2500.9, 3999.0})
      for (double z : {11500.0, 15000.2}) {
        const TVector3 p(x, -200.0, z);
        check((map.GetField(p) - region.GetField(p)).Mag() < 1e-9, "map != region away from sign border");
      }
    // On the 10 mm grid the sign flip between nodes 1860 (+1) and 1870 (-1) is linear.
    check(std::abs(map.GetField(TVector3(1865, 0, 14000)).Y()) < 1e-9, "midpoint of sign flip should be 0");
    check(map.GetField(TVector3(5000, 0, 14000)).Mag() == 0.0, "outside x range should be zero");
    check(map.GetField(TVector3(0, 0, 9000)).Mag() == 0.0, "outside z range should be zero");
    check(!map.SteelFractionApplies() && region.SteelFractionApplies(), "steel-fraction flags");
  }

  // 2. A plate-by-plate map (1 T in a plate 40 mm thick, 0 in a 40 mm gap, period 80 mm, 5 mm grid). The mean over a
  //    full period is 0.5 T wherever it starts; a 40 mm step starting 10 mm into a plate has 30 mm of steel, so 0.75 T,
  //    where a single midpoint sample (30 mm into the plate) would give 1 T.
  WriteMap(path, 0, 100, 50, 10000, 12000, 5, [](double, double z) { return std::fmod(z - 10000.0, 80.0) < 40.0 ? 1.0 : 0.0; });
  {
    const MapFieldModel map(path);
    for (double z0 : {10000.0, 10010.0, 10033.0, 10080.0, 10555.0})
      check(std::abs(map.GetMeanField(TVector3(50, 0, z0), TVector3(50, 0, z0 + 80.0)).Y() - 0.5) < 0.03,
            "mean field over a full plate+gap period should be ~0.5 T");
    const double partial = map.GetMeanField(TVector3(50, 0, 10010.0), TVector3(50, 0, 10050.0)).Y();
    check(std::abs(partial - 0.75) < 0.1, "mean over a step with 30 of 40 mm in steel should be ~0.75 T (grid-edge smearing allows ~0.69)");
    check(std::abs(map.GetField(TVector3(50, 0, 10030.0)).Y() - 1.0) < 1e-9, "the midpoint sample alone would say 1 T");
  }

  // 3. Malformed files are rejected.
  auto expectThrow = [&](const char *content, const char *what) {
    { std::ofstream out(path); out << content; }
    bool threw = false;
    try { MapFieldModel map(path); } catch (const std::runtime_error &) { threw = true; }
    check(threw, what);
  };
  expectThrow("0 0 0 10 10 10\n0 0 0 1 2 3\n", "row with 6 columns instead of 7");
  expectThrow("0 0 0 10 10 10\n0 0 0 0 1 0 1\n10 0 0 0 1 0 1\n20 0 0 0 1 0 1\n0 5 0 0 1 0 1\n", "incomplete grid");
  expectThrow("0 0 0 10 10 10\n0 0 0 0 1 0 1\n20 0 0 0 1 0 1\n", "grid spacing disagreeing with the header");
  expectThrow("# only a comment\n", "no header");
  bool threwMissing = false;
  try { MapFieldModel map("no_such_map_file.txt"); } catch (const std::runtime_error &) { threwMissing = true; }
  check(threwMissing, "missing file");

  std::remove(path.c_str());
  std::cout << "MapFieldModel_test: " << failures << " failure(s)" << std::endl;
  return failures == 0 ? 0 : 1;
}
