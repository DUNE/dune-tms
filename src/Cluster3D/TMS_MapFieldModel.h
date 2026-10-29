#ifndef _TMS_MAPFIELDMODEL_H_SEEN_
#define _TMS_MAPFIELDMODEL_H_SEEN_

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "TVector3.h"

#include "TMS_FieldModel.h"

// Field lookup from a gridded field-map text file, trilinearly interpolated.
//
// File format: the edep-sim ArbBField text map, as written by the realistic
// field mapper (mapper_integration/tms_mapper, PR #350). '#' comment lines and
// blank lines are skipped. The first other line is the grid header
// "offset_x offset_y offset_z dx dy dz" (mm); every later line is a grid point
// "x y z Bx By Bz Bmag" (mm, Tesla, global frame -- the same frame as the
// positions passed to GetField()). The points must fill a complete x-y-z grid;
// row order is not relied on (edep-sim itself wants z fastest, then y, then x).
//
// That map is the TRUE field at every point: zero in the gaps between steel
// plates, nonzero only inside them. So SteelFractionApplies() is false and the
// follower uses the value as is, instead of scaling it by the steel fraction
// of the step (which would apply the gap suppression twice). Outside the grid
// the field is zero.
class MapFieldModel : public IFieldModel {
  public:
    explicit MapFieldModel(const std::string &filename) { Load(filename); }

    // Build from points already in memory (x,y,z in mm; b in Tesla). Used by tests.
    MapFieldModel(const std::vector<TVector3> &positions, const std::vector<TVector3> &fields) {
      Build(positions, fields);
    }

    TVector3 GetField(const TVector3 &position) const override {
      size_t ix, iy, iz;
      double fx, fy, fz;
      if (!Locate(fX, position.X(), ix, fx) || !Locate(fY, position.Y(), iy, fy) ||
          !Locate(fZ, position.Z(), iz, fz))
        return TVector3(0.0, 0.0, 0.0);

      TVector3 result(0.0, 0.0, 0.0);
      for (int dx = 0; dx < 2; ++dx)
        for (int dy = 0; dy < 2; ++dy)
          for (int dz = 0; dz < 2; ++dz) {
            const double weight = (dx ? fx : 1.0 - fx) * (dy ? fy : 1.0 - fy) * (dz ? fz : 1.0 - fz);
            if (weight == 0.0) continue;
            result += weight * fB[Index(ix + dx, iy + dy, iz + dz)];
          }
      return result;
    }

    // A step can be longer than the steel plates and the gaps between them, so average the map along it
    // (midpoint rule, at least two samples per smallest grid spacing along the step).
    TVector3 GetMeanField(const TVector3 &start, const TVector3 &end) const override {
      const double length = (end - start).Mag();
      const double cell = std::min({MinSpacing(fX), MinSpacing(fY), MinSpacing(fZ)});
      const int n = std::max(1, static_cast<int>(std::ceil(2.0 * length / cell)));
      TVector3 sum(0.0, 0.0, 0.0);
      for (int i = 0; i < n; ++i) sum += GetField(start + (i + 0.5) / n * (end - start));
      return (1.0 / n) * sum;
    }

    bool SteelFractionApplies() const override { return false; }

    size_t NPoints() const { return fB.size(); }
    // Axis ranges (mm), for printing and sanity checks.
    double MinX() const { return fX.front(); }
    double MaxX() const { return fX.back(); }
    double MinY() const { return fY.front(); }
    double MaxY() const { return fY.back(); }
    double MinZ() const { return fZ.front(); }
    double MaxZ() const { return fZ.back(); }

  private:
    std::vector<double> fX, fY, fZ;  // sorted unique grid coordinates
    std::vector<TVector3> fB;        // x slowest, z fastest

    size_t Index(size_t ix, size_t iy, size_t iz) const { return (ix * fY.size() + iy) * fZ.size() + iz; }

    // Cell index and the fractional position in it along one axis; false if outside the axis range.
    static bool Locate(const std::vector<double> &axis, double value, size_t &index, double &fraction) {
      if (value < axis.front() || value > axis.back()) return false;
      if (axis.size() == 1) {
        index = 0;
        fraction = 0.0;
        return true;
      }
      const auto upper = std::upper_bound(axis.begin(), axis.end(), value);
      index = (upper == axis.end()) ? axis.size() - 2 : static_cast<size_t>(upper - axis.begin()) - 1;
      fraction = (value - axis[index]) / (axis[index + 1] - axis[index]);
      return true;
    }

    // Smallest gap between neighbouring nodes on an axis (infinite for a single node).
    static double MinSpacing(const std::vector<double> &axis) {
      double smallest = std::numeric_limits<double>::infinity();
      for (size_t i = 1; i < axis.size(); ++i) smallest = std::min(smallest, axis[i] - axis[i - 1]);
      return smallest;
    }

    static std::vector<double> UniqueSorted(std::vector<double> values) {
      std::sort(values.begin(), values.end());
      values.erase(std::unique(values.begin(), values.end()), values.end());
      return values;
    }

    void Load(const std::string &filename) {
      std::ifstream input(filename);
      if (!input.is_open()) throw std::runtime_error("MapFieldModel: could not open field map " + filename);
      std::vector<TVector3> positions, fields;
      std::string line;
      bool haveHeader = false;
      double spacing[3] = {0.0, 0.0, 0.0};
      while (std::getline(input, line)) {
        const size_t first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#') continue;
        std::istringstream row(line);
        if (!haveHeader) {
          double ox, oy, oz;
          if (!(row >> ox >> oy >> oz >> spacing[0] >> spacing[1] >> spacing[2]) || spacing[0] <= 0 ||
              spacing[1] <= 0 || spacing[2] <= 0)
            throw std::runtime_error("MapFieldModel: first line of " + filename +
                                     " must be the header 'offset_x offset_y offset_z dx dy dz', got: " + line);
          haveHeader = true;
          continue;
        }
        double x, y, z, bx, by, bz, bmag;
        if (!(row >> x >> y >> z >> bx >> by >> bz >> bmag))
          throw std::runtime_error("MapFieldModel: expected 'x y z Bx By Bz Bmag' in " + filename + ", got: " + line);
        positions.emplace_back(x, y, z);
        fields.emplace_back(bx, by, bz);
      }
      if (!haveHeader) throw std::runtime_error("MapFieldModel: no grid header found in " + filename);
      Build(positions, fields);
      // The declared spacing must be the grid's own, or the file is not what its header says.
      const std::vector<double> *axes[3] = {&fX, &fY, &fZ};
      for (int i = 0; i < 3; ++i)
        if (axes[i]->size() > 1 && std::abs((*axes[i])[1] - (*axes[i])[0] - spacing[i]) > 1e-6 * spacing[i])
          throw std::runtime_error("MapFieldModel: grid spacing in " + filename + " disagrees with its header");
    }

    void Build(const std::vector<TVector3> &positions, const std::vector<TVector3> &fields) {
      if (positions.empty() || positions.size() != fields.size())
        throw std::runtime_error("MapFieldModel: empty field map, or positions and fields differ in length");
      std::vector<double> xs, ys, zs;
      for (const TVector3 &p : positions) {
        xs.push_back(p.X());
        ys.push_back(p.Y());
        zs.push_back(p.Z());
      }
      fX = UniqueSorted(xs);
      fY = UniqueSorted(ys);
      fZ = UniqueSorted(zs);
      if (fX.size() * fY.size() * fZ.size() != positions.size())
        throw std::runtime_error("MapFieldModel: field map points do not fill a regular x-y-z grid");

      fB.assign(positions.size(), TVector3(0.0, 0.0, 0.0));
      std::vector<bool> filled(positions.size(), false);
      for (size_t i = 0; i < positions.size(); ++i) {
        const size_t ix = std::lower_bound(fX.begin(), fX.end(), positions[i].X()) - fX.begin();
        const size_t iy = std::lower_bound(fY.begin(), fY.end(), positions[i].Y()) - fY.begin();
        const size_t iz = std::lower_bound(fZ.begin(), fZ.end(), positions[i].Z()) - fZ.begin();
        const size_t index = Index(ix, iy, iz);
        if (filled[index]) throw std::runtime_error("MapFieldModel: duplicate grid point in field map");
        filled[index] = true;
        fB[index] = fields[i];
      }
    }
};

#endif
