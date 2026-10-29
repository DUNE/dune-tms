#ifndef _TMS_MAPFIELDMODEL_H_SEEN_
#define _TMS_MAPFIELDMODEL_H_SEEN_

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "TVector3.h"

#include "TMS_FieldModel.h"

// Field lookup from a gridded field-map text file, trilinearly interpolated.
//
// File format (ASSUMED -- check against edep-sim's ArbBField reader and the
// realistic-field-mapper's Mapper.txt before relying on it): plain text, one
// grid point per line, whitespace-separated "x y z Bx By Bz"; blank lines and
// lines starting with '#' are skipped. Positions are in mm, the field in
// Tesla, both in the same frame as the positions passed to GetField()
// (the lab/global frame). The points must fill a regular (but not
// necessarily evenly spaced) x-y-z grid, in any row order.
//
// Like the field it replaces, GetField() is what the Kalman follower calls,
// and the follower then scales it by the fraction of each step spent in
// steel. So until that is settled, this model is only consistent if the map
// holds the field INSIDE the steel (see IFieldModel::GetField). Outside the
// grid the field is zero.
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
      while (std::getline(input, line)) {
        const size_t first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#') continue;
        std::istringstream row(line);
        double x, y, z, bx, by, bz;
        if (!(row >> x >> y >> z >> bx >> by >> bz))
          throw std::runtime_error("MapFieldModel: expected 'x y z Bx By Bz' in " + filename + ", got: " + line);
        positions.emplace_back(x, y, z);
        fields.emplace_back(bx, by, bz);
      }
      Build(positions, fields);
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
