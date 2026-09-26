#include "TMS_ClusterLinker.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

#include "TMS_SpacePointCluster.h"

namespace TMS_ClusterLinker {

namespace {

// One end of a cluster: where it is, which way it points, and when.
struct End {
  double Z = 0.0;           // z of the outermost point layer at this end
  double X = 0.0, Y = 0.0;  // position at Z (on the end's line, if it has one)
  bool HasDirection = false;
  double DXDZ = 0.0, DYDZ = 0.0;
  bool Usable = true;       // false: the end is not line-like (a shower) -- never links
  double Time = 0.0;        // mean point time over the end's layers
};

struct ClusterEnds {
  double ZMin = 0.0, ZMax = 0.0;
  End Upstream, Downstream;
};

// Describe the end made of the given point layers (each a list of point
// indices), outermost layer first.
End DescribeEnd(const std::vector<TMS_SpacePoint> &points, const std::vector<const std::vector<int> *> &layers,
                const Config &config) {
  End end;
  std::vector<int> segment;
  for (const std::vector<int> *layer : layers) segment.insert(segment.end(), layer->begin(), layer->end());
  end.Z = points[layers.front()->front()].GetZ();
  double sumT = 0.0;
  for (int idx : segment) sumT += points[idx].GetTime();
  end.Time = sumT / segment.size();

  // Straight-line least-squares fit x(z), y(z) over the segment -- needs at
  // least two distinct layers and MinDirectionPoints points.
  if (segment.size() >= config.MinDirectionPoints && layers.size() >= 2) {
    double sz = 0, sx = 0, sy = 0, szz = 0, szx = 0, szy = 0;
    const double n = segment.size();
    for (int idx : segment) {
      const double z = points[idx].GetZ() - end.Z, x = points[idx].GetX(), y = points[idx].GetY();
      sz += z; sx += x; sy += y; szz += z * z; szx += z * x; szy += z * y;
    }
    const double det = n * szz - sz * sz;
    if (det > 1e-6) {
      end.HasDirection = true;
      end.DXDZ = (n * szx - sz * sx) / det;
      end.DYDZ = (n * szy - sz * sy) / det;
      end.X = (sx - end.DXDZ * sz) / n;  // intercept at z = end.Z
      end.Y = (sy - end.DYDZ * sz) / n;
    }
    // A shower-like end (points spread as much across as along) does not link.
    const TMS_SpacePointCluster shape(points, segment);
    if (shape.HasValidPCA() && shape.GetLinearity() < config.MinEndLinearity) end.Usable = false;
  }
  if (!end.HasDirection) {
    // Too small for a direction: the outermost layer's centroid.
    const std::vector<int> &outer = *layers.front();
    for (int idx : outer) {
      end.X += points[idx].GetX();
      end.Y += points[idx].GetY();
    }
    end.X /= outer.size();
    end.Y /= outer.size();
  }
  return end;
}

ClusterEnds DescribeCluster(const std::vector<TMS_SpacePoint> &points, const std::vector<int> &cluster,
                            const Config &config) {
  // Point layers by z (NearestY points of one plane pair share a z).
  std::map<long, std::vector<int>> byZ;
  for (int idx : cluster) byZ[std::lround(points[idx].GetZ())].push_back(idx);
  std::vector<const std::vector<int> *> layers;
  for (const auto &kv : byZ) layers.push_back(&kv.second);

  ClusterEnds ends;
  ends.ZMin = points[layers.front()->front()].GetZ();
  ends.ZMax = points[layers.back()->front()].GetZ();
  const std::size_t k = std::min<std::size_t>(config.EndLayers, layers.size());
  std::vector<const std::vector<int> *> up(layers.begin(), layers.begin() + k);
  std::vector<const std::vector<int> *> down(layers.rbegin(), layers.rbegin() + k);
  ends.Upstream = DescribeEnd(points, up, config);
  ends.Downstream = DescribeEnd(points, down, config);
  return ends;
}

// Where an end's line (or, without a direction, its centroid) puts the track
// at z.
void Project(const End &end, double z, double &x, double &y) {
  x = end.X + (end.HasDirection ? end.DXDZ * (z - end.Z) : 0.0);
  y = end.Y + (end.HasDirection ? end.DYDZ * (z - end.Z) : 0.0);
}

// Score of linking A's downstream end to B's upstream end; false if they fail
// any test.
bool LinkScore(const ClusterEnds &a, const ClusterEnds &b, const Config &config, double &score) {
  // Ordered in z, B carrying on past A (side-by-side clusters overlap and fail).
  if (b.ZMin < a.ZMax - config.MaxOverlapMM || b.ZMin <= a.ZMin || b.ZMax <= a.ZMax) return false;
  const double gap = std::max(0.0, b.ZMin - a.ZMax);
  if (gap > config.MaxGapMM) return false;

  const End &ea = a.Downstream, &eb = b.Upstream;
  if (!ea.Usable || !eb.Usable) return false;
  if (!ea.HasDirection && !eb.HasDirection) return false;
  if (std::abs(ea.Time - eb.Time) > config.MaxTimeDiffNs) return false;

  const double tolY = config.MissBaseMM + config.MissPerMeterMM * gap / 1000.0;
  const double tolX = tolY * config.BendScaleX;
  score = 0.0;
  // Each end with a direction is extrapolated onto the other end.
  if (ea.HasDirection) {
    double x, y;
    Project(ea, eb.Z, x, y);
    const double dx = (x - eb.X) / tolX, dy = (y - eb.Y) / tolY;
    if (std::abs(dx) > 1.0 || std::abs(dy) > 1.0) return false;
    score += dx * dx + dy * dy;
  }
  if (eb.HasDirection) {
    double x, y;
    Project(eb, ea.Z, x, y);
    const double dx = (x - ea.X) / tolX, dy = (y - ea.Y) / tolY;
    if (std::abs(dx) > 1.0 || std::abs(dy) > 1.0) return false;
    score += dx * dx + dy * dy;
  }
  if (ea.HasDirection && eb.HasDirection) {
    const double ax = ea.DXDZ, ay = ea.DYDZ, bx = eb.DXDZ, by = eb.DYDZ;
    const double cosAngle = (ax * bx + ay * by + 1.0) / std::sqrt((ax * ax + ay * ay + 1.0) * (bx * bx + by * by + 1.0));
    const double angle = std::acos(std::min(1.0, std::max(-1.0, cosAngle)));
    if (angle > config.MaxAngleRad) return false;
    score += (angle / config.MaxAngleRad) * (angle / config.MaxAngleRad);
  }
  return true;
}

}  // namespace

Result LinkClusters(const std::vector<TMS_SpacePoint> &points, const std::vector<std::vector<int>> &clusters,
                    const Config &config) {
  Result result;
  const int n = static_cast<int>(clusters.size());
  if (n < 2) return result;
  std::vector<ClusterEnds> ends;
  ends.reserve(n);
  for (const std::vector<int> &cluster : clusters) ends.push_back(DescribeCluster(points, cluster, config));

  // Each cluster's best downstream and best upstream partner.
  const double none = std::numeric_limits<double>::infinity();
  std::vector<int> bestDown(n, -1), bestUp(n, -1);
  std::vector<double> bestDownScore(n, none), bestUpScore(n, none);
  for (int a = 0; a < n; ++a)
    for (int b = 0; b < n; ++b) {
      if (a == b) continue;
      double score = 0.0;
      if (!LinkScore(ends[a], ends[b], config, score)) continue;
      if (score < bestDownScore[a]) { bestDownScore[a] = score; bestDown[a] = b; }
      if (score < bestUpScore[b]) { bestUpScore[b] = score; bestUp[b] = a; }
    }

  // Keep mutual best links only; follow them into chains.
  std::vector<int> next(n, -1), prev(n, -1);
  for (int a = 0; a < n; ++a) {
    const int b = bestDown[a];
    if (b >= 0 && bestUp[b] == a) {
      next[a] = b;
      prev[b] = a;
      result.Links.push_back({a, b, bestDownScore[a]});
    }
  }
  for (int a = 0; a < n; ++a) {
    if (prev[a] >= 0 || next[a] < 0) continue;  // not a chain's head
    std::vector<int> chain;
    for (int c = a; c >= 0 && chain.size() <= clusters.size(); c = next[c]) chain.push_back(c);
    result.Chains.push_back(std::move(chain));
  }
  return result;
}

}  // namespace TMS_ClusterLinker
