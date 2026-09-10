#include "TMS_GraphTrackFinder.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace TMS_GraphTrackFinder {
namespace {

struct Node {
  std::size_t InputIndex;
  std::size_t Layer;
  std::size_t HitMultiplicity;
};

struct Link {
  std::size_t NodeIndex;
  double Cost;
};

struct Hypothesis {
  std::vector<std::size_t> Nodes; // always in increasing-z order
  double Score;
};

bool BetterHypothesis(const Hypothesis &a, const Hypothesis &b) {
  if (a.Score != b.Score) return a.Score < b.Score;
  return a.Nodes.size() > b.Nodes.size();
}

void KeepBest(std::vector<Hypothesis> &items, std::size_t limit,
              Diagnostics &stats) {
  std::sort(items.begin(), items.end(), BetterHypothesis);
  if (items.size() > limit) {
    stats.HypothesesPruned += items.size() - limit;
    items.resize(limit);
  }
  stats.MaxLiveHypotheses = std::max(stats.MaxLiveHypotheses, items.size());
}

double NodeOverlapFraction(const Hypothesis &a, const Hypothesis &b) {
  if (a.Nodes.empty() || b.Nodes.empty()) return 0.0;
  std::unordered_set<std::size_t> nodes(a.Nodes.begin(), a.Nodes.end());
  std::size_t overlap = 0;
  for (std::size_t node : b.Nodes)
    if (nodes.count(node)) ++overlap;
  return static_cast<double>(overlap) /
         static_cast<double>(std::min(a.Nodes.size(), b.Nodes.size()));
}

void KeepDiverseSeeds(std::vector<Hypothesis> &seeds, const Config &config,
                      Diagnostics &stats) {
  std::sort(seeds.begin(), seeds.end(), BetterHypothesis);
  std::vector<Hypothesis> kept;
  for (const Hypothesis &seed : seeds) {
    bool duplicate = false;
    for (const Hypothesis &other : kept) {
      if (NodeOverlapFraction(seed, other) >= config.SeedOverlapFraction) {
        duplicate = true;
        break;
      }
    }
    if (!duplicate) kept.push_back(seed);
    else ++stats.HypothesesPruned;
    if (kept.size() >= config.MaxSeeds) break;
  }
  seeds.swap(kept);
}

bool SharesNativeHit(const TMS_SpacePoint &a, const TMS_SpacePoint &b) {
  return (a.GetXHitIndex() >= 0 &&
          a.GetXHitIndex() == b.GetXHitIndex()) ||
         (a.GetYHitIndex() >= 0 &&
          a.GetYHitIndex() == b.GetYHitIndex());
}

bool ConflictsWithPath(const std::vector<std::size_t> &path,
                       std::size_t candidate) {
  for (std::size_t nodeIndex : path) {
    if (nodeIndex == candidate) return true;
  }
  return false;
}

double TransitionCost(std::size_t first, std::size_t middle,
                      std::size_t last, const std::vector<Node> &nodes,
                      const std::vector<TMS_SpacePoint> &points,
                      const Config &config) {
  const TMS_SpacePoint &a = points[nodes[first].InputIndex];
  const TMS_SpacePoint &b = points[nodes[middle].InputIndex];
  const TMS_SpacePoint &c = points[nodes[last].InputIndex];
  const double dz1 = b.GetZ() - a.GetZ();
  const double dz2 = c.GetZ() - b.GetZ();
  const double raw_dxdz1 = (b.GetX() - a.GetX()) / dz1;
  const double raw_dydz1 = (b.GetY() - a.GetY()) / dz1;
  const double raw_dxdz2 = (c.GetX() - b.GetX()) / dz2;
  const double raw_dydz2 = (c.GetY() - b.GetY()) / dz2;

  const double raw_dslope_x = raw_dxdz2 - raw_dxdz1;
  const double raw_dslope_y = raw_dydz2 - raw_dydz1;
  const double quant_slope_x =
      config.PositionQuantizationX * (1.0 / std::abs(dz1) + 1.0 / std::abs(dz2));
  const double quant_slope_y =
      config.PositionQuantizationY * (1.0 / std::abs(dz1) + 1.0 / std::abs(dz2));

  const double dslope_x = std::abs(raw_dslope_x) > quant_slope_x
      ? raw_dslope_x - std::copysign(quant_slope_x, raw_dslope_x)
      : 0.0;
  const double dslope_y = std::abs(raw_dslope_y) > quant_slope_y
      ? raw_dslope_y - std::copysign(quant_slope_y, raw_dslope_y)
      : 0.0;

  return config.KinkXPenalty * dslope_x * dslope_x +
         config.KinkYPenalty * dslope_y * dslope_y;
}

double ProjectionCost(std::size_t first, std::size_t middle,
                      std::size_t last, std::size_t candidate,
                      const std::vector<Node> &nodes,
                      const std::vector<TMS_SpacePoint> &points,
                      const Config &config) {
  const TMS_SpacePoint &a = points[nodes[first].InputIndex];
  const TMS_SpacePoint &b = points[nodes[middle].InputIndex];
  const TMS_SpacePoint &c = points[nodes[last].InputIndex];
  const TMS_SpacePoint &d = points[nodes[candidate].InputIndex];

  const double dz1 = b.GetZ() - a.GetZ();
  const double dz2 = c.GetZ() - b.GetZ();
  const double dz3 = d.GetZ() - c.GetZ();
  if (std::abs(dz1) < 1e-9 || std::abs(dz2) < 1e-9 || std::abs(dz3) < 1e-9) {
    return std::numeric_limits<double>::infinity();
  }

  const double slope_x1 = (b.GetX() - a.GetX()) / dz1;
  const double slope_x2 = (c.GetX() - b.GetX()) / dz2;
  const double slope_y1 = (b.GetY() - a.GetY()) / dz1;
  const double slope_y2 = (c.GetY() - b.GetY()) / dz2;

  const double raw_dslope_x = slope_x2 - slope_x1;
  const double raw_dslope_y = slope_y2 - slope_y1;
  const double quant_slope_x =
      config.PositionQuantizationX * (1.0 / std::abs(dz1) + 1.0 / std::abs(dz2));
  const double quant_slope_y =
      config.PositionQuantizationY * (1.0 / std::abs(dz1) + 1.0 / std::abs(dz2));
  const double dslope_x = std::abs(raw_dslope_x) > quant_slope_x
      ? raw_dslope_x - std::copysign(quant_slope_x, raw_dslope_x)
      : 0.0;
  const double dslope_y = std::abs(raw_dslope_y) > quant_slope_y
      ? raw_dslope_y - std::copysign(quant_slope_y, raw_dslope_y)
      : 0.0;

  const double mid_dz = 0.5 * (dz1 + dz2);
  const double curv_x = dslope_x / mid_dz;
  const double curv_y = dslope_y / mid_dz;
  const double projected_x = c.GetX() + (slope_x2 + curv_x * dz3) * dz3;
  const double projected_y = c.GetY() + (slope_y2 + curv_y * dz3) * dz3;
  const double dx = d.GetX() - projected_x;
  const double dy = d.GetY() - projected_y;
  return config.CurvatureProjectionPenaltyX * dx * dx +
         config.CurvatureProjectionPenaltyY * dy * dy;
}

// used at the start of a seed. `path` must be given oldest-to-newest in the
// direction of travel toward the candidate (i.e. already reversed by the
// caller for backward growth).
double NextPointCost(const std::vector<std::size_t> &recent_in_order,
                    std::size_t candidate, const std::vector<Node> &nodes,
                    const std::vector<TMS_SpacePoint> &points,
                    const Config &config) {
  const std::size_t n = recent_in_order.size();
  if (config.UseCurvatureProjection && n >= 3) {
    return ProjectionCost(recent_in_order[n - 3], recent_in_order[n - 2],
                          recent_in_order[n - 1], candidate, nodes, points,
                          config);
  }
  if (n >= 2) {
    return TransitionCost(recent_in_order[n - 2], recent_in_order[n - 1],
                          candidate, nodes, points, config);
  }
  return 0.0;
}

std::size_t CountLayers(const Hypothesis &path,
                        const std::vector<Node> &nodes) {
  if (path.Nodes.empty()) return 0;
  std::size_t count = 1;
  for (std::size_t i = 1; i < path.Nodes.size(); ++i) {
    if (nodes[path.Nodes[i]].Layer != nodes[path.Nodes[i - 1]].Layer) ++count;
  }
  return count;
}

bool AcceptsSeed(const std::vector<std::vector<std::size_t> > &layers,
                 const std::vector<Node> &nodes,
                 const std::vector<TMS_SpacePoint> &points,
                 std::size_t sourceNode,
                 std::size_t targetNode,
                 const Config &config) {
  (void)points;
  (void)sourceNode;
  const std::size_t layerOccupancy = layers[nodes[targetNode].Layer].size();
  const std::size_t hitMultiplicity = nodes[targetNode].HitMultiplicity;
  return layerOccupancy <= config.MaxSeedLayerOccupancy &&
         hitMultiplicity <= config.MaxSeedHitMultiplicity;
}

// For an already-established hypothesis (>=2 accepted points) growing into
// a dense (shower-like) layer: search that layer's points directly instead
// of relying on the static forward/reverse graph, and rank purely by how
// well each point fits THIS hypothesis's own trajectory (the same
// field-aware kink/curvature-projection cost used everywhere else) rather
// than the static, history-blind per-edge cost that graph was built with.
// Deliberately does not apply AcceptsSeed's occupancy/multiplicity gate --
// that gate exists to keep fresh seeds out of ambiguous regions, not to
// stop an established trajectory from continuing through one. Returns at
// most Config::DenseLayerCandidates candidates, cheap enough that the
// combinatorics stay bounded even though the whole layer was searched.
std::vector<Link> SearchDenseLayer(std::size_t endpoint, std::size_t targetLayerIdx,
                                   const std::vector<std::vector<std::size_t> > &layers,
                                   const std::vector<Node> &nodes,
                                   const std::vector<TMS_SpacePoint> &points,
                                   const std::vector<std::size_t> &recent,
                                   const Config &config) {
  struct Candidate {
    Link link;
    double rankScore;
  };
  std::vector<Candidate> found;
  const TMS_SpacePoint &endpointPoint = points[nodes[endpoint].InputIndex];
  const std::size_t endpointLayer = nodes[endpoint].Layer;

  for (std::size_t candidate : layers[targetLayerIdx]) {
    const TMS_SpacePoint &candidatePoint = points[nodes[candidate].InputIndex];
    const bool candidateIsLower = candidatePoint.GetZ() < endpointPoint.GetZ();
    const TMS_SpacePoint &source = candidateIsLower ? candidatePoint : endpointPoint;
    const TMS_SpacePoint &target = candidateIsLower ? endpointPoint : candidatePoint;
    const double dz = target.GetZ() - source.GetZ();
    if (dz <= 0.0) continue;
    const double dxdz = (target.GetX() - source.GetX()) / dz;
    const double dydz = (target.GetY() - source.GetY()) / dz;
    if (std::abs(dxdz) > config.MaxAbsDXDZ || std::abs(dydz) > config.MaxAbsDYDZ) continue;
    if (config.MaxTimeDifference >= 0.0 &&
        std::abs(candidatePoint.GetTime() - endpointPoint.GetTime()) >
            config.MaxTimeDifference) {
      continue;
    }

    const double kink = NextPointCost(recent, candidate, nodes, points, config);
    if (!std::isfinite(kink) || kink > config.DenseLayerMaxKinkCost) continue;

    const std::size_t candidateLayer = nodes[candidate].Layer;
    const double gap = static_cast<double>(
        (candidateLayer > endpointLayer ? candidateLayer - endpointLayer
                                        : endpointLayer - candidateLayer) - 1);
    Link link;
    link.NodeIndex = candidate;
    // No occupancy/multiplicity term (the whole point of this path) and no
    // kink baked in either -- Extend() recomputes NextPointCost uniformly
    // for every candidate it processes, so including it here would double
    // count it. rankScore (below) is for selection only.
    link.Cost = -config.PointReward + config.GapPenalty * gap;
    found.push_back(Candidate{link, link.Cost + kink});
  }

  std::sort(found.begin(), found.end(), [](const Candidate &a, const Candidate &b) {
    return a.rankScore < b.rankScore;
  });
  const std::size_t keep = std::min<std::size_t>(config.DenseLayerCandidates, found.size());
  std::vector<Link> ranked;
  ranked.reserve(keep);
  for (std::size_t i = 0; i < keep; ++i) ranked.push_back(found[i].link);
  return ranked;
}

std::vector<Hypothesis> Extend(
    const Hypothesis &start, bool backward,
    const std::vector<std::vector<std::size_t> > &layers,
    const std::vector<std::vector<Link> > &forward,
    const std::vector<std::vector<Link> > &reverse,
    const std::vector<Node> &nodes,
    const std::vector<TMS_SpacePoint> &points,
    const Config &config, Diagnostics &stats) {
  std::vector<Hypothesis> active(1, start);
  std::vector<Hypothesis> finished;

  for (std::size_t step = 0; step < nodes.size() && !active.empty(); ++step) {
    std::vector<Hypothesis> next;
    for (const Hypothesis &path : active) {
      if (path.Nodes.size() >= config.MinPathPoints) finished.push_back(path);

      const std::size_t endpoint = backward ? path.Nodes.front() : path.Nodes.back();
      const std::vector<Link> &links = backward ? reverse[endpoint] : forward[endpoint];
      std::vector<std::size_t> recent;
      {
        const std::size_t n = path.Nodes.size();
        const std::size_t take = std::min<std::size_t>(3, n);
        if (backward) {
          for (std::size_t k = take; k-- > 0;) recent.push_back(path.Nodes[k]);
        } else {
          for (std::size_t k = n - take; k < n; ++k) recent.push_back(path.Nodes[k]);
        }
      }
      // A wide-baseline version of `recent`, used only for the dense-layer
      // search below: 3 points spread across as much of the established
      // path as available (up to the last 8), not just the nearest 3. Over
      // a short 2-3 point baseline, bar-pitch quantization noise (36mm
      // steps) is a large fraction of the actual slope -- a 13-point
      // established trajectory can do much better by averaging that noise
      // down over a longer lever arm, which is what actually lets a tight
      // acceptance window mean something in a shower core (see
      // DenseLayerMaxKinkCost). Falls back to `recent` when there isn't
      // enough history yet for the spread to matter.
      std::vector<std::size_t> wideRecent = recent;
      {
        const std::size_t n = path.Nodes.size();
        const std::size_t window = std::min<std::size_t>(8, n);
        if (window >= 5) {
          const std::size_t oldestOffset = window - 1;
          const std::size_t midOffset = window / 2;
          if (backward) {
            wideRecent = {path.Nodes[oldestOffset], path.Nodes[midOffset], path.Nodes[0]};
          } else {
            wideRecent = {path.Nodes[n - 1 - oldestOffset], path.Nodes[n - 1 - midOffset],
                          path.Nodes[n - 1]};
          }
        }
      }
      bool madeChild = false;
      bool resourceLimitHit = false;
      // The last (up to) 3 already-accepted points, ordered far-to-near
      // relative to the candidate, so NextPointCost can fit a curvature
      // trend from them regardless of growth direction. Shared by both the
      // static-graph path below and the dense-layer search (which passes
      // wideRecent instead of recent).
      auto tryLink = [&](const Link &link, bool enforceAcceptsSeed,
                         const std::vector<std::size_t> &recentForKink) {
        if (stats.HypothesesCreated >= config.MaxHypotheses) {
          resourceLimitHit = true;
          return;
        }
        if (enforceAcceptsSeed &&
            !AcceptsSeed(layers, nodes, points, endpoint, link.NodeIndex, config)) {
          return;
        }
        const double kink = NextPointCost(recentForKink, link.NodeIndex, nodes, points, config);
        if (!std::isfinite(kink)) return;
        Hypothesis child = path;
        child.Score += link.Cost + kink;
        if (backward) child.Nodes.insert(child.Nodes.begin(), link.NodeIndex);
        else child.Nodes.push_back(link.NodeIndex);
        next.push_back(std::move(child));
        ++stats.HypothesesCreated;
        madeChild = true;
      };

      for (const Link &link : links) {
        tryLink(link, /*enforceAcceptsSeed=*/true, recent);
        if (resourceLimitHit) break;
      }

      // See SearchDenseLayer / Config::DenseLayerCandidates: once a real
      // trajectory exists, reach past the static graph into any dense
      // layer within range instead of relying on AcceptsSeed to (wrongly)
      // veto it outright.
      if (config.UseDenseLayerSearch && !resourceLimitHit && wideRecent.size() >= 2) {
        const std::size_t endpointLayer = nodes[endpoint].Layer;
        const std::size_t gapLimit = static_cast<std::size_t>(config.MaxLayerGap);
        bool haveRange = true;
        std::size_t loLayer = 0, hiLayer = 0;
        if (backward) {
          if (endpointLayer == 0) {
            haveRange = false;
          } else {
            hiLayer = endpointLayer - 1;
            loLayer = (endpointLayer > gapLimit) ? (endpointLayer - gapLimit) : 0;
          }
        } else {
          if (endpointLayer + 1 >= layers.size()) {
            haveRange = false;
          } else {
            loLayer = endpointLayer + 1;
            hiLayer = std::min(layers.size() - 1, endpointLayer + gapLimit);
          }
        }
        if (haveRange) {
          for (std::size_t layerIdx = loLayer; layerIdx <= hiLayer && !resourceLimitHit;
               ++layerIdx) {
            if (layers[layerIdx].size() <= config.MaxSeedLayerOccupancy) continue;
            for (const Link &link : SearchDenseLayer(endpoint, layerIdx, layers, nodes,
                                                      points, wideRecent, config)) {
              tryLink(link, /*enforceAcceptsSeed=*/false, wideRecent);
              if (resourceLimitHit) break;
            }
          }
        }
      }
      if (resourceLimitHit) stats.ResourceLimitReached = true;

      if (!madeChild) finished.push_back(path);
      if (stats.ResourceLimitReached) break;
    }

    KeepBest(next, config.BeamWidth, stats);
    KeepBest(finished, config.BeamWidth, stats);
    active.swap(next);
    if (stats.ResourceLimitReached) break;
  }
  finished.insert(finished.end(), active.begin(), active.end());
  KeepBest(finished, config.MaxPathsPerSeed, stats);
  return finished;
}

double OverlapFraction(const Path &a, const Path &b) {
  if (a.SpacePointIndices.empty() || b.SpacePointIndices.empty()) return 0.0;
  std::unordered_set<std::size_t> indices(a.SpacePointIndices.begin(),
                                          a.SpacePointIndices.end());
  std::size_t overlap = 0;
  for (std::size_t index : b.SpacePointIndices) {
    if (indices.count(index)) ++overlap;
  }
  return static_cast<double>(overlap) /
         static_cast<double>(std::min(a.SpacePointIndices.size(),
                                      b.SpacePointIndices.size()));
}

} // namespace

Finder::Finder(const Config &config) : fConfig(config) {
  if (fConfig.LayerZTolerance <= 0.0 || fConfig.MaxLayerGap < 1 ||
      fConfig.SeedLength < 2 || fConfig.MaxLinksPerTargetLayer == 0 ||
      fConfig.BeamWidth == 0 || fConfig.MaxSeeds == 0 ||
      fConfig.SeedOverlapFraction < 0.0 ||
      fConfig.SeedOverlapFraction > 1.0 ||
      fConfig.DuplicateOverlapFraction < 0.0 ||
      fConfig.DuplicateOverlapFraction > 1.0) {
    throw std::invalid_argument("Invalid Graph Track Finder configuration");
  }
}

Result Finder::Find(const std::vector<TMS_SpacePoint> &spacePoints) const {
  Result result;
  Diagnostics &stats = result.Stats;
  stats.InputPoints = spacePoints.size();
  if (spacePoints.size() < fConfig.SeedLength) return result;

  // Stable z ordering makes all subsequent graph operations deterministic.
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

  std::unordered_map<int, std::size_t> xMultiplicity;
  std::unordered_map<int, std::size_t> yMultiplicity;
  for (const TMS_SpacePoint &point : spacePoints) {
    if (point.GetXHitIndex() >= 0) ++xMultiplicity[point.GetXHitIndex()];
    if (point.GetYHitIndex() >= 0) ++yMultiplicity[point.GetYHitIndex()];
  }

  std::vector<Node> nodes;
  std::vector<std::vector<std::size_t> > layers;
  for (std::size_t inputIndex : order) {
    const TMS_SpacePoint &point = spacePoints[inputIndex];
    if (layers.empty() ||
        std::abs(point.GetZ() -
                 spacePoints[nodes[layers.back().front()].InputIndex].GetZ()) >
            fConfig.LayerZTolerance) {
      layers.push_back(std::vector<std::size_t>());
    }
    Node node;
    node.InputIndex = inputIndex;
    node.Layer = layers.size() - 1;
    const std::size_t xCount = point.GetXHitIndex() >= 0
        ? xMultiplicity[point.GetXHitIndex()] : 0;
    const std::size_t yCount = point.GetYHitIndex() >= 0
        ? yMultiplicity[point.GetYHitIndex()] : 0;
    node.HitMultiplicity = std::max(xCount, yCount);
    nodes.push_back(node);
    layers.back().push_back(nodes.size() - 1);
  }
  stats.Layers = layers.size();

  std::vector<std::vector<Link> > forward(nodes.size());
  std::vector<std::vector<Link> > reverse(nodes.size());

  for (std::size_t sourceLayer = 0; sourceLayer < layers.size(); ++sourceLayer) {
    const std::size_t lastTarget = std::min(
        layers.size() - 1, sourceLayer + static_cast<std::size_t>(fConfig.MaxLayerGap));
    for (std::size_t targetLayer = sourceLayer + 1;
         targetLayer <= lastTarget; ++targetLayer) {
      for (std::size_t sourceNode : layers[sourceLayer]) {
        std::vector<Link> candidates;
        const TMS_SpacePoint &source = spacePoints[nodes[sourceNode].InputIndex];
        for (std::size_t targetNode : layers[targetLayer]) {
          ++stats.LinksTested;
          const TMS_SpacePoint &target = spacePoints[nodes[targetNode].InputIndex];
          const double dz = target.GetZ() - source.GetZ();
          if (dz <= 0.0) continue;
          const double dxdz = (target.GetX() - source.GetX()) / dz;
          const double dydz = (target.GetY() - source.GetY()) / dz;
          if (std::abs(dxdz) > fConfig.MaxAbsDXDZ ||
              std::abs(dydz) > fConfig.MaxAbsDYDZ) continue;
          if (fConfig.MaxTimeDifference >= 0.0 &&
              std::abs(target.GetTime() - source.GetTime()) >
                  fConfig.MaxTimeDifference) continue;

          const double gap = static_cast<double>(targetLayer - sourceLayer - 1);
          const double occupancy = static_cast<double>(layers[targetLayer].size() - 1);
          const double multiplicity = static_cast<double>(
              nodes[targetNode].HitMultiplicity > 0
                  ? nodes[targetNode].HitMultiplicity - 1
                  : 0);
          Link link;
          link.NodeIndex = targetNode;
          link.Cost = -fConfig.PointReward + fConfig.GapPenalty * gap +
                      fConfig.SlopePenalty * (dxdz * dxdz + dydz * dydz) +
                      fConfig.OccupancyPenalty * occupancy +
                      fConfig.HitMultiplicityPenalty * multiplicity;
          candidates.push_back(link);
        }
        std::sort(candidates.begin(), candidates.end(), [](const Link &a,
                                                           const Link &b) {
          return a.Cost < b.Cost;
        });
        if (candidates.size() > fConfig.MaxLinksPerTargetLayer)
          candidates.resize(fConfig.MaxLinksPerTargetLayer);
        for (const Link &link : candidates) {
          forward[sourceNode].push_back(link);
          reverse[link.NodeIndex].push_back(Link{sourceNode, link.Cost});
          ++stats.LinksAccepted;
        }
      }
    }
  }

  // Form short, clean stubs.  This is intentionally global: no first-plane or
  // entry-point assumption is made.
  std::vector<Hypothesis> seeds;
  for (std::size_t source = 0; source < nodes.size(); ++source) {
    for (const Link &link : forward[source]) {
      if (!AcceptsSeed(layers, nodes, spacePoints, source, link.NodeIndex, fConfig))
        continue;
      seeds.push_back(Hypothesis{{source, link.NodeIndex}, link.Cost});
      ++stats.HypothesesCreated;
    }
  }
  KeepBest(seeds, fConfig.MaxSeedFrontier, stats);

  while (!seeds.empty() && seeds.front().Nodes.size() < fConfig.SeedLength) {
    std::vector<Hypothesis> grown;
    for (const Hypothesis &seed : seeds) {
      const std::size_t endpoint = seed.Nodes.back();
      for (const Link &link : forward[endpoint]) {
        if (stats.HypothesesCreated >= fConfig.MaxHypotheses) {
          stats.ResourceLimitReached = true;
          break;
        }
        if (!AcceptsSeed(layers, nodes, spacePoints, endpoint, link.NodeIndex, fConfig) ||
          ConflictsWithPath(seed.Nodes, link.NodeIndex))
          continue;
        Hypothesis child = seed;
        const std::size_t n = seed.Nodes.size();
        const std::size_t take = std::min<std::size_t>(3, n);
        std::vector<std::size_t> recent(seed.Nodes.begin() + (n - take), seed.Nodes.end());
        child.Score += link.Cost +
            NextPointCost(recent, link.NodeIndex, nodes, spacePoints, fConfig);
        child.Nodes.push_back(link.NodeIndex);
        grown.push_back(std::move(child));
        ++stats.HypothesesCreated;
      }
      if (stats.ResourceLimitReached) break;
    }
    seeds.swap(grown);
    KeepBest(seeds, fConfig.MaxSeedFrontier, stats);
    if (stats.ResourceLimitReached) break;
  }
  stats.SeedsGenerated = seeds.size();
  KeepDiverseSeeds(seeds, fConfig, stats);
  stats.SeedsRetained = seeds.size();

  std::vector<Path> candidatePaths;
  for (const Hypothesis &seed : seeds) {
    std::vector<Hypothesis> lowZ = Extend(seed, true, layers, forward, reverse,
                                          nodes, spacePoints, fConfig, stats);
    std::vector<Hypothesis> complete;
    for (const Hypothesis &partial : lowZ) {
      std::vector<Hypothesis> highZ = Extend(partial, false, layers, forward, reverse,
                                             nodes, spacePoints, fConfig, stats);
      complete.insert(complete.end(), highZ.begin(), highZ.end());
    }
    KeepBest(complete, fConfig.MaxPathsPerSeed, stats);
    for (const Hypothesis &path : complete) {
      if (path.Nodes.size() < fConfig.MinPathPoints) continue;
      Path output;
      output.Score = path.Score;
      output.DistinctLayers = CountLayers(path, nodes);
      for (std::size_t nodeIndex : path.Nodes)
        output.SpacePointIndices.push_back(nodes[nodeIndex].InputIndex);
      candidatePaths.push_back(std::move(output));
    }
    if (stats.ResourceLimitReached) break;
  }

  std::sort(candidatePaths.begin(), candidatePaths.end(), [](const Path &a,
                                                              const Path &b) {
    if (a.Score != b.Score) return a.Score < b.Score;
    return a.DistinctLayers > b.DistinctLayers;
  });
  stats.PathsBeforeDeduplication = candidatePaths.size();
  for (const Path &candidate : candidatePaths) {
    bool duplicate = false;
    for (const Path &kept : result.Paths) {
      if (OverlapFraction(candidate, kept) >= fConfig.DuplicateOverlapFraction) {
        duplicate = true;
        break;
      }
    }
    if (!duplicate) result.Paths.push_back(candidate);
    if (result.Paths.size() >= fConfig.MaxOutputPaths) break;
  }
  stats.PathsAfterDeduplication = result.Paths.size();
  return result;
}

} // namespace TMS_GraphTrackFinder
