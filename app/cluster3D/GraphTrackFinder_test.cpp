#include "TMS_GraphTrackFinder.h"

#include <algorithm>
#include <iostream>
#include <random>
#include <set>
#include <vector>

namespace {

struct Sample {
  std::vector<TMS_SpacePoint> Points;
  std::set<std::size_t> MuonIndices;
  std::set<std::size_t> SecondaryIndices;
};

Sample MakeMessyEntranceEvent() {
  Sample sample;
  std::mt19937 random(73021);
  std::normal_distribution<double> resolution(0.0, 2.5);
  std::uniform_real_distribution<double> showerX(-900.0, 1300.0);
  std::uniform_real_distribution<double> showerY(-700.0, 800.0);
  std::uniform_real_distribution<double> showerTime(85.0, 115.0);
  int nextHit = 1000;

  // A through-going, gently bending muon. Curvature is in x, as expected for
  // the final Y-Y-X layout (Y bars measure x).
  std::vector<int> muonXHit(36);
  std::vector<int> muonYHit(36);
  for (int layer = 0; layer < 36; ++layer) {
    const double z = 100.0 * layer;
    const double x = 50.0 + 0.055 * z + 0.000008 * z * z + resolution(random);
    const double y = -35.0 + 0.020 * z + resolution(random);
    muonXHit[layer] = nextHit++;
    muonYHit[layer] = nextHit++;
    sample.Points.push_back(TMS_SpacePoint(
        x, y, z, muonXHit[layer], muonYHit[layer], 100.0 + 0.05 * layer));
  }

  // A shorter real prong in the same DBSCAN cluster. The finder may retain it
  // as a second path rather than forcing the cluster into one track.
  const int firstSecondaryHit = nextHit;
  for (int layer = 2; layer < 16; ++layer) {
    const double z = 100.0 * layer;
    const int xHit = nextHit++;
    const int yHit = nextHit++;
    sample.Points.push_back(TMS_SpacePoint(
        110.0 - 0.18 * z + resolution(random),
        -20.0 + 0.11 * z + resolution(random), z,
        xHit, yHit, 101.0 + 0.04 * layer));
  }
  const int lastSecondaryHit = nextHit;

  // Shower activity makes the upstream layers unsuitable for seeding. Some
  // points reuse one muon native hit, mimicking space-point ghosts.
  for (int layer = 0; layer < 10; ++layer) {
    const double z = 100.0 * layer;
    for (int i = 0; i < 28; ++i) {
      int xHit = nextHit++;
      int yHit = nextHit++;
      if (i < 5) xHit = muonXHit[layer];
      if (i >= 5 && i < 10) yHit = muonYHit[layer];
      sample.Points.push_back(TMS_SpacePoint(
          showerX(random), showerY(random), z, xHit, yHit,
          showerTime(random)));
    }
  }

  std::shuffle(sample.Points.begin(), sample.Points.end(), random);

  // Recover truth labels after shuffling. This is test bookkeeping only.
  for (std::size_t i = 0; i < sample.Points.size(); ++i) {
    const int xHit = sample.Points[i].GetXHitIndex();
    const int yHit = sample.Points[i].GetYHitIndex();
    bool isMuon = false;
    for (int layer = 0; layer < 36; ++layer) {
      if (xHit == muonXHit[layer] && yHit == muonYHit[layer]) {
        isMuon = true;
        break;
      }
    }
    if (isMuon) sample.MuonIndices.insert(i);
    else if (xHit >= firstSecondaryHit && xHit < lastSecondaryHit &&
             yHit >= firstSecondaryHit && yHit < lastSecondaryHit)
      sample.SecondaryIndices.insert(i);
  }
  return sample;
}

Sample MakeDenseEntranceEvent() {
  Sample sample;
  std::mt19937 random(98765);
  std::normal_distribution<double> resolution(0.0, 2.5);
  std::uniform_real_distribution<double> showerX(-1000.0, 1500.0);
  std::uniform_real_distribution<double> showerY(-800.0, 900.0);
  std::uniform_real_distribution<double> showerTime(80.0, 120.0);
  int nextHit = 2000;

  std::vector<int> muonXHit(24);
  std::vector<int> muonYHit(24);
  for (int layer = 0; layer < 24; ++layer) {
    const double z = 100.0 * layer;
    const double x = 80.0 + 0.06 * z + 0.00001 * z * z + resolution(random);
    const double y = -42.0 + 0.018 * z + resolution(random);
    muonXHit[layer] = nextHit++;
    muonYHit[layer] = nextHit++;
    sample.Points.push_back(TMS_SpacePoint(
        x, y, z, muonXHit[layer], muonYHit[layer], 95.0 + 0.04 * layer));
  }

  for (int layer = 0; layer < 6; ++layer) {
    const double z = 100.0 * layer;
    for (int i = 0; i < 42; ++i) {
      int xHit = nextHit++;
      int yHit = nextHit++;
      if (i < 10) xHit = muonXHit[layer];
      if (i >= 10 && i < 20) yHit = muonYHit[layer];
      sample.Points.push_back(TMS_SpacePoint(
          showerX(random), showerY(random), z, xHit, yHit, showerTime(random)));
    }
  }

  for (int layer = 6; layer < 24; ++layer) {
    const double z = 100.0 * layer;
    for (int i = 0; i < 8; ++i) {
      int xHit = nextHit++;
      int yHit = nextHit++;
      sample.Points.push_back(TMS_SpacePoint(
          120.0 + 0.25 * z + resolution(random),
          -20.0 + 0.07 * z + resolution(random), z,
          xHit, yHit, 95.0 + 0.04 * layer));
    }
  }

  std::shuffle(sample.Points.begin(), sample.Points.end(), random);
  for (std::size_t i = 0; i < sample.Points.size(); ++i) {
    const int xHit = sample.Points[i].GetXHitIndex();
    const int yHit = sample.Points[i].GetYHitIndex();
    bool isMuon = false;
    for (int layer = 0; layer < 24; ++layer) {
      if (xHit == muonXHit[layer] && yHit == muonYHit[layer]) {
        isMuon = true;
        break;
      }
    }
    if (isMuon) sample.MuonIndices.insert(i);
  }
  return sample;
}

std::size_t CountTruth(const TMS_GraphTrackFinder::Path &path,
                       const std::set<std::size_t> &truth) {
  std::size_t count = 0;
  for (std::size_t index : path.SpacePointIndices)
    if (truth.count(index)) ++count;
  return count;
}

} // namespace

int main() {
  const Sample sample = MakeMessyEntranceEvent();
  TMS_GraphTrackFinder::Config config;
  config.MaxSeedLayerOccupancy = 12;
  config.MaxSeedHitMultiplicity = 4;
  config.MaxAbsDXDZ = 0.8;
  config.MaxAbsDYDZ = 0.8;

  const TMS_GraphTrackFinder::Result result =
      TMS_GraphTrackFinder::Finder(config).Find(sample.Points);

  std::cout << "Graph Track Finder synthetic trial\n"
            << "  points: " << result.Stats.InputPoints << '\n'
            << "  z layers: " << result.Stats.Layers << '\n'
            << "  links tested/kept: " << result.Stats.LinksTested << "/"
            << result.Stats.LinksAccepted << '\n'
            << "  seeds generated/retained: " << result.Stats.SeedsGenerated
            << "/" << result.Stats.SeedsRetained << '\n'
            << "  hypotheses created/pruned: "
            << result.Stats.HypothesesCreated << "/"
            << result.Stats.HypothesesPruned << '\n'
            << "  paths before/after deduplication: "
            << result.Stats.PathsBeforeDeduplication << "/"
            << result.Stats.PathsAfterDeduplication << '\n';

  std::size_t bestMuon = 0;
  std::size_t bestSecondary = 0;
  for (std::size_t i = 0; i < result.Paths.size(); ++i) {
    const std::size_t muon = CountTruth(result.Paths[i], sample.MuonIndices);
    const std::size_t secondary =
        CountTruth(result.Paths[i], sample.SecondaryIndices);
    bestMuon = std::max(bestMuon, muon);
    bestSecondary = std::max(bestSecondary, secondary);
    std::cout << "  path " << i << ": points="
              << result.Paths[i].SpacePointIndices.size()
              << ", layers=" << result.Paths[i].DistinctLayers
              << ", score=" << result.Paths[i].Score
              << ", muon truth=" << muon
              << ", secondary truth=" << secondary << '\n';
  }

  // Success means a clean downstream seed recovered the unambiguous muon tail
  // without relying on a unique first-plane entry point. It is intentionally
  // allowed to stop rather than absorb dubious points in the shower core.
  if (bestMuon < 24) {
    std::cerr << "FAIL: best path contains only " << bestMuon
              << " of 36 muon points\n";
    return 1;
  }
  if (result.Stats.ResourceLimitReached) {
    std::cerr << "FAIL: configured hypothesis limit was reached\n";
    return 2;
  }

  const Sample dense = MakeDenseEntranceEvent();
  const TMS_GraphTrackFinder::Result denseResult =
      TMS_GraphTrackFinder::Finder(config).Find(dense.Points);
  std::size_t denseMuon = 0;
  for (const TMS_GraphTrackFinder::Path &path : denseResult.Paths) {
    denseMuon = std::max(denseMuon, CountTruth(path, dense.MuonIndices));
  }
  if (denseMuon < 10) {
    std::cerr << "FAIL: dense-entrance case recovered only " << denseMuon
              << "/24 muon points before the path breaks\n";
    return 3;
  }

  std::cout << "PASS: recovered " << bestMuon
            << "/36 muon points; best secondary path contains "
            << bestSecondary << "/14 points\n";
  std::cout << "PASS: dense-entry case recovered " << denseMuon
            << "/24 early muon points\n";
  return 0;
}
