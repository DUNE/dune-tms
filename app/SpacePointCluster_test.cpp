// Stage 1 standalone test: verifies TMS_KDTree's radius query against a
// brute-force O(N^2) check, then runs TMS_SpacePointDBScan on synthetic 3D
// data (two Gaussian blobs + two lines) and plots the resulting clusters.
// Mirrors app/DBSCAN_test.cpp's structure and plotting style.

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "TMS_KDTree.h"
#include "TMS_SpacePoint.h"
#include "TMS_SpacePointDBScan.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TString.h"

// Returns true if TMS_KDTree::RadiusQuery agrees exactly (same index set)
// with a brute-force O(N^2) distance check, for every point in `pts`.
bool CheckKDTreeAgainstBruteForce(const std::vector<std::array<double, 3>> &pts, double radius) {
  TMS_KDTree tree(pts);
  bool all_ok = true;
  for (size_t i = 0; i < pts.size(); ++i) {
    std::vector<int> tree_result = tree.RadiusQuery(static_cast<int>(i), radius);
    std::vector<int> brute_result;
    for (size_t j = 0; j < pts.size(); ++j) {
      double dx = pts[i][0] - pts[j][0];
      double dy = pts[i][1] - pts[j][1];
      double dz = pts[i][2] - pts[j][2];
      if (std::sqrt(dx * dx + dy * dy + dz * dz) <= radius) brute_result.push_back(static_cast<int>(j));
    }
    std::sort(tree_result.begin(), tree_result.end());
    std::sort(brute_result.begin(), brute_result.end());
    if (tree_result != brute_result) {
      std::cout << "MISMATCH at query point " << i << ": tree found " << tree_result.size()
                << " neighbours, brute force found " << brute_result.size() << std::endl;
      all_ok = false;
    }
  }
  return all_ok;
}

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cout << "Need 2 arguments: epsilon and nmin" << std::endl;
    return -1;
  }
  double eps = std::atof(argv[1]);
  int nMin = std::atoi(argv[2]);

  // Synthetic 3D data: two well-separated Gaussian blobs (should NOT cluster
  // together, and should have low PCA linearity when that's tested in Stage
  // 2) plus two lines with different orientations (should cluster as tight,
  // elongated groups).
  TRandom3 rnd(42); // fixed non-zero seed -- TRandom3(0) auto-seeds from system entropy, giving non-reproducible test data
  std::vector<std::array<double, 3>> pts;

  for (int i = 0; i < 30; ++i) pts.push_back({rnd.Gaus(0, 3), rnd.Gaus(0, 3), rnd.Gaus(0, 3)});
  for (int i = 0; i < 30; ++i) pts.push_back({rnd.Gaus(100, 3), rnd.Gaus(100, 3), rnd.Gaus(100, 3)});
  for (int i = 0; i < 40; ++i) {
    double z = i * 2.0;
    pts.push_back({rnd.Gaus(50, 0.3), rnd.Gaus(50, 0.3), z});
  }
  for (int i = 0; i < 40; ++i) {
    double t = i * 2.0;
    pts.push_back({-50 + 0.5 * t, -50 + 0.3 * t, -50 + 0.8 * t + rnd.Gaus(0, 0.3)});
  }

  // ---- Correctness gate: KD-tree radius query vs brute force ----
  std::vector<double> test_epsilons = {1.0, 5.0, 10.0, 20.0, 50.0};
  bool all_ok = true;
  for (double e : test_epsilons) {
    bool ok = CheckKDTreeAgainstBruteForce(pts, e);
    std::cout << "KD-tree vs brute force @ epsilon=" << e << ": " << (ok ? "MATCH" : "MISMATCH") << std::endl;
    all_ok = all_ok && ok;
  }
  if (!all_ok) {
    std::cerr << "KD-tree radius query does not match brute force! Aborting." << std::endl;
    return -1;
  }
  std::cout << "KD-tree radius query verified against brute force for all test epsilons." << std::endl;

  // ---- Cluster the same points via TMS_SpacePointDBScan ----
  std::vector<TMS_SpacePoint> space_points;
  for (const auto &p : pts) space_points.emplace_back(p[0], p[1], p[2], -1, -1, 0.0);

  TMS_SpacePointDBScan dbscan(space_points, static_cast<unsigned int>(nMin), eps);
  std::vector<std::vector<int>> clusters = dbscan.RunAndGetClusterIndices();

  std::vector<bool> in_cluster(space_points.size(), false);
  for (auto &c : clusters)
    for (int idx : c) in_cluster[idx] = true;
  int n_noise = 0;
  for (bool b : in_cluster)
    if (!b) ++n_noise;

  std::cout << "Number of points: " << space_points.size() << std::endl;
  std::cout << "Number of clusters: " << clusters.size() << std::endl;
  std::cout << "Number of noise points: " << n_noise << std::endl;

  // ---- Plot: X-Z and Y-Z projections, colored by cluster ----
  TCanvas canv("canv", "canv", 1600, 800);
  canv.Divide(2, 1);
  TString canvname = Form("spacepoint_clusters_%2.2f_%i.pdf", eps, nMin);
  canv.Print(canvname + "[");

  const int nClusters = static_cast<int>(clusters.size());
  std::vector<TGraph *> graphs_xz(nClusters), graphs_yz(nClusters);
  for (int i = 0; i < nClusters; ++i) {
    graphs_xz[i] = new TGraph(static_cast<int>(clusters[i].size()));
    graphs_yz[i] = new TGraph(static_cast<int>(clusters[i].size()));
    for (TGraph *g : {graphs_xz[i], graphs_yz[i]}) {
      g->SetMarkerSize(1.5);
      g->SetMarkerStyle(kFullCircle);
      g->SetMarkerColor(i + 2);
    }
    for (size_t k = 0; k < clusters[i].size(); ++k) {
      int idx = clusters[i][k];
      graphs_xz[i]->SetPoint(static_cast<int>(k), pts[idx][2], pts[idx][0]);
      graphs_yz[i]->SetPoint(static_cast<int>(k), pts[idx][2], pts[idx][1]);
    }
  }

  TGraph noise_xz(n_noise), noise_yz(n_noise);
  noise_xz.SetMarkerStyle(kFullSquare);
  noise_xz.SetMarkerColor(kBlack);
  noise_xz.SetMarkerSize(1.2);
  noise_yz.SetMarkerStyle(kFullSquare);
  noise_yz.SetMarkerColor(kBlack);
  noise_yz.SetMarkerSize(1.2);
  int npoint = 0;
  for (size_t idx = 0; idx < pts.size(); ++idx) {
    if (!in_cluster[idx]) {
      noise_xz.SetPoint(npoint, pts[idx][2], pts[idx][0]);
      noise_yz.SetPoint(npoint, pts[idx][2], pts[idx][1]);
      ++npoint;
    }
  }

  // Invisible full-range graphs so the axes cover all points regardless of draw order.
  TGraph frame_xz(static_cast<int>(pts.size()));
  TGraph frame_yz(static_cast<int>(pts.size()));
  for (size_t i = 0; i < pts.size(); ++i) {
    frame_xz.SetPoint(static_cast<int>(i), pts[i][2], pts[i][0]);
    frame_yz.SetPoint(static_cast<int>(i), pts[i][2], pts[i][1]);
  }
  frame_xz.SetMarkerColorAlpha(kWhite, 0);
  frame_yz.SetMarkerColorAlpha(kWhite, 0);
  frame_xz.SetTitle("X-Z projection;Z;X");
  frame_yz.SetTitle("Y-Z projection;Z;Y");

  TLegend leg(0.12, 0.7, 0.5, 0.88);
  leg.SetFillStyle(0);
  leg.SetLineWidth(0);
  leg.SetBorderSize(0);
  for (int i = 0; i < nClusters; ++i) leg.AddEntry(graphs_xz[i], Form("Cluster %i (n=%zu)", i, clusters[i].size()), "p");
  leg.AddEntry(&noise_xz, Form("Noise (n=%i)", n_noise), "p");

  canv.cd(1);
  frame_xz.Draw("AP");
  for (int i = 0; i < nClusters; ++i) graphs_xz[i]->Draw("P,same");
  noise_xz.Draw("P,same");
  leg.Draw("same");

  canv.cd(2);
  frame_yz.Draw("AP");
  for (int i = 0; i < nClusters; ++i) graphs_yz[i]->Draw("P,same");
  noise_yz.Draw("P,same");

  canv.Print(canvname);
  canv.Print(canvname + "]");

  std::cout << "Saved plot to " << canvname << std::endl;
  return 0;
}
