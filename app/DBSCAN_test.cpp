#include <iostream>
#include "TMS_DBScan.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TCanvas.h"

int main(int argc, char **argv) {

  if (argc != 3) {
    std::cout << "Need 2 arguments, epsilon and nmin" << std::endl;
    return -1;
  }

  double eps = std::atof(argv[1]);
  int nMin = std::atoi(argv[2]);

  std::vector<TMS_DBScan_Point> Points;
  const int nPoints = 100;
  TRandom3 rnd(0);
  //int factor = 0;
  for (int i = 0; i < nPoints; ++i) {
    // Make some points around somewhere
    double x,y = -1;
    if (i < nPoints-20) {
      x = rnd.Gaus(4,1);
      y = rnd.Gaus(4,1);
    } else {
      x = 4+(nPoints-i);
      y = 4+(nPoints-i);
    }

    TMS_DBScan_Point point(x,y,0);
    Points.push_back(point);
  }

  // Draw as an example
  TCanvas canv("canv", "canv", 1024, 1024);
  canv.Print("points.pdf[");
  TGraph g(nPoints);
  for (int i = 0; i < nPoints; ++i) {
    g.SetPoint(i, Points[i].x, Points[i].y);
  }
  g.Draw("AP");
  g.SetMarkerSize(3);
  g.SetMarkerStyle(kCircle);
  canv.Print("points.pdf");

  TMS_DBScan cluster(nMin, eps, Points);

  cluster.RunDBScan();

  std::vector<TMS_DBScan_Point> Scanned = cluster.GetPoints();
  int nCluster = cluster.GetNClusters();
  int nNoise = cluster.GetNNoise();
  if (nCluster == 0) {
    std::cout << "Found no clusters" << std::endl;
    return -1;
  }

  TGraph **Graphs = new TGraph*[nCluster];
  for (int i = 0; i < nCluster; ++i) {
    Graphs[i] = new TGraph(nPoints);
    Graphs[i]->SetMarkerSize(3);
    Graphs[i]->SetMarkerStyle(kCircle);
    Graphs[i]->SetMarkerColor(i+2);
  }
  Graphs[nCluster-1]->SetMarkerStyle(20);
  Graphs[nCluster-1]->SetMarkerColor(kBlack);

  // Now make some TGraphs showing each cluster
  int npoint = 0;
  for (auto &i: Scanned) {
    if (i.ClusterID != kNoise) Graphs[i.ClusterID-1]->SetPoint(npoint, i.x, i.y);
    else Graphs[nCluster-1]->SetPoint(npoint, i.x, i.y);
    npoint++;
  }

  std::cout << "Number of points: " << Scanned.size() << std::endl;
  std::cout << "Number of clusters: " << nCluster << std::endl;
  std::cout << "Number of noise hits: " << nNoise << std::endl;
  std::cout << "Number of unclassified: " << cluster.GetNUnclassified() << std::endl;

  Graphs[nCluster-1]->Draw("P,same");
  for (int i = 0; i < nCluster-1; ++i) {
    Graphs[i]->Draw("P,same");
  }
  canv.Print("points.pdf");

  canv.Print("points.pdf]");
}
