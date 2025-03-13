#include <iostream>
#include "TMS_DBScan.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"

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
  const double mean = 4;
  for (int i = 0; i < nPoints; ++i) {
    // Make some points around somewhere
    double x,y = -1;
    if (i < nPoints-40) {
      x = rnd.Gaus(mean,1);
      y = rnd.Gaus(mean,1);
    } else if (i < nPoints - 10) {
      x = mean+(i-(nPoints-40));
      y = mean+1.5*(i-(nPoints-40));
    } else {
      x = mean+(nPoints-i);
      y = mean-1.5*(nPoints-i);
    }

    TMS_DBScan_Point point(x,y,0);
    Points.push_back(point);
  }

  // Draw as an example
  TCanvas canv("canv", "canv", 1024, 1024);
  TString canvname = Form("points_%2.2f_%i.pdf", eps, nMin);
  canv.Print(canvname+"[");
  TGraph g(nPoints);
  for (int i = 0; i < nPoints; ++i) {
    g.SetPoint(i, Points[i].x, Points[i].y);
  }
  g.SetMarkerSize(3);
  g.SetMarkerStyle(kCircle);

  TMS_DBScan cluster(nMin, eps, Points);

  cluster.RunDBScan();

  std::vector<TMS_DBScan_Point> Scanned = cluster.GetPoints();
  int nCluster = cluster.GetNClusters();
  int nNoise = cluster.GetNNoise();
  if (nCluster == 0) {
    std::cout << "Found no clusters" << std::endl;
    return -1;
  }

  // Make a graph for each of the clusters
  TGraph **Graphs = new TGraph*[nCluster];
  for (int i = 0; i < nCluster; ++i) {
    Graphs[i] = new TGraph(nPoints);
    Graphs[i]->SetMarkerSize(3);
    Graphs[i]->SetMarkerStyle(kCircle);
    Graphs[i]->SetMarkerColor(i+2);
  }
  // And one graph for the noise hits
  TGraph *Noise = new TGraph(nPoints);
  Noise->SetMarkerSize(3);
  Noise->SetMarkerStyle(20);
  Noise->SetMarkerColor(kBlack);
  Noise->SetTitle("Noise");

  // Now make some TGraphs showing each cluster
  int npoint = 0;
  for (auto &i: Scanned) {
    if (i.ClusterID != int(PointClassifier::kNoise)) Graphs[i.ClusterID-1]->SetPoint(npoint, i.x, i.y);
    else Noise->SetPoint(npoint, i.x, i.y);
    npoint++;
  }

  std::cout << "Number of points: " << Scanned.size() << std::endl;
  std::cout << "Number of clusters: " << nCluster << std::endl;
  std::cout << "Number of noise hits: " << nNoise << std::endl;
  std::cout << "Number of unclassified: " << cluster.GetNUnclassified() << std::endl;

  TLegend leg(0.1, 0.7, 0.65, 0.85);
  leg.AddEntry(&g, "Graph", "p");
  leg.SetFillStyle(0);
  leg.SetLineWidth(0);
  leg.SetBorderSize(0);
  g.Draw("AP");
  for (int i = 0; i < nCluster; ++i) {
    Graphs[i]->Draw("P,same");
    leg.AddEntry(Graphs[i], Form("Cluster %i", i), "p");
  }
  Noise->Draw("P,same");
  leg.AddEntry(Noise, "Noise", "p");

  leg.Draw("same");
  canv.Print(canvname);

  canv.Print(canvname+"]");
}
