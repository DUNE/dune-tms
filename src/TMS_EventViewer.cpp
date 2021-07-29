#include "TMS_EventViewer.h"

TMS_EventViewer::TMS_EventViewer() :
DrawTrackFinding(false) {

  nDraws = 0;

  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(255);

  const double zmin = TMS_Const::TMS_Start[2];
  const double zmax = TMS_Const::TMS_End[2];
  const double xmin = TMS_Const::TMS_Start[0];
  const double xmax = TMS_Const::TMS_End[0];
  const double ymin = TMS_Const::TMS_Start[1];
  const double ymax = TMS_Const::TMS_End[1];
  // Scint bars are 4 by 1 cm
  //const int nbinsz = ((zmax-zmin)/10)/5;
  const int nbinsz = ((zmax-zmin)/100)*2;
  const int nbinsx = ((xmax-xmin)/100)*3;
  const int nbinsy = ((ymax-ymin)/100)*3;

  // The 2D views
  xz_view = new TH2D("TMS_Viewer_xz", "TMS viewer xz;z (mm); x (mm); Energy Deposit (MeV)", nbinsz, zmin, zmax, nbinsx, xmin, xmax);
  yz_view = new TH2D("TMS_Viewer_yz", "TMS viewer yz;z (mm); y (mm); Energy Deposit (MeV)", nbinsz, zmin, zmax, nbinsy, ymin, ymax);

  xz_view->SetMinimum(-0.01);
  yz_view->SetMinimum(-0.01);
  xz_view->SetMaximum(5);
  yz_view->SetMaximum(5);

  yz_view->GetZaxis()->SetTitleOffset(yz_view->GetZaxis()->GetTitleOffset()*1.2);
  xz_view->GetZaxis()->SetTitleOffset(xz_view->GetZaxis()->GetTitleOffset()*1.2);

  yz_view->GetYaxis()->SetTitleOffset(1.8);
  xz_view->GetYaxis()->SetTitleOffset(1.8);

  yz_view->GetYaxis()->SetMaxDigits(3);
  yz_view->GetXaxis()->SetMaxDigits(3);
  xz_view->GetYaxis()->SetMaxDigits(3);
  xz_view->GetXaxis()->SetMaxDigits(3);

  TGaxis::SetExponentOffset(0, -0.031, "x");

  // The canvas
  Canvas = new TCanvas("TMS_EventViewer", "TMS_EventViewer", 1024, 1024);
  Canvas->Divide(3);
  //Canvas->cd(1)->SetLeftMargin(Canvas->GetLeftMargin()*1.2);
  //Canvas->cd(1)->SetRightMargin(Canvas->GetRightMargin()*1.5);

  // Full view from inspecting all hits
  xz_box_Full = new TBox(TMS_Const::TMS_Thin_Start,
      -3485,
      TMS_Const::TMS_Thick_End,
      3485);
  xz_box_Full->SetLineColor(kGreen);
  xz_box_Full->SetFillStyle(0);

  yz_box_Full = new TBox(TMS_Const::TMS_Thin_Start,
      -2340,
      TMS_Const::TMS_Thick_End,
      870);
  yz_box_Full->SetLineColor(kGreen);
  yz_box_Full->SetFillStyle(0);

  // FV just taking 50 cm in from the full
  xz_box_FV = new TBox(xz_box_Full->GetX1(),
      xz_box_Full->GetY1()+500,
      xz_box_Full->GetX2()-500,
      xz_box_Full->GetY2()-500);
  xz_box_FV->SetLineColor(kRed);
  xz_box_FV->SetLineStyle(kDashed);
  xz_box_FV->SetFillStyle(0);

  yz_box_FV = new TBox(yz_box_Full->GetX1(),
      yz_box_Full->GetY1()+500,
      yz_box_Full->GetX2()-500,
      yz_box_Full->GetY2()-500);
  yz_box_FV->SetLineColor(kRed);
  yz_box_FV->SetLineStyle(kDashed);
  yz_box_FV->SetFillStyle(0);

  // Include the dead region boxes
  xz_dead_top = new TBox(TMS_Const::TMS_Thin_Start,
      TMS_Const::TMS_Dead_Top[0],
      TMS_Const::TMS_Thick_End,
      TMS_Const::TMS_Dead_Top[1]);
  xz_dead_center = new TBox(TMS_Const::TMS_Thin_Start,
      TMS_Const::TMS_Dead_Center[0],
      TMS_Const::TMS_Thick_End,
      TMS_Const::TMS_Dead_Center[0]);
  xz_dead_bottom = new TBox(TMS_Const::TMS_Thin_Start,
      TMS_Const::TMS_Dead_Bottom[0],
      TMS_Const::TMS_Thick_End,
      TMS_Const::TMS_Dead_Bottom[1]);
  xz_dead_top->SetFillStyle(3003);
  xz_dead_center->SetFillStyle(3003);
  xz_dead_bottom->SetFillStyle(3003);
  xz_dead_top->SetFillColor(kGray);
  xz_dead_center->SetFillColor(kGray);
  xz_dead_bottom->SetFillColor(kGray);
  
  // And a line at the thin/thick divide
  xz_Thin_Thick = new TLine(TMS_Const::TMS_Thick_Start,
      TMS_Const::TMS_Start_Exact[0],
      TMS_Const::TMS_Thick_Start,
      TMS_Const::TMS_End_Exact[0]);
  xz_Thin_Thick->SetLineColor(kGray);
  xz_Thin_Thick->SetLineStyle(kDashed);

  yz_Thin_Thick = new TLine(TMS_Const::TMS_Thick_Start,
      TMS_Const::TMS_Start_Exact[1],
      TMS_Const::TMS_Thick_Start,
      TMS_Const::TMS_End_Exact[1]);
  yz_Thin_Thick->SetLineColor(kGray);
  yz_Thin_Thick->SetLineStyle(kDashed);

  // Make the output file
  std::string filename = TMS_Manager::GetInstance().GetFileName();

  while (filename.find("/") != std::string::npos) {
    filename = filename.substr(filename.find("/")+1, filename.size());
  }
  CanvasName = filename.c_str();
  CanvasName.ReplaceAll(".root", "_TMS_EventViewer");
  gErrorIgnoreLevel = kWarning;
  CanvasName += "_"+TMS_Manager::GetInstance().Get_Reco_TrackMethod();
  CanvasName += Form("_Cluster%i", TMS_Manager::GetInstance().Get_Reco_Clustering());
  Canvas->Print(CanvasName+".pdf[");
  gErrorIgnoreLevel = kInfo;
}

// Draw the finished event
void TMS_EventViewer::Draw(TMS_Event &event) {

  xz_view->Reset();
  yz_view->Reset();

  int EventNumber = event.GetEventNumber();
  xz_view->SetTitle(Form("TMS viewer xz, Event %i", EventNumber));
  yz_view->SetTitle(Form("TMS viewer yz, Event %i", EventNumber));

  std::vector<TMS_Hit> TMS_Hits = event.GetHits();

  // Check that there are hits; overlaps with nMinHits in TMS_Reco
  // Make this a constant
  if (TMS_Hits.size() < 10 || TMS_TrackFinder::GetFinder().GetCleanedHits().size() < 10) {
    //std::cout << "Trying to draw an event that has no hits in the TMS, returning..." << std::endl;
    return;
  }

  // Loop over the hits and add them
  for (std::vector<TMS_Hit>::iterator it = TMS_Hits.begin(); it != TMS_Hits.end(); ++it) {
    TMS_Bar bar = (*it).GetBar();  
    double x = bar.GetX();
    double y = bar.GetY();
    double z = bar.GetZ();
    int BarType = bar.GetBarType();
    double e = (*it).GetE();

    // Bar along y (no x info)
    if (BarType == TMS_Bar::BarType::kYBar) {
      xz_view->Fill(z, x, e);
    }
    // Bar along x (no y info)
    else if (BarType == TMS_Bar::BarType::kXBar) {
      yz_view->Fill(z, y, e);
    }
  }

  // Get the true particle's trajectories
  std::vector<TMS_TrueParticle> traj = event.GetTrueParticles();
  int ntraj = traj.size();
  // Make a TGraph for each trajectory
  std::vector<TGraph*> trajgraphs(ntraj);
  // Loop over the trajectories
  int it = 0;
  for (auto i = traj.begin(); i != traj.end(); ++i,it++) {
    TGraph *tempgraph = new TGraph((*i).GetPositionPoints().size());
    int npoints = int(((*i).GetPositionPoints()).size());
    for (int j = 0; j < npoints; ++j) {
      tempgraph->SetPoint(j, (*i).GetPositionPoints()[j].Z(), (*i).GetPositionPoints()[j].X());
    }

    // Set a specific marker from primary particles
    tempgraph->SetMarkerStyle(24);
    tempgraph->SetMarkerSize(0.4);
    if ((*i).GetParent() == -1) {
      tempgraph->SetMarkerStyle(25);
      tempgraph->SetMarkerSize(1.0);
    } 

    if (abs((*i).GetPDG()) == 13) {
      tempgraph->SetMarkerColor(kYellow-3);
    } else if (abs((*i).GetPDG()) == 11) {
      tempgraph->SetMarkerColor(kGreen-2);
    } else if (abs((*i).GetPDG()) == 211) {
      tempgraph->SetMarkerColor(kRed-7);
    } else if (abs((*i).GetPDG()) == 2212) {
      tempgraph->SetMarkerColor(kMagenta-7);
    } else if (abs((*i).GetPDG()) == 2112) {
      tempgraph->SetMarkerColor(kCyan-3);
      tempgraph->SetMarkerSize(0.1);
      tempgraph->SetMarkerStyle(31);
    } else if (abs((*i).GetPDG()) == 22) {
      tempgraph->SetMarkerColor(kGray);
      tempgraph->SetMarkerSize(0.1);
      tempgraph->SetMarkerStyle(31);
    }

    trajgraphs[it] = tempgraph;
  }

  // Loop over the reconstructed tracks to overlay with hits
  //std::vector<TMS_Hit> Candidates = TMS_TrackFinder::GetFinder().GetCandidates();
  //std::vector<std::vector<TMS_Hit> > TotalCandidates = TMS_TrackFinder::GetFinder().GetTotalCandidates();

  // Get the Hough candidates
  std::vector<std::vector<TMS_Hit> > HoughCandidates = TMS_TrackFinder::GetFinder().GetHoughCandidates();
  // Now loop over the cluster candidates, make them TGraphs
  int nLines = HoughCandidates.size();
  std::vector<TGraph*> HoughGraphVect(nLines);
  for (int i = 0; i < nLines; ++i) {
    HoughGraphVect[i] = new TGraph(HoughCandidates[i].size());
    HoughGraphVect[i]->SetLineColor(i);
    HoughGraphVect[i]->SetLineWidth(2);
    HoughGraphVect[i]->SetMarkerColor(i);
    HoughGraphVect[i]->SetMarkerSize(1.5);
    HoughGraphVect[i]->SetMarkerStyle(25);
  }
  int LineIt = 0;
  for (auto Line: HoughCandidates) {
    int HitIt = 0;
    for (auto HoughHit: Line) {
      HoughGraphVect[LineIt]->SetPoint(HitIt, HoughHit.GetZ(), HoughHit.GetNotZ());
      HitIt++;
    }
    LineIt++;
  }
  
  // Get the cluster candidates
  std::vector<std::vector<TMS_Hit> > ClusterCandidates = TMS_TrackFinder::GetFinder().GetClusterCandidates();
  // Now loop over the cluster candidates, make them TGraphs
  int nClusters = ClusterCandidates.size();
  std::vector<TGraph*> GraphVect(nClusters);
  for (int i = 0; i < nClusters; ++i) {
    GraphVect[i] = new TGraph(ClusterCandidates[i].size());
    GraphVect[i]->SetLineColor(nLines+i);
    GraphVect[i]->SetLineWidth(2);
    GraphVect[i]->SetMarkerColor(nLines+i);
    GraphVect[i]->SetMarkerSize(1);
    GraphVect[i]->SetMarkerStyle(4);
  }
  int ClusterIt = 0;
  for (auto Cluster: ClusterCandidates) {
    int HitIt = 0;
    for (auto ClusterHit: Cluster) {
      GraphVect[ClusterIt]->SetPoint(HitIt, ClusterHit.GetZ(), ClusterHit.GetNotZ());
      HitIt++;
    }
    ClusterIt++;
  }

  // Get all the hough lines
  std::vector<std::pair<bool, TF1*> > HoughLines = TMS_TrackFinder::GetFinder().GetHoughLines();

  gStyle->SetPalette(kBird);

  // Draw the raw hits
  Canvas->cd(1); 
  xz_view->Draw("colz");
  xz_box_FV->Draw("same");
  xz_box_Full->Draw("same");
  xz_dead_top->Draw("same");
  xz_dead_center->Draw("same");
  xz_dead_bottom->Draw("same");
  xz_Thin_Thick->Draw("same");

  // Draw the true trajectories
  Canvas->cd(2); 
  xz_view->Draw("colz");
  xz_box_FV->Draw("same");
  xz_box_Full->Draw("same");
  xz_dead_top->Draw("same");
  xz_dead_center->Draw("same");
  xz_dead_bottom->Draw("same");
  xz_Thin_Thick->Draw("same");
  for (int i = 0; i < int(trajgraphs.size()); ++i) {
    if (trajgraphs[i] != NULL) trajgraphs[i]->Draw("P, same");
  }

  // Draw the reconstructed tracks 
  Canvas->cd(3);
  xz_view->Draw("colz");
  xz_box_FV->Draw("same");
  xz_box_Full->Draw("same");
  xz_dead_top->Draw("same");
  xz_dead_center->Draw("same");
  xz_dead_bottom->Draw("same");
  xz_Thin_Thick->Draw("same");
  for (auto graph: HoughGraphVect) graph->Draw("P,same");
  it = 0;
  for (auto i : HoughLines) {
    if (i.first == true) {
      i.second->SetLineColor(it);
      i.second->Draw("same");
      it++;
    }
  }
  for (auto graph: GraphVect) graph->Draw("P,same");

  xz_view->SetTitle(Form("#splitline{%s}{nLines: %i, nCluster: %i}", xz_view->GetTitle(), nLines, nClusters));

  /*
  Canvas->cd(2); 
  yz_view->Draw("colz");
  yz_box_FV->Draw("same");
  yz_box_Full->Draw("same");
  yz_Thin_Thick->Draw("same");
  for (auto i : HoughLines) if (i.first == false) i.second->Draw("same");
  */

  // ROOT is very verbose when printing, turn this off
  gErrorIgnoreLevel = kWarning;
  Canvas->Print(CanvasName+".pdf");

  if (DrawTrackFinding) {
    gStyle->SetPalette(87);
    Canvas->cd(1);
    TH2D *accumulator_xz = TMS_TrackFinder::GetFinder().AccumulatorToTH2D(false);
    accumulator_xz->Draw("colz");

    Canvas->cd(2);
    TH2D *accumulator_yz = TMS_TrackFinder::GetFinder().AccumulatorToTH2D(true);
    accumulator_yz->Draw("colz");
    Canvas->Print(CanvasName+".pdf");

    delete accumulator_xz;
    delete accumulator_yz;
  }

  for (int i = 0; i < nClusters; ++i) {
    delete GraphVect[i];
  }
  for (int i = 0; i < nLines; ++i) {
    delete HoughGraphVect[i];
  }

  for (int i = 0; i < int(trajgraphs.size()); ++i) {
    delete trajgraphs[i];
  }

  gErrorIgnoreLevel = kInfo;

  nDraws++;
}

