#include "TMS_EventViewer.h"

TMS_EventViewer::TMS_EventViewer() :
DrawTrackFinding(false) {

  nDraws = 0;

  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(255);

  const double zmin = TMS_Const::TMS_Start[2]/1.E3;
  const double zmax = TMS_Const::TMS_End[2]/1.E3;
  const double xmin = TMS_Const::TMS_Start[0]/1.E3;
  const double xmax = TMS_Const::TMS_End[0]/1.E3;
  const double ymin = TMS_Const::TMS_Start[1]/1.E3;
  const double ymax = TMS_Const::TMS_End[1]/1.E3;
  // Scint bars are 4 by 1 cm
  //const int nbinsz = ((zmax-zmin)/10)/5;
  const int nbinsz = ((zmax-zmin)/0.1)*2;
  const int nbinsx = ((xmax-xmin)/0.1)*3;
  const int nbinsy = ((ymax-ymin)/0.1)*3;

  // The 2D views
  xz_view = new TH2D("TMS_Viewer_xz", "TMS viewer xz;z (m); x (m); Energy Deposited (MeV)", nbinsz, zmin, zmax, nbinsx, xmin, xmax);
  yz_view = new TH2D("TMS_Viewer_yz", "TMS viewer yz;z (m); y (m); Energy Deposited (MeV)", nbinsz, zmin, zmax, nbinsy, ymin, ymax);

  xz_view->SetMinimum(0);
  yz_view->SetMinimum(0);
  xz_view->SetMaximum(3);
  yz_view->SetMaximum(3);

  yz_view->GetZaxis()->SetTitleOffset(yz_view->GetZaxis()->GetTitleOffset()*1.2);
  xz_view->GetZaxis()->SetTitleOffset(xz_view->GetZaxis()->GetTitleOffset()*1.2);

  yz_view->GetYaxis()->SetTitleOffset(1.1);
  xz_view->GetYaxis()->SetTitleOffset(1.1);

  yz_view->GetYaxis()->SetMaxDigits(3);
  yz_view->GetXaxis()->SetMaxDigits(3);
  xz_view->GetYaxis()->SetMaxDigits(3);
  xz_view->GetXaxis()->SetMaxDigits(3);

  TGaxis::SetExponentOffset(0, -0.031, "x");

  // The canvas
  Canvas = new TCanvas("TMS_EventViewer", "TMS_EventViewer", 1024, 1024);
  Canvas->SetTopMargin(Canvas->GetTopMargin()*1.4);
  Canvas->SetRightMargin(Canvas->GetRightMargin()*1.4);

  // Full view from inspecting all hits
  xz_box_Full = new TBox(TMS_Const::TMS_Thin_Start/1E3,
      -3485/1E3,
      TMS_Const::TMS_Thick_End/1E3,
      3485/1E3);
  xz_box_Full->SetLineColor(kGreen);
  xz_box_Full->SetFillStyle(0);

  yz_box_Full = new TBox(TMS_Const::TMS_Thin_Start/1E3,
      -2340/1E3,
      TMS_Const::TMS_Thick_End/1E3,
      870/1E3);
  yz_box_Full->SetLineColor(kGreen);
  yz_box_Full->SetFillStyle(0);

  // FV just taking 50 cm in from the full
  xz_box_FV = new TBox(xz_box_Full->GetX1(),
      xz_box_Full->GetY1()+500/1E3,
      xz_box_Full->GetX2()-500/1E3,
      xz_box_Full->GetY2()-500/1E3);
  xz_box_FV->SetLineColor(kRed);
  xz_box_FV->SetLineStyle(kDashed);
  xz_box_FV->SetFillStyle(0);

  yz_box_FV = new TBox(yz_box_Full->GetX1(),
      yz_box_Full->GetY1()+500/1E3,
      yz_box_Full->GetX2()-500/1E3,
      yz_box_Full->GetY2()-500/1E3);
  yz_box_FV->SetLineColor(kRed);
  yz_box_FV->SetLineStyle(kDashed);
  yz_box_FV->SetFillStyle(0);

  // Include the dead region boxes
  xz_dead_top = new TBox(TMS_Const::TMS_Thin_Start/1E3,
      TMS_Const::TMS_Dead_Top[0]/1E3,
      TMS_Const::TMS_Thick_End/1E3,
      TMS_Const::TMS_Dead_Top[1]/1E3);
  xz_dead_center = new TBox(TMS_Const::TMS_Thin_Start/1E3,
      TMS_Const::TMS_Dead_Center[0]/1E3,
      TMS_Const::TMS_Thick_End/1E3,
      TMS_Const::TMS_Dead_Center[0]/1E3);
  xz_dead_bottom = new TBox(TMS_Const::TMS_Thin_Start/1E3,
      TMS_Const::TMS_Dead_Bottom[0]/1E3,
      TMS_Const::TMS_Thick_End/1E3,
      TMS_Const::TMS_Dead_Bottom[1]/1E3);
  xz_dead_top->SetFillStyle(3003);
  xz_dead_center->SetFillStyle(3003);
  xz_dead_bottom->SetFillStyle(3003);
  xz_dead_top->SetFillColor(kGray);
  xz_dead_center->SetFillColor(kGray);
  xz_dead_bottom->SetFillColor(kGray);
  
  // And a line at the thin/thick divide
  xz_Thin_Thick = new TLine(TMS_Const::TMS_Thick_Start/1E3,
      TMS_Const::TMS_Start_Exact[0]/1E3,
      TMS_Const::TMS_Thick_Start/1E3,
      TMS_Const::TMS_End_Exact[0]/1E3);
  xz_Thin_Thick->SetLineColor(kGray);
  xz_Thin_Thick->SetLineStyle(kDashed);

  yz_Thin_Thick = new TLine(TMS_Const::TMS_Thick_Start/1E3,
      TMS_Const::TMS_Start_Exact[1]/1E3,
      TMS_Const::TMS_Thick_Start/1E3,
      TMS_Const::TMS_End_Exact[1]/1E3);
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
  xz_view->SetTitle(Form("Event %i", EventNumber));
  yz_view->SetTitle(Form("Event %i", EventNumber));

  std::vector<TMS_Hit> TMS_Hits = TMS_TrackFinder::GetFinder().GetCleanedHits();

  // Check that there are hits; overlaps with nMinHits in TMS_Reco
  // Make this a constant
  if (TMS_Hits.size() < 1 || TMS_TrackFinder::GetFinder().GetCleanedHits().size() < 1) {
    //std::cout << "Trying to draw an event that has no hits in the TMS, returning..." << std::endl;
    return;
  }

  // Loop over the hits and add them
  for (std::vector<TMS_Hit>::iterator it = TMS_Hits.begin(); it != TMS_Hits.end(); ++it) {
    TMS_Bar bar = (*it).GetBar();  
    double x = bar.GetX()/1E3;
    double y = bar.GetY()/1E3;
    double z = bar.GetZ()/1E3;
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
  int truemuon_traj = -1;
  for (auto i = traj.begin(); i != traj.end(); ++i,it++) {
    TGraph *tempgraph = new TGraph((*i).GetPositionPoints().size());
    int npoints = int(((*i).GetPositionPoints()).size());
    for (int j = 0; j < npoints; ++j) {
      tempgraph->SetPoint(j, (*i).GetPositionPoints()[j].Z()/1E3, (*i).GetPositionPoints()[j].X()/1E3);
    }

    // Set a specific marker from primary particles
    tempgraph->SetMarkerStyle(24);
    tempgraph->SetMarkerSize(0.4);
    if ((*i).GetParent() == -1) {
      tempgraph->SetMarkerStyle(25);
      tempgraph->SetMarkerSize(1.5);
    } 

    if (abs((*i).GetPDG()) == 13) {
      truemuon_traj = it;
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

  // Get the Hough candidates
  std::vector<std::vector<TMS_Hit> > HoughCandidatesOne = TMS_TrackFinder::GetFinder().GetHoughCandidatesOne();
  std::vector<std::vector<TMS_Hit> > HoughCandidatesOther = TMS_TrackFinder::GetFinder().GetHoughCandidatesOther();
  // Now loop over the cluster candidates, make them TGraphs
  int nLinesOne = HoughCandidatesOne.size();
  int nLinesOther = HoughCandidatesOther.size();
  std::vector<TGraph*> HoughGraphVectOne(nLinesOne);
  std::vector<TGraph*> HoughGraphVectOther(nLinesOther);
  for (int i = 0; i < nLinesOne; ++i) {
    HoughGraphVectOne[i] = new TGraph(HoughCandidatesOne[i].size());
    HoughGraphVectOne[i]->SetLineColor(i);
    HoughGraphVectOne[i]->SetLineWidth(2);
    HoughGraphVectOne[i]->SetMarkerColor(i+1);
    HoughGraphVectOne[i]->SetMarkerSize(1.5);
    HoughGraphVectOne[i]->SetMarkerStyle(25);
  }
  for (int i = 0; i < nLinesOther; ++i) {
    HoughGraphVectOther[i] = new TGraph(HoughCandidatesOther[i].size());
    HoughGraphVectOther[i]->SetLineColor(i);
    HoughGraphVectOther[i]->SetLineWidth(2);
    HoughGraphVectOther[i]->SetMarkerColor(i+1);
    HoughGraphVectOther[i]->SetMarkerSize(1.5);
    HoughGraphVectOther[i]->SetMarkerStyle(25);
  }
  int LineIt = 0;
  for (auto &Line: HoughCandidatesOne) {
    int HitIt = 0;
    for (auto HoughHit: Line) {
      HoughGraphVectOne[LineIt]->SetPoint(HitIt, HoughHit.GetZ()/1E3, HoughHit.GetNotZ()/1E3);
      HitIt++;
    }
    LineIt++;
  }
  LineIt = 0;
  for (auto &Line: HoughCandidatesOther) {
    int HitIt = 0;
    for (auto HoughHit: Line) {
      HoughGraphVectOther[LineIt]->SetPoint(HitIt, HoughHit.GetZ()/1E3, HoughHit.GetNotZ()/1E3);
      HitIt++;
    }
    LineIt++;
  }
  
  // Get the cluster candidates
  std::vector<std::vector<TMS_Hit> > ClusterCandidatesOne = TMS_TrackFinder::GetFinder().GetClusterCandidatesOne();
  std::vector<std::vector<TMS_Hit> > ClusterCandidatesOther = TMS_TrackFinder::GetFinder().GetClusterCandidatesOther();
  // Now loop over the cluster candidates, make them TGraphs
  int nClustersOne = ClusterCandidatesOne.size();
  int nClustersOther = ClusterCandidatesOther.size();
  std::vector<TGraph*> GraphVectOne(nClustersOne);
  std::vector<TGraph*> GraphVectOther(nClustersOther);
  for (int i = 0; i < nClustersOne; ++i) {
    GraphVectOne[i] = new TGraph(ClusterCandidatesOne[i].size());
    GraphVectOne[i]->SetLineColor(nLinesOne+i);
    GraphVectOne[i]->SetLineWidth(2);
    GraphVectOne[i]->SetMarkerColor(nLinesOne+i);
    GraphVectOne[i]->SetMarkerSize(1.2);
    GraphVectOne[i]->SetMarkerStyle(4);
  }
  for (int i = 0; i < nClustersOther; ++i) {
    GraphVectOther[i] = new TGraph(ClusterCandidatesOther[i].size());
    GraphVectOther[i]->SetLineColor(nLinesOther+i);
    GraphVectOther[i]->SetLineWidth(2);
    GraphVectOther[i]->SetMarkerColor(nLinesOther+i);
    GraphVectOther[i]->SetMarkerSize(1.2);
    GraphVectOther[i]->SetMarkerStyle(4);
  }
  int ClusterIt = 0;
  for (auto &Cluster: ClusterCandidatesOne) {
    int HitIt = 0;
    for (auto ClusterHit: Cluster) {
      GraphVectOne[ClusterIt]->SetPoint(HitIt, ClusterHit.GetZ()/1E3, ClusterHit.GetNotZ()/1E3);
      HitIt++;
    }
    ClusterIt++;
  }
  ClusterIt = 0;
  for (auto &Cluster: ClusterCandidatesOther) {
    int HitIt = 0;
    for (auto ClusterHit: Cluster) {
      GraphVectOther[ClusterIt]->SetPoint(HitIt, ClusterHit.GetZ()/1E3, ClusterHit.GetNotZ()/1E3);
      HitIt++;
    }
    ClusterIt++;
  }
  

  // Get all the hough lines
  std::vector<std::pair<bool, TF1*> > HoughLinesOne = TMS_TrackFinder::GetFinder().GetHoughLinesOne();
  std::vector<std::pair<bool, TF1*> > HoughLinesOther = TMS_TrackFinder::GetFinder().GetHoughLinesOther();

  int pdg = event.GetNeutrinoPDG();
  double enu = event.GetNeutrinoP4().E();
  TString reaction = event.GetReaction().c_str();
  // Get the true muon info
  double muonke = event.GetMuonTrueKE();
  double Emu = muonke + 106.5;
  Emu *= 1E-3;

  reaction.ReplaceAll(";",",");
  xz_view->SetTitle(Form("#splitline{Event %i, all reco hits, #nu^{true} PDG: %i, E^{true}_{#nu}=%.2f GeV, E^{true}_{#mu}=%.2f GeV}{%s}", EventNumber, pdg, enu, Emu, reaction.Data()));

  gStyle->SetPalette(89);
  xz_view->Draw("colz");
  xz_box_FV->Draw("same");
  xz_box_Full->Draw("same");
  xz_dead_top->Draw("same");
  xz_dead_center->Draw("same");
  xz_dead_bottom->Draw("same");
  xz_Thin_Thick->Draw("same");
  gErrorIgnoreLevel = kWarning;
  Canvas->Print(CanvasName+".pdf");

  // Draw the true trajectories
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
  // Draw the muon trajectory on top
  if (truemuon_traj >= 0) trajgraphs[truemuon_traj]->Draw("P,same");
  Canvas->Print(CanvasName+".pdf");

  // Draw the reconstructed tracks 
  xz_view->SetTitle(Form("#splitline{Event %i reconstructed}{nLinesOne: %i, nLinesOther: %i, nClustersOne: %i, nClustersOther: %i}", EventNumber, nLinesOne, nLinesOther, nClustersOne, nClustersOther));
  xz_view->Draw("colz");
  xz_box_FV->Draw("same");
  xz_box_Full->Draw("same");
  xz_dead_top->Draw("same");
  xz_dead_center->Draw("same");
  xz_dead_bottom->Draw("same");
  xz_Thin_Thick->Draw("same");
  // Draw the clusters first
  for (auto &graph: GraphVectOne) graph->Draw("P,same");
  for (auto &graph: GraphVectOther) graph->Draw("P,same");
  // Then the track candidates
  for (auto &graph: HoughGraphVectOne) graph->Draw("P,same");
  for (auto &graph: HoughGraphVectOther) graph->Draw("P,same");
  it = 0;
  for (auto &i : HoughLinesOne) {
    if (i.first == true) {
      i.second->SetLineColor(it+1);
      i.second->Draw("same");
      it++;
    }
  }
  for (auto &i : HoughLinesOther) {
    if (i.first == true) {
      i.second->SetLineColor(it+1);
      i.second->Draw("same");
      it++;
    }
  }
  Canvas->Print(CanvasName+".pdf");

  if (DrawTrackFinding) {
    gStyle->SetPalette(87);
    TH2D *accumulator_xz = TMS_TrackFinder::GetFinder().AccumulatorToTH2D(false);
    accumulator_xz->RebinX(100);
    accumulator_xz->RebinY(100);
    accumulator_xz->SetMinimum(-0.01);
    accumulator_xz->Draw("colz");
    Canvas->Print(CanvasName+".pdf");

    //Canvas->cd(2);
    //TH2D *accumulator_yz = TMS_TrackFinder::GetFinder().AccumulatorToTH2D(true);
    //accumulator_yz->Draw("colz");
    //Canvas->Print(CanvasName+".pdf");

    delete accumulator_xz;
    //delete accumulator_yz;
  }

  for (int i = 0; i < nClustersOne; ++i) {
    delete GraphVectOne[i];
  }
  for (int i = 0; i < nClustersOther; ++i) {
    delete GraphVectOther[i];
  }
  for (int i = 0; i < nLinesOne; ++i) {
    delete HoughGraphVectOne[i];
  }
  for (int i = 0; i < nLinesOther; ++i) {
    delete HoughGraphVectOther[i];
  }

  for (int i = 0; i < int(trajgraphs.size()); ++i) {
    delete trajgraphs[i];
  }

  gErrorIgnoreLevel = kInfo;

  nDraws++;
}

