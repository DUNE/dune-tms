// python code was very slow...
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TGeoManager.h"
#include "TStopwatch.h"
#include "TLorentzVector.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGeoManager.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TBox.h"

// EDepSim includes
#include "EDepSim/TG4Event.h"

#include "TMS_Event.h"
#include "TMS_Constants.h"
#include "TMS_Reco.h" // Needed for hit cleaning

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cout << "Need one argument, the input filename" << std::endl;
    return -1;
  }

  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(255);
  gStyle->SetPalette(89);

  bool DrawTraj = true;
  bool DrawRecoHits = true;

  std::string filename = std::string(argv[1]);

  std::cout << "Got " << filename << std::endl;

  // The input file
  TFile *input = new TFile(filename.c_str(), "open");
  // The EDepSim events
  TTree *events = (TTree*)input->Get("EDepSimEvents");
  // Get the detector geometry
  TGeoManager *geom = (TGeoManager*)input->Get("EDepSimGeometry");

  // The global manager
  TMS_Manager::GetInstance().SetFileName(filename);

  // Load up the geometry
  TMS_Geom::GetInstance().SetGeometry(geom);
  // The generator pass-through information
  TTree *gRoo = (TTree*)input->Get("DetSimPassThru/gRooTracker");
  // Get the true neutrino vector from the gRooTracker object
  int StdHepPdg[__EDEP_SIM_MAX_PART__];
  double StdHepP4[__EDEP_SIM_MAX_PART__][4];
  gRoo->SetBranchStatus("*", false);
  gRoo->SetBranchStatus("StdHepPdg", true);
  gRoo->SetBranchStatus("StdHepP4", true);
  gRoo->SetBranchAddress("StdHepPdg", StdHepPdg);
  gRoo->SetBranchAddress("StdHepP4", StdHepP4);

  // Get the event
  TG4Event *event = NULL;
  events->SetBranchAddress("Event", &event);

  int N_entries = events->GetEntries();
  N_entries = 100;
  const double zmin = 3000;
  const double zmax = 19000;
  const double xmin = -4000;
  const double xmax = 4000;

  // Make some boxes for the fiducial volumes
  // Full view from inspecting all hits
  TBox *xz_box_Full = new TBox(TMS_Const::TMS_Thin_Start/1E3,
      -3485/1E3,
      TMS_Const::TMS_Thick_End/1E3,
      3485/1E3);
  xz_box_Full->SetLineColor(kMagenta+2);
  xz_box_Full->SetFillStyle(0);

  TBox *xz_box_LAr_Full = new TBox(TMS_Const::LAr_Start_Exact[2]/1E3,
      TMS_Const::LAr_Start_Exact[0]/1E3,
      TMS_Const::LAr_End_Exact[2]/1E3,
      TMS_Const::LAr_End_Exact[0]/1E3);
  xz_box_LAr_Full->SetLineColor(kMagenta+2);
  xz_box_LAr_Full->SetFillStyle(0);

  // FV just taking 50 cm in from the full
  TBox *xz_box_FV = new TBox(xz_box_Full->GetX1(),
      (xz_box_Full->GetY1()+500/1E3),
      (xz_box_Full->GetX2()-500/1E3),
      (xz_box_Full->GetY2()-500/1E3));
  xz_box_FV->SetLineColor(kRed);
  xz_box_FV->SetLineStyle(kDashed);
  xz_box_FV->SetFillStyle(0);

  TBox *xz_box_LAr_FV = new TBox(xz_box_LAr_Full->GetX1()+500/1E3,
      xz_box_LAr_Full->GetY1()+500/1E3,
      xz_box_LAr_Full->GetX2()-500/1E3,
      xz_box_LAr_Full->GetY2()-500/1E3);
  xz_box_LAr_FV->SetLineColor(kRed);
  xz_box_LAr_FV->SetLineStyle(kDashed);
  xz_box_LAr_FV->SetFillStyle(0);


  // Include the dead region boxes
  TBox *xz_dead_top = new TBox(TMS_Const::TMS_Thin_Start/1E3,
      TMS_Const::TMS_Dead_Top[0]/1E3,
      TMS_Const::TMS_Thick_End/1E3,
      TMS_Const::TMS_Dead_Top[1]/1E3);
  TBox *xz_dead_center = new TBox(TMS_Const::TMS_Thin_Start/1E3,
      TMS_Const::TMS_Dead_Center[0]/1E3,
      TMS_Const::TMS_Thick_End/1E3,
      TMS_Const::TMS_Dead_Center[0]/1E3);
  TBox *xz_dead_bottom = new TBox(TMS_Const::TMS_Thin_Start/1E3,
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
  TLine *xz_Thin_Thick = new TLine(TMS_Const::TMS_Thick_Start/1E3,
      TMS_Const::TMS_Start_Exact[0]/1E3,
      TMS_Const::TMS_Thick_Start/1E3,
      TMS_Const::TMS_End_Exact[0]/1E3);
  xz_Thin_Thick->SetLineColor(kGray);
  xz_Thin_Thick->SetLineStyle(kDashed);

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->SetLeftMargin(canv->GetLeftMargin()*1.4);
  canv->SetRightMargin(canv->GetRightMargin()*1.5);
  canv->SetTopMargin(canv->GetTopMargin()*1.2);
  // Make the output name
  std::string outputname = filename;
  while (outputname.find("/") != std::string::npos) {
    outputname = outputname.substr(outputname.find("/")+1, outputname.size());
  }
  TString Outputname = outputname.c_str();
  Outputname.ReplaceAll(".root", "_Deposits.pdf");
  canv->Print(Outputname+"[");

  TString Outputname_root = outputname.c_str();
  Outputname_root.ReplaceAll(".root", "_Deposits.root");
  TFile *outfile = new TFile(Outputname_root, "recreate");
  TTree *outtree = new TTree("energy_dep", "energy_dep");
  float pos[4];
  float energy;
  outtree->Branch("Position", pos, "Position[4]/F");
  outtree->Branch("Energy", &energy, "Energy/F");

  TPaveText *text = new TPaveText(0,0,1,1,"NDC");
  text->AddText(Form("Exec: %s", argv[0]));
  text->AddText(Form("Input: %s", filename.c_str()));
  text->AddText(Form("Output: %s", Outputname.Data()));
  text->SetBorderSize(0);
  text->SetFillColor(0);
  text->Draw();
  canv->Print(Outputname);

  gErrorIgnoreLevel = kWarning;
  int i = 0;
  TH2D *plot = new TH2D("xz", "xz-view;z (m); x (m); Energy deposit (MeV)", 300, zmin/1E3, zmax/1E3, 200, xmin/1E3, xmax/1E3);
  plot->GetZaxis()->SetTitleOffset(plot->GetZaxis()->GetTitleOffset()*1.3);

  TString zaxistitle = plot->GetZaxis()->GetTitle();
  for (; i < N_entries; ++i) {
    events->GetEntry(i);
    gRoo->GetEntry(i);

    plot->Reset();
    if (i % (N_entries/10) == 0) {
      std::cout << "Processed " << i << "/" << N_entries << " (" << double(i)*100./N_entries << "%)" << std::endl;
    }

    TMS_Event tms_event = TMS_Event(*event, true);
    tms_event.FillTruthFromGRooTracker(StdHepPdg, StdHepP4);

    int pdg = tms_event.GetNeutrinoPDG();
    double enu = tms_event.GetNeutrinoP4().E();
    TString reaction = tms_event.GetReaction().c_str();
    reaction.ReplaceAll(";",",");
    plot->SetTitle(Form("#splitline{Event %i, all true hits, #nu PDG: %i, E_{#nu}=%.2f GeV}{%s}", i, pdg, enu, reaction.Data()));

    std::vector<TMS_Hit> TMS_Hits = tms_event.GetHits();
    // First draw the hits
    for (auto i = TMS_Hits.begin(); i != TMS_Hits.end(); ++i) {
      double x = (*i).GetTrueHit().GetX();
      double y = (*i).GetTrueHit().GetY();
      double z = (*i).GetTrueHit().GetZ();
      double t = (*i).GetTrueHit().GetT();
      double e = (*i).GetTrueHit().GetE();
      plot->Fill(z/1E3, x/1E3, e);
      pos[0] = x;
      pos[1] = y;
      pos[2] = z;
      pos[3] = t;
      energy = e;
      outtree->Fill();
    }

    // Get the true particle's trajectories
    std::vector<TMS_TrueParticle> traj = tms_event.GetTrueParticles();
    int ntraj = traj.size();
    // Make a TGraph for each trajectory
    std::vector<TGraph*> trajgraphs(ntraj);
    // Loop over the trajectories
    int it = 0;
    int muonindex = 0;
    for (auto i = traj.begin(); i != traj.end(); ++i,it++) {
      TGraph *tempgraph = new TGraph((*i).GetPositionPoints().size());
      int npoints = int(((*i).GetPositionPoints()).size());
      for (int j = 0; j < npoints; ++j) {
        tempgraph->SetPoint(j, (*i).GetPositionPoints()[j].Z()/1E3, (*i).GetPositionPoints()[j].X()/1E3);
      }

      // Set a specific marker from primary particles
      tempgraph->SetMarkerStyle(24);
      tempgraph->SetMarkerSize(1.0);
      if ((*i).GetParent() == -1) {
        tempgraph->SetMarkerStyle(25);
        tempgraph->SetMarkerSize(1.0);
      } 

      if (abs((*i).GetPDG()) == 13) {
        tempgraph->SetMarkerColor(kYellow-3);
        muonindex = it;
      } else if (abs((*i).GetPDG()) == 11) {
        tempgraph->SetMarkerColor(kGreen-2);
      } else if (abs((*i).GetPDG()) == 211) {
        tempgraph->SetMarkerColor(kRed-7);
      } else if (abs((*i).GetPDG()) == 2212) {
        tempgraph->SetMarkerColor(kMagenta-7);
      } else if (abs((*i).GetPDG()) == 2112) {
        tempgraph->SetMarkerColor(kCyan-3);
        tempgraph->SetMarkerSize(0.0);
        tempgraph->SetMarkerStyle(26);
      } else if (abs((*i).GetPDG()) == 22) {
        tempgraph->SetMarkerColor(kGray);
        tempgraph->SetMarkerSize(0.0);
        tempgraph->SetMarkerStyle(27);
      }

      trajgraphs[it] = tempgraph;
    }

    // Now make the TGraphs for the true trajectory points

    if (plot->Integral() == 0) continue;

    plot->SetMinimum(0);
    plot->SetMaximum(2);

    plot->GetZaxis()->SetTitle(zaxistitle);
    plot->Draw("colz");
    xz_box_FV->Draw("same");
    xz_box_LAr_FV->Draw("same");
    xz_box_Full->Draw("same");
    xz_box_LAr_Full->Draw("same");
    xz_dead_top->Draw("same");
    xz_dead_center->Draw("same");
    xz_dead_bottom->Draw("same");
    xz_Thin_Thick->Draw("same");

    canv->Print(Outputname);

    // Go from back to front to have muon draw on top
    if (DrawTraj) {
      for (int i = trajgraphs.size()-1; i >= 0; --i) {
        if (trajgraphs[i] != NULL) trajgraphs[i]->Draw("P, same");
      }
      // Draw muons on top
      trajgraphs[muonindex]->Draw("P,same");
      TString plottitle = plot->GetTitle();
      plottitle.ReplaceAll("all true hits", "all true hits and traj.");
      plot->SetTitle(plottitle);
      canv->Print(Outputname);
    }

    if (DrawRecoHits) {
      // Clean up hits
      std::vector<TMS_Hit> TMS_Hits_Cleaned = TMS_TrackFinder::GetFinder().CleanHits(TMS_Hits);
      plot->Reset();
      // First draw the hits
      for (auto i = TMS_Hits_Cleaned.begin(); i != TMS_Hits_Cleaned.end(); ++i) {
        double x = (*i).GetNotZ();
        double z = (*i).GetZ();
        plot->Fill(z/1E3, x/1E3);
      }
      plot->Draw("colz");
      xz_box_FV->Draw("same");
      xz_box_LAr_FV->Draw("same");
      xz_box_Full->Draw("same");
      xz_box_LAr_Full->Draw("same");
      xz_dead_top->Draw("same");
      xz_dead_center->Draw("same");
      xz_dead_bottom->Draw("same");
      xz_Thin_Thick->Draw("same");

      TString plottitle = plot->GetTitle();
      plottitle.ReplaceAll("all true hits and traj.", "all binary reco hits");
      plot->SetTitle(plottitle);
      // Also change the z-axis
      plot->GetZaxis()->SetTitle("Binary hit");
      plot->SetMaximum(1);

      canv->Print(Outputname);

      if (DrawTraj) {
      for (int i = trajgraphs.size()-1; i >= 0; --i) {
        if (trajgraphs[i] != NULL) trajgraphs[i]->Draw("P, same");
      }
      // Draw muons on top
      trajgraphs[muonindex]->Draw("P,same");
      TString plottitle = plot->GetTitle();
      plottitle.ReplaceAll("all binary reco hits", "all binary reco hits and true traj.");
      plot->SetTitle(plottitle);
      canv->Print(Outputname);
      }
    }

    for (int i = 0; i < int(trajgraphs.size()); ++i) {
      delete trajgraphs[i];
    }
    outfile->cd();
    canv->Write(Form("Event_%i", i));
  }

  outfile->cd();
  outtree->Write();
  outfile->Close();

  input->Close();
  std::cout << "Wrote pdf file to " << Outputname << std::endl;
  std::cout << "Wrote root file to " << Outputname_root << std::endl;

  canv->Print(Outputname+"]");
}
