#include <iostream>
#include <vector>
#include <algorithm>

#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraph.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TBranchRef.h"

#include "TApplication.h"

void draw_1d(std::string filename) {
  gStyle->SetPalette(kBird);
  gErrorIgnoreLevel = kError;

  TFile *file = new TFile(filename.c_str());
  TTree *tree = (TTree*)file->Get("tree");

  int nEntries = tree->GetEntries();

  // plot the xpt and zpt
  TH2D *xzTraj = new TH2D("xz", "xz;z (cm);x (cm)", 200, -50, 1500, 200, -400, 400);
  xzTraj->SetMinimum(-0.001);
  xzTraj->SetMaximum(20);

  std::vector<float> *xpt = NULL;
  std::vector<float> *zpt = NULL;
  //std::vector<float> *xptLep = NULL;
  //std::vector<float> *zptLep = NULL;
  char *reac = NULL;
  float vtx[3];
  int ievt;
  int muonReco;
  float Enu;
  float lepE;
  int lepPdg;
  int PDGnu;
  float muonDeath[3];
  float muonBirth[3];
  float muonExitPt[3];

  tree->SetBranchStatus("*", false);
  tree->SetBranchStatus("vtx", true);
  tree->SetBranchAddress("vtx", vtx);
  tree->SetBranchStatus("lepDeath", true);
  tree->SetBranchAddress("lepDeath", muonDeath);
  tree->SetBranchStatus("muonBirth", true);
  tree->SetBranchAddress("muonBirth", muonBirth);
  tree->SetBranchStatus("muonExitPt", true);
  tree->SetBranchAddress("muonExitPt", muonExitPt);
  tree->SetBranchStatus("ievt", true);
  tree->SetBranchAddress("ievt", &ievt);
  tree->SetBranchStatus("xpt", true);
  tree->SetBranchAddress("xpt", &xpt);
  tree->SetBranchStatus("zpt", true);
  tree->SetBranchAddress("zpt", &zpt);
  tree->SetBranchStatus("reac", true);
  tree->SetBranchAddress("reac", &reac);
  tree->SetBranchStatus("Ev", true);
  tree->SetBranchAddress("Ev", &Enu);
  tree->SetBranchStatus("lepE", true);
  tree->SetBranchAddress("lepE", &lepE);
  tree->SetBranchStatus("lepPdg", true);
  tree->SetBranchAddress("lepPdg", &lepPdg);
  tree->SetBranchStatus("PDGv", true);
  tree->SetBranchAddress("PDGv", &PDGnu);

  tree->SetBranchStatus("muonReco", true);
  tree->SetBranchAddress("muonReco", &muonReco);

  /*
  tree->SetBranchStatus("xptLep", true);
  tree->SetBranchAddress("xptLep", &xptLep);
  tree->SetBranchStatus("zptLep", true);
  tree->SetBranchAddress("zptLep", &zptLep);
  */

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  //canv->Divide(2,2);
  TPad *p1 = new TPad("p1", "p1", 0.0, 0.0, 1.0, 0.9);
  p1->Draw();

  p1->SetBottomMargin(0.10);
  p1->SetTopMargin(0.08);
  p1->SetRightMargin(0.13);

  canv->SetRightMargin(canv->GetRightMargin()*1.2);
  TString canvname = filename.c_str();
  canvname.ReplaceAll(".root", "");
  canvname+="_evtdisplay_full_xzonly_allhits";
  canv->Print(canvname+".pdf[");

  // Make a TGraph in xz and yz of the vertex location
  TGraph *xzVert = new TGraph(1);
  xzVert->SetMarkerStyle(24);
  xzVert->SetMarkerSize(2);
  xzVert->SetMarkerColor(kGreen);

  TGraph *xzDeath = new TGraph(1);
  xzDeath->SetMarkerStyle(25);
  xzDeath->SetMarkerSize(2);
  xzDeath->SetMarkerColor(kRed);

  TGraph *xzBirth = new TGraph(1);
  xzBirth->SetMarkerStyle(25);
  xzBirth->SetMarkerSize(2);
  xzBirth->SetMarkerColor(kGreen);

  TGraph *xzExit = new TGraph(1);
  xzExit->SetMarkerStyle(25);
  xzExit->SetMarkerSize(2);
  xzExit->SetMarkerColor(kRed);

  bool reset = false;
  int nBad = 0;
  double gev_cut = 5;
  int nPass = 0;
  for (int i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);

    //if (i > 1000) break;
    if (i%1000==0) {
      std::cout << "Entry " << i << std::endl;
    }

    // Only 100 events
    //if (nPass > 1000) break;
    if (muonReco != 2) continue;

    if (abs(vtx[0]) > 300 || abs(vtx[1]) > 100 || vtx[2] < 50 || vtx[2] > 350) continue;

    if (fabs(Enu-lepE/1E3) > gev_cut) nBad++;
    nPass++;

    int nxpoints = xpt->size();
    //int nxpointsLep = xptLep->size();

    xzTraj->Reset();

    for (int j = 0; j < nxpoints; ++j) xzTraj->Fill((*zpt)[j], (*xpt)[j]);
    //for (int j = 0; j < nxpointsLep; ++j) xzTraj->Fill((*zptLep)[j], (*xptLep)[j]);

    xzVert->SetPoint(0, vtx[2], vtx[0]);

    xzDeath->SetPoint(0, muonDeath[2], muonDeath[0]);

    xzBirth->SetPoint(0, muonBirth[2], muonBirth[0]);

    xzExit->SetPoint(0, muonExitPt[2], muonExitPt[0]);

    // Build a string depending on where the muon is reconstructed
    std::string Where;
    if (muonReco == 0) Where = "Not cont.";
    if (muonReco == 1) Where = "LAr cont.";
    if (muonReco == 2) Where = "TMS cont.";

    std::string flav;
    if (lepPdg == 11) flav = "e^{-}";
    if (lepPdg == -11) flav = "e^{+}";
    if (lepPdg == 13) flav = "#mu^{-}";
    if (lepPdg == -13) flav = "#mu^{+}";

    std::string nuflav;
    if (PDGnu == 12) nuflav = "#nu_{e}";
    if (PDGnu == -12) nuflav = "#bar{#nu_{e}}";
    if (PDGnu == 14) nuflav = "#nu_{#mu}";
    if (PDGnu == -14) nuflav = "#bar{#nu_{#mu}}";

    // Draw lines around fiducial volume
    TBox *box3 = new TBox(50, -300, 350, 300);
    box3->SetLineColor(kRed);
    box3->SetLineStyle(kDashed);
    TBox *box3_full = new TBox(0, -350, 510, 350);
    box3_full->SetLineColor(kGreen);
    TBox *box3_TMS = new TBox(730, -300, 1365, 300);
    box3_TMS->SetLineColor(kRed);
    box3_TMS->SetLineStyle(kDashed);
    TBox *box3_TMS_full = new TBox(730, -350, 1415, 350);
    box3_TMS_full->SetLineColor(kGreen);

    p1->cd();
    box3->Draw("same");
    box3_TMS->Draw("same");
    box3_full->Draw("same");
    box3_TMS_full->Draw("same");
    xzVert->Draw("P,same");
    xzDeath->Draw("P,same");
    xzBirth->Draw("P,same");
    xzExit->Draw("P,same");

    xzTraj->SetTitle("xz, all particles");

    xzTraj->Draw("colz");

    box3->Draw("same");
    box3_TMS->Draw("same");
    box3_full->Draw("same");
    box3_TMS_full->Draw("same");
    xzVert->Draw("P,same");
    xzDeath->Draw("P,same");
    xzBirth->Draw("P,same");
    xzExit->Draw("P,same");

    lepE /= 1.E3;

    // Add a little box with the info
    TPaveText *text = new TPaveText(0, 0.9, 1, 0.99, "NDC");
    canv->cd(0);
    text->SetBorderSize(0);
    text->AddText(Form("Event %i, E_{%s}=%2.2f GeV, E_{%s}=%2.2f GeV, %s", ievt, nuflav.c_str(), Enu, flav.c_str(), lepE, Where.c_str()));
    std::cout << reac << std::endl;
    //text->AddText((reac));
    text->Draw("same");

    canv->Print(canvname+".pdf");

    delete text;
    delete box3;
    delete box3_full;
    delete box3_TMS;
    delete box3_TMS_full;
  }
  std::cout << nPass << "/" << nEntries << " (" << nPass*100./nEntries << "%) passed" << std::endl;
  std::cout << nBad << "/" << nEntries << " (" << nBad*100./nEntries << "%) events with hadronic energy > " << gev_cut << " GeV" << std::endl;
  canv->Print(canvname+".pdf]");
}

