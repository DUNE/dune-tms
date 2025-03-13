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

struct SpaceTime {
  float x;
  float y;
  float z;
  float t;
};

bool operator< (SpaceTime S1, SpaceTime S2) {
    return (S1.t < S2.t);
}

void draw(std::string filename) {
  gStyle->SetPalette(55);
  gErrorIgnoreLevel = kError;

  TFile *file = new TFile(filename.c_str());
  TTree *tree = (TTree*)file->Get("tree");

  int nEntries = tree->GetEntries();

  // plot the xpt and zpt
  TH2D *xzTraj = new TH2D("xz", "xz;z (cm);x (cm)", 100, -50, 1500, 100, -400, 400);
  TH2D *yzTraj = new TH2D("yz", "yz;z (cm);y (cm)", 100, -50, 1500, 100, -300, 200);
  TH2D *xzTrajLep = new TH2D("xzLep", "xzLep;z (cm);x (cm)", 100, -50, 1500, 100, -400, 400);
  TH2D *yzTrajLep = new TH2D("yzLep", "yzLep;z (cm);y (cm)", 100, -50, 1500, 100, -300, 200);
  xzTraj->SetMinimum(-0.001);
  xzTraj->SetMaximum(20);
  yzTraj->SetMinimum(-0.001);
  yzTraj->SetMaximum(20);
  xzTrajLep->SetMinimum(-0.001);
  xzTrajLep->SetMaximum(20);
  yzTrajLep->SetMinimum(-0.001);
  yzTrajLep->SetMaximum(20);

  std::vector<float> *xpt = NULL;
  std::vector<float> *ypt = NULL;
  std::vector<float> *zpt = NULL;
  std::vector<float> *xptLep = NULL;
  std::vector<float> *yptLep = NULL;
  std::vector<float> *zptLep = NULL;
  std::vector<float> *tptLep = NULL;
  std::string *reac = NULL;
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
  tree->SetBranchStatus("ypt", true);
  tree->SetBranchAddress("ypt", &ypt);
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

  tree->SetBranchStatus("xptLep", true);
  tree->SetBranchAddress("xptLep", &xptLep);
  tree->SetBranchStatus("yptLep", true);
  tree->SetBranchAddress("yptLep", &yptLep);
  tree->SetBranchStatus("zptLep", true);
  tree->SetBranchAddress("zptLep", &zptLep);
  //tree->SetBranchStatus("tptLep", true);
  //tree->SetBranchAddress("tptLep", &tptLep);

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  //canv->Divide(2,2);
  TPad *p1 = new TPad("p1", "p1", 0.0, 0.0, 0.5, 0.45);
  p1->Draw();
  TPad *p2 = new TPad("p2", "p2", 0.5, 0.0, 1.0, 0.45);
  p2->Draw();
  TPad *p3 = new TPad("p3", "p3", 0.0, 0.45, 0.5, 0.90);
  p3->Draw();
  TPad *p4 = new TPad("p4", "p4", 0.5, 0.45, 1.0, 0.90);
  p4->Draw();

  p1->SetBottomMargin(0.10);
  p2->SetBottomMargin(0.10);
  p3->SetBottomMargin(0.10);
  p4->SetBottomMargin(0.10);
  p1->SetTopMargin(0.08);
  p2->SetTopMargin(0.08);
  p3->SetTopMargin(0.08);
  p4->SetTopMargin(0.08);
  p1->SetRightMargin(0.13);
  p2->SetRightMargin(0.13);
  p3->SetRightMargin(0.13);
  p4->SetRightMargin(0.13);

  canv->SetRightMargin(canv->GetRightMargin()*1.2);
  TString canvname = filename.c_str();
  canvname.ReplaceAll(".root", "");
  canvname+="_processed";
  canv->Print(canvname+".pdf[");

  // Make a TGraph in xz and yz of the vertex location
  TGraph *xzVert = new TGraph(1);
  xzVert->SetMarkerStyle(24);
  xzVert->SetMarkerSize(2);
  xzVert->SetMarkerColor(kGreen);
  TGraph *yzVert = new TGraph(1);
  yzVert->SetMarkerStyle(24);
  yzVert->SetMarkerSize(2);
  yzVert->SetMarkerColor(kGreen);

  TGraph *xzDeath = new TGraph(1);
  xzDeath->SetMarkerStyle(25);
  xzDeath->SetMarkerSize(2);
  xzDeath->SetMarkerColor(kRed);
  TGraph *yzDeath = new TGraph(1);
  yzDeath->SetMarkerStyle(25);
  yzDeath->SetMarkerSize(2);
  yzDeath->SetMarkerColor(kRed);

  TGraph *xzBirth = new TGraph(1);
  xzBirth->SetMarkerStyle(25);
  xzBirth->SetMarkerSize(2);
  xzBirth->SetMarkerColor(kGreen);
  TGraph *yzBirth = new TGraph(1);
  yzBirth->SetMarkerStyle(25);
  yzBirth->SetMarkerSize(2);
  yzBirth->SetMarkerColor(kGreen);

  TGraph *xzExit = new TGraph(1);
  xzExit->SetMarkerStyle(25);
  xzExit->SetMarkerSize(2);
  xzExit->SetMarkerColor(kRed);
  TGraph *yzExit = new TGraph(1);
  yzExit->SetMarkerStyle(25);
  yzExit->SetMarkerSize(2);
  yzExit->SetMarkerColor(kRed);

  bool reset = false;
  int nBad = 0;
  double gev_cut = 5;
  for (int i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);

    if (i%1000==0) {
      std::cout << "Entry " << i << std::endl;
    }

    //if (i > 100) break;

    if (abs(vtx[0]) > 300 || abs(vtx[1]) > 100 || vtx[2] < 50 || vtx[2] > 350) continue;

    if (fabs(Enu-lepE/1E3) > gev_cut) nBad++;

    int nxpoints = xpt->size();
    int nypoints = ypt->size();

    int nxpointsLep = xptLep->size();
    int nypointsLep = yptLep->size();

    xzTraj->Reset();
    yzTraj->Reset();
    xzTrajLep->Reset();
    yzTrajLep->Reset();

    for (int j = 0; j < nxpoints; ++j) xzTraj->Fill((*zpt)[j], (*xpt)[j]);
    for (int j = 0; j < nypoints; ++j) yzTraj->Fill((*zpt)[j], (*ypt)[j]);
    for (int j = 0; j < nxpointsLep; ++j) xzTrajLep->Fill((*zptLep)[j], (*xptLep)[j]);
    for (int j = 0; j < nypointsLep; ++j) yzTrajLep->Fill((*zptLep)[j], (*yptLep)[j]);

    xzVert->SetPoint(0, vtx[2], vtx[0]);
    yzVert->SetPoint(0, vtx[2], vtx[1]);

    xzDeath->SetPoint(0, muonDeath[2], muonDeath[0]);
    yzDeath->SetPoint(0, muonDeath[2], muonDeath[1]);

    xzBirth->SetPoint(0, muonBirth[2], muonBirth[0]);
    yzBirth->SetPoint(0, muonBirth[2], muonBirth[1]);

    xzExit->SetPoint(0, muonExitPt[2], muonExitPt[0]);
    yzExit->SetPoint(0, muonExitPt[2], muonExitPt[1]);

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

    p3->cd();
    // Draw lines around fiducial volume
    TBox *box3 = new TBox(50, -300, 350, 300);
    box3->SetLineColor(kRed);
    TBox *box3_TMS = new TBox(730, -300, 1300, 300);
    box3_TMS->SetLineColor(kRed);
    xzTrajLep->SetTitle("xz, lepton");
    xzTrajLep->Draw("colz");
    box3->Draw("same");
    box3_TMS->Draw("same");
    xzVert->Draw("P,same");
    xzDeath->Draw("P,same");
    xzBirth->Draw("P,same");
    xzExit->Draw("P,same");

    p4->cd();
    TBox *box4 = new TBox(50, -100, 350, 100);
    box4->SetLineColor(kRed);
    TBox *box4_TMS = new TBox(730, -260, 1300, 25);
    box4_TMS->SetLineColor(kRed);
    yzTrajLep->SetTitle("yz, lepton");
    yzTrajLep->Draw("colz");
    box4->Draw("same");
    box4_TMS->Draw("same");
    yzVert->Draw("P,same");
    yzDeath->Draw("P,same");
    yzBirth->Draw("P,same");
    yzExit->Draw("P,same");


    p1->cd();
    xzTraj->SetTitle("xz, oth");
    xzTraj->Draw("colz");
    box3->Draw("same");
    box3_TMS->Draw("same");
    xzVert->Draw("P,same");
    xzDeath->Draw("P,same");
    xzBirth->Draw("P,same");
    xzExit->Draw("P,same");

    p2->cd();
    yzTraj->SetTitle("yz, oth");
    yzTraj->Draw("colz");
    box4->Draw("same");
    box4_TMS->Draw("same");
    yzVert->Draw("P,same");
    yzDeath->Draw("P,same");
    yzBirth->Draw("P,same");
    yzExit->Draw("P,same");

    lepE /= 1.E3;

    // Add a little box with the info
    canv->cd(0);
    TPaveText *text = new TPaveText(0, 0.9, 1, 0.99, "NDC");
    text->SetBorderSize(0);
    text->AddText(Form("Event %i, E_{%s}=%2.2f GeV, E_{%s}=%2.2f GeV, %s", ievt, nuflav.c_str(), Enu, flav.c_str(), lepE, Where.c_str()));
    text->AddText((*reac).c_str());
    text->Draw("same");

    canv->Print(canvname+".pdf");
    delete text;
    delete box3;
    delete box4;
  }
  std::cout << nBad << "/" << nEntries << " (" << nBad*100./nEntries << "%) events with hadronic energy > " << gev_cut << " GeV" << std::endl;
  canv->Print(canvname+".pdf]");
}

/*
int main(int argc, char** argv) {
  TApplication a("a", 0, 0);  // just to make sure that the autoloading of ROOT libraries works
  if (argc != 2) {
    std::cerr << "Need arg" << std::endl;
    return -1;
  }
  draw(std::string(argv[1]));
  return 0;
}
*/
