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

void second_deposits(std::string filename) {
  gStyle->SetPalette(55);
  gErrorIgnoreLevel = kError;

  TFile *file = new TFile(filename.c_str());
  TTree *tree = (TTree*)file->Get("tree");

  int nEntries = tree->GetEntries();

  std::vector<float> *xpt = NULL;
  std::vector<float> *ypt = NULL;
  std::vector<float> *zpt = NULL;
  std::vector<float> *deposits = NULL;

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
  float muonExitKE;

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
  tree->SetBranchStatus("sec_deposit", true);
  tree->SetBranchAddress("sec_deposit", &deposits);
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

  tree->SetBranchStatus("muonExitKE", true);
  tree->SetBranchAddress("muonExitKE", &muonExitKE);

  tree->SetBranchStatus("muonReco", true);
  tree->SetBranchAddress("muonReco", &muonReco);

  TH2D *DepvsDistToEnd = new TH2D("DepvsDist", "DepvsDist;Distance to Death (cm); Energy Deposited", 200, 0, 700, 100, 0, 7);
  TH2D *DepvsDistToEndZ = new TH2D("DepvsDistZ", "DepvsDist;Distance to Death Z (cm); Energy Deposited", 200, 0, 700, 100, 0, 7);

  TH2D *DepvsDistToEnd_zoom = new TH2D("DepvsDist_zoom", "DepvsDist;Distance to Death (cm); Energy Deposited", 200, 0, 20, 100, 0, 7);
  TH2D *DepvsDistToEndZ_zoom = new TH2D("DepvsDistZ_zoom", "DepvsDist;Distance to Death Z (cm); Energy Deposited", 200, 0, 20, 100, 0, 7);

  for (int i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);

    if (i%10000==0) {
      std::cout << "Entry " << i << std::endl;
    }

    if (abs(vtx[0]) > 300 || abs(vtx[1]) > 100 || vtx[2] < 50 || vtx[2] > 350) continue;

    // Fiducial cut in TMS in z
    if (muonBirth[2] > 740. || muonDeath[2] > 1375. || 
        muonDeath[1] > 25.+50. || muonDeath[1] < -260.+50. ||
        abs(muonDeath[0]) > 300. ||
        muonExitKE <= 0 ||
        abs(muonDeath[0]) > 165. || abs(muonDeath[0]) < 10. ) continue;

    if (muonDeath[2] > 939) continue;

    if (muonReco != 2) continue;

    // Only contains deposits in the LAr
    int nDeposits = deposits->size();
    // Contains the hits in the LAr and the TMS
    int nHits = xpt->size();
    int start = nHits-nDeposits;
    double olddist = 1E5;
    double sum = 0;
    double oldz = 1E5;
    for (int j = 0; j < nDeposits; ++j) {
      // Calculate the distance to the finish point
      double xdist = muonDeath[0] - (*xpt)[j+start];
      double ydist = muonDeath[1] - (*ypt)[j+start];
      double zdist = muonDeath[2] - (*zpt)[j+start];
      double tempdeposit = (*deposits)[j];

      double dist = sqrt(xdist*xdist+ydist*ydist+zdist*zdist);
      // Skip hits that are further away than previous -> indicates secondary delayed deposits
      //if (dist > olddist || dist == 0) continue;
      if (dist > olddist) continue;
      if (tempdeposit < 0.05) continue;
      sum += tempdeposit;
      if ((*zpt)[j+start] - oldz > 4) {
      //if (abs(dist-olddist) > 4) {
      /*
        DepvsDistToEnd->Fill(dist, tempdeposit);
        DepvsDistToEnd_zoom->Fill(dist, tempdeposit);
        DepvsDistToEndZ->Fill(zdist, tempdeposit);
        DepvsDistToEndZ_zoom->Fill(zdist, tempdeposit);
        */
        DepvsDistToEnd->Fill(dist, sum);
        DepvsDistToEnd_zoom->Fill(dist, sum);
        DepvsDistToEndZ->Fill(zdist, sum);
        DepvsDistToEndZ_zoom->Fill(zdist, sum);
        sum = 0;
      } 
      oldz = (*zpt)[j+start];
      olddist = dist;
    }

  } // end for loop
  std::cout << "Done" << std::endl;

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->SetRightMargin(canv->GetRightMargin()*1.3);
  canv->SetLeftMargin(canv->GetLeftMargin()*1.2);
  TString canvname = filename.c_str();
  canvname.ReplaceAll(".root", "_SecDepEnergy_Distance_Thin");
  canv->Print(canvname+".pdf[");

  DepvsDistToEnd->Draw("colz");
  canv->Print(canvname+".pdf");

  DepvsDistToEnd_zoom->Draw("colz");
  canv->Print(canvname+".pdf");

  DepvsDistToEndZ->Draw("colz");
  canv->Print(canvname+".pdf");

  DepvsDistToEndZ_zoom->Draw("colz");
  canv->Print(canvname+".pdf");

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
