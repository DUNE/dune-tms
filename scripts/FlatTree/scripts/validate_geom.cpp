#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
void validate_geom(std::string filename) {
  gStyle->SetPalette(55);
  gErrorIgnoreLevel = kError;
  TFile *file = new TFile(filename.c_str());
  TTree *tree = (TTree*)file->Get("tree");
  int nEntries = tree->GetEntries();

  tree->SetBranchStatus("*", false);

  std::vector<float> *xpt = NULL;
  std::vector<float> *ypt = NULL;
  std::vector<float> *zpt = NULL;

  float muonBirth[3];
  float muonDeath[3];

  tree->SetBranchStatus("TMSBirth", true);
  tree->SetBranchAddress("TMSBirth", muonBirth);

  tree->SetBranchStatus("TMSDeath", true);
  tree->SetBranchAddress("TMSDeath", muonDeath);

  double zmin = 9999;
  double zmax = 0;
  for (int i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    if (muonDeath[2] > zmax) zmax = muonDeath[2];
    if (muonBirth[2] < zmin) zmin = muonBirth[2];
  }
  std::cout << "zmin: " << zmin << std::endl;
  std::cout << "zmax: " << zmax << std::endl;

  TH1D *zhits = new TH1D("zhits", "zhits", 10000, zmin-10, zmax+10);
  TH1D *gap = new TH1D("gap", "gap", 10000, 5, 40);

  tree->SetBranchStatus("xpt", true);
  tree->SetBranchAddress("xpt", &xpt);
  tree->SetBranchStatus("ypt", true);
  tree->SetBranchAddress("ypt", &ypt);
  tree->SetBranchStatus("zpt", true);
  tree->SetBranchAddress("zpt", &zpt);
  std::cout << "Starting main loop" << std::endl;

  int npass = 0;
  for (int i = 0; i < nEntries; ++i) {
    if (i > 10000) break;
    tree->GetEntry(i);
    int nDep = zpt->size();
    // Calculate track length
    double length = 0;
    double prevx = 0;
    double prevy = 0;
    double prevz = 0;
    double prevdist = 0;
    // Check the last entry
    if (zpt->back() < zmax-1) continue;
    for (int j = 0; j < nDep; ++j) {
      double x = (*xpt)[j];
      double y = (*ypt)[j];
      double z = (*zpt)[j];
      // Sometimes the first TMS hit is registered as in the LAr!
      if (z < zmin-1) continue;
      // If we're in the first TMS hit
      if (z < zmin+1) {
        prevx = x;
        prevy = y;
        prevz = z;
        continue;
      }

      double distx = x-prevx;
      double disty = y-prevy;
      double distz = z-prevz;

      double dist = sqrt(distx*distx+disty*disty+distz*distz);
      // Only want forward going
      // Sometimes get multiple hits in a bar
      if (distz < 0 || abs(distz) < 4.5) continue;

      std::cout << z << " " << prevz << std::endl;
      std::cout << distz << std::endl;

      // If hits are separated far the particle has gone out of the active material, skip these
      if (abs(z-prevz) > 10) continue;
      length += dist;
      zhits->Fill(z);
      gap->Fill(z-prevz);

      prevx = x;
      prevy = y;
      prevz = z;
    }
    npass++;
  }
  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->Print("geom.pdf[");

  zhits->SetFillStyle(1001);
  zhits->SetFillColor(kRed);
  zhits->SetLineColor(kRed);

  double previous = 0;
  for (int i = 0; i < zhits->GetXaxis()->GetNbins(); ++i) {
    if (zhits->GetBinContent(i+1) == 0) continue;
    double center = zhits->GetBinCenter(i+1);

    zhits->GetXaxis()->SetBinLabel(i+1, Form("+%2.2f", center-previous));
    previous = center;
    zhits->SetBinContent(i+1, 1);
  }
  zhits->Draw("hist");
  canv->SetLogy();
  canv->Print("geom.pdf");

  zhits->GetXaxis()->SetRange(0,zhits->GetXaxis()->FindBin(850));
  zhits->Draw("hist");
  canv->SetLogy();
  canv->Print("geom.pdf");

  zhits->GetXaxis()->SetRange(zhits->GetXaxis()->FindBin(850), zhits->GetXaxis()->FindBin(910));
  zhits->Draw("hist");
  canv->SetLogy();
  canv->Print("geom.pdf");

  zhits->GetXaxis()->SetRange(zhits->GetXaxis()->FindBin(910), zhits->GetXaxis()->FindBin(1000));
  zhits->Draw("hist");
  canv->SetLogy();
  canv->Print("geom.pdf");

  zhits->GetXaxis()->SetRange(zhits->GetXaxis()->FindBin(1000), zhits->GetXaxis()->FindBin(1100));
  zhits->Draw("hist");
  canv->SetLogy();
  canv->Print("geom.pdf");

  zhits->GetXaxis()->SetRange(zhits->GetXaxis()->FindBin(1100), zhits->GetXaxis()->FindBin(1200));
  zhits->Draw("hist");
  canv->SetLogy();
  canv->Print("geom.pdf");

  zhits->GetXaxis()->SetRange(zhits->GetXaxis()->FindBin(1200), zhits->GetXaxis()->FindBin(1300));
  zhits->Draw("hist");
  canv->SetLogy();
  canv->Print("geom.pdf");

  zhits->GetXaxis()->SetRange(zhits->GetXaxis()->FindBin(1300), zhits->GetXaxis()->FindBin(1400));
  zhits->Draw("hist");
  canv->SetLogy();
  canv->Print("geom.pdf");

  for (int i = 0; i < gap->GetXaxis()->GetNbins(); ++i) {
    if (gap->GetBinContent(i+1) == 0) continue;
    gap->GetXaxis()->SetBinLabel(i+1, Form("%2.2f", gap->GetBinCenter(i+1)));
    gap->SetBinContent(i+1, 1);
  }
  gap->Draw("hist");
  canv->SetLogy();
  canv->Print("geom.pdf");

  canv->Print("geom.pdf]");

}

