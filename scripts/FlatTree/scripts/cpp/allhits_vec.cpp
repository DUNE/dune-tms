#include <iostream>
#include <vector>
#include <algorithm>

#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TError.h"

void allhits_vec(std::string filename) {
  gStyle->SetPalette(55);
  gErrorIgnoreLevel = kError;

  TFile *file = new TFile(filename.c_str());
  TTree *tree = (TTree*)file->Get("tree");

  int nEntries = tree->GetEntries();

  std::vector<float> *xpt = NULL;
  std::vector<float> *ypt = NULL;
  std::vector<float> *zpt = NULL;

  tree->SetBranchStatus("*", false);
  //tree->SetBranchStatus("xpt", true);
  //tree->SetBranchAddress("xpt", &xpt);
  //tree->SetBranchStatus("ypt", true);
  //tree->SetBranchAddress("ypt", &ypt);
  tree->SetBranchStatus("zpt", true);
  tree->SetBranchAddress("zpt", &zpt);

  // Keep some arrays of all hits
  //std::vector<float> xhits;
  //std::vector<float> yhits;
  std::vector<float> zhits;

  int nhits = 0;
  for (int i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    /*
    if (i%1000==0) {
      std::cout << "Entry " << i << "/" << nEntries << std::endl;
      std::cout << "nhits = " << nhits << std::endl;
    }
    */
    int nxpoints = zpt->size();

    if (nhits > 5E7) break;
    for (int j = 0; j < nxpoints; ++j) {
      //xhits.push_back((*xpt)[j]);
      //yhits.push_back((*ypt)[j]);
      zhits.push_back((*zpt)[j]);
      nhits++;
    }
  }

  /*
  TGraph *zy = new TGraph(zhits.size(), &zhits[0], &yhits[0]);
  zy->SetMarkerStyle(20);
  zy->SetMarkerSize(0.2);

  TGraph *zx = new TGraph(zhits.size(), &zhits[0], &xhits[0]);
  zx->SetMarkerStyle(20);
  zx->SetMarkerSize(0.2);

  TString canvname = filename.c_str();
  canvname.ReplaceAll(".root", "");
  canvname+="_allhits_graph.root";
  TFile *out = new TFile(canvname, "recreate");
  zx->Write("zx");
  zy->Write("zx");
  out->Close();
  */

  // Now sort the vector in z
  std::sort(zhits.begin(), zhits.end());
  // Loop over and find when gaps are greater than 1 cm
  float prevval = 0;
  float width = 0;
  float thin = 0;
  float thick = 0;
  float trans = 0;
  float prevdist = 0;
  int nthin = 0;
  int nthick = 0;
  int ntrans = 0;
  for (std::vector<float>::iterator it = zhits.begin(); it != zhits.end(); ++it) {
    float val = (*it);
    if (val < 730) continue;
    if (val - prevval < 1) continue;

    // Calculate the width
    if (val-prevval < 4.5) width = val-prevval;

    // Skip the width
    if (val-prevval == 1) continue;

    if (fabs((val-prevval)-prevdist) > 0.01 && prevdist != 0) {
      std::cout << "Found change at " << prevval << "-" << val << std::endl;
      std::cout << "Previous distance: " << prevdist << std::endl;
      std::cout << "New distance: " << val-prevval << std::endl;
      std::cout << std::endl;
    }

    if (val < 940 ) thin = val-prevval;
    if (val > 957) thick = val-prevval;
    if (val > 940 && val < 957) trans = val-prevval;

    // Count up
    if (fabs(val-prevval - 5.5) < 1E-2) nthin++;
    else if (fabs(val-prevval - 8) < 1E-2) nthick++;
    else if (fabs(val-prevval - 9.5) < 1E-2) ntrans++;

    //std::cout << val << ", " << prevval << std::endl;
    //std::cout << val-prevval << std::endl;
    prevdist = val-prevval;
    prevval = val;
  }

  std::cout << std::endl;
  std::cout << "Bar width: " << width << std::endl;
  std::cout << "Thin: " << thin << std::endl;
  std::cout << "Trans: " << trans << std::endl;
  std::cout << "Thick: " << thick << std::endl;

  std::cout << "nthin-2: " << nthin << std::endl;
  std::cout << "nthick-2: " << nthick << std::endl;
  std::cout << "ntrans: " << ntrans << std::endl;

}

int main(int argc, char** argv) {
  allhits_vec(std::string(argv[1]));
  return 0;
}

