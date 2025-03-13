double disc(double val, int dim) {
  double mini,maxi,widthi;
  if (dim == 1) {
    mini = -348.5;
    maxi = 348.5;
    widthi = 4; // 4 cm wide

    int bin = (val-mini)/widthi;

    return mini+bin*widthi;
  }
  else if (dim == 2) {
    mini = 730.8;
    maxi = 1415;
    // Account for thin and thick regions
    if (val > 949) {
      mini = 949;
      widthi = 8;
      int bin = (val-mini)/widthi;
      return mini + bin*widthi;
    }
    else if (val < 939) {
      mini = 730.8;
      widthi = 5.5;
      int bin = (val-mini)/widthi;
      return mini + bin*widthi;
    } else {
      return (949-939)/2;
    }


  }

  //std::cout << val << " " << mini << " " << widthi << " " << bin <<  std::endl;

}

void magbend(std::string filename) {
  gStyle->SetPalette(55);

  TFile *file = new TFile(filename.c_str());
  //TTree *tree = (TTree*)file->Get("tree");
  TChain *tree = new TChain("tree");
  TString fileset = filename.c_str();
  fileset = fileset(0,fileset.Last('/')+1);

  // Get directory, add all root files
  tree->Add(fileset+"*FlatTree.root");

  int nEntries = tree->GetEntries();
  std::vector<float> *xptLep = NULL;
  std::vector<float> *yptLep = NULL;
  std::vector<float> *zptLep = NULL;
  float vtx[3];
  float Enu;
  float muonDeath[3];
  float muonBirth[3];
  float muonExitPt[3];
  float muonExitKE;
  int lepPdg;
  float rmmsKE;

  //tree->Print();
  //return;

  tree->SetBranchStatus("*", false);
  tree->SetBranchStatus("Ev", true);
  tree->SetBranchAddress("Ev", &Enu);

  tree->SetBranchStatus("xpt", true);
  tree->SetBranchAddress("xpt", &xptLep);
  tree->SetBranchStatus("ypt", true);
  tree->SetBranchAddress("ypt", &yptLep);
  tree->SetBranchStatus("zpt", true);
  tree->SetBranchAddress("zpt", &zptLep);

  tree->SetBranchStatus("vtx", true);
  tree->SetBranchAddress("vtx", vtx);

  tree->SetBranchStatus("muonDeath", true);
  tree->SetBranchAddress("muonDeath", muonDeath);
  tree->SetBranchStatus("muonBirth", true);
  tree->SetBranchAddress("muonBirth", muonBirth);
  tree->SetBranchStatus("muonExitPt", true);
  tree->SetBranchAddress("muonExitPt", muonExitPt);

  tree->SetBranchStatus("muonExitKE", true);
  tree->SetBranchAddress("muonExitKE", &muonExitKE);

  tree->SetBranchStatus("lepPdg", true);
  tree->SetBranchAddress("lepPdg", &lepPdg);

  tree->SetBranchStatus("rmmsKE", true);
  tree->SetBranchAddress("rmmsKE", &rmmsKE);

  int muonReco;
  tree->SetBranchStatus("muonReco", true);
  tree->SetBranchAddress("muonReco", &muonReco);


  TH1D *LepAll = new TH1D("Lep_All", "Lep_All", 100, -200, 200);
  TH1D *LepCorrect = new TH1D("Lep_Correct", "Lep_Correct",  100, -200, 200);
  LepCorrect->SetLineColor(kBlue);
  TH1D *LepWrong = new TH1D("Lep_Wrong", "Lep_Wrong",  100, -200, 200);
  LepWrong->SetLineColor(kRed);

  TString canvname = filename.c_str();
  canvname.ReplaceAll(".root", "_signdistupd6_both_onlycent");
  TString outputname = canvname+".root";
  TFile *output = new TFile(outputname, "recreate");
  // Write some of the branches out

  float signed_dist;
  float LArKE;
  float TMSKE;
  float Enu_out;
  int pdg;

  TTree *outtree = new TTree("mag_sum", "mag_sum");
  outtree->Branch("signed_dist", &signed_dist);
  outtree->Branch("LArKE", &LArKE);
  outtree->Branch("TMSKE", &TMSKE);
  outtree->Branch("Enu", &Enu_out);
  outtree->Branch("pdg", &pdg);

  for (int i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);

    // Update the variables
    signed_dist = 0;
    LArKE = muonExitKE;
    pdg = lepPdg;
    Enu_out = Enu;
    TMSKE = rmmsKE;

    //if (lepPdg > 0) continue;
    if (muonReco != 2) continue;

    if (abs(vtx[0]) > 300. || abs(vtx[1]) > 100. || vtx[2] < 50. || vtx[2] > 350.) continue;

    if (muonBirth[2] > 735. || muonDeath[2] > 1365. || 
        muonDeath[1] > 87-50. || muonDeath[1] < -234+50. ||
        abs(muonDeath[0]) > 300. ||
        muonExitKE <= 0 ) continue; //||
        //abs(muonDeath[0]) > 165. || abs(muonDeath[0]) < 10. ) continue;

    double x1 = muonExitPt[0];
    double z1 = muonExitPt[2];

    double x2 = muonBirth[0];
    double z2 = muonBirth[2];
    x2 = disc(x2,1);
    z2 = disc(z2, 2);
    //double x2disc = disc(x2,1);
    //std::cout << x2 << "->" << x2disc << std::endl;

    double x3 = muonDeath[0];
    double z3 = muonDeath[2];
    x3 = disc(x3,1);
    z3 = disc(z3, 2);
    //double x3disc = disc(x3,1);

    // Check that the whole event is contained
    int nhits = (*xptLep).size();
    bool bad = false;
    for (int j = 0; j < nhits; ++j) {
      double x = (*xptLep)[j];
      if (fabs(x) > 160) {
        bad=true;
        break;
      }
    }
    if (bad) continue;

    //signed_dist = (-1*(z2-z1)*x3 + (x2-x1)*z3 + x1*z2 - z1*x2)/sqrt((x2-x1)*(x2-x1)+(z2-z1)*(z2-z1));
    //signed_dist = sqrt((z3-z2)*(z3-z2)+(x3-x2)*(x3-x2))/sqrt((z2-z1)*(z2-z1)+(x2-x1)*(x2-x1))*(z2-z1);
    //signed_dist = (x3-x1)+(z3-z1)*(x2-x1)/(z2-z1);
    // projected x position
    double x4 = x1-(x1-x2)*(z3-z1)/(z2-z1);
    //if (x3 > x4) signed_dist = x3-x4;
    //else signed_dist = x4-x3;
    signed_dist = x4-x3;

    LepAll->Fill(signed_dist);
    if ((lepPdg > 0 && signed_dist > 0) || (lepPdg < 0 && signed_dist < 0)) LepCorrect->Fill(signed_dist);
    else if ((lepPdg < 0 && signed_dist > 0) || (lepPdg > 0 && signed_dist < 0)) LepWrong->Fill(signed_dist);

    outtree->Fill();

  }

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->Print(canvname+".pdf[");

  LepAll->Draw();
  canv->Print(canvname+".pdf");

  LepCorrect->Draw();
  canv->Print(canvname+".pdf");

  LepWrong->Draw();
  canv->Print(canvname+".pdf");

  LepAll->Draw();
  LepCorrect->Draw("same");
  LepWrong->Draw("same");

  LepCorrect->SetTitle(Form("%2.2f/%2.2f=%2.2f", LepCorrect->Integral(), LepAll->Integral(), double(LepCorrect->Integral())/LepAll->Integral()));
  LepWrong->SetTitle(Form("%2.2f/%2.2f=%2.2f", LepWrong->Integral(), LepAll->Integral(), double(LepWrong->Integral())/LepAll->Integral()));
  canv->BuildLegend();

  canv->Print(canvname+".pdf");
  canv->Print(canvname+".pdf]");

  output->cd();
  outtree->Write();
  LepAll->Write();
  LepCorrect->Write();
  LepWrong->Write();
  output->Close();
}
