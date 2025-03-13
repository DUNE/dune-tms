void eff(std::string filename) {
  TFile *file = new TFile(filename.c_str());
  TTree *tree = (TTree*)file->Get("tree");
  //TChain *tree = new TChain("tree");
  //TString fileset = filename.c_str();
  //fileset = fileset(0,fileset.Last('/')+1);
  // Get directory, add all root files
  //tree->Add(fileset+"*FlatTree.root");

  int muonReco;
  int muonStart;
  float lep[3];
  float vtx[3];
  float lepKE;
  float Enu;
  int leppdg;

  tree->SetBranchStatus("*", false);

  tree->SetBranchStatus("muonReco", true);
  tree->SetBranchAddress("muonReco", &muonReco);

  tree->SetBranchStatus("muonStart", true);
  tree->SetBranchAddress("muonStart", &muonStart);

  tree->SetBranchStatus("lepKE", true);
  tree->SetBranchAddress("lepKE", &lepKE);

  tree->SetBranchStatus("lepPdg", true);
  tree->SetBranchAddress("lepPdg", &leppdg);

  tree->SetBranchStatus("p3lep", true);
  tree->SetBranchAddress("p3lep", lep);

  tree->SetBranchStatus("Ev", true);
  tree->SetBranchAddress("Ev", &Enu);

  tree->SetBranchStatus("vtx", true);
  tree->SetBranchAddress("vtx", vtx);

  int nEntries = tree->GetEntries();
  float theta;

  TString canvname = filename.c_str();
  canvname.ReplaceAll(".root", "_efficiency_wvtx_wmode.root");
  TFile *output = new TFile(canvname, "recreate");
  TTree *outtree = new TTree("out", "out");
  outtree->Branch("theta", &theta);
  outtree->Branch("muonReco", &muonReco);
  outtree->Branch("muonStart", &muonStart);
  outtree->Branch("lep[3]", lep);
  outtree->Branch("leppdg", &leppdg);
  outtree->Branch("Enu", &Enu);
  outtree->Branch("vtx[3]", vtx);

  TVector3 beam_angle(0,0,1);
  beam_angle.RotateX(0.101);

  for (int i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    theta = -999;

    TVector3 mu(lep[0], lep[1], lep[2]);
    theta = mu.Angle(beam_angle)*180./3.1415;

    //std::cout << vtx[0] << " " << vtx[1] << " " << vtx[2] << std::endl;

    //if (leppdg != 13) continue;
    //if (theta > 30) continue;

    outtree->Fill();
  }

  output->cd();
  outtree->Write();
  output->Close();

}
