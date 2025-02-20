void draw_dep(std::string filename) {
  gStyle->SetPalette(55);
  TFile *f = new TFile(filename.c_str());
  TTree *tree = (TTree*)f->Get("tree");

  std::vector<double> *xpt = NULL;
  std::vector<double> *ypt = NULL;
  std::vector<double> *zpt = NULL;
  tree->SetBranchAddress("xpt", &xpt);
  tree->SetBranchAddress("ypt", &ypt);
  tree->SetBranchAddress("zpt", &zpt);
  int muonReco;
  tree->SetBranchAddress("muonReco", &muonReco);
  int muonStart;
  tree->SetBranchAddress("muonStart", &muonStart);
  int ievt;
  tree->SetBranchAddress("ievt", &ievt);

  int nentries = tree->GetEntries();

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->Print("hits.pdf[");
  canv->Divide(2);

  TH2D *xz = new TH2D("xz", "xz", 500, 11E3, 19E3, 500, -4E3, 4E3);
  TH2D *yz = new TH2D("yz", "yz", 500, 11E3, 19E3, 500, -4E3, 1.5E3);
  for (int i = 0; i < nentries; ++i) {
    tree->GetEntry(i);
    if (i > 10000) break;
    if (muonReco != 2) continue;

    xz->Reset();
    yz->Reset();
    int nhits = xpt->size();
    for (int j = 0; j < nhits; ++j) {
      xz->Fill((*zpt)[j], (*xpt)[j]);
      yz->Fill((*zpt)[j], (*ypt)[j]);
    }
    canv->cd(1);
    xz->SetTitle(Form("Event %i xz", ievt));
    xz->Draw("colz");

    canv->cd(2);
    yz->SetTitle(Form("Event %i yz", ievt));
    yz->Draw("colz");
    canv->Print("hits.pdf");

  }
  canv->Print("hits.pdf]");
}
