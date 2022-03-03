void muonke(std::string filename) {
  int nLinesCut = 1; // How many lines can our events have?
  int nClustersCut = 0; // Only allow for less clusters than this
  float OccupancyCut = 0.9; // Only allow for higher occupancy tracks than this
  bool CCmuOnly = false; // Only include events with a true CC muon (not necessarily selected as the track, but present in the event)
  float ClusterEnergyCut = 5; // How many MeV of energy in all clusters do we cut on (greater than this number gets excluded)
  bool AtLeastOneLine = true; // Do we require at least one line? (necessary for track length measurement)
  bool AllDet = true; // Muon starts in the whole detector? Or just thin region

  TFile *f = new TFile(filename.c_str());
  TTree *truth = (TTree*)f->Get("Truth_Info");
  TTree *reco = (TTree*)f->Get("Line_Candidates");

  // Reco info
  const int MAX_LINES = 20;
  const int MAX_HITS = 1500;
  const int MAX_CLUSTERS = 25;
  int nLines;
  float TotalTrackEnergy[MAX_LINES];
  float TrackLength[MAX_LINES];
  float Occupancy[MAX_LINES];
  float FirstHoughHit[MAX_LINES][2];
  float LastHoughHit[MAX_LINES][2];
  bool TMSStart;

  float RecoHitEnergy[MAX_HITS];
  float ClusterEnergy[MAX_CLUSTERS];
  int nClusters;

  reco->SetBranchStatus("*", false);

  reco->SetBranchStatus("nLines", true);
  reco->SetBranchAddress("nLines", &nLines);
  reco->SetBranchStatus("TotalTrackEnergy", true);
  reco->SetBranchAddress("TotalTrackEnergy", TotalTrackEnergy);
  reco->SetBranchStatus("TrackLength", true);
  reco->SetBranchAddress("TrackLength", TrackLength);
  reco->SetBranchStatus("Occupancy", true);
  reco->SetBranchAddress("Occupancy", Occupancy);
  reco->SetBranchStatus("FirstHoughHit", true);
  reco->SetBranchAddress("FirstHoughHit", FirstHoughHit);
  reco->SetBranchStatus("LastHoughHit", true);
  reco->SetBranchAddress("LastHoughHit", LastHoughHit);
  reco->SetBranchStatus("RecoHitEnergy", true);
  reco->SetBranchAddress("RecoHitEnergy", RecoHitEnergy);

  reco->SetBranchStatus("ClusterEnergy", true);
  reco->SetBranchAddress("ClusterEnergy", ClusterEnergy);
  reco->SetBranchStatus("nClusters", true);
  reco->SetBranchAddress("nClusters", &nClusters);

  // Truth info
  float MuonP4[4];
  float Muon_Vertex[4];
  float Muon_Death[4];
  float Muon_TrueKE;
  int nParticles;
  int NeutrinoPDG;
  float NeutrinoP4[4];

  truth->SetBranchStatus("*", false);

  truth->SetBranchStatus("MuonP4", true);
  truth->SetBranchAddress("MuonP4", MuonP4);
  truth->SetBranchStatus("Muon_Vertex", true);
  truth->SetBranchAddress("Muon_Vertex", Muon_Vertex);
  truth->SetBranchStatus("Muon_Death", true);
  truth->SetBranchAddress("Muon_Death", Muon_Death);
  truth->SetBranchStatus("Muon_TrueKE", true);
  truth->SetBranchAddress("Muon_TrueKE", &Muon_TrueKE);
  //truth->SetBranchStatus("nParticles", true);
  //truth->SetBranchAddress("nParticles", &nParticles);
  //truth->SetBranchStatus("NeutrinoPDG", true);
  //truth->SetBranchAddress("NeutrinoPDG", &NeutrinoPDG);
  //truth->SetBranchStatus("NeutrinoP4", true);
  //truth->SetBranchAddress("NeutrinoP4", NeutrinoP4);

  TH2D *KE = new TH2D("KE", "KE;True muon KE (MeV);Track length of best track (g/cm^{2})", 100, 0, 5000, 50, 0, 2500);
  TH1D *h_Occupancy = new TH1D("Occ", "Occupancy; Occupancy of longest track; Number of events", 100, 0, 1);

  TH2D *KEest = new TH2D("KEest", "KE estimator; True muon KE (MeV); KE estimate", 50, 0, 5000, 50, 0, 5000);
  KE->GetYaxis()->SetTitleOffset(KE->GetYaxis()->GetTitleOffset()*1.4);
  KE->GetZaxis()->SetTitleOffset(KE->GetZaxis()->GetTitleOffset()*1.4);

  KEest->GetYaxis()->SetTitleOffset(KEest->GetYaxis()->GetTitleOffset()*1.4);
  KEest->GetZaxis()->SetTitleOffset(KEest->GetZaxis()->GetTitleOffset()*1.4);

  int ngood = 0;
  int nentries = truth->GetEntries();
  int trklen_counter = 0;
  std::cout << nentries << " events..." << std::endl;
  for (int i = 0; i < nentries; ++i) {
    truth->GetEntry(i);
    reco->GetEntry(i);

    if (i % int(nentries/100.) == 0) std::cout << "Event " << i << std::endl;

    // The event has to have a muon
    if (Muon_TrueKE < 0) continue;

    // Ask for a true muon in the event if requested
    if (CCmuOnly && Muon_Vertex[2] < 0) continue;

    // Only include events with lines
    if (AtLeastOneLine && nLines < 1) continue;

    // Run the cuts
    if (nLines > nLinesCut) continue;
    if (nClusters > nClustersCut) continue;
    // Sum up the total cluster energy
    float cluster_en = 0;
    for (int j = 0; j < nClusters; ++j) {
      cluster_en += ClusterEnergy[j];
    }
    if (cluster_en > ClusterEnergyCut) continue;

    // Find the best track
    int besttrack = 0;
    for (int j = 0; j < nLines; ++j) {
      if (Occupancy[j] > Occupancy[besttrack]) besttrack = j;
    }

    // Also check the track with the longest track length
    int lon_trklen = 0;
    for (int j = 0; j < nLines; ++j) {
      if (TrackLength[j] > TrackLength[lon_trklen]) lon_trklen = j;
    }

    // And also check longest track
    float longest = 0;
    int longtrack = 0;
    for (int j = 0; j < nLines; ++j) {
      float xdist = LastHoughHit[j][0]-FirstHoughHit[j][0];
      float ydist = LastHoughHit[j][1]-FirstHoughHit[j][1];
      float dist = sqrt(xdist*xdist+ydist*ydist);
      if (dist > longest) longtrack = j;
    }

    if (longtrack != lon_trklen) {
      trklen_counter++;
    }
    longtrack = lon_trklen;

    // Look only at events with true muons that die inside the detector
    if (Muon_Death[1] > 1159 || Muon_Death[1] < -3864) continue;

    // Check that the longest track stops in the detector, and starts in the detector FV
    if (AllDet) {
      if (FirstHoughHit[longtrack][0] < 11362+55*2) continue;
    } else {
      if (FirstHoughHit[longtrack][0] < 11362+55*2 || FirstHoughHit[longtrack][0] > 13600) continue;
    }
    if (LastHoughHit[longtrack][0] > 18294-80*2) continue;
    //if (LastHoughHit[longtrack][0] > 13600) continue;

    // 10 cm inwards
    if (fabs(FirstHoughHit[longtrack][1]) > 3520-200) continue;
    if (fabs(LastHoughHit[longtrack][1]) > 3520-200) continue;

    h_Occupancy->Fill(Occupancy[longtrack]);

    // Ask for only small amount of other energy deposits
    if (Occupancy[longtrack] < OccupancyCut) continue;

    float best_tracklength = TrackLength[longtrack];
    KEest->Fill(Muon_TrueKE, 80+1.75*best_tracklength);
    KE->Fill(Muon_TrueKE, best_tracklength);
    ngood++;
  }
  std::cout << ngood << "/" << nentries << std::endl;
  std::cout << trklen_counter << std::endl;

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->SetLeftMargin(canv->GetLeftMargin()*1.5);
  canv->SetRightMargin(canv->GetRightMargin()*1.4);

  // Fix filename

  while (filename.find("/") != std::string::npos) {
    filename = filename.substr(filename.find("/")+1, filename.size());
  }

  TString canvname = Form("MuonKE_%s", filename.c_str());
  canvname += Form("_nLinesCut%i_nClusterCut%i_OccupancyCut%.2f_CCmuOnly%o_ClusterEnCut%.2f_AtLeastOneLine%o_AllDet%o", nLinesCut, nClustersCut, OccupancyCut, CCmuOnly, ClusterEnergyCut, AtLeastOneLine, AllDet);
  canvname += ".pdf";
  canv->Print(canvname+"[");

  gStyle->SetPalette(55);
  KE->Draw("colz");
  canv->Print(canvname);
  h_Occupancy->Draw();
  canv->Print(canvname);

  TH1D *arith = new TH1D("arith", "arith", KE->GetXaxis()->GetNbins(), KE->GetXaxis()->GetBinLowEdge(1), KE->GetXaxis()->GetBinLowEdge(KE->GetXaxis()->GetNbins()+1));
  // Now make the muon KE
  gStyle->SetOptStat(1111);
  for (int i = 0; i < KE->GetXaxis()->GetNbins(); ++i) {
    double center = KE->GetXaxis()->GetBinCenter(i);
    TH1D *proj = KE->ProjectionY(Form("KE %.2f", center), i, i);
    double mean = proj->GetMean();
    double mode = proj->GetBinCenter(proj->GetMaximumBin());
    double rms = proj->GetRMS();
    arith->SetBinContent(i+1, mode);
    arith->SetBinError(i+1, rms);
    proj->Draw();
    proj->SetStats(1);
    canv->Print(canvname);
  }
  TF1 *fit = new TF1("fit", "[0]+[1]*x", KE->GetYaxis()->GetBinLowEdge(1), KE->GetYaxis()->GetBinLowEdge(KE->GetYaxis()->GetNbins()+1));
  gStyle->SetOptStat(0);
  KE->Draw("colz");
  arith->Draw("same");
  arith->Fit(fit, "S", "", 700, 2500);
  TLegend *leg = new TLegend(0.2, 0.5, 0.6, 0.9);
  leg->AddEntry(arith, "Artihmetic mean and RMS", "le");
  leg->AddEntry(fit, Form("y=%.2f x + %.2f", fit->GetParameter(1), fit->GetParameter(0)), "l");
  fit->SetLineColor(kRed);
  arith->GetXaxis()->SetTitle("Muon True KE");
  leg->Draw("same");

  canv->Print(canvname);

  KEest->Draw("colz");
  canv->Print(canvname);

  // Now make a simple neutrino energy estimator assuming a QE event

  canv->Print(canvname+"]");
}
