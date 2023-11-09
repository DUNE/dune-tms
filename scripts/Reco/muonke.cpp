// Simple script to run on TMS reco output
// Uses the track length and some basic selection to get a KE estimate
//
void muonke(std::string filename) {
  int nLinesCut = 1; // How many lines can our events have?
  int nClustersCut = 0; // Only allow for less clusters than this
  float OccupancyCut = 0.9; // Only allow for higher occupancy tracks than this
  bool CCmuOnly = false; // Only include events with a true CC muon (not necessarily selected as the track, but present in the event)
  float ClusterEnergyCut = 0; // How many MeV of energy in all clusters do we cut on (greater than this number gets excluded)
  bool AtLeastOneLine = true; // Do we require at least one line? (necessary for track length measurement)
  bool AllDet = true; // Muon starts in the whole detector? Or just thin region

  TFile *f = new TFile(filename.c_str());
  TTree *truth = (TTree*)f->Get("Truth_Info");
  TTree *reco = (TTree*)f->Get("Line_Candidates");

  // Reco info
  const int MAX_LINES = 20;
  const int MAX_HITS = 1500;
  const int MAX_CLUSTERS = 25;
  int nLinesOne;
  int nLinesOther;
  float TotalTrackEnergyOne[MAX_LINES];
  float TotalTrackEnergyOther[MAX_LINES];
  float TrackLengthOne[MAX_LINES];
  float TrackLengthOther[MAX_LINES];
  float OccupancyOne[MAX_LINES];
  float OccupancyOther[MAX_LINES];
  float FirstHoughHitOne[MAX_LINES][2];
  float FirstHoughHitOther[MAX_LINES][2];
  float LastHoughHitOne[MAX_LINES][2];
  float LastHoughHitOther[MAX_LINES][2];
  bool TMSStart;
  int EventNum_reco;

  float RecoHitEnergy[MAX_HITS];
  float ClusterEnergyOne[MAX_CLUSTERS];
  float ClusterEnergyOther[MAX_CLUSTERS];
  int nClustersOne;
  int nClustersOther;

  reco->SetBranchStatus("*", false);

  reco->SetBranchStatus("nLinesOne", true);
  reco->SetBranchStatus("nLinesOther", true);
  reco->SetBranchAddress("nLinesOne", &nLinesOne);
  reco->SetBranchAddress("nLinesOther", &nLinesOther);
  reco->SetBranchStatus("TotalTrackEnergyOne", true);
  reco->SetBranchStatus("TotalTrackEnergyOther", true);
  reco->SetBranchAddress("TotalTrackEnergyOne", TotalTrackEnergyOne);
  reco->SetBranchAddress("TotalTrackEnergyOther", TotalTrackEnergyOther);
  reco->SetBranchStatus("TrackLengthOne", true);
  reco->SetBranchStatus("TrackLengthOther", true);
  reco->SetBranchAddress("TrackLengthOne", TrackLengthOne);
  reco->SetBranchAddress("TrackLengthOther", TrackLengthOther);
  reco->SetBranchStatus("OccupancyOne", true);
  reco->SetBranchStatus("OccupancyOther", true);
  reco->SetBranchAddress("OccupancyOne", OccupancyOne);
  reco->SetBranchAddress("OccupancyOther", OccupancyOther);
  reco->SetBranchStatus("FirstHoughHitOne", true);
  reco->SetBranchStatus("FirstHoughHitOther", true);
  reco->SetBranchAddress("FirstHoughHitOne", FirstHoughHitOne);
  reco->SetBranchAddress("FirstHoughHitOther", FirstHoughHitOther);
  reco->SetBranchStatus("LastHoughHitOne", true);
  reco->SetBranchStatus("LastHoughHitOther", true);
  reco->SetBranchAddress("LastHoughHitOne", LastHoughHitOne);
  reco->SetBranchAddress("LastHoughHitOther", LastHoughHitOther);
  reco->SetBranchStatus("RecoHitEnergy", true);
  reco->SetBranchAddress("RecoHitEnergy", RecoHitEnergy);

  reco->SetBranchStatus("ClusterEnergyOne", true);
  reco->SetBranchStatus("ClusterEnergyOther", true);
  reco->SetBranchAddress("ClusterEnergyOne", ClusterEnergyOne);
  reco->SetBranchAddress("ClusterEnergyOther", ClusterEnergyOther);
  reco->SetBranchStatus("nClustersOne", true);
  reco->SetBranchStatus("nClustersOther", true);
  reco->SetBranchAddress("nClustersOne", &nClustersOne);
  reco->SetBranchAddress("nClustersOther", &nClustersOther);

  reco->SetBranchStatus("EventNo", true);
  reco->SetBranchAddress("EventNo", &EventNum_reco);

  // Truth info
  float MuonP4[4];
  float Muon_Vertex[4];
  float Muon_Death[4];
  float Muon_TrueKE;
  int nParticles;
  int NeutrinoPDG;
  float NeutrinoP4[4];
  int EventNum_true;

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
  truth->SetBranchStatus("NeutrinoP4", true);
  truth->SetBranchAddress("NeutrinoP4", NeutrinoP4);
  truth->SetBranchStatus("EventNo", true);
  truth->SetBranchAddress("EventNo", &EventNum_true);

  TH2D *KEOne = new TH2D("KEOne", "KEOne;True muon KE (MeV);Track length of best track (g/cm^{2})", 100, 0, 5000, 50, 0, 2500);
  TH2D *KEOther = new TH2D("KEOther", "KEOther;True muon KE (MeV);Track length of best track (g/cm^{2})", 100, 0, 5000, 50, 0, 2500);
  TH1D *h_OccupancyOne = new TH1D("OccOne", "OccupancyOne; Occupancy of longest track; Number of events", 110, 0, 1.1);
  TH1D *h_OccupancyOther = new TH1D("OccOther", "OccupancyOther; Occupancy of longest track; Number of events", 110, 0, 1.1);

  TH2D *KEestOne = new TH2D("KEestOne", "KEOne estimator; True muon KE (MeV); KE estimate", 50, 0, 5000, 50, 0, 5000);
  KEOne->GetYaxis()->SetTitleOffset(KEOne->GetYaxis()->GetTitleOffset()*1.4);
  KEOne->GetZaxis()->SetTitleOffset(KEOne->GetZaxis()->GetTitleOffset()*1.4);

  TH2D *KEestOther = new TH2D("KEestOther", "KEOther estimator; True muon KE (MeV); KE estimate", 50, 0, 5000, 50, 0, 5000);
  KEOther->GetYaxis()->SetTitleOffset(KEOther->GetYaxis()->GetTitleOffset()*1.4);
  KEOther->GetZaxis()->SetTitleOffset(KEOther->GetZaxis()->GetTitleOffset()*1.4);

  KEestOne->GetYaxis()->SetTitleOffset(KEestOne->GetYaxis()->GetTitleOffset()*1.4);
  KEestOne->GetZaxis()->SetTitleOffset(KEestOne->GetZaxis()->GetTitleOffset()*1.4);

  KEestOther->GetYaxis()->SetTitleOffset(KEestOther->GetYaxis()->GetTitleOffset()*1.4);
  KEestOther->GetZaxis()->SetTitleOffset(KEestOther->GetZaxis()->GetTitleOffset()*1.4);

  int ngood = 0;
  int nentries = truth->GetEntries();
  int trklenOne_counter = 0;
  int trklenOther_counter = 0;
  std::cout << nentries << " events..." << std::endl;
  int true_entry = 0;
  int reco_entry = 0;

  std::cout << "nLinesCut: " << nLinesCut << std::endl;

  for (int i = 0; i < nentries; ++i, ++true_entry, ++reco_entry) {
    truth->GetEntry(true_entry);
    reco->GetEntry(reco_entry);

    if (EventNum_reco != EventNum_true) {
      std::cout << "Event " << i << " " << EventNum_reco << " " << EventNum_true << std::endl;
      std::cout << true_entry << " " << reco_entry << std::endl;
      true_entry++;
      continue;
    }

    if (i % int(nentries/100.) == 0) std::cout << "Event " << i << std::endl;

//    if (nLinesOne) std::cout << "nLinesOne: " << nLinesOne << std::endl;
//    if (nClustersOne) std::cout << "nClustersOne: " << nClustersOne << std::endl;
//    if (nLinesOther) std::cout << "nLinesOther: " << nLinesOther << std::endl;
//    if (nClustersOther) std::cout << "nClustersOther: " << nClustersOther << std::endl;

    // The event has to have a muon
    if (Muon_TrueKE < 0) continue;

    // Ask for a true muon in the event if requested
    if (CCmuOnly && Muon_Vertex[2] < 0) continue;

    // Only include events with lines
    if (AtLeastOneLine && nLinesOne < 1 && nLinesOther < 1) continue;

    // Run the cuts
    if ((nLinesOne > nLinesCut) && (nLinesOther > nLinesCut)) continue;
    if ((nClustersOne > nClustersCut) && (nClustersOther > nClustersCut)) continue;
    // Sum up the total cluster energy
    float clusterOne_en = 0;
    float clusterOther_en = 0;
    for (int j = 0; j < nClustersOne; ++j) {
      clusterOne_en += ClusterEnergyOne[j];
    }

    for (int j = 0; j < nClustersOther; ++j) {
      clusterOther_en += ClusterEnergyOther[j];
    }
    if ((clusterOne_en > ClusterEnergyCut) || (clusterOther_en > ClusterEnergyCut)) continue;

    // Find the best track
    int besttrackOne = 0;
    for (int j = 0; j < nLinesOne; ++j) {
      if (OccupancyOne[j] > OccupancyOne[besttrackOne]) besttrackOne = j;
    }
    
    if (nLinesOne) std::cout << "besttrackOne: " << besttrackOne << std::endl;

    int besttrackOther = 0;
    for (int j = 0; j < nLinesOther; ++j) {
      if (OccupancyOther[j] > OccupancyOther[besttrackOther]) besttrackOther = j;
    }

    if (nLinesOther) std::cout << "besttrackOther: " << besttrackOther << std::endl;

    // Also check the track with the longest track length
    int lon_trklenOne = 0;
    for (int j = 0; j < nLinesOne; ++j) {
      if (TrackLengthOne[j] > TrackLengthOne[lon_trklenOne]) lon_trklenOne = j;
    }

    if (nLinesOne) std::cout << "lon_trklenOne: " << lon_trklenOne << std::endl;

    int lon_trklenOther = 0;
    for (int j = 0; j < nLinesOther; ++j) {
      if (TrackLengthOther[j] > TrackLengthOther[lon_trklenOther]) lon_trklenOther = j;
    }

    if (nLinesOther) std::cout << "lon_trklenOther: " << lon_trklenOther << std::endl;

    // And also check longest track
    float longestOne = 0;
    int longtrackOne = 0;
    for (int j = 0; j < nLinesOne; ++j) {
      float xdist = LastHoughHitOne[j][0]-FirstHoughHitOne[j][0];
      float ydist = LastHoughHitOne[j][1]-FirstHoughHitOne[j][1];
      float dist = sqrt(xdist*xdist+ydist*ydist);
      if (dist > longestOne) longtrackOne = j;
    }

    if (longtrackOne != lon_trklenOne) {
      trklenOne_counter++;
      if (nLinesOne) std::cout << "longtrackOne: " << longtrackOne << " | lon_trklenOne: " << lon_trklenOne << std::endl;
    }
    longtrackOne = lon_trklenOne;

    if (nLinesOne) std::cout << "longtrackOne: " << longtrackOne << std::endl;

    float longestOther = 0;
    int longtrackOther = 0;
    for (int j = 0; j < nLinesOther; ++j) {
      float xdist = LastHoughHitOther[j][0]-FirstHoughHitOther[j][0];
      float ydist = LastHoughHitOther[j][1]-FirstHoughHitOther[j][1];
      float dist = sqrt(xdist*xdist+ydist*ydist);
      if (dist > longestOther) longtrackOther = j;
    }

    if (longtrackOther != lon_trklenOther) {
      trklenOther_counter++;
      if (nLinesOther) std::cout << "longtrackOther: " << longtrackOther << " | lon_trklenOther: " << lon_trklenOther << std::endl;
    }
    longtrackOther = lon_trklenOther;

    if (nLinesOther) std::cout << "longtrackOther: " << longtrackOther << std::endl;

    // Look only at events with true muons that die inside the detector
    if (Muon_Death[1] > 1159 || Muon_Death[1] < -3864) continue;

    std::cout << "Muon didn't die outside" << std::endl;

    // Check that the longest track stops in the detector, and starts in the detector FV
    if (AllDet) {
      if (FirstHoughHitOne[longtrackOne][0] < 11362+55*2) continue;
      if (FirstHoughHitOther[longtrackOther][0] < 11362+55*2) continue;
    } else {
      if (FirstHoughHitOne[longtrackOne][0] < 11362+55*2 || FirstHoughHitOne[longtrackOne][0] > 13600) continue;
      if (FirstHoughHitOther[longtrackOther][0] < 11362+55*2 || FirstHoughHitOther[longtrackOther][0] > 13600) continue;
    }
    if (LastHoughHitOne[longtrackOne][0] > 18294-80*2) continue;
    if (LastHoughHitOther[longtrackOther][0] > 18294-80*2) continue;
    //if (LastHoughHit[longtrack][0] > 13600) continue;

    std::cout << "Longest track stops in detector" << std::endl;

    // 20 cm inwards
    if ((fabs(FirstHoughHitOne[longtrackOne][1]) > 3520-200) || (fabs(FirstHoughHitOther[longtrackOther][1]) > 3520-200)) continue;
    if ((fabs(LastHoughHitOne[longtrackOne][1]) > 3520-200) || (fabs(LastHoughHitOther[longtrackOther][1]) > 3520-200)) continue;

    h_OccupancyOne->Fill(OccupancyOne[longtrackOne]);
    h_OccupancyOther->Fill(OccupancyOther[longtrackOther]);

    std::cout << "20 cm inwards" << std::endl;
    std::cout << "Occupancy Cut: " << OccupancyCut << std::endl;
    
    // Ask for only small amount of other energy deposits
    //if (OccupancyOne[longtrackOne] < OccupancyCut) continue;
    //if (OccupancyOther[longtrackOther] < OccupancyCut) continue;

    std::cout << "only small amount of other energy deposited" << std::endl;

    float best_tracklengthOne = TrackLengthOne[longtrackOne];
    KEestOne->Fill(Muon_TrueKE, 82+1.75*best_tracklengthOne);
    KEOne->Fill(Muon_TrueKE, best_tracklengthOne);
    
    float best_tracklengthOther = TrackLengthOther[longtrackOther];
    KEestOther->Fill(Muon_TrueKE, 82+1.75*best_tracklengthOther);
    KEOther->Fill(Muon_TrueKE, best_tracklengthOther);

    if (nLinesOne) std::cout << "Best track length One: " << best_tracklengthOne << std::endl;
    if (nLinesOther) std::cout << "Best track length Other: " << best_tracklengthOther << std::endl;

    ngood++;
  }
  std::cout << ngood << "/" << nentries << std::endl;
  std::cout << trklenOne_counter << std::endl;
  std::cout << trklenOther_counter << std::endl;

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
  KEOne->Draw("colz");
  canv->Print(canvname);
  KEOther->Draw("colz");
  canv->Print(canvname);
  h_OccupancyOne->Draw();
  canv->Print(canvname);
  h_OccupancyOther->Draw();
  canv->Print(canvname);

  TH1D *arithOne = new TH1D("arith", "arith", KEOne->GetXaxis()->GetNbins(), KEOne->GetXaxis()->GetBinLowEdge(1), KEOne->GetXaxis()->GetBinLowEdge(KEOne->GetXaxis()->GetNbins()+1));
  gStyle->SetOptStat(1111);
  TH1D *arithOther = new TH1D("arith", "arith", KEOther->GetXaxis()->GetNbins(), KEOther->GetXaxis()->GetBinLowEdge(1), KEOther->GetXaxis()->GetBinLowEdge(KEOther->GetXaxis()->GetNbins()+1));
  // Now make the muon KE
  gStyle->SetOptStat(1111);
  /*for (int i = 0; i < KE->GetXaxis()->GetNbins(); ++i) {
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
  }*/
  TF1 *fitOne = new TF1("fit", "[0]+[1]*x", KEOne->GetYaxis()->GetBinLowEdge(1), KEOne->GetYaxis()->GetBinLowEdge(KEOne->GetYaxis()->GetNbins()+1));
  gStyle->SetOptStat(0);
  KEOne->Draw("colz");
  arithOne->Draw("same");
  arithOne->Fit(fitOne, "S", "", 700, 2500);
  TLegend *legOne = new TLegend(0.2, 0.5, 0.6, 0.9);
  legOne->AddEntry(arithOne, "Artihmetic mean and RMS", "le");
  legOne->AddEntry(fitOne, Form("y=%.2f x + %.2f", fitOne->GetParameter(1), fitOne->GetParameter(0)), "l");
  fitOne->SetLineColor(kRed);
  arithOne->GetXaxis()->SetTitle("Muon True KE");
  legOne->Draw("same");

  canv->Print(canvname);

  TF1 *fitOther = new TF1("fit", "[0]+[1]*x", KEOther->GetYaxis()->GetBinLowEdge(1), KEOther->GetYaxis()->GetBinLowEdge(KEOther->GetYaxis()->GetNbins()+1));
  gStyle->SetOptStat(0);
  KEOther->Draw("colz");
  arithOther->Draw("same");
  arithOther->Fit(fitOther, "S", "", 700, 2500);
  TLegend *legOther = new TLegend(0.2, 0.5, 0.6, 0.9);
  legOther->AddEntry(arithOther, "Arithmetic mean and RMS", "le");
  legOther->AddEntry(fitOther, Form("x=%.2f x + %.2f", fitOther->GetParameter(1), fitOther->GetParameter(0)), "l");
  fitOther->SetLineColor(kRed);
  arithOther->GetXaxis()->SetTitle("Muon True KE");
  legOther->Draw("same");

  canv->Print(canvname);

  KEestOne->Draw("colz");
  canv->Print(canvname);

  KEestOther->Draw("colz");
  canv->Print(canvname);

  canv->Print(canvname+"]");
}
