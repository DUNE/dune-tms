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
  int nLinesU;
  int nLinesV;
  float TotalTrackEnergyU[MAX_LINES];
  float TotalTrackEnergyV[MAX_LINES];
  float TrackLengthU[MAX_LINES];
  float TrackLengthV[MAX_LINES];
  float OccupancyU[MAX_LINES];
  float OccupancyV[MAX_LINES];
  float FirstHoughHitU[MAX_LINES][2];
  float FirstHoughHitV[MAX_LINES][2];
  float LastHoughHitU[MAX_LINES][2];
  float LastHoughHitV[MAX_LINES][2];
  bool TMSStart;
  int EventNum_reco;

  float RecoHitEnergy[MAX_HITS];
  float ClusterEnergyU[MAX_CLUSTERS];
  float ClusterEnergyV[MAX_CLUSTERS];
  int nClustersU;
  int nClustersV;

  reco->SetBranchStatus("*", false);

  reco->SetBranchStatus("nLinesU", true);
  reco->SetBranchStatus("nLinesV", true);
  reco->SetBranchAddress("nLinesU", &nLinesU);
  reco->SetBranchAddress("nLinesV", &nLinesV);
  reco->SetBranchStatus("TotalTrackEnergyU", true);
  reco->SetBranchStatus("TotalTrackEnergyV", true);
  reco->SetBranchAddress("TotalTrackEnergyU", TotalTrackEnergyU);
  reco->SetBranchAddress("TotalTrackEnergyV", TotalTrackEnergyV);
  reco->SetBranchStatus("TrackLengthU", true);
  reco->SetBranchStatus("TrackLengthV", true);
  reco->SetBranchAddress("TrackLengthU", TrackLengthU);
  reco->SetBranchAddress("TrackLengthV", TrackLengthV);
  reco->SetBranchStatus("OccupancyU", true);
  reco->SetBranchStatus("OccupancyV", true);
  reco->SetBranchAddress("OccupancyU", OccupancyU);
  reco->SetBranchAddress("OccupancyV", OccupancyV);
  reco->SetBranchStatus("FirstHoughHitU", true);
  reco->SetBranchStatus("FirstHoughHitV", true);
  reco->SetBranchAddress("FirstHoughHitU", FirstHoughHitU);
  reco->SetBranchAddress("FirstHoughHitV", FirstHoughHitOV);
  reco->SetBranchStatus("LastHoughHitU", true);
  reco->SetBranchStatus("LastHoughHitV", true);
  reco->SetBranchAddress("LastHoughHitU", LastHoughHitU);
  reco->SetBranchAddress("LastHoughHitV", LastHoughHitV);
  reco->SetBranchStatus("RecoHitEnergy", true);
  reco->SetBranchAddress("RecoHitEnergy", RecoHitEnergy);

  reco->SetBranchStatus("ClusterEnergyU", true);
  reco->SetBranchStatus("ClusterEnergyV", true);
  reco->SetBranchAddress("ClusterEnergyU", ClusterEnergyU);
  reco->SetBranchAddress("ClusterEnergyV", ClusterEnergyV);
  reco->SetBranchStatus("nClustersU", true);
  reco->SetBranchStatus("nClustersV", true);
  reco->SetBranchAddress("nClustersU", &nClustersU);
  reco->SetBranchAddress("nClustersV", &nClustersV);

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

  TH2D *KEU = new TH2D("KEU", "KEU;True muon KE (MeV);Track length of best track (g/cm^{2})", 100, 0, 5000, 50, 0, 2500);
  TH2D *KEV = new TH2D("KEV", "KEV;True muon KE (MeV);Track length of best track (g/cm^{2})", 100, 0, 5000, 50, 0, 2500);
  TH1D *h_OccupancyU = new TH1D("OccU", "OccupancyU; Occupancy of longest track; Number of events", 110, 0, 1.1);
  TH1D *h_OccupancyV = new TH1D("OccV", "OccupancyV; Occupancy of longest track; Number of events", 110, 0, 1.1);

  TH2D *KEestU = new TH2D("KEestU", "KEU estimator; True muon KE (MeV); KE estimate", 50, 0, 5000, 50, 0, 5000);
  KEU->GetYaxis()->SetTitleOffset(KEU->GetYaxis()->GetTitleOffset()*1.4);
  KEU->GetZaxis()->SetTitleOffset(KEU->GetZaxis()->GetTitleOffset()*1.4);

  TH2D *KEestV = new TH2D("KEestV", "KEV estimator; True muon KE (MeV); KE estimate", 50, 0, 5000, 50, 0, 5000);
  KEV->GetYaxis()->SetTitleOffset(KEV->GetYaxis()->GetTitleOffset()*1.4);
  KEV->GetZaxis()->SetTitleOffset(KEV->GetZaxis()->GetTitleOffset()*1.4);

  KEestU->GetYaxis()->SetTitleOffset(KEestU->GetYaxis()->GetTitleOffset()*1.4);
  KEestU->GetZaxis()->SetTitleOffset(KEestU->GetZaxis()->GetTitleOffset()*1.4);

  KEestV->GetYaxis()->SetTitleOffset(KEestV->GetYaxis()->GetTitleOffset()*1.4);
  KEestV->GetZaxis()->SetTitleOffset(KEestV->GetZaxis()->GetTitleOffset()*1.4);

  int ngood = 0;
  int nentries = truth->GetEntries();
  int trklenU_counter = 0;
  int trklenV_counter = 0;
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
    if (AtLeastOneLine && nLinesU < 1 && nLinesV < 1) continue;

    // Run the cuts
    if ((nLinesU > nLinesCut) && (nLinesV > nLinesCut)) continue;
    if ((nClustersU > nClustersCut) && (nClustersV > nClustersCut)) continue;
    // Sum up the total cluster energy
    float clusterU_en = 0;
    float clusterV_en = 0;
    for (int j = 0; j < nClustersU; ++j) {
      clusterU_en += ClusterEnergyU[j];
    }

    for (int j = 0; j < nClustersV; ++j) {
      clusterV_en += ClusterEnergyV[j];
    }
    if ((clusterU_en > ClusterEnergyCut) || (clusterV_en > ClusterEnergyCut)) continue;

    // Find the best track
    int besttrackU = 0;
    for (int j = 0; j < nLinesU; ++j) {
      if (OccupancyU[j] > OccupancyU[besttrackU]) besttrackU = j;
    }
    
    if (nLinesU) std::cout << "besttrackU: " << besttrackU << std::endl;

    int besttrackV = 0;
    for (int j = 0; j < nLinesV; ++j) {
      if (OccupancyV[j] > OccupancyV[besttrackV]) besttrackV = j;
    }

    if (nLinesV) std::cout << "besttrackV: " << besttrackV << std::endl;

    // Also check the track with the longest track length
    int lon_trklenU = 0;
    for (int j = 0; j < nLinesU; ++j) {
      if (TrackLengthU[j] > TrackLengthU[lon_trklenU]) lon_trklenU = j;
    }

    if (nLinesU) std::cout << "lon_trklenU: " << lon_trklenU << std::endl;

    int lon_trklenV = 0;
    for (int j = 0; j < nLinesV; ++j) {
      if (TrackLengthV[j] > TrackLengthV[lon_trklenV]) lon_trklenV = j;
    }

    if (nLinesV) std::cout << "lon_trklenV: " << lon_trklenV << std::endl;

    // And also check longest track
    float longestU = 0;
    int longtrackU = 0;
    for (int j = 0; j < nLinesU; ++j) {
      float xdist = LastHoughHitU[j][0]-FirstHoughHitU[j][0];
      float ydist = LastHoughHitU[j][1]-FirstHoughHitU[j][1];
      float dist = sqrt(xdist*xdist+ydist*ydist);
      if (dist > longestU) longtrackU = j;
    }

    if (longtrackU != lon_trklenU) {
      trklenU_counter++;
      if (nLinesU) std::cout << "longtrackU: " << longtrackU << " | lon_trklenU: " << lon_trklenU << std::endl;
    }
    longtrackU = lon_trklenU;

    if (nLinesU) std::cout << "longtrackU: " << longtrackU << std::endl;

    float longestV = 0;
    int longtrackV = 0;
    for (int j = 0; j < nLinesV; ++j) {
      float xdist = LastHoughHitV[j][0]-FirstHoughHitV[j][0];
      float ydist = LastHoughHitV[j][1]-FirstHoughHitV[j][1];
      float dist = sqrt(xdist*xdist+ydist*ydist);
      if (dist > longestV) longtrackV = j;
    }

    if (longtrackV != lon_trklenV) {
      trklenV_counter++;
      if (nLinesV) std::cout << "longtrackV: " << longtrackV << " | lon_trklenV: " << lon_trklenV << std::endl;
    }
    longtrackV = lon_trklenV;

    if (nLinesV) std::cout << "longtrackV: " << longtrackV << std::endl;

    // Look only at events with true muons that die inside the detector
    if (Muon_Death[1] > 1159 || Muon_Death[1] < -3864) continue;

    std::cout << "Muon didn't die outside" << std::endl;

    // Check that the longest track stops in the detector, and starts in the detector FV
    if (AllDet) {
      if (FirstHoughHitU[longtrackU][0] < 11185+55*2) continue;
      if (FirstHoughHitV[longtrackV][0] < 11185+55*2) continue;
    } else {
      if (FirstHoughHitU[longtrackU][0] < 11185+55*2 || FirstHoughHitU[longtrackU][0] > 14435) continue;  // Changed here to the new TMS_Thick_Start
      if (FirstHoughHitV[longtrackV][0] < 11185+55*2 || FirstHoughHitV[longtrackV][0] > 14435) continue;  // Same here
    }
    if (LastHoughHitU[longtrackU][0] > 18535-80*2) continue;
    if (LastHoughHitV[longtrackV][0] > 18535-80*2) continue;
    //if (LastHoughHit[longtrack][0] > 14435) continue; // Same here

    std::cout << "Longest track stops in detector" << std::endl;

    // 20 cm inwards
    if ((fabs(FirstHoughHitU[longtrackU][1]) > 3520-200) || (fabs(FirstHoughHitV[longtrackV][1]) > 3520-200)) continue;
    if ((fabs(LastHoughHitU[longtrackU][1]) > 3520-200) || (fabs(LastHoughHitV[longtrackV][1]) > 3520-200)) continue;

    h_OccupancyU->Fill(OccupancyU[longtrackU]);
    h_OccupancyV->Fill(OccupancyV[longtrackV]);

    std::cout << "20 cm inwards" << std::endl;
    std::cout << "Occupancy Cut: " << OccupancyCut << std::endl;
    
    // Ask for only small amount of other energy deposits
    //if (OccupancyOne[longtrackOne] < OccupancyCut) continue;
    //if (OccupancyOther[longtrackOther] < OccupancyCut) continue;

    std::cout << "only small amount of other energy deposited" << std::endl;

    float best_tracklengthU = TrackLengthU[longtrackU];
    KEestU->Fill(Muon_TrueKE, 82+1.75*best_tracklengthU);
    KEU->Fill(Muon_TrueKE, best_tracklengthU);
    
    float best_tracklengthV = TrackLengthV[longtrackV];
    KEestV->Fill(Muon_TrueKE, 82+1.75*best_tracklengthV);
    KEV->Fill(Muon_TrueKE, best_tracklengthV);

    if (nLinesU) std::cout << "Best track length U: " << best_tracklengthU << std::endl;
    if (nLinesV) std::cout << "Best track length V: " << best_tracklengthV << std::endl;

    ngood++;
  }
  std::cout << ngood << "/" << nentries << std::endl;
  std::cout << trklenU_counter << std::endl;
  std::cout << trklenV_counter << std::endl;

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
  KEU->Draw("colz");
  canv->Print(canvname);
  KEV->Draw("colz");
  canv->Print(canvname);
  h_OccupancyU->Draw();
  canv->Print(canvname);
  h_OccupancyV->Draw();
  canv->Print(canvname);

  TH1D *arithU = new TH1D("arith", "arith", KEU->GetXaxis()->GetNbins(), KEU->GetXaxis()->GetBinLowEdge(1), KEU->GetXaxis()->GetBinLowEdge(KEU->GetXaxis()->GetNbins()+1));
  gStyle->SetOptStat(1111);
  TH1D *arithV = new TH1D("arith", "arith", KEV->GetXaxis()->GetNbins(), KEV->GetXaxis()->GetBinLowEdge(1), KEV->GetXaxis()->GetBinLowEdge(KEV->GetXaxis()->GetNbins()+1));
  // Now make the muon KE
  gStyle->SetOptStat(1111);
  /*for (int i = 0; i < KE->GetXaxis()->GetNbins(); ++i) {
    double center = KE->GetXaxis()->GetBinCenter(i);
    TH1D *proj = KE->ProjectionY(Form("KE %.2f", center), i, i);
    double mean = proj->GetMean() // Same here;
    double mode = proj->GetBinCenter(proj->GetMaximumBin());
    double rms = proj->GetRMS();
    arith->SetBinContent(i+1, mode);
    arith->SetBinError(i+1, rms);
    proj->Draw();
    proj->SetStats(1);
    canv->Print(canvname);
  }*/
  TF1 *fitU = new TF1("fit", "[0]+[1]*x", KEU->GetYaxis()->GetBinLowEdge(1), KEU->GetYaxis()->GetBinLowEdge(KEU->GetYaxis()->GetNbins()+1));
  gStyle->SetOptStat(0);
  KEU->Draw("colz");
  arithU->Draw("same");
  arithU->Fit(fitU, "S", "", 700, 2500);
  TLegend *legU = new TLegend(0.2, 0.5, 0.6, 0.9);
  legU->AddEntry(arithU, "Artihmetic mean and RMS", "le");
  legU->AddEntry(fitU, Form("y=%.2f x + %.2f", fitU->GetParameter(1), fitU->GetParameter(0)), "l");
  fitU->SetLineColor(kRed);
  arithU->GetXaxis()->SetTitle("Muon True KE");
  legU->Draw("same");

  canv->Print(canvname);

  TF1 *fitV = new TF1("fit", "[0]+[1]*x", KEV->GetYaxis()->GetBinLowEdge(1), KEV->GetYaxis()->GetBinLowEdge(KEV->GetYaxis()->GetNbins()+1));
  gStyle->SetOptStat(0);
  KEV->Draw("colz");
  arithV->Draw("same");
  arithV->Fit(fitV, "S", "", 700, 2500);
  TLegend *legV = new TLegend(0.2, 0.5, 0.6, 0.9);
  legV->AddEntry(arithV, "Arithmetic mean and RMS", "le");
  legV->AddEntry(fitV, Form("x=%.2f x + %.2f", fitV->GetParameter(1), fitV->GetParameter(0)), "l");
  fitV->SetLineColor(kRed);
  arithV->GetXaxis()->SetTitle("Muon True KE");
  legV->Draw("same");

  canv->Print(canvname);

  KEestU->Draw("colz");
  canv->Print(canvname);

  KEestV->Draw("colz");
  canv->Print(canvname);

  canv->Print(canvname+"]");
}
