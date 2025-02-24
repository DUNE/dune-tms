void angle(std::string filename) {
  int nLinesCut = 1; // How many lines can our events have?
  int nClustersCut = 100; // Only allow for this many clusters
  float OccupancyCut = 0.80; // Only allow for higher occupancy tracks than this
  bool CCmuOnly = false; // Only include events with a true CC muon (not necessarily selected as the track, but present in the event)
  float ClusterEnergyCut = 100; // How many MeV of energy in all clusters do we cut on (greater than this number gets excluded)
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
  int nHitsInTrack[MAX_LINES];
  float TrackHitPos[MAX_LINES][200][2];

  float RecoHitEnergy[MAX_HITS];
  float ClusterEnergy[MAX_CLUSTERS];
  int nClusters;

  int EventNum_reco;

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

  reco->SetBranchStatus("TrackHitPos", true);
  reco->SetBranchAddress("TrackHitPos", TrackHitPos);
  reco->SetBranchStatus("nHitsInTrack", true);
  reco->SetBranchAddress("nHitsInTrack", nHitsInTrack);


  reco->SetBranchStatus("ClusterEnergy", true);
  reco->SetBranchAddress("ClusterEnergy", ClusterEnergy);
  reco->SetBranchStatus("nClusters", true);
  reco->SetBranchAddress("nClusters", &nClusters);

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

  TH2D *AnglevsTrue = new TH2D("AnglevsTrue", "AnglevsTrue; muon z true angle (degrees); muon z reco angle (degrees)", 360, 0, 90, 360, 0, 90);
  AnglevsTrue->GetYaxis()->SetTitleOffset(AnglevsTrue->GetYaxis()->GetTitleOffset()*1.4);
  AnglevsTrue->GetZaxis()->SetTitleOffset(AnglevsTrue->GetZaxis()->GetTitleOffset()*1.4);

  // Prepare the neutrino
  TVector3 nuvect(0, 0, 1);
  nuvect.RotateX(0.101);

  int ngood = 0;
  int nentries = truth->GetEntries();
  int trklen_counter = 0;
  //nentries = 1E5;
  std::cout << nentries << " events..." << std::endl;
  int true_entry = 0;
  int reco_entry = 0;
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
    //std::cout << "Event " << i << std::endl;

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
    //if (LastHoughHit[longtrack][0] > 18294-80*2) continue;
    //if (LastHoughHit[longtrack][0] > 13600) continue;

    // 10 cm inwards
    if (fabs(FirstHoughHit[longtrack][1]) > 3520-200) continue;
    if (fabs(LastHoughHit[longtrack][1]) > 3520-200) continue;

    // Ask for only small amount of other energy deposits
    if (Occupancy[longtrack] < OccupancyCut) continue;

    // Calculate the angle in the x-z plane relative the neutrino
    // Recalculate from the track start positions
    // Find where the kink in the track happens (if at all)
    double firstx = TrackHitPos[longtrack][0][1];
    int kink = 0;
    //std::cout << "***" << std::endl;
    //std::cout << "Event " << i << std::endl;
    for (int j = 0; j < nHitsInTrack[longtrack]; ++j) {
      double z = TrackHitPos[longtrack][j][0];
      double x = TrackHitPos[longtrack][j][1];
      //std::cout << "hit " << j << " " << z << " " << x << std::endl;
      if (fabs(firstx - x) > 30.54*2 && kink == 0) {
        //std::cout << "kink" << std::endl;
        kink = j;
        break;
      }
    }
    if (kink == 0) {
      //std::cout << "found no kink" << std::endl;
      kink = nHitsInTrack[longtrack]-1;
      //std::cout << TrackHitPos[longtrack][kink][1] << " to " << firstx << std::endl;
      //std::cout << TrackHitPos[longtrack][kink][0] << " to " << TrackHitPos[longtrack][0][0] << std::endl;
    }
    // Calculate the slope until the kink
    double dx = TrackHitPos[longtrack][kink][1] - firstx;
    double dz = TrackHitPos[longtrack][kink][0] - TrackHitPos[longtrack][0][0];
    //std::cout << dz << std::endl;
    //double theta = atan(dx/dx);
    TVector3 mu(dx, 0, dz);
    mu = mu.Unit();
    //double recoangle = nuvect.Angle(mu)*180./3.1415;
    double recoangle = acos(mu.Z())*180./3.1415;

    // check true angle
    //TVector3 truemu(MuonP4[0], MuonP4[1], MuonP4[2]);
    TVector3 truemu(MuonP4[0], 0, MuonP4[2]);
    truemu = truemu.Unit();
    //double trueangle = nuvect.Angle(truemu)*180./3.1415;
    double trueangle = acos(truemu.Z())*180./3.1415;

    //std::cout << "true: " << trueangle << " reco: " << recoangle << std::endl;
    
    AnglevsTrue->Fill(trueangle, recoangle);
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

  TString canvname = Form("Angle_%s", filename.c_str());
  canvname += Form("_nLinesCut%i_nClusterCut%i_OccupancyCut%.2f_CCmuOnly%o_ClusterEnCut%.2f_AtLeastOneLine%o_AllDet%o", nLinesCut, nClustersCut, OccupancyCut, CCmuOnly, ClusterEnergyCut, AtLeastOneLine, AllDet);
  canvname += ".pdf";
  canv->Print(canvname+"[");

  gStyle->SetPalette(55);
  AnglevsTrue->Draw("colz");
  canv->Print(canvname);

  TH1D *arith = new TH1D("arith", "arith", AnglevsTrue->GetXaxis()->GetNbins(), AnglevsTrue->GetXaxis()->GetBinLowEdge(1), AnglevsTrue->GetXaxis()->GetBinLowEdge(AnglevsTrue->GetXaxis()->GetNbins()+1));
  // Now make the muon AnglevsTrue
  gStyle->SetOptStat(1111);
  for (int i = 0; i < AnglevsTrue->GetXaxis()->GetNbins(); ++i) {
    double center = AnglevsTrue->GetXaxis()->GetBinCenter(i);
    TH1D *proj = AnglevsTrue->ProjectionY(Form("AnglevsTrue %.2f", center), i, i);
    double mean = proj->GetMean();
    double mode = proj->GetBinCenter(proj->GetMaximumBin());
    double rms = proj->GetRMS();
    arith->SetBinContent(i+1, mode);
    arith->SetBinError(i+1, rms);
    proj->Draw();
    proj->SetStats(1);
    canv->Print(canvname);
  }

  TF1 *fit = new TF1("fit", "[0]+[1]*x", AnglevsTrue->GetYaxis()->GetBinLowEdge(1), AnglevsTrue->GetYaxis()->GetBinLowEdge(AnglevsTrue->GetYaxis()->GetNbins()+1));
  gStyle->SetOptStat(0);
  AnglevsTrue->Draw("colz");
  arith->Draw("same");
  arith->Fit(fit, "S", "", 700, 2500);
  TLegend *leg = new TLegend(0.2, 0.5, 0.6, 0.9);
  leg->AddEntry(arith, "Artihmetic mean and RMS", "le");
  leg->AddEntry(fit, Form("y=%.2f x + %.2f", fit->GetParameter(1), fit->GetParameter(0)), "l");
  fit->SetLineColor(kRed);
  arith->GetXaxis()->SetTitle("Muon True KE");
  leg->Draw("same");
  canv->Print(canvname);

  TString rootfilename = canvname;
  rootfilename = rootfilename.ReplaceAll(".pdf",".root");
  TFile *output = new TFile(rootfilename, "recreate");
  AnglevsTrue->Write();
  output->Close();

  canv->Print(canvname+"]");
}
