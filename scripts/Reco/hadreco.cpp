void hadreco() {
  TFile *f = new TFile("neutrino_merge_0p1Ecut_updtrue.root");
  TTree *truth = (TTree*)f->Get("Truth_Info");
  TTree *reco = (TTree*)f->Get("Line_Candidates");

  // Reco info
  const int MAX_LINES = 20;
  const int MAX_HITS = 1500;
  const int MAX_CLUSTERS = 20;
  int nLines;
  float TrackLength[MAX_LINES];
  float Occupancy[MAX_LINES];
  float FirstHoughHit[MAX_LINES][2];
  float LastHoughHit[MAX_LINES][2];
  bool TMSStart;
  float TotalTrackEnergy[MAX_LINES];

  float RecoHitEnergy[MAX_HITS];
  float ClusterEnergy[MAX_CLUSTERS];
  int nClusters;
  reco->SetBranchAddress("nLines", &nLines);
  reco->SetBranchAddress("TrackLength", TrackLength);
  reco->SetBranchAddress("Occupancy", Occupancy);
  reco->SetBranchAddress("FirstHoughHit", FirstHoughHit);
  reco->SetBranchAddress("LastHoughHit", LastHoughHit);
  reco->SetBranchAddress("TMSStart", TMSStart);
  reco->SetBranchAddress("RecoHitEnergy", RecoHitEnergy);
  reco->SetBranchAddress("ClusterEnergy", ClusterEnergy);
  reco->SetBranchAddress("nClusters", &nClusters);
  reco->SetBranchAddress("TotalTrackEnergy", TotalTrackEnergy);

  // Truth info
  float MuonP4[4];
  float Muon_Vertex[4];
  float Muon_Death[4];
  float Muon_TrueKE;
  int nParticles;
  int NeutrinoPDG;
  float NeutrinoP4[4];

  truth->SetBranchAddress("MuonP4", MuonP4);
  truth->SetBranchAddress("Muon_Vertex", Muon_Vertex);
  truth->SetBranchAddress("Muon_Death", Muon_Death);
  truth->SetBranchAddress("Muon_TrueKE", &Muon_TrueKE);
  truth->SetBranchAddress("nParticles", &nParticles);
  truth->SetBranchAddress("NeutrinoPDG", &NeutrinoPDG);
  truth->SetBranchAddress("NeutrinoP4", NeutrinoP4);

  TH2D *HadvsReco = new TH2D("Had", "E_{#nu}-E_{#mu} true vs cluster energy; true E_{#nu}-E_{#mu}; Total Cluster Energy", 10, 0, 2, 10, 0, 50);

  int nentries = truth->GetEntries();
  for (int i = 0; i < nentries; ++i) {
    truth->GetEntry(i);
    reco->GetEntry(i);

    // The event has to have a muon
    if (Muon_TrueKE < 0) continue;

    // Only include events with lines
    if (nLines < 1) continue;

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

    // And also check longest track
    float longest = 0;
    int longtrack = 0;
    for (int j = 0; j < nLines; ++j) {
      float xdist = LastHoughHit[j][0]-FirstHoughHit[j][0];
      float ydist = LastHoughHit[j][1]-FirstHoughHit[j][1];
      float dist = sqrt(xdist*xdist+ydist*ydist);
      if (dist > longest) longtrack = j;
    }
    longtrack = lon_trklen;

    // Look only at events with true muons that die inside the detector
    if (Muon_Death[1] > 1159 || Muon_Death[1] < -3864) continue;

    // Check that the longest track stops in the detector, and starts in the detector FV
    //if (FirstHoughHit[longtrack][0] < 11362+55*2) continue;
    if (FirstHoughHit[longtrack][0] < 11362+55*2 || FirstHoughHit[longtrack][0] > 13600) continue;
    if (LastHoughHit[longtrack][0] > 18294-80*2) continue;
    //if (LastHoughHit[longtrack][0] > 13600) continue;

    // 10 cm inwards
    if (fabs(FirstHoughHit[longtrack][1]) > 3520-200) continue;
    if (fabs(LastHoughHit[longtrack][1]) > 3520-200) continue;

    if (Occupancy[longtrack] < 0.5) continue;

    // Sum up the total cluster energy
    float cluster_en = 0;
    for (int j = 0; j < nClusters; ++j) {
      cluster_en += ClusterEnergy[j];
    }
    // Also sum up all the energy from tracks that aren't the muon candidates
    double trk_en = 0;
    for (int j = 0; j < nLines; ++j) {
      if (j == longtrack) continue;
      trk_en += TotalTrackEnergy[j];
    }
    cluster_en += trk_en;
    if (cluster_en < 2) continue;
    //std::cout << trk_en << " " << cluster_en << std::endl;
    //std::cout << NeutrinoP4[3] << " " << MuonP4[3]/1000. << std::endl;

    //std::cout << cluster_en << std::endl;
    HadvsReco->Fill(NeutrinoP4[3]-MuonP4[3]/1000., cluster_en);
  }
  gStyle->SetPalette(55);
  HadvsReco->Draw("colz");
}
