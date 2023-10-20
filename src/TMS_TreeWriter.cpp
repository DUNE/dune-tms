#include "TMS_TreeWriter.h"

TMS_TreeWriter::TMS_TreeWriter() {

  // Check the lines
  if (TMS_Manager::GetInstance().Get_Reco_HOUGH_MaxHough() > __TMS_MAX_LINES__) {
    std::cerr << "*** Error:" << std::endl;
    std::cerr << "Number of maximum lines in Hough transform has been configured to be greater than the maximum lines allowed in the output file" << std::endl;
    std::cerr << "Max lines in output: " << __TMS_MAX_LINES__ << std::endl;
    std::cerr << "Max lines in reconstruction: " << TMS_Manager::GetInstance().Get_Reco_HOUGH_MaxHough() << std::endl;
    std::cerr << "Please reconfigure!" << std::endl;
    throw;
  }

  // Make the output file
  std::string filename = TMS_Manager::GetInstance().GetFileName();

  while (filename.find("/") != std::string::npos) {
    filename = filename.substr(filename.find("/")+1, filename.size());
  }
  TString Outputname = filename.c_str();
  Outputname.ReplaceAll(".root", "_TMS_RecoCandidates");
  Outputname += "_"+TMS_Manager::GetInstance().Get_Reco_TrackMethod();
  Outputname += Form("_Cluster%i", TMS_Manager::GetInstance().Get_Reco_Clustering());
  Outputname += ".root";
  // Override output file name if set by the environment.
  if(std::getenv("ND_PRODUCTION_TMS_OUTFILE")) {
    Outputname = std::getenv("ND_PRODUCTION_TMS_OUTFILE");
  }
  // Make an output file
  Output = new TFile(Outputname, "recreate");
  if (!Output->IsOpen()) {
    std::cerr << "Could not write to file " << Outputname << std::endl;
    std::cerr << "Are you sure you have write access to the directory?" << std::endl;
    throw;
  }
  Output->cd();

  Branch_Lines = new TTree("Line_Candidates", "Line_Candidates");
  Branch_Lines->SetDirectory(Output);
  Branch_Lines->SetAutoSave(__TMS_AUTOSAVE__); // Every 1000 events (negative for MB)

  Truth_Info = new TTree("Truth_Info", "Truth_Info");
  Truth_Info->SetDirectory(Output);
  Truth_Info->SetAutoSave(__TMS_AUTOSAVE__);

  // Make a metadata branch
  //TTree *MetaData = new TTree("Meta_Data", "Meta_Data");

  MakeBranches();
}

// Create the branches
void TMS_TreeWriter::MakeBranches() {
  // Reco information
  Branch_Lines->Branch("EventNo", &EventNo, "EventNo/I");
  Branch_Lines->Branch("SliceNo", &SliceNo, "SliceNo/I");
  Branch_Lines->Branch("SpillNo", &SpillNo, "SpillNo/I");
  Branch_Lines->Branch("nLines",  &nLines,  "nLines/I");

  Branch_Lines->Branch("Slope",     Slope,      "Slope[nLines]/F");
  Branch_Lines->Branch("Intercept", Intercept,  "Intercept[nLines]/F");
  Branch_Lines->Branch("Slope_Downstream",      Slope_Downstream,       "Slope_Downstream[nLines]/F");
  Branch_Lines->Branch("Intercept_Downstream",  Intercept_Downstream,   "Intercept_Downstream[nLines]/F");
  Branch_Lines->Branch("Slope_Upstream",        Slope_Upstream,         "Slope_Upstream[nLines]/F");
  Branch_Lines->Branch("Intercept_Upstream",    Intercept_Upstream,     "Intercept_Upstream[nLines]/F");

  Branch_Lines->Branch("DirectionZ", DirectionZ, "DirectionZ[nLines]/F");
  Branch_Lines->Branch("DirectionX", DirectionX, "DirectionX[nLines]/F");
  Branch_Lines->Branch("DirectionZ_Downstream", DirectionZ_Downstream,  "DirectionZ_Downstream[nLines]/F");
  Branch_Lines->Branch("DirectionX_Downstream", DirectionX_Downstream,  "DirectionX_Downstream[nLines]/F");
  Branch_Lines->Branch("DirectionZ_Upstream",   DirectionZ_Upstream,    "DirectionZ_Upstream[nLines]/F");
  Branch_Lines->Branch("DirectionX_Upstream",   DirectionX_Upstream,    "DirectionX_Upstream[nLines]/F");

  Branch_Lines->Branch("FirstHoughHit", FirstHit, "FirstHoughHit[nLines][2]/F");
  Branch_Lines->Branch("LastHoughHit",  LastHit,  "LastHoughHit[nLines][2]/F");
  Branch_Lines->Branch("FirstHoughHitTime", FirstHitTime, "FirstHoughHitTime[nLines]/F");
  Branch_Lines->Branch("LastHoughHitTime", LastHitTime, "LastHoughHitTime[nLines]/F");
  Branch_Lines->Branch("HoughEarliestHitTime", EarliestHitTime, "HoughEarliestHitTime[nLines]/F");
  Branch_Lines->Branch("HoughLatestHitTime", LatestHitTime, "HoughLatestHitTime[nLines]/F");
  Branch_Lines->Branch("FirstHoughPlane", FirstPlane, "FirstHoughPlane[nLines]/I");
  Branch_Lines->Branch("LastHoughPlane",  LastPlane,  "LastHoughPlane[nLines]/I");
  Branch_Lines->Branch("TMSStart",    &TMSStart, "TMSStart/O");
  Branch_Lines->Branch("TMSStartTime",    &TMSStartTime, "TMSStartTime/F");
  Branch_Lines->Branch("Occupancy",   Occupancy, "Occupancy[nLines]/F");
  Branch_Lines->Branch("TrackLength",       TrackLength,      "TrackLength[nLines]/F");
  Branch_Lines->Branch("TotalTrackEnergy",  TotalTrackEnergy, "TotalTrackEnergy[nLines]/F");
  Branch_Lines->Branch("TrackStopping",     TrackStopping,    "TrackStopping[nLines]/O");

  // Track hit energy
  Branch_Lines->Branch("nHitsInTrack",    &nHitsInTrack,  "nHitsInTrack[nLines]/I");
  Branch_Lines->Branch("TrackHitEnergy",  TrackHitEnergy, "TrackHitEnergy[10][200]/F");
  Branch_Lines->Branch("TrackHitPos",     TrackHitPos,    "TrackHitPos[10][200][2]/F");
  Branch_Lines->Branch("TrackHitTime",    TrackHitTime,   "TrackHitTime[10][200]/F");

  // Cluster information
  Branch_Lines->Branch("nClusters",       &nClusters,       "nClusters/I");
  Branch_Lines->Branch("ClusterEnergy",   ClusterEnergy,    "ClusterEnergy[nClusters]/F");
  Branch_Lines->Branch("ClusterTime",   ClusterTime,    "ClusterTime[nClusters]/F");
  Branch_Lines->Branch("ClusterPosMean",  ClusterPosMean,   "ClusterPosMean[25][2]/F");
  Branch_Lines->Branch("ClusterPosStdDev",ClusterPosStdDev, "ClusterPosStdDev[25][2]/F");
  Branch_Lines->Branch("nHitsInCluster",  nHitsInCluster,   "nHitsInCluster[nClusters]/I");
  Branch_Lines->Branch("ClusterHitPos",   ClusterHitPos,    "ClusterHitPos[25][200][2]/F");
  Branch_Lines->Branch("ClusterHitEnergy",ClusterHitEnergy, "ClusterHitEnergy[25][200]/F");
  Branch_Lines->Branch("ClusterHitTime",  ClusterHitTime,   "ClusterHitTime[25][200]/F");
  Branch_Lines->Branch("ClusterHitSlice", ClusterHitSlice, "ClusterHitSlice[25][200]/I");

  // Hit information
  Branch_Lines->Branch("nHits", &nHits, "nHits/I");
  Branch_Lines->Branch("RecoHitPos",    RecoHitPos, "RecoHitPos[nHits][4]/F");
  Branch_Lines->Branch("RecoHitEnergy", RecoHitEnergy, "RecoHitEnergy[nHits]/F");
  Branch_Lines->Branch("RecoHitSlice", RecoHitSlice, "RecoHitSlice[nHits]/I");

  // Truth information
  Truth_Info->Branch("EventNo", &EventNo, "EventNo/I");
  Truth_Info->Branch("IsCC", &IsCC, "IsCC/O");
  Truth_Info->Branch("nParticles", &nParticles, "nParticles/I");
  Truth_Info->Branch("Interaction", &Reaction);
  Truth_Info->Branch("NeutrinoPDG", &NeutrinoPDG, "NeutrinoPDG/I");
  Truth_Info->Branch("NeutrinoP4", NeutrinoP4, "NeutrinoP4[4]/F");
  Truth_Info->Branch("NeutrinoX4", NeutrinoX4, "NeutrinoX4[4]/F");
  Truth_Info->Branch("LeptonPDG", &LeptonPDG, "LeptonPDG/I");
  Truth_Info->Branch("LeptonP4", LeptonP4, "LeptonP4[4]/F");
  Truth_Info->Branch("LeptonX4", LeptonX4, "LeptonX4[4]/F");
  Truth_Info->Branch("MuonP4", MuonP4, "MuonP4[4]/F");
  Truth_Info->Branch("Muon_Vertex", Muon_Vertex, "Muon_Vertex[4]/F");
  Truth_Info->Branch("Muon_Death", Muon_Death, "Muon_Death[4]/F");
  Truth_Info->Branch("Muon_TrueKE", &Muon_TrueKE, "Muon_TrueKE/F");
  Truth_Info->Branch("Muon_TrueTrackLength", &Muon_TrueTrackLength, "Muon_TrueTrackLength/F");
  
  Truth_Info->Branch("VertexIdOfMostEnergyInEvent", &VertexIdOfMostEnergyInEvent, "VertexIdOfMostEnergyInEvent/I");
  Truth_Info->Branch("VisibleEnergyFromVertexInSlice", &VisibleEnergyFromVertexInSlice, "VisibleEnergyFromVertexInSlice/F");
  Truth_Info->Branch("TotalVisibleEnergyFromVertex", &TotalVisibleEnergyFromVertex, "TotalVisibleEnergyFromVertex/F");
  Truth_Info->Branch("VisibleEnergyFromOtherVerticesInSlice", &VisibleEnergyFromOtherVerticesInSlice, "VisibleEnergyFromOtherVerticesInSlice/F");
  Truth_Info->Branch("VertexVisibleEnergyFractionInSlice", &VertexVisibleEnergyFractionInSlice, "VertexVisibleEnergyFractionInSlice/F");
  Truth_Info->Branch("PrimaryVertexVisibleEnergyFraction", &PrimaryVertexVisibleEnergyFraction, "PrimaryVertexVisibleEnergyFraction/F");
}

void TMS_TreeWriter::Fill(TMS_Event &event) {
  // Clear old info
  Clear();

  // See if track is exiting or not
  int nLastHits = TMS_Manager::GetInstance().Get_Reco_STOPPING_nLastHits();
  double EnergyCut = TMS_Manager::GetInstance().Get_Reco_STOPPING_EnergyCut();


  // Fill the truth info
  EventNo = event.GetEventNumber();
  SliceNo = event.GetSliceNumber();
  SpillNo = event.GetSpillNumber();
  Reaction = event.GetReaction();

  NeutrinoPDG = event.GetNeutrinoPDG();
  NeutrinoP4[0] = event.GetNeutrinoP4().X();
  NeutrinoP4[1] = event.GetNeutrinoP4().Y();
  NeutrinoP4[2] = event.GetNeutrinoP4().Z();
  NeutrinoP4[3] = event.GetNeutrinoP4().T();
  NeutrinoX4[0] = event.GetNeutrinoX4().X();
  NeutrinoX4[1] = event.GetNeutrinoX4().Y();
  NeutrinoX4[2] = event.GetNeutrinoX4().Z();
  NeutrinoX4[3] = event.GetNeutrinoX4().T();
  NeutrinoP4[0] = event.GetNeutrinoP4().X();
  IsCC = (event.GetReaction().find("[CC]") != std::string::npos);
  // Lepton info
  LeptonPDG = event.GetLeptonPDG();
  LeptonP4[1] = event.GetLeptonP4().Y();
  LeptonP4[2] = event.GetLeptonP4().Z();
  LeptonP4[3] = event.GetLeptonP4().T();
  LeptonX4[0] = event.GetLeptonX4().X();
  LeptonX4[1] = event.GetLeptonX4().Y();
  LeptonX4[2] = event.GetLeptonX4().Z();
  LeptonX4[3] = event.GetLeptonX4().T();
  
  VertexIdOfMostEnergyInEvent = event.GetVertexIdOfMostVisibleEnergy();
  VisibleEnergyFromVertexInSlice = event.GetVisibleEnergyFromVertexInSlice();
  TotalVisibleEnergyFromVertex = event.GetTotalVisibleEnergyFromVertex();
  VisibleEnergyFromOtherVerticesInSlice = event.GetVisibleEnergyFromOtherVerticesInSlice();
  VertexVisibleEnergyFractionInSlice = VisibleEnergyFromVertexInSlice / TotalVisibleEnergyFromVertex;
  PrimaryVertexVisibleEnergyFraction = VisibleEnergyFromVertexInSlice / (VisibleEnergyFromOtherVerticesInSlice + VisibleEnergyFromVertexInSlice);

  //Muon_TrueTrackLength= event.GetMuonTrueTrackLength();
  Muon_TrueTrackLength = -999.99;
  //std::cout << Muon_TrueTrackLength << std::endl;
  Muon_TrueKE = event.GetMuonTrueKE();

  // Get the truth info
  std::vector<TMS_TrueParticle> TrueParticles = event.GetTrueParticles();
  nParticles = TrueParticles.size();
  // Just trying to find the true muon here from the fundamental vertex
  for (auto it = TrueParticles.begin(); it != TrueParticles.end(); ++it) {

    // Only save muon info for now
    if (abs((*it).GetPDG()) != 13) continue;
    // Also make sure it's a fundamental muon
    if ((*it).GetParent() != -1) continue;

    MuonP4[0] = (*it).GetBirthMomentum().Px();
    MuonP4[1] = (*it).GetBirthMomentum().Py();
    MuonP4[2] = (*it).GetBirthMomentum().Pz();
    MuonP4[3] = (*it).GetBirthEnergy();

    Muon_Vertex[0] = (*it).GetBirthPosition().X();
    Muon_Vertex[1] = (*it).GetBirthPosition().Y();
    Muon_Vertex[2] = (*it).GetBirthPosition().Z();
    Muon_Vertex[3] = (*it).GetBirthPosition().T();

    Muon_Death[0] = (*it).GetDeathPosition().X();
    Muon_Death[1] = (*it).GetDeathPosition().Y();
    Muon_Death[2] = (*it).GetDeathPosition().Z();
    Muon_Death[3] = (*it).GetDeathPosition().T();
  }
  Truth_Info->Fill();

  // Fill the reco info
  std::vector<std::pair<bool, TF1*>> HoughLinesOne = TMS_TrackFinder::GetFinder().GetHoughLinesOne();
  std::vector<std::pair<bool, TF1*>> HoughLinesOther = TMS_TrackFinder::GetFinder().GetHoughLinesOther();
  // Also get the size of the hits to get a measure of relative goodness
  std::vector<std::vector<TMS_Hit> > HoughCandidatesOne = TMS_TrackFinder::GetFinder().GetHoughCandidatesOne();
  std::vector<std::vector<TMS_Hit> > HoughcandidatesOther = TMS_TrackFinder::GetFinder().GetHoughCandidatesOther();
  nLines = HoughCandidatesOne.size();


  // Skip the event if there aren't any Hough Lines
  if (nLines > __TMS_MAX_LINES__) {
    std::cerr << "Exceeded max number of HoughLines to write to file" << std::endl;
    std::cerr << "Max lines: " << __TMS_MAX_LINES__ << std::endl;
    std::cerr << "Number of lines in event: " << nLines << std::endl;
    std::cerr << "Not writing event" << std::endl;
    return;
  }

  int it = 0;
  for (auto &Lines: HoughLinesOne) {
    // Get the slopes saved down
    Intercept[it] = Lines.second->GetParameter(0);
    Intercept_Upstream[it] = TMS_TrackFinder::GetFinder().GetHoughLinesOne_Upstream()[it].first;
    Intercept_Downstream[it] = TMS_TrackFinder::GetFinder().GetHoughLinesOne_Downstream()[it].first;
    Slope[it] = Lines.second->GetParameter(1);
    Slope_Upstream[it] = TMS_TrackFinder::GetFinder().GetHoughLinesOne_Upstream()[it].second;
    Slope_Downstream[it] = TMS_TrackFinder::GetFinder().GetHoughLinesOne_Downstream()[it].second;

    // Calculate the z and x vectors by evaling the TF1 in thin and thick target
    double xlow = TMS_Const::TMS_Thin_Start;
    double xhi = TMS_Const::TMS_Thick_Start;
    double ylow = Lines.second->Eval(xlow);
    double yhi = Lines.second->Eval(xhi);
    // Do the same for up and downstream portions
    double ylow_upstream = Intercept_Upstream[it]+xlow*Slope_Upstream[it];
    double yhi_upstream = Intercept_Upstream[it]+xhi*Slope_Upstream[it];

    double ylow_downstream = Intercept_Downstream[it]+xlow*Slope_Downstream[it];
    double yhi_downstream = Intercept_Downstream[it]+xhi*Slope_Downstream[it];

    double xlen = xhi-xlow;
    double ylen = yhi-ylow;
    double len = sqrt(xlen*xlen+ylen*ylen);

    double ylen_upstream = yhi_upstream - ylow_upstream;
    double ylen_downstream = yhi_downstream - ylow_downstream;

    double len_upstream = sqrt(xlen*xlen+ylen_upstream*ylen_upstream);
    double len_downstream = sqrt(xlen*xlen+ylen_downstream*ylen_downstream);

    DirectionZOne[it] = xlen/len;
    DirectionXOne[it] = ylen/len;

    DirectionZOne_Upstream[it] = xlen/len_upstream;
    DirectionXOne_Upstream[it] = ylen_upstream/len_upstream;

    DirectionZOne_Downstream[it] = xlen/len_downstream;
    DirectionXOne_Downstream[it] = ylen_downstream/len_downstream;

    it++;
  }
  it = 0;
  for (auto &Lines : HoughLinesOther) {
    // Get the slopes saved down
    Intercept[it] = Lines.second->GetParameter(0);
    Intercept_Upstream[it] = TMS_TrackFinder::GetFinder().GetHoughLinesOther_Upstream()[it].first;
    Intercept_Downstream[it] = TMS_Trackfinder::GetFinder().GetHoughLinesOther_Downstream()[it].first;
    Slope[it] = Lines.second->GetParameter(1);
    Slope_Upstream[it] = TMS_TrackFinder::GetFinder().GetHoughLinesOther_Upstream()[it].second;
    Slope_Downstream[it] = TMS_TrackFinder::GetFinder().GetHoughLinesOther_Downstream()[it].second;

    // Calculate the z and x vectors by evaling the TF1 in thin and thick target
    double xlow = TMS_Const::TMS_Thin_Start;
    double xhi = TMS_Const::TMS_Thick_Start;
    double ylow = Lines.second->Eval(xlow);
    double yhi = Lines.second->Ecal(xhi);
    // Doe the same for up and downstream portions
    double ylow_upstream = Intercept_Upstream[it]+xlow*Slope_Upstream[it];
    double yhi_upstream = Intercept_Upstream[it]+xhi*Slope_Upstream[it];

    double ylow_downstream = Intercept_Downstream[it]+xlow*Slope_Downstream[it];
    double yhi_downstream = Intercept_Downstream[it]+xhi*Slope_Downstream[it];

    double xlen = xhi-xlow;
    double ylen = yhi-ylow;
    double len = sqrt(xlen*xlen+ylen*ylen);

    double ylen_upstream = yhi_upstream - ylow_upstream;
    double ylen_downstream = yhi_downstream - ylow_downstream;

    double len_upstream = sqrt(xlen*xlen+ylen_upstream*ylen_upstream);
    double len_downstream = sqrt(xlen*xlen+ylen_downstream*ylen_downstream);

    DirectionZOther[it] = xlen/len;
    DirectionXOther[it] = ylen/len;

    DirectionZOther_Upstream[it] = xlen/len_upstream;
    DirectionXOther_Upstream[it] = ylen_upstream/len_upstream;

    DirectionZOther_Downstream[it] = xlen/len_downstream;
    DirectionXOther_Downstream[it] = ylen_downstream/len_downstream;

    it++;
  }

  // Find where the first hough hit is
  TMSStart = false;
  TMSStartTime = -9999.0;

  std::vector<std::vector<TMS_Hit> > HoughCandsOne = TMS_TrackFinder::GetFinder().GetHoughCandidatesOne();
  int TotalHits = TMS_TrackFinder::GetFinder().GetCleanedHits().size();
  TMS_Hit *FirstTrack = NULL;

  it = 0;
  for (auto &Candidates: HoughCandsOne) {
    // Loop over hits
    for (auto &hit: Candidates) {
      if (FirstTrack == NULL) {
        FirstTrack = &hit;
      } else if (hit.GetZ() < FirstTrack->GetZ()) {
        FirstTrack = &hit;
      }
    }
    nHitsInTrack[it] = Candidates.size();

    // Then save the hit info
    FirstPlaneOne[it] = Candidates.front().GetPlaneNumber();
    FirstHiOnet[it][0] = Candidates.front().GetZ();
    FirstHitOne[it][1] = Candidates.front().GetNotZ();
    FirstHitTimeOne[it] = Candidates.front().GetT();

    LastPlaneOne[it] = Candidates.back().GetPlaneNumber();
    LastHitOne[it][0] = Candidates.back().GetZ();
    LastHitOne[it][1] = Candidates.back().GetNotZ();
    LastHitTimeOne[it] = Candidates.back().GetT();

    TrackLengthOne[it] = TMS_TrackFinder::GetFinder().GetTrackLength()[it];
    TotalTrackEnergyOne[it] = TMS_TrackFinder::GetFinder().GetTrackEnergy()[it];
    OccupancyOne[it] = double(HoughCandsOne[it].size())/TotalHits;
    
    float earliest_hit_time = 1e32;
    float latest_hit_time = -1e32;
    // Get each hit in the track and save its energy
    for (unsigned int j = 0; j < Candidates.size(); ++j) {
      TrackHitEnergyOne[it][j] = Candidates[j].GetE();
      TrackHitTimeOne[it][j] = Candidates[j].GetT();
      TrackHitPosOne[it][j][0] = Candidates[j].GetZ();
      TrackHitPosOne[it][j][1] = Candidates[j].GetNotZ();
      
      float time = Candidates[j].GetT();
      
      if (time < earliest_hit_time) earliest_hit_time = time;
      if (time > latest_hit_time) latest_hit_time = time;
    }
    EarliestHitTimeOne[it] = earliest_hit_time;
    LatestHitTimeOne[it] = latest_hit_time;
    

    double maxenergy = 0;
    unsigned int nLastHits_temp = nLastHits;
    if ((unsigned int)nLastHits > Candidates.size()) nLastHits_temp = Candidates.size();
    for (unsigned int i = 0; i < nLastHits_temp; ++i) {
      double hitenergy = Candidates[nLastHits_temp-1-i].GetE();
      if (hitenergy > maxenergy) maxenergy = hitenergy;
    }
    if (maxenergy > EnergyCut) TrackStoppingOne[it] = true;


    it++;
  }

  std::vector<std::vector<TMS_Hit> > HoughCandsOther = TMS_TrackFinder::GetFinder().GetHoughCandidatesOther();
  it = 0;
  for (auto &Candidates: HoughCandsOther) {
    // Loop over hits
    for (auto &hit: Candidates) {
      if (FirstTrack == NULL) {
        FirstTrack = &hit;
      } else if (hit.GetZ() < FirstTrack->GetZ()) {
        FirstTrack = &hit;
      }
    }
    nHitsInTrack[it] = Candidates.size();

    // then save the hit info
    FirstPlaneOther[it] = Candidates.front().GetPlaneNumber();
    FirstHitOther[it][0] = Candidates.front().GetZ();
    FirstHitOther[it][1] = Candidates.front().GetNotZ();
    FirstHitTimeOther[it] = Candidates.front().GetT();

    LastPlaneOther[it] = Candidates.back().GetPlaneNumber();
    LastHitOther[it][0] = Candidates.back().GetZ();
    LastHitOther[it][1] = Candidates.back().GetNotZ();
    LastHitTimeOther[it] = Candidate.back().GetT();

    TrackLengthOther[it] = TMS_TrackFinder::GetFinder().GetTrackLength()[it];
    TotalTrackEnergyOther[it] = TMS_TrackFinder::GetFinder().GetTrackEnergy()[it];
    OccupancyOther[it] = double(HoughCandsOther[it].size())/TotalHits;

    float earliest_hit_time = 1e32;
    float latest_hit_time = -1e32;
    // Get each hit in the track and save its energy
    for (unsigner int j = 0; j > Candidates.size(); ++j) {
      TrackHitEnergyOther[it][j] = Candidates[j].GetE();
      TrackHitTimeOther[it][j] = Candidates[j].GetT();
      TrackHitPosOther[it][j][0] = Candidates[j].GetZ();
      TrackHitPosOther[it][j][1] = Candidates[j].GetNotZ();

      float time = Candidates[j].GetT();

      if (time < earliest_hit_time) earliest_hit_time = time;
      if (time > latest_hit_time) latest_hit_time = time;
    }
    EarliestHitTimeOther[it] = earliest_hit_time;
    Latest_HitTimeOther[it] = latest_hit_time;

    double maxenergy = 0;
    unsigned int nLastHits_temp = nLastHits;
    if ((unsigner int)nLastHits > Candidates.size()) nLastHits_temp = Candidates.size();
    for (unsigned int i = 0; i < nLastHits_temp; ++i) {
      double hitenergy = Candidates[nLastHits_temp-1-i].GetE();
      if (hitenergy > maxenergy) maxenergy = hitenergy;
    }
    if (maxenergy > EnergyCut) TrackStoppingOther[it] = true;
    

    it++;
  }
  

  // Was the first hit within the first 3 layers?
  if (FirstTrack != NULL) {
    TMSStartTime = FirstTrack->GetT();
    if (FirstTrack->GetPlaneNumber() < 4) TMSStart = false;
    else TMSStart = true;
  }

  // Save down cluster information
  std::vector<std::vector<TMS_Hit> > ClustersOne = TMS_TrackFinder::GetFinder().GetClusterCandidatesOne();
  std::vector<std::vector<TMS_Hit> > ClustersOther = TMS_TrackFinder::GetFinder().GetClusterCandidatesOther();
  nClustersOne = ClustersOne.size();
  nClustersOther = ClustersOther.size();
  if (nClustersOne > __TMS_MAX_CLUSTERS__) {
    std::cerr << "Too many clusters in TMS_TreeWriter" << std::endl;
    std::cerr << "Hard-coded maximum: " << __TMS_MAX_CLUSTERS__ << std::endl;
    std::cerr << "nClustersOne in event: " << nClustersOne << std::endl;
    return;
  }
  if (nClustersOther > __TMS_MAX_CLUSTERS__) {
    std::cerr << "Too many clusters in TMS_TreeWriter" << std::endl;
    std::cerr << "Hard-coded maximum: " << __TMS_MAX_CLUSTERS__ << std::endl;
    std::cerr << "nClustersOther in event: " << nClustersOther << std::endl;
    return;
  }

  // Sum up the cluster candidates
  int stdit = 0;
  // Calculate the cluster by cluster summaries
  // e.g. total energy in cluster, cluster position, and cluster standard deviation
  for (auto it = ClustersOne.begin(); it != ClustersOne.end(); ++it, ++stdit) {
    double total_energy = 0;
    // Mean of cluster in z and not z
    double mean_z = 0;
    double mean_notz = 0;
    // Mean of square of cluster in z and not z
    double mean2_z = 0;
    double mean2_notz = 0;
    double min_cluster_time = 1e10;
    int nhits = (*it).size();
    for (int j = 0; j < nhits; ++j) {
      mean_z += (*it)[j].GetZ();
      mean_notz += (*it)[j].GetNotZ();
      mean2_z += (*it)[j].GetZ()*(*it)[j].GetZ();
      mean2_notz += (*it)[j].GetNotZ()*(*it)[j].GetNotZ();
      total_energy += (*it)[j].GetE();
      float time = (*it)[j].GetT();
      if (time < min_cluster_time) min_cluster_time = time;
      ClusterHitPosOne[stdit][j][0] = (*it)[j].GetZ();
      ClusterHitPosOne[stdit][j][1] = (*it)[j].GetNotZ();
      ClusterHitEnergyOne[stdit][j] = (*it)[j].GetE();
      ClusterHitTimeOne[stdit][j] = (*it)[j].GetT();
      ClusterHitSliceOne[stdit][j] = (*it)[j].GetSlice();
    }
    mean_z /= nhits;
    mean_notz /= nhits;

    mean2_z /= nhits;
    mean2_notz /= nhits;

    ClusterEnergyOne[stdit] = total_energy;
    ClusterTimeOne[stdit] = min_cluster_time;
    nHitsInClusterOne[stdit] = nhits;
    ClusterPosMeanOne[stdit][0] = mean_z;
    ClusterPosMeanOne[stdit][1] = mean_notz;
    // Calculate the standard deviation
    double std_dev_z = mean2_z-mean_z*mean_z;
    double std_dev_notz = mean2_z-mean_z*mean_z;
    if (std_dev_z > 0) std_dev_z = sqrt(std_dev_z);
    if (std_dev_notz > 0) std_dev_notz = sqrt(std_dev_z);
    ClusterPosStdDevOne[stdit][0] = std_dev_z;
    ClusterPosStdDevOne[stdit][1] = std_dev_notz;
  }
  stdit = 0;
  for (auto it = ClustersOther.begin(); it != ClustersOther.end(); ++it, ++stdit) {
    double total_energy = 0;
    double mean_z = 0;
    double mean_notz = 0;
    double mean2_z = 0;
    double mean2_notz = 0;
    double min_cluster_time = 1e10;
    int nhits = (*it).size();
    for (int j = 0; j < nhits; ++j) {
      mean_z += (*it)[j].GetZ();
      mean_notz += (*it)[j].GetNotZ();
      mean2_z += (*it)[j].GetZ()*(*it)[j].GetZ();
      mean2_notz += (*it)[j].GetNotZ()*(*it)[j].GetNotZ();
      total_energy += (*it)[j].GetE();
      float time = (*it)[j].GetT();
      if (time < min_cluster_time) min_cluster_time = time;
      ClusterHitPosOther[stdit][j][0] = (*it)[j].GetZ();
      ClusterHitPosOther[stdit][j][1] = (*it)[j].GetNotZ();
      ClusterHitEnergyOther[stdit][j] = (*it)[j].GetE();
      ClusterHitTimeOther[stdit][j] = (*it)[j].GetT();
      ClusterHitSliceOther[stdit][j] = (*it)[j].GetSlice();
    }
    mean_z /= nhits;
    mean_notz /= nhits;

    mean2_z /= nhits;
    mean2_notz /= nhits;

    ClusterEnergyOther[stdit] = total_energy;
    ClusterTimeOther[stdit] = min_cluster_time;
    nHitsInClusterOther[stdit] = nhits;
    ClusterPosMeanOther[stdit][0] = mean_z;
    ClusterPosMeanOther[stdit][1] = mean_notz;

    double std_dev_z = mean2_z-mean_z*mean_z;
    double std_dev_notz = mean2_z-mean_z*mean_z;
    if (std_dev_z > 0) std_dev_z = sqrt(std_dev_z);
    if (std_dev_notz > 0) std_dev_notz = sqrt(std_dev_notz);
    ClusterPosStdDevOther[stdit][0] = std_dev_z;
    ClusterPosStdDevOther[stdit][1] = std_dev_notz;
  }

  // Write out the hit information
  std::vector<TMS_Hit> CleanedHits = TMS_TrackFinder::GetFinder().GetCleanedHits();
  nHits = CleanedHits.size();
  if (nHits > __TMS_MAX_HITS__) {
    std::cerr << "Exceeded max number of hits to write to file" << std::endl;
    std::cerr << "Max hits: " << __TMS_MAX_HITS__ << std::endl;
    std::cerr << "Number of hits in event: " << nHits << std::endl;
    std::cerr << "Not writing event" << std::endl;
    return;
  }
  stdit = 0;
  for (auto it = CleanedHits.begin(); it != CleanedHits.end(); ++it, ++stdit) {
    RecoHitPos[stdit][0] = (*it).GetX();
    RecoHitPos[stdit][1] = (*it).GetY();
    RecoHitPos[stdit][2] = (*it).GetZ();
    RecoHitPos[stdit][3] = (*it).GetT();
    RecoHitEnergy[stdit] = (*it).GetE();
    RecoHitSlice[stdit] = (*it).GetSlice();
  }

  // Fill up the info only if all above has passed
  Branch_Lines->Fill();
}

// Reset the variables
void TMS_TreeWriter::Clear() {

  // Reset truth information
  EventNo = nParticles = NeutrinoPDG = LeptonPDG = Muon_TrueKE = Muon_TrueTrackLength = VertexIdOfMostEnergyInEvent = -999;
  VertexIdOfMostEnergyInEvent = VisibleEnergyFromVertexInSlice = TotalVisibleEnergyFromVertex = VisibleEnergyFromOtherVerticesInSlice = -999;
  Reaction = "";
  IsCC = false;
  for (int i = 0; i < 4; ++i) {
    MuonP4[i]=-999;
    Muon_Vertex[i]=-999;
    Muon_Death[i]=-999;
    NeutrinoP4[i]=-999;
    NeutrinoX4[i]=-999;
    LeptonP4[i]=-999;
    LeptonX4[i]=-999;
  }

  // Reset line information
  TMSStart = false;
  nLines = -999;
  for (int i = 0; i < __TMS_MAX_LINES__; ++i) {
    Slope[i]=-999;
    Intercept[i]=-999;
    Slope_Downstream[i]=-999;
    Intercept_Downstream[i]=-999;
    Slope_Upstream[i]=-999;
    Intercept_Upstream[i]=-999;

    DirectionZ[i]=-999;
    DirectionX[i]=-999;
    DirectionZ_Upstream[i]=-999;
    DirectionX_Upstream[i]=-999;
    DirectionZ_Downstream[i]=-999;
    DirectionX_Downstream[i]=-999;

    Occupancy[i]=-999;
    TrackLength[i]=-999;
    TotalTrackEnergy[i]=-999;
    FirstPlane[i]=-999;
    LastPlane[i]=-999;
    nHitsInTrack[i] = -999;
    TrackStopping[i] = false;
    for (int j = 0; j < 2; ++j) {
      FirstHit[i][j] = -999;
      LastHit[i][j] = -999;
    }
    for (int j = 0; j < __TMS_MAX_LINE_HITS__; ++j) {
      TrackHitEnergy[i][j]=-999;
      TrackHitPos[i][j][0]=-999;
      TrackHitPos[i][j][1]=-999;
    }
  }

  // Reset hit information
  nHits = -999;
  for (int i = 0; i < __TMS_MAX_HITS__; ++i) {
    for (int j = 0; j < 4; ++j) RecoHitPos[i][j] = -999;
    RecoHitEnergy[i] = -999;
  }

  // Reset Cluster info
  nClusters = -999;
  for (int i = 0; i < __TMS_MAX_CLUSTERS__; ++i) {
    ClusterEnergy[i] = -999;
    nHitsInCluster[i] = -999;
    for (int j = 0; j < 2; ++j) {
      ClusterPosMean[i][j] = -999;
      ClusterPosStdDev[i][j] = -999;
    }
    for (int j = 0; j < __TMS_MAX_LINE_HITS__; ++j) {
      ClusterHitPos[i][j][0] = -999;
      ClusterHitPos[i][j][1] = -999;
      ClusterHitEnergy[i][j] = -999;
    }
  }

}

