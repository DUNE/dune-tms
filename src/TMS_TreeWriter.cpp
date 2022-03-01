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
  Branch_Lines->Branch("nLines", &nLines, "nLines/I");
  Branch_Lines->Branch("Slope", Slope, "Slope[nLines]/F");
  Branch_Lines->Branch("Intercept", Intercept, "Intercept[nLines]/F");
  Branch_Lines->Branch("DirectionZ", DirectionZ, "DirectionZ[nLines]/F");
  Branch_Lines->Branch("DirectionX", DirectionX, "DirectionX[nLines]/F");
  Branch_Lines->Branch("FirstHoughHit", FirstHit, "FirstHoughHit[nLines][2]/F");
  Branch_Lines->Branch("LastHoughHit", LastHit, "LastHoughHit[nLines][2]/F");
  Branch_Lines->Branch("FirstHoughPlane", FirstPlane, "FirstHoughPlane[nLines]/I");
  Branch_Lines->Branch("LastHoughPlane", LastPlane, "LastHoughPlane[nLines]/I");
  Branch_Lines->Branch("TMSStart", &TMSStart, "TMSStart/O");
  Branch_Lines->Branch("Occupancy", Occupancy, "Occupancy[nLines]/F");
  Branch_Lines->Branch("TrackLength", TrackLength, "TrackLength[nLines]/F");
  Branch_Lines->Branch("TotalTrackEnergy", TotalTrackEnergy, "TotalTrackEnergy[nLines]/F");

  // Track hit energy
  Branch_Lines->Branch("nHitsInTrack", &nHitsInTrack, "nHitsInTrack[nLines]/I");
  Branch_Lines->Branch("TrackHitEnergy", TrackHitEnergy, "TrackHitEnergy[10][200]/F");
  Branch_Lines->Branch("TrackHitPos", TrackHitPos, "TrackHitPos[10][200][2]/F");

  // Cluster information
  Branch_Lines->Branch("nClusters", &nClusters, "nClusters/I");
  Branch_Lines->Branch("ClusterEnergy", ClusterEnergy, "ClusterEnergy[nClusters]/F");
  Branch_Lines->Branch("ClusterPosMean", ClusterPosMean, "ClusterPosMean[25][2]/F");
  Branch_Lines->Branch("ClusterPosStdDev", ClusterPosStdDev, "ClusterPosStdDev[25][2]/F");
  Branch_Lines->Branch("nHitsInCluster", nHitsInCluster, "nHitsInCluster[nClusters]/I");
  Branch_Lines->Branch("ClusterHitPos", ClusterHitPos, "ClusterHitPos[25][200][2]/F");
  Branch_Lines->Branch("ClusterHitEnergy", ClusterHitEnergy, "ClusterHitEnergy[25][200]/F");

  // Hit information
  Branch_Lines->Branch("nHits", &nHits, "nHits/I");
  Branch_Lines->Branch("RecoHitPos", RecoHitPos, "RecoHitPos[1000][4]/F");
  Branch_Lines->Branch("RecoHitEnergy", RecoHitEnergy, "RecoHitEnergy[1000]/F");

  // Truth information
  Truth_Info->Branch("EventNo", &EventNo, "EventNo/I");
  Truth_Info->Branch("IsCC", &IsCC, "IsCC/O");
  Truth_Info->Branch("nParticles", &nParticles, "nParticles/I");
  Truth_Info->Branch("Interaction", &Reaction);
  Truth_Info->Branch("NeutrinoPDG", &NeutrinoPDG, "NeutrinoPDG/I");
  Truth_Info->Branch("NeutrinoP4", NeutrinoP4, "NeutrinoP4[4]/F");
  Truth_Info->Branch("MuonP4", MuonP4, "MuonP4[4]/F");
  Truth_Info->Branch("Muon_Vertex", Muon_Vertex, "Muon_Vertex[4]/F");
  Truth_Info->Branch("Muon_Death", Muon_Death, "Muon_Death[4]/F");
  Truth_Info->Branch("Muon_TrueKE", &Muon_TrueKE, "Muon_TrueKE/F");
  Truth_Info->Branch("Muon_TrueTrackLength", &Muon_TrueTrackLength, "Muon_TrueTrackLength/F");
}

void TMS_TreeWriter::Fill(TMS_Event &event) {
  // Clear old info
  Clear();

  // Fill the truth info
  EventNo = event.GetEventNumber();
  Reaction = event.GetReaction();

  NeutrinoPDG = event.GetNeutrinoPDG();
  NeutrinoP4[0] = event.GetNeutrinoP4().X();
  NeutrinoP4[1] = event.GetNeutrinoP4().Y();
  NeutrinoP4[2] = event.GetNeutrinoP4().Z();
  NeutrinoP4[3] = event.GetNeutrinoP4().T();
  IsCC = (event.GetReaction().find("[CC]") != std::string::npos);

  Muon_TrueTrackLength= event.GetMuonTrueTrackLength();
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
  std::vector<std::pair<bool, TF1*>> HoughLines = TMS_TrackFinder::GetFinder().GetHoughLines();
  // Also get the size of the hits to get a measure of relative goodness
  std::vector<std::vector<TMS_Hit> > HoughCandidates = TMS_TrackFinder::GetFinder().GetHoughCandidates();
  nLines = HoughCandidates.size();
  // Get the cleaned hits for a reference of the "total"

  // Skip the event if there aren't any Hough Lines
  if (nLines > __TMS_MAX_LINES__) {
    std::cerr << "Exceeded max number of HoughLines to write to file" << std::endl;
    std::cerr << "Max lines: " << __TMS_MAX_LINES__ << std::endl;
    std::cerr << "Number of lines in event: " << nLines << std::endl;
    std::cerr << "Not writing event" << std::endl;
    return;
  }

  int it = 0;
  for (auto &Lines: HoughLines) {
    Intercept[it] = Lines.second->GetParameter(0);
    Slope[it] = Lines.second->GetParameter(1);

    // Calculate the z and x vectors by evaling the TF1 in thin and thick target
    double xlow = TMS_Const::TMS_Thin_Start;
    double xhi = TMS_Const::TMS_Thick_Start;
    double ylow = Lines.second->Eval(xlow);
    double yhi = Lines.second->Eval(xhi);

    double xlen = xhi-xlow;
    double ylen = yhi-ylow;
    double len = sqrt(xlen*xlen+ylen*ylen);
    xlen = xlen/len;
    ylen = ylen/len;
    DirectionZ[it] = xlen;
    DirectionX[it] = ylen;

    it++;
  }

  // Find where the first hough hit is
  TMSStart = false;

  std::vector<std::vector<TMS_Hit> > HoughCands = TMS_TrackFinder::GetFinder().GetHoughCandidates();
  int TotalHits = TMS_TrackFinder::GetFinder().GetCleanedHits().size();
  TMS_Hit *FirstTrack = NULL;

  it = 0;
  for (auto &Candidates: HoughCands) {
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
    FirstPlane[it] = Candidates.front().GetPlaneNumber();
    FirstHit[it][0] = Candidates.front().GetZ();
    FirstHit[it][1] = Candidates.front().GetNotZ();

    LastPlane[it] = Candidates.back().GetPlaneNumber();
    LastHit[it][0] = Candidates.back().GetZ();
    LastHit[it][1] = Candidates.back().GetNotZ();

    TrackLength[it] = TMS_TrackFinder::GetFinder().GetTrackLength()[it];
    TotalTrackEnergy[it] = TMS_TrackFinder::GetFinder().GetTrackEnergy()[it];
    Occupancy[it] = double(HoughCands[it].size())/TotalHits;
    // Get each hit in the track and save its energy
    for (unsigned int j = 0; j < Candidates.size(); ++j) {
      TrackHitEnergy[it][j] = Candidates[j].GetE();
      TrackHitPos[it][j][0] = Candidates[j].GetZ();
      TrackHitPos[it][j][1] = Candidates[j].GetNotZ();
    }

    it++;
  }

  // Was the first hit within the first 3 layers?
  if (FirstTrack != NULL) {
    if (FirstTrack->GetPlaneNumber() < 4) TMSStart = false;
    else TMSStart = true;
  }

  // Save down cluster information
  std::vector<std::vector<TMS_Hit> > Clusters = TMS_TrackFinder::GetFinder().GetClusterCandidates();
  nClusters = Clusters.size();
  if (nClusters > __TMS_MAX_CLUSTERS__) {
    std::cerr << "Too many clusters in TMS_TreeWriter" << std::endl;
    std::cerr << "Hard-coded maximum: " << __TMS_MAX_CLUSTERS__ << std::endl;
    std::cerr << "nClusters in event: " << nClusters << std::endl;
    return;
  }

  // Sum up the cluster candidates
  int stdit = 0;
  // Calculate the cluster by cluster summaries
  // e.g. total energy in cluster, cluster position, and cluster standard deviation
  for (auto it = Clusters.begin(); it != Clusters.end(); ++it, ++stdit) {
    double total_energy = 0;
    // Mean of cluster in z and not z
    double mean_z = 0;
    double mean_notz = 0;
    // Mean of square of cluster in z and not z
    double mean2_z = 0;
    double mean2_notz = 0;
    int nhits = (*it).size();
    for (int j = 0; j < nhits; ++j) {
      mean_z += (*it)[j].GetZ();
      mean_notz += (*it)[j].GetNotZ();
      mean2_z += (*it)[j].GetZ()*(*it)[j].GetZ();
      mean2_notz += (*it)[j].GetNotZ()*(*it)[j].GetNotZ();
      total_energy += (*it)[j].GetE();
      ClusterHitPos[stdit][j][0] = (*it)[j].GetZ();
      ClusterHitPos[stdit][j][1] = (*it)[j].GetNotZ();
      ClusterHitEnergy[stdit][j] = (*it)[j].GetE();
    }
    mean_z /= nhits;
    mean_notz /= nhits;

    mean2_z /= nhits;
    mean2_notz /= nhits;

    ClusterEnergy[stdit] = total_energy;
    nHitsInCluster[stdit] = nhits;
    ClusterPosMean[stdit][0] = mean_z;
    ClusterPosMean[stdit][1] = mean_notz;
    // Calculate the standard deviation
    double std_dev_z = mean2_z-mean_z*mean_z;
    double std_dev_notz = mean2_z-mean_z*mean_z;
    if (std_dev_z > 0) std_dev_z = sqrt(std_dev_z);
    if (std_dev_notz > 0) std_dev_notz = sqrt(std_dev_z);
    ClusterPosStdDev[stdit][0] = std_dev_z;
    ClusterPosStdDev[stdit][1] = std_dev_notz;
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
  }

  // Fill up the info only if all above has passed
  Branch_Lines->Fill();
}

// Reset the variables
void TMS_TreeWriter::Clear() {

  // Reset truth information
  EventNo = nParticles = NeutrinoPDG = Muon_TrueKE = Muon_TrueTrackLength = -999;
  Reaction = "";
  IsCC = false;
  for (int i = 0; i < 4; ++i) {
    MuonP4[i]=-999;
    Muon_Vertex[i]=-999;
    Muon_Death[i]=-999;
    NeutrinoP4[i]=-999;
  }

  // Reset line information
  TMSStart = false;
  nLines = -999;
  for (int i = 0; i < __TMS_MAX_LINES__; ++i) {
    Slope[i]=-999;
    Intercept[i]=-999;
    DirectionZ[i]=-999;
    DirectionX[i]=-999;
    Occupancy[i]=-999;
    TrackLength[i]=-999;
    TotalTrackEnergy[i]=-999;
    FirstPlane[i]=-999;
    LastPlane[i]=-999;
    nHitsInTrack[i] = -999;
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

