#include "TMS_TreeWriter.h"
#include "TMS_Reco.h"

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

  Reco_Tree = new TTree("Reco_Tree", "Reco_Tree");
  Reco_Tree->SetDirectory(Output);
  Reco_Tree->SetAutoSave(__TMS_AUTOSAVE__); // Every 1000 events (negative for MB)

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
  Branch_Lines->Branch("nLinesOne",   &nLinesOne,   "nLinesOne/I");
  Branch_Lines->Branch("nLinesOther", &nLinesOther, "nLinesOther/I");
  Branch_Lines->Branch("nLines3D",    &nLines3D,    "nLines3D/I");

  Branch_Lines->Branch("SlopeOne",     SlopeOne,      "SlopeOne[nLinesOne]/F");
  Branch_Lines->Branch("InterceptOne", InterceptOne,  "InterceptOne[nLinesOne]/F");
  Branch_Lines->Branch("Slope_DownstreamOne",      Slope_DownstreamOne,       "Slope_DownstreamOne[nLinesOne]/F");
  Branch_Lines->Branch("Intercept_DownstreamOne",  Intercept_DownstreamOne,   "Intercept_DownstreamOne[nLinesOne]/F");
  Branch_Lines->Branch("Slope_UpstreamOne",        Slope_UpstreamOne,         "Slope_UpstreamOne[nLinesOne]/F");
  Branch_Lines->Branch("Intercept_UpstreamOne",    Intercept_UpstreamOne,     "Intercept_UpstreamOne[nLinesOne]/F");

  Branch_Lines->Branch("SlopeOther",     SlopeOther,     "SlopeOther[nLinesOther]/F");
  Branch_Lines->Branch("InterceptOther", InterceptOther, "InterceptOther[nLinesOther]/F");
  Branch_Lines->Branch("Slope_DownstreamOther",     Slope_DownstreamOther,     "Slope_DownstreamOther[nLinesOther]/F");
  Branch_Lines->Branch("Intercept_DownstreamOther", Intercept_DownstreamOther, "Intercept_DownstreamOther[nLinesOther]/F");
  Branch_Lines->Branch("Slope_UpstreamOther",       Slope_UpstreamOther,       "Slope_UpstreamOther[nLinesOther]/F");
  Branch_Lines->Branch("Intercept_UpstreamOther",   Intercept_UpstreamOther,   "Intercept_UpstreamOther[nLinesOther]/F");

  Branch_Lines->Branch("DirectionZOne",   DirectionZOne,   "DirectionZOne[nLinesOne]/F");
  Branch_Lines->Branch("DirectionZOther", DirectionZOther, "DirectionZOther[nLinesOther]/F");
  Branch_Lines->Branch("DirectionXOne",   DirectionXOne,   "DirectionXOne[nLinesOne]/F");
  Branch_Lines->Branch("DirectionXOther", DirectionXOther, "DirectionXOther[nLinesOther]/F");
  Branch_Lines->Branch("DirectionZOne_Downstream",   DirectionZOne_Downstream,   "DirectionZOne_Downstream[nLinesOne]/F");
  Branch_Lines->Branch("DirectionXOne_Downstream",   DirectionXOne_Downstream,   "DirectionXOne_Downstream[nLinesOne]/F");
  Branch_Lines->Branch("DirectionZOne_Upstream",     DirectionZOne_Upstream,     "DirectionZOne_Upstream[nLinesOne]/F");
  Branch_Lines->Branch("DirectionXOne_Upstream",     DirectionXOne_Upstream,     "DirectionXOne_Upstream[nLinesOne]/F");
  Branch_Lines->Branch("DirectionZOther_Downstream", DirectionZOther_Downstream, "DirectionZOther_Downstream[nLinesOther]/F");
  Branch_Lines->Branch("DirectionXOther_Downstream", DirectionXOther_Downstream, "DirectionXOther_Downstream[nLinesOther]/F");
  Branch_Lines->Branch("DirectionZOther_Upstream",   DirectionZOther_Upstream,   "DirectionZOther_Upstream[nLinesOther]/F");
  Branch_Lines->Branch("DirectionXOther_Upstream",   DirectionXOther_Upstream,   "DirectionXOther_Upstream[nLinesOther]/F");

  Branch_Lines->Branch("FirstHoughHitOne",       FirstHitOne,       "FirstHoughHitOne[nLinesOne][2]/F");
  Branch_Lines->Branch("FirstHoughHitOther",     FirstHitOther,     "FirstHoughHitOther[nLinesOther][2]/F");
  Branch_Lines->Branch("LastHoughHitOne",        LastHitOne,        "LastHoughHitOne[nLinesOne][2]/F");
  Branch_Lines->Branch("LastHoughHitOther",      LastHitOther,      "LastHoughHitOther[nLinesOther][2]/F");
  Branch_Lines->Branch("FirstHoughHitTimeOne",   FirstHitTimeOne,   "FirstHoughHitTimeOne[nLinesOne]/F");
  Branch_Lines->Branch("FirstHoughHitTimeOther", FirstHitTimeOther, "FirstHoughHitTimeOther[nLinesOther]/F"); 
  Branch_Lines->Branch("LastHoughHitTimeOne",    LastHitTimeOne,    "LastHoughHitTimeOne[nLinesOne]/F");
  Branch_Lines->Branch("LastHoughHitTimeOther",  LastHitTimeOther,  "LastHoughHitTimeOther[nLinesOther]/F");
  Branch_Lines->Branch("HoughEarliestHitTimeOne",   EarliestHitTimeOne,   "HoughEarliestHitTimeOne[nLinesOne]/F");
  Branch_Lines->Branch("HoughEarliestHitTimeOther", EarliestHitTimeOther, "HoughEarliestHitTimeOther[nLinesOther]/F");
  Branch_Lines->Branch("HoughLatestHitTimeOne",     LatestHitTimeOne,     "HoughLatestHitTimeOne[nLinesOne]/F");
  Branch_Lines->Branch("HoughLatestHitTimeOther",   LatestHitTimeOther,   "HoughLatestHitTimeOther[nLinesOther]/F");
  Branch_Lines->Branch("FirstHoughPlaneOne",   FirstPlaneOne,   "FirstHoughPlaneOne[nLinesOne]/I");
  Branch_Lines->Branch("FirstHoughPlaneOther", FirstPlaneOther, "FirstHoughPlaneOther[nLinesOther]/I");
  Branch_Lines->Branch("LastHoughPlaneOne",    LastPlaneOne,    "LastHoughPlaneOne[nLinesOne]/I");
  Branch_Lines->Branch("LastHoughPlaneOther",  LastPlaneOther,  "LastHoughPlaneOther[nLinesOther]/I");
  Branch_Lines->Branch("TMSStart",          &TMSStart,        "TMSStart/O");
  Branch_Lines->Branch("TMSStartTime",      &TMSStartTime,    "TMSStartTime/F");
  Branch_Lines->Branch("OccupancyOne",      OccupancyOne,     "OccupancyOne[nLinesOne]/F");
  Branch_Lines->Branch("OccupancyOther",    OccupancyOther,   "OccupancyOther[nLinesOther]/F");
  Branch_Lines->Branch("TrackLengthOne",    TrackLengthOne,   "TrackLengthOne[nLinesOne]/F");
  Branch_Lines->Branch("TrackLengthOther",  TrackLengthOther, "TrackLengthOther[nLinesOther]/F");
  Branch_Lines->Branch("TotalTrackEnergyOne",   TotalTrackEnergyOne,   "TotalTrackEnergyOne[nLinesOne]/F");
  Branch_Lines->Branch("TotalTrackEnergyOther", TotalTrackEnergyOther, "TotalTrackEnergyOther[nLinesOther]/F");
  Branch_Lines->Branch("TrackStoppingOne",      TrackStoppingOne,      "TrackStoppingOne[nLinesOne]/O");
  Branch_Lines->Branch("TrackStoppingOther",    TrackStoppingOther,    "TrackStoppingOther[nLinesOther]/O");

  // 3D Track information
  Branch_Lines->Branch("FirstHoughHit3D",        FirstHit3D,        "FirstHoughHit3D[nLines3D][3]/F");
  Branch_Lines->Branch("LastHoughHit3D",         LastHit3D,         "LastHoughHit3D[nLines3D][3]/F");
  Branch_Lines->Branch("FirstHoughHitTime3D",    FirstHitTime3D,    "FirstHitTime3D[nLines3D]/F");
  Branch_Lines->Branch("LastHoughHitTime3D",     LastHitTime3D,     "LastHitTime3D[nLines3D]/F");
  Branch_Lines->Branch("HoughEarliestHitTime3D", EarliestHitTime3D, "HoughEarliestHitTime3D[nLines3D]/F");
  Branch_Lines->Branch("HoughLatestHitTime3D",   LatestHitTime3D,   "HoughLatestHitTime3D[nLines3D]/F");
  Branch_Lines->Branch("FirstHoughPlane3D",      FirstPlane3D,      "FirstHoughPlane3D[nLines3D]/I");
  Branch_Lines->Branch("LastHoughPlane3D",       LastPlane3D,       "LastHoughPlane3D[nLines3D]/I");
  Branch_Lines->Branch("Occupancy3D",            Occupancy3D,       "Occupancy3D[nLines3D]/F");
  Branch_Lines->Branch("TrackLength3D",          TrackLength3D,     "TrackLength3D[nLines3D]/F");
  Branch_Lines->Branch("TotalTrackEnergy3D",     TotalTrackEnergy3D,"TotalTrackEnergy3D[nLines3D]/F");
  Branch_Lines->Branch("TrackStopping3D",        TrackStopping3D,   "TrackStopping3D[nLines3D]/O");

  Branch_Lines->Branch("nHitsInTrack3D",      &nHitsInTrack3D,    "nHitsInTrack3D[nLines3D]/I");
  Branch_Lines->Branch("TrackHitEnergy3D",    TrackHitEnergy3D,   "TrackHitEnergy3D[10][200]/F");
  Branch_Lines->Branch("TrackHitPos3D",       TrackHitPos3D,      "TrackHitPos3D[10][200][3]/F");
  Branch_Lines->Branch("TrackHitTime3D",      TrackHitTime3D,     "TrackHitTime3D[10][200]/F");

  // Track hit energy
  Branch_Lines->Branch("nHitsInTrackOne",     &nHitsInTrackOne,    "nHitsInTrackOne[nLinesOne]/I");
  Branch_Lines->Branch("nHitsInTrackOther",   &nHitsInTrackOther,  "nHitsInTrackOther[nLinesOther]/I");
  Branch_Lines->Branch("TrackHitEnergyOne",   TrackHitEnergyOne,   "TrackHitEnergyOne[10][200]/F");
  Branch_Lines->Branch("TrackHitEnergyOther", TrackHitEnergyOther, "TrackHitEnergyOther[10][200]/F");
  Branch_Lines->Branch("TrackHitPosOne",      TrackHitPosOne,      "TrackHitPosOne[10][200][2]/F");
  Branch_Lines->Branch("TrackHitPosOther",    TrackHitPosOther,    "TrackHitPosOther[10][200][2]/F");
  Branch_Lines->Branch("TrackHitTimeOne",     TrackHitTimeOne,     "TrackHitTimeOne[10][200]/F");
  Branch_Lines->Branch("TrackHitTimeOther",   TrackHitTimeOther,   "TrackHitTimeOther[10][200]/F");

  // Cluster information
  Branch_Lines->Branch("nClustersOne",          &nClustersOne,         "nClustersOne/I");
  Branch_Lines->Branch("nClustersOther",        &nClustersOther,       "nClustersOther/I");
  Branch_Lines->Branch("ClusterEnergyOne",      ClusterEnergyOne,      "ClusterEnergyOne[nClustersOne]/F");
  Branch_Lines->Branch("ClusterEnergyOther",    ClusterEnergyOther,    "ClusterEnergyOther[nClustersOther]/F");
  Branch_Lines->Branch("ClusterTimeOne",        ClusterTimeOne,        "ClusterTimeOne[nClustersOne]/F");
  Branch_Lines->Branch("ClusterTimeOther",      ClusterTimeOther,      "ClusterTimeOther[nClustersOther]/F");
  Branch_Lines->Branch("ClusterPosMeanOne",     ClusterPosMeanOne,     "ClusterPosMeanOne[25][2]/F");
  Branch_Lines->Branch("ClusterPosMeanOther",   ClusterPosMeanOther,   "ClusterPosMeanOther[25][2]/F");
  Branch_Lines->Branch("ClusterPosStdDevOne",   ClusterPosStdDevOne,   "ClusterPosStdDevOne[25][2]/F");
  Branch_Lines->Branch("ClusterPosStdDevOther", ClusterPosStdDevOther, "ClusterPosStdDevOther[25][2]/F");
  Branch_Lines->Branch("nHitsInClusterOne",     nHitsInClusterOne,     "nHitsInClusterOne[nClustersOne]/I");
  Branch_Lines->Branch("nHitsInClusterOther",   nHitsInClusterOther,   "nHitsInClusterOther[nClustersOther]/I");
  Branch_Lines->Branch("ClusterHitPosOne",      ClusterHitPosOne,      "ClusterHitPosOne[25][200][2]/F");
  Branch_Lines->Branch("ClusterHitPosOther",    ClusterHitPosOther,    "ClusterHitPosOther[25][200][2]/F");
  Branch_Lines->Branch("ClusterHitEnergyOne",   ClusterHitEnergyOne,   "ClusterHitEnergyOne[25][200]/F");
  Branch_Lines->Branch("ClusterHitEnergyOther", ClusterHitEnergyOther, "ClusterHitEnergyOther[25][200]/F");
  Branch_Lines->Branch("ClusterHitTimeOne",     ClusterHitTimeOne,     "ClusterHitTimeOne[25][200]/F");
  Branch_Lines->Branch("ClusterHitTimeOther",   ClusterHitTimeOther,   "ClusterHitTimeOther[25][200]/F");
  Branch_Lines->Branch("ClusterHitSliceOne",    ClusterHitSliceOne,    "ClusterHitSliceOne[25][200]/I");
  Branch_Lines->Branch("ClusterHitSliceOther",  ClusterHitSliceOther,  "ClusterHitSliceOther[25][200]/I");

  // Hit information
  Branch_Lines->Branch("nHits",         &nHits,        "nHits/I");
  Branch_Lines->Branch("RecoHitPos",    RecoHitPos,    "RecoHitPos[nHits][4]/F");
  Branch_Lines->Branch("RecoHitEnergy", RecoHitEnergy, "RecoHitEnergy[nHits]/F");
  Branch_Lines->Branch("RecoHitSlice",  RecoHitSlice,  "RecoHitSlice[nHits]/I");

  // TODO: Fill these properly
  Reco_Tree->Branch("EventNo", &EventNo, "EventNo/I");
  Reco_Tree->Branch("SliceNo", &SliceNo, "SliceNo/I");
  Reco_Tree->Branch("SpillNo", &SpillNo, "SpillNo/I");

  Reco_Tree->Branch("nTracks",      &nTracks,               "nTracks/I");
  Reco_Tree->Branch("nHits",        &nHitsIn3DTrack,        "nHits[nTracks]/I");
  Reco_Tree->Branch("StartPos",     RecoTrackStartPos,      "StartPos[nTracks][3]/F");
  Reco_Tree->Branch("Direction",    RecoTrackDirection,     "Direction[nTracks][3]/F");
  Reco_Tree->Branch("EndPos",       RecoTrackEndPos,        "EndPos[nTracks][3]/F");
  Reco_Tree->Branch("Energy",       RecoTrackEnergy,        "Energy[nTracks]/F");
  Reco_Tree->Branch("EnergyDeposit",RecoTrackEnergyDeposit, "EnergyDeposit[nTracks]/F");
  Reco_Tree->Branch("Length",       RecoTrackLength,        "Length[nTracks]/F");


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
  std::vector<std::vector<TMS_Hit> > HoughCandidatesOther = TMS_TrackFinder::GetFinder().GetHoughCandidatesOther();
  nLinesOne = HoughCandidatesOne.size();
  nLinesOther = HoughCandidatesOther.size();


  // Skip the event if there aren't any Hough Lines
  if (nLinesOne > __TMS_MAX_LINES__) {
    std::cerr << "Exceeded max number of HoughLines to write to file" << std::endl;
    std::cerr << "Max lines: " << __TMS_MAX_LINES__ << std::endl;
    std::cerr << "Number of lines in event (one): " << nLinesOne << std::endl;
    std::cerr << "Not writing event" << std::endl;
    return;
  }

  if (nLinesOther > __TMS_MAX_LINES__) {
    std::cerr << "Exceeded max number of HoughLines to write to file" << std::endl;
    std::cerr << "Max lines: " << __TMS_MAX_LINES__ << std::endl;
    std::cerr << "Number of lines in event (other): " << nLinesOther << std::endl;
    std::cerr << "Not writing evet" << std::endl;
    return;
  }

  int it = 0;
  for (auto &Lines: HoughLinesOne) {
    // Get the slopes saved down
    InterceptOne[it] = Lines.second->GetParameter(0);
    Intercept_UpstreamOne[it] = TMS_TrackFinder::GetFinder().GetHoughLinesOne_Upstream()[it].first;
    Intercept_DownstreamOne[it] = TMS_TrackFinder::GetFinder().GetHoughLinesOne_Downstream()[it].first;
    SlopeOne[it] = Lines.second->GetParameter(1);
    Slope_UpstreamOne[it] = TMS_TrackFinder::GetFinder().GetHoughLinesOne_Upstream()[it].second;
    Slope_DownstreamOne[it] = TMS_TrackFinder::GetFinder().GetHoughLinesOne_Downstream()[it].second;

    // Calculate the z and x vectors by evaling the TF1 in thin and thick target
    double xlow = TMS_Const::TMS_Thin_Start;
    double xhi = TMS_Const::TMS_Thick_Start;
    double ylow = Lines.second->Eval(xlow);
    double yhi = Lines.second->Eval(xhi);
    // Do the same for up and downstream portions
    double ylow_upstream = Intercept_UpstreamOne[it]+xlow*Slope_UpstreamOne[it];
    double yhi_upstream = Intercept_UpstreamOne[it]+xhi*Slope_UpstreamOne[it];

    double ylow_downstream = Intercept_DownstreamOne[it]+xlow*Slope_DownstreamOne[it];
    double yhi_downstream = Intercept_DownstreamOne[it]+xhi*Slope_DownstreamOne[it];

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
    InterceptOther[it] = Lines.second->GetParameter(0);
    Intercept_UpstreamOther[it] = TMS_TrackFinder::GetFinder().GetHoughLinesOther_Upstream()[it].first;
    Intercept_DownstreamOther[it] = TMS_TrackFinder::GetFinder().GetHoughLinesOther_Downstream()[it].first;
    SlopeOther[it] = Lines.second->GetParameter(1);
    Slope_UpstreamOther[it] = TMS_TrackFinder::GetFinder().GetHoughLinesOther_Upstream()[it].second;
    Slope_DownstreamOther[it] = TMS_TrackFinder::GetFinder().GetHoughLinesOther_Downstream()[it].second;

    // Calculate the z and x vectors by evaling the TF1 in thin and thick target
    double xlow = TMS_Const::TMS_Thin_Start;
    double xhi = TMS_Const::TMS_Thick_Start;
    double ylow = Lines.second->Eval(xlow);
    double yhi = Lines.second->Eval(xhi);
    // Doe the same for up and downstream portions
    double ylow_upstream = Intercept_UpstreamOther[it]+xlow*Slope_UpstreamOther[it];
    double yhi_upstream = Intercept_UpstreamOther[it]+xhi*Slope_UpstreamOther[it];

    double ylow_downstream = Intercept_DownstreamOther[it]+xlow*Slope_DownstreamOther[it];
    double yhi_downstream = Intercept_DownstreamOther[it]+xhi*Slope_DownstreamOther[it];

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
    nHitsInTrackOne[it] = Candidates.size();

    // Then save the hit info
    FirstPlaneOne[it] = Candidates.front().GetPlaneNumber();
    FirstHitOne[it][0] = Candidates.front().GetZ();
    FirstHitOne[it][1] = Candidates.front().GetNotZ();
    FirstHitTimeOne[it] = Candidates.front().GetT();

    LastPlaneOne[it] = Candidates.back().GetPlaneNumber();
    LastHitOne[it][0] = Candidates.back().GetZ();
    LastHitOne[it][1] = Candidates.back().GetNotZ();
    LastHitTimeOne[it] = Candidates.back().GetT();

    TrackLengthOne[it] = TMS_TrackFinder::GetFinder().GetTrackLengthOne()[it];
    TotalTrackEnergyOne[it] = TMS_TrackFinder::GetFinder().GetTrackEnergyOne()[it];
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
    nHitsInTrackOther[it] = Candidates.size();

    // then save the hit info
    FirstPlaneOther[it] = Candidates.front().GetPlaneNumber();
    FirstHitOther[it][0] = Candidates.front().GetZ();
    FirstHitOther[it][1] = Candidates.front().GetNotZ();
    FirstHitTimeOther[it] = Candidates.front().GetT();

    LastPlaneOther[it] = Candidates.back().GetPlaneNumber();
    LastHitOther[it][0] = Candidates.back().GetZ();
    LastHitOther[it][1] = Candidates.back().GetNotZ();
    LastHitTimeOther[it] = Candidates.back().GetT();

    TrackLengthOther[it] = TMS_TrackFinder::GetFinder().GetTrackLengthOther()[it];
    TotalTrackEnergyOther[it] = TMS_TrackFinder::GetFinder().GetTrackEnergyOther()[it];
    OccupancyOther[it] = double(HoughCandsOther[it].size())/TotalHits;

    float earliest_hit_time = 1e32;
    float latest_hit_time = -1e32;
    // Get each hit in the track and save its energy
    for (unsigned int j = 0; j < Candidates.size(); ++j) {
      TrackHitEnergyOther[it][j] = Candidates[j].GetE();
      TrackHitTimeOther[it][j] = Candidates[j].GetT();
      TrackHitPosOther[it][j][0] = Candidates[j].GetZ();
      TrackHitPosOther[it][j][1] = Candidates[j].GetNotZ();

      float time = Candidates[j].GetT();

      if (time < earliest_hit_time) earliest_hit_time = time;
      if (time > latest_hit_time) latest_hit_time = time;
    }
    EarliestHitTimeOther[it] = earliest_hit_time;
    LatestHitTimeOther[it] = latest_hit_time;

    double maxenergy = 0;
    unsigned int nLastHits_temp = nLastHits;
    if ((unsigned int)nLastHits > Candidates.size()) nLastHits_temp = Candidates.size();
    for (unsigned int i = 0; i < nLastHits_temp; ++i) {
      double hitenergy = Candidates[nLastHits_temp-1-i].GetE();
      if (hitenergy > maxenergy) maxenergy = hitenergy;
    }
    if (maxenergy > EnergyCut) TrackStoppingOther[it] = true;
    

    it++;
  }

/*  std::vector<TMS_Track> HoughCands3D = TMS_TrackFinder::GetFinder().GetHoughTrack3D(); //TODO this should be done by Liam with Reco_Tree
  it = 0;
  for (auto &Candidates: HoughCands3D) {
    // Loop over hits
    for (auto &hit: Candidates) {
      if (FirstTrack == NULL) {
        FirstTrack = &hit;
      } else if (hit.GetZ() < FirstTrack->GetZ()) {
        FirstTrack = &hit;
      }
    }
    nHitsInTrack3D[it] = Candidates.size();

    // then save the hit info
    FirstPlane3D[it] = Candidates.front().GetPlaneNumber();
    FirstHit3D[it][0] = Candidates.front().GetZ();
    FirstHit3D[it][1] = Candidates.front().GetNotZ();
    FirstHitTime3D[it] = Candidates.front().GetT();

    LastPlane3D[it] = Candidates.back().GetPlaneNumber();
    LastHit3D[it][0] = Candidates.back().GetZ();
    LastHit3D[it][1] = Candidates.back().GetNotZ();
    LastHitTime3D[it] = Candidates.back().GetT();

    TrackLength3D[it] = TMS_TrackFinder::GetFinder().GetTrackLength3D()[it];
    TotalTrackEnergy3D[it] = TMS_TrackFinder::GetFinder().GetTrackEnergy3D()[it];
    Occupancy3D[it] = double(HoughCands3D[it].size())/TotalHits;

    float earliest_hit_time = 1e32;
    float latest_hit_time = -1e32;
    // Get each hit in the track and save its energy
    for (unsigned int j = 0; j < Candidates.size(); ++j) {
      TrackHitEnergy3D[it][j] = Candidates[j].GetE();
      TrackHitTime3D[it][j] = Candidates[j].GetT();
      TrackHitPos3D[it][j][0] = Candidates[j].GetZ();
      TrackHitPos3D[it][j][1] = Candidates[j].GetNotZ();

      float time = Candidates[j].GetT();

      if (time < earliest_hit_time) earliest_hit_time = time;
      if (time > latest_hit_time) latest_hit_time = time;
    }
    EarliestHitTime3D[it] = earliest_hit_time;
    LatestHitTime3D[it] = latest_hit_time;

    double maxenergy = 0;
    unsigned int nLastHits_temp = nLastHits;
    if ((unsigned int)nLastHits > Candidates.size()) nLastHits_temp = Candidates.size();
    for (unsigned int i = 0; i < nLastHits_temp; ++i) {
      double hitenergy = Candidates[nLastHits_temp-1-i].GetE();
      if (hitenergy > maxenergy) maxenergy = hitenergy;
    }
    if (maxenergy > EnergyCut) TrackStopping3D[it] = true;

    it++;
  }*/
  

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
    std::cout << "Cluster One nHits: " << nhits << std::endl;
    for (int j = 0; j < nhits; ++j) {
      mean_z += (*it)[j].GetZ();
      mean_notz += (*it)[j].GetNotZ();
      mean2_z += (*it)[j].GetZ()*(*it)[j].GetZ();
      mean2_notz += (*it)[j].GetNotZ()*(*it)[j].GetNotZ();
      total_energy += (*it)[j].GetE();
      float time = (*it)[j].GetT();
      if (time < min_cluster_time) min_cluster_time = time;
      std::cout << (*it)[j].GetZ() << " " << (*it)[j].GetNotZ() << std::endl;
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
    if (std_dev_notz > 0) std_dev_notz = sqrt(std_dev_notz);
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
    std::cout << "Cluster Other nHits: " << nhits << std::endl;
    for (int j = 0; j < nhits; ++j) {
      mean_z += (*it)[j].GetZ();
      mean_notz += (*it)[j].GetNotZ();
      mean2_z += (*it)[j].GetZ()*(*it)[j].GetZ();
      mean2_notz += (*it)[j].GetNotZ()*(*it)[j].GetNotZ();
      total_energy += (*it)[j].GetE();
      float time = (*it)[j].GetT();
      if (time < min_cluster_time) min_cluster_time = time;
      std::cout << (*it)[j].GetZ() << " " << (*it)[j].GetNotZ() << std::endl;
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


  // Fill the 3D Tracks
  // First get the tracks for this event:
  //TODO: Function here that uses the info from ^^^^^ to fill the 3DTrack objects

  int itTrack= 0;
  std::vector<TMS_Track> Reco_Tracks = TMS_TrackFinder::GetFinder().GetHoughTracks3D();
  nTracks = Reco_Tracks.size();

  for (auto RecoTrack = Reco_Tracks.begin(); RecoTrack != Reco_Tracks.end(); ++RecoTrack, ++itTrack) {
    nHitsIn3DTrack[itTrack]         = (int) RecoTrack->Hits.size(); // Do we need to cast it? idk
    RecoTrackEnergy[itTrack]        =       RecoTrack->GetEnergyRange();
    RecoTrackLength[itTrack]        =       RecoTrack->Length;
    RecoTrackEnergyDeposit[itTrack] =       RecoTrack->GetEnergyDeposit();
    for (int j = 0; j < 3; j++) {
      RecoTrackStartPos[itTrack][j]  = RecoTrack->Start[j];
      RecoTrackEndPos[itTrack][j]    = RecoTrack->End[j];
      RecoTrackDirection[itTrack][j] = RecoTrack->Direction[j];
    }

    // Can manually compute direction if it hasn't been set
    if ( (RecoTrackDirection[itTrack][0] == 0) && (RecoTrackDirection[itTrack][1] == 0) && (RecoTrackDirection[itTrack][2] == 0) )
    { // If true^ it seems the direction hasn't been set
      for (int j = 0; j < 3; j++)
      { // Right now no need to make sure this is a unit vector
        RecoTrackDirection[itTrack][j] = RecoTrack->End[j] - RecoTrack->Start[j];
      }
    }
  }

  Reco_Tree->Fill();
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
  nLinesOne = -999;
  nLinesOther = -999;
  for (int i = 0; i < __TMS_MAX_LINES__; ++i) {
    SlopeOne[i]=-999;
    SlopeOther[i]=-999;
    InterceptOne[i]=-999;
    InterceptOther[i]=-999;
    Slope_DownstreamOne[i]=-999;
    Slope_DownstreamOther[i]=-999;
    Intercept_DownstreamOne[i]=-999;
    Intercept_DownstreamOther[i]=-999;
    Slope_UpstreamOne[i]=-999;
    Slope_UpstreamOther[i]=-999;
    Intercept_UpstreamOne[i]=-999;
    Intercept_UpstreamOther[i]=-999;

    DirectionZOne[i]=-999;
    DirectionXOne[i]=-999;
    DirectionZOne_Upstream[i]=-999;
    DirectionXOne_Upstream[i]=-999;
    DirectionZOne_Downstream[i]=-999;
    DirectionXOne_Downstream[i]=-999;

    DirectionZOther[i]=-999;
    DirectionXOther[i]=-999;
    DirectionZOther_Upstream[i]=-999;
    DirectionXOther_Upstream[i]=-999;
    DirectionZOther_Downstream[i]=-999;
    DirectionXOther_Downstream[i]=-999;

    OccupancyOne[i]=-999;
    OccupancyOther[i]=-999;
    TrackLengthOne[i]=-999;
    TrackLengthOther[i]=-999;
    TotalTrackEnergyOne[i]=-999;
    TotalTrackEnergyOther[i]=-999;
    FirstPlaneOne[i]=-999;
    FirstPlaneOther[i]=-999;
    LastPlaneOne[i]=-999;
    LastPlaneOther[i]=-999;
    nHitsInTrackOne[i] = -999;
    nHitsInTrackOther[i] = -999;
    TrackStoppingOne[i] = false;
    TrackStoppingOther[i] = false;

    Occupancy3D[i] = -999;
    TrackLength3D[i] = -999;
    TotalTrackEnergy3D[i] = -999;
    FirstPlane3D[i] = -999;
    LastPlane3D[i] = -999;
    nHitsInTrack3D[i] = -999;
    TrackStopping3D[i] = false;

    for (int j = 0; j < 2; ++j) {
      FirstHitOne[i][j] = -999;
      FirstHitOther[i][j] = -999;
      LastHitOne[i][j] = -999;
      LastHitOther[i][j] = -999;
    }
    for (int j = 0; j < 3; ++j) {
      FirstHit3D[i][j] = -999;
      LastHit3D[i][j] = -999;
    }
    for (int j = 0; j < __TMS_MAX_LINE_HITS__; ++j) {
      TrackHitEnergyOne[i][j]=-999;
      TrackHitEnergyOther[i][j]=-999;
      TrackHitPosOne[i][j][0]=-999;
      TrackHitPosOne[i][j][1]=-999;
      TrackHitPosOther[i][j][0]=-999;
      TrackHitPosOther[i][j][1]=-999;

      TrackHitEnergy3D[i][j] = -999;
      TrackHitPos3D[i][j][0] = -999;
      TrackHitPos3D[i][j][1] = -999;
      TrackHitPos3D[i][j][2] = -999;
    }
  }

  // Reset hit information
  nHits = -999;
  for (int i = 0; i < __TMS_MAX_HITS__; ++i) {
    for (int j = 0; j < 4; ++j) RecoHitPos[i][j] = -999;
    RecoHitEnergy[i] = -999;
  }

  // Reset Cluster info
  nClustersOne = -999;
  nClustersOther = -999;
  for (int i = 0; i < __TMS_MAX_CLUSTERS__; ++i) {
    ClusterEnergyOne[i] = -999;
    ClusterEnergyOther[i] = -999;
    nHitsInClusterOne[i] = -999;
    nHitsInClusterOther[i] = -999;
    for (int j = 0; j < 2; ++j) {
      ClusterPosMeanOne[i][j] = -999;
      ClusterPosStdDevOne[i][j] = -999;
      
      ClusterPosMeanOther[i][j] = -999;
      ClusterPosStdDevOther[i][j] = -999;
    }
    for (int j = 0; j < __TMS_MAX_LINE_HITS__; ++j) {
      ClusterHitPosOne[i][j][0] = -999;
      ClusterHitPosOne[i][j][1] = -999;
      ClusterHitEnergyOne[i][j] = -999;

      ClusterHitPosOther[i][j][0] = -999;
      ClusterHitPosOther[i][j][1] = -999;
      ClusterHitEnergyOther[i][j] = -999;
    }
  }

}

