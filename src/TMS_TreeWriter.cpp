#include "TMS_TreeWriter.h"
#include "TMS_Reco.h"
#include "TMS_Utils.h"

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
  if(std::getenv("ND_PRODUCTION_TMSRECO_OUTFILE")) {
    Outputname = std::getenv("ND_PRODUCTION_TMSRECO_OUTFILE");
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
  Branch_Lines->Branch("nLinesU", &nLinesU, "nLinesU/I");
  Branch_Lines->Branch("nLinesV", &nLinesV, "nLinesV/I");
  Branch_Lines->Branch("nLinesX", &nLinesX, "nLinesX/I");
  Branch_Lines->Branch("nLines3D",    &nLines3D,    "nLines3D/I");

  Branch_Lines->Branch("SlopeU",     SlopeU,      "SlopeU[nLinesU]/F");
  Branch_Lines->Branch("InterceptU", InterceptU,  "InterceptU[nLinesU]/F");
  Branch_Lines->Branch("Slope_DownstreamU",      Slope_DownstreamU,       "Slope_DownstreamU[nLinesU]/F");
  Branch_Lines->Branch("Intercept_DownstreamU",  Intercept_DownstreamU,   "Intercept_DownstreamU[nLinesU]/F");
  Branch_Lines->Branch("Slope_UpstreamU",        Slope_UpstreamU,         "Slope_UpstreamU[nLinesU]/F");
  Branch_Lines->Branch("Intercept_UpstreamU",    Intercept_UpstreamU,     "Intercept_UpstreamU[nLinesU]/F");

  Branch_Lines->Branch("SlopeV",     SlopeV,     "SlopeV[nLinesV]/F");
  Branch_Lines->Branch("InterceptV", InterceptV, "InterceptV[nLinesV]/F");
  Branch_Lines->Branch("Slope_DownstreamV",     Slope_DownstreamV,      "Slope_DownstreamV[nLinesV]/F");
  Branch_Lines->Branch("Intercept_DownstreamV", Intercept_DownstreamV,  "Intercept_DownstreamV[nLinesV]/F");
  Branch_Lines->Branch("Slope_UpstreamV",       Slope_UpstreamV,        "Slope_UpstreamV[nLinesV]/F");
  Branch_Lines->Branch("Intercept_UpstreamV",   Intercept_UpstreamV,    "Intercept_UpstreamV[nLinesV]/F");

  Branch_Lines->Branch("SlopeX",                SlopeX,                 "SlopeX[nLinesX]/F");
  Branch_Lines->Branch("InterceptX",            InterceptX,             "InterceptX[nLinesX]/F");
  Branch_Lines->Branch("Slope_DownstreamX",     Slope_DownstreamX,      "Slope_DownstreamX[nLinesX]/F");
  Branch_Lines->Branch("Intercept_DownstreamX", Intercept_DownstreamX,  "Intercept_DownstreamX[nLinesX]/F");
  Branch_Lines->Branch("Slope_UpstreamX",       Slope_UpstreamX,        "Slope_UpstreamX[nLinesX]/F");
  Branch_Lines->Branch("Intercept_UpstreamX",   Intercept_UpstreamX,    "Intercept_UpstreamX[nLinesX]/F");

  Branch_Lines->Branch("DirectionZU", DirectionZU, "DirectionZU[nLinesU]/F");
  Branch_Lines->Branch("DirectionZV", DirectionZV, "DirectionZV[nLinesV]/F");
  Branch_Lines->Branch("DirectionZX", DirectionZX, "DirectionZX[nLinesX]/F");
  Branch_Lines->Branch("DirectionXU", DirectionXU, "DirectionXU[nLinesU]/F");
  Branch_Lines->Branch("DirectionXV", DirectionXV, "DirectionXV[nLinesV]/F");
  Branch_Lines->Branch("DirectionYX", DirectionYX, "DirectionYX[nLinesX]/F");
  Branch_Lines->Branch("DirectionZU_Downstream",   DirectionZU_Downstream,   "DirectionZU_Downstream[nLinesU]/F");
  Branch_Lines->Branch("DirectionXU_Downstream",   DirectionXU_Downstream,   "DirectionXU_Downstream[nLinesU]/F");
  Branch_Lines->Branch("DirectionZU_Upstream",     DirectionZU_Upstream,     "DirectionZU_Upstream[nLinesU]/F");
  Branch_Lines->Branch("DirectionXU_Upstream",     DirectionXU_Upstream,     "DirectionXU_Upstream[nLinesU]/F");
  Branch_Lines->Branch("DirectionZV_Downstream", DirectionZV_Downstream, "DirectionZV_Downstream[nLinesV]/F");
  Branch_Lines->Branch("DirectionXV_Downstream", DirectionXV_Downstream, "DirectionXV_Downstream[nLinesV]/F");
  Branch_Lines->Branch("DirectionZV_Upstream",   DirectionZV_Upstream,   "DirectionZV_Upstream[nLinesV]/F");
  Branch_Lines->Branch("DirectionXV_Upstream",   DirectionXV_Upstream,   "DirectionXV_Upstream[nLinesV]/F");
  Branch_Lines->Branch("DirectionZX_Downstream",  DirectionZX_Downstream, "DirectionZX_Downstream[nLinesX]/F");
  Branch_Lines->Branch("DirectionYX_Downstream",  DirectionYX_Downstream, "DirectionYX_Downstream[nLinesX]/F");
  Branch_Lines->Branch("DirectionZX_Upstream",    DirectionZX_Upstream,   "DirectionZX_Upstream[nLinesX]/F");
  Branch_Lines->Branch("DirectionYX_Upstream",    DirectionYX_Upstream,   "DirectionYX_Upstream[nLinesX]/F");

  Branch_Lines->Branch("FirstHoughHitU",     FirstHitU,     "FirstHoughHitU[nLinesU][2]/F");
  Branch_Lines->Branch("FirstHoughHitV",     FirstHitV,     "FirstHoughHitV[nLinesV][2]/F");
  Branch_Lines->Branch("FirstHoughHitX",     FirstHitX,     "FirstHoughHitX[nLinesX][2]/F");
  Branch_Lines->Branch("LastHoughHitU",      LastHitU,      "LastHoughHitU[nLinesU][2]/F");
  Branch_Lines->Branch("LastHoughHitV",      LastHitV,      "LastHoughHitV[nLinesV][2]/F");
  Branch_Lines->Branch("LastHoughHitX",      LastHitX,      "LastHoughHitX[nLinesX][2]/F");
  Branch_Lines->Branch("FirstHoughHitTimeU", FirstHitTimeU, "FirstHoughHitTimeU[nLinesU]/F");
  Branch_Lines->Branch("FirstHoughHitTimeV", FirstHitTimeV, "FirstHoughHitTimeV[nLinesV]/F"); 
  Branch_Lines->Branch("FirstHoughHitTimeX", FirstHitTimeX, "FirstHoughHitTimeX[nLinesX]/F");
  Branch_Lines->Branch("LastHoughHitTimeU",  LastHitTimeU,  "LastHoughHitTimeU[nLinesU]/F");
  Branch_Lines->Branch("LastHoughHitTimeV",  LastHitTimeV,  "LastHoughHitTimeV[nLinesV]/F");
  Branch_Lines->Branch("LastHoughHitTImeX",  LastHitTimeX,  "LastHoughHitTimeX[nLinesX]/F");
  Branch_Lines->Branch("HoughEarliestHitTimeU", EarliestHitTimeU, "HoughEarliestHitTimeU[nLinesU]/F");
  Branch_Lines->Branch("HoughEarliestHitTimeV", EarliestHitTimeV, "HoughEarliestHitTimeV[nLinesV]/F");
  Branch_Lines->Branch("HoughEarliestHitTimeX", EarliestHitTimeX, "HoughEarliestHitTimeX[nLinesX]/F");
  Branch_Lines->Branch("HoughLatestHitTimeU",   LatestHitTimeU,   "HoughLatestHitTimeU[nLinesU]/F");
  Branch_Lines->Branch("HoughLatestHitTimeV",   LatestHitTimeV,   "HoughLatestHitTimeV[nLinesV]/F");
  Branch_Lines->Branch("HoughLatestHitTimeX",   LatestHitTimeX,   "HoughLatestHitTimeX[nLinesX]/F");
  Branch_Lines->Branch("FirstHoughPlaneU", FirstPlaneU,   "FirstHoughPlaneU[nLinesU]/I");
  Branch_Lines->Branch("FirstHoughPlaneV", FirstPlaneV,   "FirstHoughPlaneV[nLinesV]/I");
  Branch_Lines->Branch("FirstHoughPlaneX", FirstPlaneX,   "FirstHoughPlaneX[nLinesX]/I");
  Branch_Lines->Branch("LastHoughPlaneU",  LastPlaneU,    "LastHoughPlaneU[nLinesU]/I");
  Branch_Lines->Branch("LastHoughPlaneV",  LastPlaneV,    "LastHoughPlaneV[nLinesV]/I");
  Branch_Lines->Branch("LastHoughPlaneX",  LastPlaneX,    "LastHoughPlaneX[nLinesX]/I");
  Branch_Lines->Branch("TMSStart",         &TMSStart,     "TMSStart/O");
  Branch_Lines->Branch("TMSStartTime",     &TMSStartTime, "TMSStartTime/F");
  Branch_Lines->Branch("OccupancyU",       OccupancyU,    "OccupancyU[nLinesU]/F");
  Branch_Lines->Branch("OccupancyV",       OccupancyV,    "OccupancyV[nLinesV]/F");
  Branch_Lines->Branch("OccupancyX",       OccupancyX,    "OccupancyX[nLinesX]/F");
  Branch_Lines->Branch("TrackLengthU",     TrackLengthU,  "TrackLengthU[nLinesU]/F");
  Branch_Lines->Branch("TrackLengthV",     TrackLengthV,  "TrackLengthV[nLinesV]/F");
  Branch_Lines->Branch("TrackLengthX",     TrackLengthX,  "TrackLengthX[nLinesX]/F");
  Branch_Lines->Branch("TotalTrackEnergyU", TotalTrackEnergyU, "TotalTrackEnergyU[nLinesU]/F");
  Branch_Lines->Branch("TotalTrackEnergyV", TotalTrackEnergyV, "TotalTrackEnergyV[nLinesV]/F");
  Branch_Lines->Branch("TotalTrackEnergyX", TotalTrackEnergyX, "TotalTrackEnergyX[nLinesX]/F");
  Branch_Lines->Branch("TrackStoppingU",    TrackStoppingU,    "TrackStoppingU[nLinesU]/O");
  Branch_Lines->Branch("TrackStoppingV",    TrackStoppingV,    "TrackStoppingV[nLinesV]/O");
  Branch_Lines->Branch("TrackStoppingX",    TrackStoppingX,    "TrackStoppingX[nLinesX]/O");

/*  // 3D Track information
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
  Branch_Lines->Branch("TrackHitTime3D",      TrackHitTime3D,     "TrackHitTime3D[10][200]/F");*/

  // Track hit energy
  Branch_Lines->Branch("nHitsInTrackU",   &nHitsInTrackU,  "nHitsInTrackU[nLinesU]/I");
  Branch_Lines->Branch("nHitsInTrackV",   &nHitsInTrackV,  "nHitsInTrackV[nLinesV]/I");
  Branch_Lines->Branch("nHitsInTrackX",   &nHitsInTrackX,  "nHitsInTrackX[nLinesX]/I");
  Branch_Lines->Branch("TrackHitEnergyU", TrackHitEnergyU, "TrackHitEnergyU[10][200]/F");
  Branch_Lines->Branch("TrackHitEnergyV", TrackHitEnergyV, "TrackHitEnergyV[10][200]/F");
  Branch_Lines->Branch("TrackHitEnergyX", TrackHitEnergyX, "TrackHitEnergyX[10][200]/F");
  Branch_Lines->Branch("TrackHitPosU",    TrackHitPosU,    "TrackHitPosU[10][200][2]/F");
  Branch_Lines->Branch("TrackHitPosV",    TrackHitPosV,    "TrackHitPosV[10][200][2]/F");
  Branch_Lines->Branch("TrackHitPosX",    TrackHitPosX,    "TrackHitPosX[10][200][2]/F");
  Branch_Lines->Branch("TrackHitTimeU",   TrackHitTimeU,   "TrackHitTimeU[10][200]/F");
  Branch_Lines->Branch("TrackHitTimeV",   TrackHitTimeV,   "TrackHitTimeV[10][200]/F");
  Branch_Lines->Branch("TrackHitTimeX",   TrackHitTimeX,   "TrackHitTimeX[10][200]/F");

  // Cluster information
  Branch_Lines->Branch("nClustersU",        &nClustersU,       "nClustersU/I");
  Branch_Lines->Branch("nClustersV",        &nClustersV,       "nClustersV/I");
  Branch_Lines->Branch("nClusterX",         &nClustersX,       "nClustersX/I");
  Branch_Lines->Branch("ClusterEnergyU",    ClusterEnergyU,    "ClusterEnergyU[nClustersU]/F");
  Branch_Lines->Branch("ClusterEnergyV",    ClusterEnergyV,    "ClusterEnergyV[nClustersV]/F");
  Branch_Lines->Branch("ClusterEnergyX",    ClusterEnergyX,    "ClusterEnergyX[nClustersX]/F");
  Branch_Lines->Branch("ClusterTimeU",      ClusterTimeU,      "ClusterTimeU[nClustersU]/F");
  Branch_Lines->Branch("ClusterTimeV",      ClusterTimeV,      "ClusterTimeV[nClustersV]/F");
  Branch_Lines->Branch("ClusterTimeX",      ClusterTimeX,      "ClusterTimeX[nClustersX]/F");
  Branch_Lines->Branch("ClusterPosMeanU",   ClusterPosMeanU,   "ClusterPosMeanU[25][2]/F");
  Branch_Lines->Branch("ClusterPosMeanV",   ClusterPosMeanV,   "ClusterPosMeanV[25][2]/F");
  Branch_Lines->Branch("ClusterPosMeanX",   ClusterPosMeanX,   "ClusterPosMeanX[25][2]/F");
  Branch_Lines->Branch("ClusterPosStdDevU", ClusterPosStdDevU, "ClusterPosStdDevU[25][2]/F");
  Branch_Lines->Branch("ClusterPosStdDevV", ClusterPosStdDevV, "ClusterPosStdDevV[25][2]/F");
  Branch_Lines->Branch("ClusterPosStdDevX", ClusterPosStdDevX, "ClusterPosStdDevX[25][2]/F");
  Branch_Lines->Branch("nHitsInClusterU",   nHitsInClusterU,   "nHitsInClusterU[nClustersU]/I");
  Branch_Lines->Branch("nHitsInClusterV",   nHitsInClusterV,   "nHitsInClusterV[nClustersV]/I");
  Branch_Lines->Branch("nHitsInClusterX",   nHitsInClusterX,   "nHitsInClusterX[nClustersX]/I");
  Branch_Lines->Branch("ClusterHitPosU",    ClusterHitPosU,    "ClusterHitPosU[25][200][2]/F");
  Branch_Lines->Branch("ClusterHitPosV",    ClusterHitPosV,    "ClusterHitPosV[25][200][2]/F");
  Branch_Lines->Branch("ClusterHitPosX",    ClusterHitPosX,    "ClusterHitPosX[25][200][2]/F");
  Branch_Lines->Branch("ClusterHitEnergyU", ClusterHitEnergyU, "ClusterHitEnergyU[25][200]/F");
  Branch_Lines->Branch("ClusterHitEnergyV", ClusterHitEnergyV, "ClusterHitEnergyV[25][200]/F");
  Branch_Lines->Branch("ClusterHitEnergyX", ClusterHitEnergyX, "ClusterHitEnergyX[25][200]/F");
  Branch_Lines->Branch("ClusterHitTimeU",   ClusterHitTimeU,   "ClusterHitTimeU[25][200]/F");
  Branch_Lines->Branch("ClusterHitTimeV",   ClusterHitTimeV,   "ClusterHitTimeV[25][200]/F");
  Branch_Lines->Branch("ClusterHitTimeX",   ClusterHitTimeX,   "ClusterHitTimeX[25][200]/F");
  Branch_Lines->Branch("ClusterHitSliceU",  ClusterHitSliceU,  "ClusterHitSliceU[25][200]/I");
  Branch_Lines->Branch("ClusterHitSliceV",  ClusterHitSliceV,  "ClusterHitSliceV[25][200]/I");
  Branch_Lines->Branch("ClusterHitSliceX",  ClusterHitSliceX,  "ClusterHitSliceX[25][200]/I");

  // Hit information
  Branch_Lines->Branch("nHits",         &nHits,        "nHits/I");
  Branch_Lines->Branch("RecoHitPos",    RecoHitPos,    "RecoHitPos[nHits][4]/F");
  Branch_Lines->Branch("RecoHitEnergy", RecoHitEnergy, "RecoHitEnergy[nHits]/F");
  Branch_Lines->Branch("RecoHitSlice",  RecoHitSlice,  "RecoHitSlice[nHits]/I");

  // Track information
  // TODO: Fill these properly
  Reco_Tree->Branch("EventNo", &EventNo, "EventNo/I");
  Reco_Tree->Branch("SliceNo", &SliceNo, "SliceNo/I");
  Reco_Tree->Branch("SpillNo", &SpillNo, "SpillNo/I");

  Reco_Tree->Branch("nTracks",      &nTracks,               "nTracks/I");
  Reco_Tree->Branch("nHits",        nHitsIn3DTrack,         "nHits[nTracks]/I");
  Reco_Tree->Branch("TrackHitPos",  RecoTrackHitPos,        "TrackHitPos[nTracks][200][3]/F");
  Reco_Tree->Branch("StartPos",     RecoTrackStartPos,      "StartPos[nTracks][3]/F");
  Reco_Tree->Branch("Direction",    RecoTrackDirection,     "Direction[nTracks][3]/F");
  Reco_Tree->Branch("EndPos",       RecoTrackEndPos,        "EndPos[nTracks][3]/F");
  Reco_Tree->Branch("EnergyRange",  RecoTrackEnergyRange,   "EnergyRange[nTracks]/F");
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
  Truth_Info->Branch("VisibleEnergyFromUVertexInSlice", &VisibleEnergyFromUVertexInSlice, "VisibleEnergyFromUVertexInSlice/F");
  Truth_Info->Branch("TotalVisibleEnergyFromVertex", &TotalVisibleEnergyFromVertex, "TotalVisibleEnergyFromVertex/F");
  Truth_Info->Branch("VisibleEnergyFromVVerticesInSlice", &VisibleEnergyFromVVerticesInSlice, "VisibleEnergyFromVVerticesInSlice/F");
  Truth_Info->Branch("VertexVisibleEnergyFractionInSlice", &VertexVisibleEnergyFractionInSlice, "VertexVisibleEnergyFractionInSlice/F");
  Truth_Info->Branch("PrimaryVertexVisibleEnergyFraction", &PrimaryVertexVisibleEnergyFraction, "PrimaryVertexVisibleEnergyFraction/F");
  
  Truth_Info->Branch("RecoTrackN", &RecoTrackN, "RecoTrackN/I");
  Truth_Info->Branch("RecoTrackTrueVisibleEnergy", RecoTrackTrueVisibleEnergy,
                     "RecoTrackTrueVisibleEnergy[RecoTrackN]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleIndex", RecoTrackPrimaryParticleIndex, "RecoTrackPrimaryParticleIndex[RecoTrackN]/I");
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueVisibleEnergy", RecoTrackPrimaryParticleTrueVisibleEnergy,
                     "RecoTrackPrimaryParticleTrueVisibleEnergy[RecoTrackN]/F");
  Truth_Info->Branch("RecoTrackSecondaryParticleIndex", RecoTrackSecondaryParticleIndex, "RecoTrackSecondaryParticleIndex[RecoTrackN]/I");
  Truth_Info->Branch("RecoTrackSecondaryParticleTrueVisibleEnergy", RecoTrackSecondaryParticleTrueVisibleEnergy,
                     "RecoTrackSecondaryParticleTrueVisibleEnergy[RecoTrackN]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueMomentumTrackStart", RecoTrackPrimaryParticleTrueMomentumTrackStart,
                     "RecoTrackPrimaryParticleTrueMomentumTrackStart[RecoTrackN][4]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleTruePositionTrackStart", RecoTrackPrimaryParticleTruePositionTrackStart,
                     "RecoTrackPrimaryParticleTruePositionTrackStart[RecoTrackN][4]/F"); 
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueMomentumTrackEnd", RecoTrackPrimaryParticleTrueMomentumTrackEnd,
                     "RecoTrackPrimaryParticleTrueMomentumTrackEnd[RecoTrackN][4]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleTruePositionTrackEnd", RecoTrackPrimaryParticleTruePositionTrackEnd,
                     "RecoTrackPrimaryParticleTruePositionTrackEnd[RecoTrackN][4]/F"); 
                 
  Truth_Info->Branch("nTrueParticles", &nTrueParticles, "nTrueParticles/I");
  Truth_Info->Branch("VertexID", VertexID, "VertexID[nTrueParticles]/I");
  Truth_Info->Branch("Parent", Parent, "Parent[nTrueParticles]/I");
  Truth_Info->Branch("TrackId", TrackId, "TrackId[nTrueParticles]/I");
  Truth_Info->Branch("PDG", PDG, "PDG[nTrueParticles]/I");
  Truth_Info->Branch("TrueVisibleEnergy", TrueVisibleEnergy, "TrueVisibleEnergy[nTrueParticles]/F");
  
  Truth_Info->Branch("BirthMomentum", BirthMomentum, "BirthMomentum[nTrueParticles][4]/F");
  Truth_Info->Branch("BirthPosition", BirthPosition, "BirthPosition[nTrueParticles][4]/F");
  Truth_Info->Branch("DeathMomentum", DeathMomentum, "DeathMomentum[nTrueParticles][4]/F");
  Truth_Info->Branch("DeathPosition", DeathPosition, "DeathPosition[nTrueParticles][4]/F");
  
  // IsInside-based start/end
  Truth_Info->Branch("MomentumLArStart", MomentumLArStart, "MomentumLArStart[nTrueParticles][4]/F");
  Truth_Info->Branch("PositionLArStart", PositionLArStart, "PositionLArStart[nTrueParticles][4]/F");
  Truth_Info->Branch("MomentumLArEnd", MomentumLArEnd, "MomentumLArEnd[nTrueParticles][4]/F");
  Truth_Info->Branch("PositionLArEnd", PositionLArEnd, "PositionLArEnd[nTrueParticles][4]/F");
  Truth_Info->Branch("MomentumTMSStart", MomentumTMSStart, "MomentumTMSStart[nTrueParticles][4]/F");
  Truth_Info->Branch("PositionTMSStart", PositionTMSStart, "PositionTMSStart[nTrueParticles][4]/F");
  Truth_Info->Branch("MomentumTMSFirstTwoModulesEnd", MomentumTMSFirstTwoModulesEnd, "MomentumTMSFirstTwoModulesEnd[nTrueParticles][4]/F");
  Truth_Info->Branch("PositionTMSFirstTwoModulesEnd", PositionTMSFirstTwoModulesEnd, "PositionTMSFirstTwoModulesEnd[nTrueParticles][4]/F"); 
  Truth_Info->Branch("MomentumTMSThinEnd", MomentumTMSThinEnd, "MomentumTMSThinEnd[nTrueParticles][4]/F");
  Truth_Info->Branch("PositionTMSThinEnd", PositionTMSThinEnd, "PositionTMSThinEnd[nTrueParticles][4]/F"); 
  Truth_Info->Branch("MomentumTMSEnd", MomentumTMSEnd, "MomentumTMSEnd[nTrueParticles][4]/F");
  Truth_Info->Branch("PositionTMSEnd", PositionTMSEnd, "PositionTMSEnd[nTrueParticles][4]/F"); 
  
  // Z-based start/end
  Truth_Info->Branch("MomentumZIsLArEnd", MomentumZIsLArEnd, "MomentumZIsLArEnd[nTrueParticles][4]/F");
  Truth_Info->Branch("PositionZIsLArEnd", PositionZIsLArEnd, "PositionZIsLArEnd[nTrueParticles][4]/F");
  Truth_Info->Branch("MomentumZIsTMSStart", MomentumZIsTMSStart, "MomentumZIsTMSStart[nTrueParticles][4]/F");
  Truth_Info->Branch("PositionZIsTMSStart", PositionZIsTMSStart, "PositionZIsTMSStart[nTrueParticles][4]/F");
  Truth_Info->Branch("MomentumZIsTMSEnd", MomentumZIsTMSEnd, "MomentumZIsTMSEnd[nTrueParticles][4]/F");
  Truth_Info->Branch("PositionZIsTMSEnd", PositionZIsTMSEnd, "PositionZIsTMSEnd[nTrueParticles][4]/F"); 
}

static void setMomentum(float *branch, TVector3 momentum, double energy = -9999) {
    branch[0] = momentum.X();
    branch[1] = momentum.Y();
    branch[2] = momentum.Z();
    branch[3] = energy;
}

static void setPosition(float *branch, TLorentzVector position) {
    branch[0] = position.X();
    branch[1] = position.Y();
    branch[2] = position.Z();
    branch[3] = position.T();
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
  VisibleEnergyFromUVertexInSlice = event.GetVisibleEnergyFromUVertexInSlice();
  TotalVisibleEnergyFromVertex = event.GetTotalVisibleEnergyFromVertex();
  VisibleEnergyFromVVerticesInSlice = event.GetVisibleEnergyFromVVerticesInSlice();
  VertexVisibleEnergyFractionInSlice = VisibleEnergyFromUVertexInSlice / TotalVisibleEnergyFromVertex;
  PrimaryVertexVisibleEnergyFraction = VisibleEnergyFromUVertexInSlice / (VisibleEnergyFromVVerticesInSlice + VisibleEnergyFromUVertexInSlice);

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
  
  
    
  nTrueParticles = TrueParticles.size();
  for (auto it = TrueParticles.begin(); it != TrueParticles.end(); ++it) {
    int index = it - TrueParticles.begin();
    
    if (index >= __TMS_MAX_TRUE_PARTICLES__) {
      std::cerr<<"WARNING: Found more particles than __TMS_MAX_TRUE_PARTICLES__. Stopping loop early. If this happens often, increase the max"<<std::endl;
      break;
    }
  
    VertexID[index] = (*it).GetVertexID();
    Parent[index] = (*it).GetParent();
    TrackId[index] = (*it).GetTrackId();
    PDG[index] = (*it).GetPDG();
    TrueVisibleEnergy[index] = (*it).GetTrueVisibleEnergy();
    
    setMomentum(BirthMomentum[index], (*it).GetBirthMomentum(), (*it).GetBirthEnergy());
    setPosition(BirthPosition[index], (*it).GetBirthPosition());
    
    setMomentum(DeathMomentum[index], (*it).GetDeathMomentum());
    setPosition(DeathPosition[index], (*it).GetDeathPosition());
    
    setMomentum(MomentumZIsLArEnd[index], (*it).GetMomentumZIsLArEnd());
    setPosition(PositionZIsLArEnd[index], (*it).GetPositionZIsLArEnd());
    
    setMomentum(MomentumZIsTMSStart[index], (*it).GetMomentumZIsTMSStart());
    setPosition(PositionZIsTMSStart[index], (*it).GetPositionZIsTMSStart());
    
    setMomentum(MomentumZIsTMSEnd[index], (*it).GetMomentumZIsTMSEnd());
    setPosition(PositionZIsTMSEnd[index], (*it).GetPositionZIsTMSEnd());
    
    setMomentum(MomentumLArStart[index], (*it).GetMomentumEnteringLAr());
    setPosition(PositionLArStart[index], (*it).GetPositionEnteringLAr());
    
    setMomentum(MomentumLArEnd[index], (*it).GetMomentumLeavingLAr());
    setPosition(PositionLArEnd[index], (*it).GetPositionLeavingLAr());
    
    setMomentum(MomentumTMSStart[index], (*it).GetMomentumEnteringTMS());
    setPosition(PositionTMSStart[index], (*it).GetPositionEnteringTMS());
    
    setMomentum(MomentumTMSEnd[index], (*it).GetMomentumLeavingTMS());
    setPosition(PositionTMSEnd[index], (*it).GetPositionLeavingTMS());
    
    setMomentum(MomentumTMSThinEnd[index], (*it).GetMomentumLeavingTMSThin());
    setPosition(PositionTMSThinEnd[index], (*it).GetPositionLeavingTMSThin());
    
    setMomentum(MomentumTMSFirstTwoModulesEnd[index], (*it).GetMomentumLeavingTMSFirstTwoModules());
    setPosition(PositionTMSFirstTwoModulesEnd[index], (*it).GetPositionLeavingTMSFirstTwoModules());
    
  }

  // Fill the reco info
  std::vector<std::pair<bool, TF1*>> HoughLinesU = TMS_TrackFinder::GetFinder().GetHoughLinesU();
  std::vector<std::pair<bool, TF1*>> HoughLinesV = TMS_TrackFinder::GetFinder().GetHoughLinesV();
  std::vector<std::pair<bool, TF1*>> HoughLinesX = TMS_TrackFinder::GetFinder().GetHoughLinesX();
  // Also get the size of the hits to get a measure of relative goodness
  std::vector<std::vector<TMS_Hit> > HoughCandidatesU = TMS_TrackFinder::GetFinder().GetHoughCandidatesU();
  std::vector<std::vector<TMS_Hit> > HoughCandidatesV = TMS_TrackFinder::GetFinder().GetHoughCandidatesV();
  std::vector<std::vector<TMS_Hit> > HoughCandidatesX = TMS_TrackFinder::GetFinder().GetHoughCandidatesX();
  nLinesU = HoughCandidatesU.size();
  nLinesV = HoughCandidatesV.size();
  nLinesX = HoughCandidatesX.size();


  // Skip the event if there aren't any Hough Lines
  if (nLinesU > __TMS_MAX_LINES__) {
    std::cerr << "Exceeded max number of HoughLines to write to file" << std::endl;
    std::cerr << "Max lines: " << __TMS_MAX_LINES__ << std::endl;
    std::cerr << "Number of lines in event (u): " << nLinesU << std::endl;
    std::cerr << "Not writing event" << std::endl;
    return;
  }

  if (nLinesV > __TMS_MAX_LINES__) {
    std::cerr << "Exceeded max number of HoughLines to write to file" << std::endl;
    std::cerr << "Max lines: " << __TMS_MAX_LINES__ << std::endl;
    std::cerr << "Number of lines in event (v): " << nLinesV << std::endl;
    std::cerr << "Not writing evet" << std::endl;
    return;
  }

  if (nLinesX > __TMS_MAX_LINES__) {
    std::cerr << "Exceeded max number of HoughLines to write to file" << std::endl;
    std::cerr << "Max lines: " << __TMS_MAX_LINES__ << std::endl;
    std::cerr << "Number of lines in event (x): " << nLinesX << std::endl;
    std::cerr << "Not writing evet" << std::endl;
    return;
  }

  int it = 0;
  for (auto &Lines: HoughLinesU) {
    // Get the slopes saved down
    InterceptU[it] = Lines.second->GetParameter(0);
    Intercept_UpstreamU[it] = TMS_TrackFinder::GetFinder().GetHoughLinesU_Upstream()[it].first;
    Intercept_DownstreamU[it] = TMS_TrackFinder::GetFinder().GetHoughLinesU_Downstream()[it].first;
    SlopeU[it] = Lines.second->GetParameter(1);
    Slope_UpstreamU[it] = TMS_TrackFinder::GetFinder().GetHoughLinesU_Upstream()[it].second;
    Slope_DownstreamU[it] = TMS_TrackFinder::GetFinder().GetHoughLinesU_Downstream()[it].second;

    // Calculate the z and x vectors by evaling the TF1 in thin and thick target
    double xlow = TMS_Const::TMS_Thin_Start;
    double xhi = TMS_Const::TMS_Thick_Start;
    double ylow = Lines.second->Eval(xlow);
    double yhi = Lines.second->Eval(xhi);
    // Do the same for up and downstream portions
    double ylow_upstream = Intercept_UpstreamU[it]+xlow*Slope_UpstreamU[it];
    double yhi_upstream = Intercept_UpstreamU[it]+xhi*Slope_UpstreamU[it];

    double ylow_downstream = Intercept_DownstreamU[it]+xlow*Slope_DownstreamU[it];
    double yhi_downstream = Intercept_DownstreamU[it]+xhi*Slope_DownstreamU[it];

    double xlen = xhi-xlow;
    double ylen = yhi-ylow;
    double len = sqrt(xlen*xlen+ylen*ylen);

    double ylen_upstream = yhi_upstream - ylow_upstream;
    double ylen_downstream = yhi_downstream - ylow_downstream;

    double len_upstream = sqrt(xlen*xlen+ylen_upstream*ylen_upstream);
    double len_downstream = sqrt(xlen*xlen+ylen_downstream*ylen_downstream);

    DirectionZU[it] = xlen/len;
    DirectionXU[it] = ylen/len;

    DirectionZU_Upstream[it] = xlen/len_upstream;
    DirectionXU_Upstream[it] = ylen_upstream/len_upstream;

    DirectionZU_Downstream[it] = xlen/len_downstream;
    DirectionXU_Downstream[it] = ylen_downstream/len_downstream;

    it++;
  }
  it = 0;
  for (auto &Lines : HoughLinesV) {
    // Get the slopes saved down
    InterceptV[it] = Lines.second->GetParameter(0);
    Intercept_UpstreamV[it] = TMS_TrackFinder::GetFinder().GetHoughLinesV_Upstream()[it].first;
    Intercept_DownstreamV[it] = TMS_TrackFinder::GetFinder().GetHoughLinesV_Downstream()[it].first;
    SlopeV[it] = Lines.second->GetParameter(1);
    Slope_UpstreamV[it] = TMS_TrackFinder::GetFinder().GetHoughLinesV_Upstream()[it].second;
    Slope_DownstreamV[it] = TMS_TrackFinder::GetFinder().GetHoughLinesV_Downstream()[it].second;

    // Calculate the z and x vectors by evaling the TF1 in thin and thick target
    double xlow = TMS_Const::TMS_Thin_Start;
    double xhi = TMS_Const::TMS_Thick_Start;
    double ylow = Lines.second->Eval(xlow);
    double yhi = Lines.second->Eval(xhi);
    // Doe the same for up and downstream portions
    double ylow_upstream = Intercept_UpstreamV[it]+xlow*Slope_UpstreamV[it];
    double yhi_upstream = Intercept_UpstreamV[it]+xhi*Slope_UpstreamV[it];

    double ylow_downstream = Intercept_DownstreamV[it]+xlow*Slope_DownstreamV[it];
    double yhi_downstream = Intercept_DownstreamV[it]+xhi*Slope_DownstreamV[it];

    double xlen = xhi-xlow;
    double ylen = yhi-ylow;
    double len = sqrt(xlen*xlen+ylen*ylen);

    double ylen_upstream = yhi_upstream - ylow_upstream;
    double ylen_downstream = yhi_downstream - ylow_downstream;

    double len_upstream = sqrt(xlen*xlen+ylen_upstream*ylen_upstream);
    double len_downstream = sqrt(xlen*xlen+ylen_downstream*ylen_downstream);

    DirectionZV[it] = xlen/len;
    DirectionXV[it] = ylen/len;

    DirectionZV_Upstream[it] = xlen/len_upstream;
    DirectionXV_Upstream[it] = ylen_upstream/len_upstream;

    DirectionZV_Downstream[it] = xlen/len_downstream;
    DirectionXV_Downstream[it] = ylen_downstream/len_downstream;

    it++;
  }

  it = 0;
  for (auto &Lines: HoughLinesX) {
    // Get the slopes saved down
    InterceptX[it] = Lines.second->GetParameter(0);
    Intercept_UpstreamX[it] = TMS_TrackFinder::GetFinder().GetHoughLinesX_Upstream()[it].first;
    Intercept_DownstreamX[it] = TMS_TrackFinder::GetFinder().GetHoughLinesX_Downstream()[it].first;
    SlopeX[it] = Lines.second->GetParameter(1);
    Slope_UpstreamX[it] = TMS_TrackFinder::GetFinder().GetHoughLinesX_Upstream()[it].second;
    Slope_DownstreamX[it] = TMS_TrackFinder::GetFinder().GetHoughLinesX_Downstream()[it].second;

    // Calculate the z and x vectors by evaling the TF1 in thin and thick target
    double xlow = TMS_Const::TMS_Thin_Start;
    double xhi = TMS_Const::TMS_Thick_Start;
    double ylow = Lines.second->Eval(xlow);
    double yhi = Lines.second->Eval(xhi);
    // Do the same for up and downstream portions
    double ylow_upstream = Intercept_UpstreamX[it]+xlow*Slope_UpstreamX[it];
    double yhi_upstream = Intercept_UpstreamX[it]+xhi*Slope_UpstreamX[it];

    double ylow_downstream = Intercept_DownstreamX[it]+xlow*Slope_DownstreamX[it];
    double yhi_downstream = Intercept_DownstreamX[it]+xhi*Slope_DownstreamX[it];

    double xlen = xhi-xlow;
    double ylen = yhi-ylow;
    double len = sqrt(xlen*xlen+ylen*ylen);

    double ylen_upstream = yhi_upstream - ylow_upstream;
    double ylen_downstream = yhi_downstream - ylow_downstream;

    double len_upstream = sqrt(xlen*xlen+ylen_upstream*ylen_upstream);
    double len_downstream = sqrt(xlen*xlen+ylen_downstream*ylen_downstream);

    DirectionZX[it] = xlen/len;
    DirectionYX[it] = ylen/len;

    DirectionZX_Upstream[it] = xlen/len_upstream;
    DirectionYX_Upstream[it] = ylen_upstream/len_upstream;

    DirectionZX_Downstream[it] = xlen/len_downstream;
    DirectionYX_Downstream[it] = ylen_downstream/len_downstream;

    it++;
  }


  // Find where the first hough hit is
  TMSStart = false;
  TMSStartTime = -9999.0;

  std::vector<std::vector<TMS_Hit> > HoughCandsU = TMS_TrackFinder::GetFinder().GetHoughCandidatesU();
  int TotalHits = TMS_TrackFinder::GetFinder().GetCleanedHits().size();
  TMS_Hit *FirstTrack = NULL;

  it = 0;
  for (auto &Candidates: HoughCandsU) {
    // Loop over hits
    for (auto &hit: Candidates) {
      if (FirstTrack == NULL) {
        FirstTrack = &hit;
      } else if (hit.GetZ() < FirstTrack->GetZ()) {
        FirstTrack = &hit;
      }
    }
    nHitsInTrackU[it] = Candidates.size();

    // Then save the hit info
    FirstPlaneU[it] = Candidates.front().GetPlaneNumber();
    FirstHitU[it][0] = Candidates.front().GetZ();
    FirstHitU[it][1] = Candidates.front().GetNotZ();
    FirstHitTimeU[it] = Candidates.front().GetT();

    LastPlaneU[it] = Candidates.back().GetPlaneNumber();
    LastHitU[it][0] = Candidates.back().GetZ();
    LastHitU[it][1] = Candidates.back().GetNotZ();
    LastHitTimeU[it] = Candidates.back().GetT();

    TrackLengthU[it] = TMS_TrackFinder::GetFinder().GetTrackLengthU()[it];
    TotalTrackEnergyU[it] = TMS_TrackFinder::GetFinder().GetTrackEnergyU()[it];
    OccupancyU[it] = double(HoughCandsU[it].size())/TotalHits;
    
    float earliest_hit_time = 1e32;
    float latest_hit_time = -1e32;
    // Get each hit in the track and save its energy
    for (unsigned int j = 0; j < Candidates.size(); ++j) {
      TrackHitEnergyU[it][j] = Candidates[j].GetE();
      TrackHitTimeU[it][j] = Candidates[j].GetT();
      TrackHitPosU[it][j][0] = Candidates[j].GetZ();
      TrackHitPosU[it][j][1] = Candidates[j].GetNotZ();
      
      float time = Candidates[j].GetT();
      
      if (time < earliest_hit_time) earliest_hit_time = time;
      if (time > latest_hit_time) latest_hit_time = time;
    }
    EarliestHitTimeU[it] = earliest_hit_time;
    LatestHitTimeU[it] = latest_hit_time;
    

    double maxenergy = 0;
    unsigned int nLastHits_temp = nLastHits;
    if ((unsigned int)nLastHits > Candidates.size()) nLastHits_temp = Candidates.size();
    for (unsigned int i = 0; i < nLastHits_temp; ++i) {
      double hitenergy = Candidates[nLastHits_temp-1-i].GetE();
      if (hitenergy > maxenergy) maxenergy = hitenergy;
    }
    if (maxenergy > EnergyCut) TrackStoppingU[it] = true;


    it++;
  }

  std::vector<std::vector<TMS_Hit> > HoughCandsV = TMS_TrackFinder::GetFinder().GetHoughCandidatesV();
  it = 0;
  for (auto &Candidates: HoughCandsV) {
    // Loop over hits
    for (auto &hit: Candidates) {
      if (FirstTrack == NULL) {
        FirstTrack = &hit;
      } else if (hit.GetZ() < FirstTrack->GetZ()) {
        FirstTrack = &hit;
      }
    }
    nHitsInTrackV[it] = Candidates.size();

    // then save the hit info
    FirstPlaneV[it] = Candidates.front().GetPlaneNumber();
    FirstHitV[it][0] = Candidates.front().GetZ();
    FirstHitV[it][1] = Candidates.front().GetNotZ();
    FirstHitTimeV[it] = Candidates.front().GetT();

    LastPlaneV[it] = Candidates.back().GetPlaneNumber();
    LastHitV[it][0] = Candidates.back().GetZ();
    LastHitV[it][1] = Candidates.back().GetNotZ();
    LastHitTimeV[it] = Candidates.back().GetT();

    TrackLengthV[it] = TMS_TrackFinder::GetFinder().GetTrackLengthV()[it];
    TotalTrackEnergyV[it] = TMS_TrackFinder::GetFinder().GetTrackEnergyV()[it];
    OccupancyV[it] = double(HoughCandsV[it].size())/TotalHits;

    float earliest_hit_time = 1e32;
    float latest_hit_time = -1e32;
    // Get each hit in the track and save its energy
    for (unsigned int j = 0; j < Candidates.size(); ++j) {
      TrackHitEnergyV[it][j] = Candidates[j].GetE();
      TrackHitTimeV[it][j] = Candidates[j].GetT();
      TrackHitPosV[it][j][0] = Candidates[j].GetZ();
      TrackHitPosV[it][j][1] = Candidates[j].GetNotZ();

      float time = Candidates[j].GetT();

      if (time < earliest_hit_time) earliest_hit_time = time;
      if (time > latest_hit_time) latest_hit_time = time;
    }
    EarliestHitTimeV[it] = earliest_hit_time;
    LatestHitTimeV[it] = latest_hit_time;

    double maxenergy = 0;
    unsigned int nLastHits_temp = nLastHits;
    if ((unsigned int)nLastHits > Candidates.size()) nLastHits_temp = Candidates.size();
    for (unsigned int i = 0; i < nLastHits_temp; ++i) {
      double hitenergy = Candidates[nLastHits_temp-1-i].GetE();
      if (hitenergy > maxenergy) maxenergy = hitenergy;
    }
    if (maxenergy > EnergyCut) TrackStoppingV[it] = true;
    

    it++;
  }

  std::vector<std::vector<TMS_Hit> > HoughCandsX = TMS_TrackFinder::GetFinder().GetHoughCandidatesX();
  it = 0;
  for (auto &Candidates: HoughCandsX) {
    // Loop over hits
    for (auto &hit: Candidates) {
      if (FirstTrack == NULL) {
        FirstTrack = &hit;
      } else if (hit.GetZ() < FirstTrack->GetZ()) {
        FirstTrack = &hit;
      }
    }
    nHitsInTrackX[it] = Candidates.size();

    // Then save the hit info
    FirstPlaneX[it] = Candidates.front().GetPlaneNumber();
    FirstHitX[it][0] = Candidates.front().GetZ();
    FirstHitX[it][1] = Candidates.front().GetNotZ();
    FirstHitTimeX[it] = Candidates.front().GetT();

    LastPlaneX[it] = Candidates.back().GetPlaneNumber();
    LastHitX[it][0] = Candidates.back().GetZ();
    LastHitX[it][1] = Candidates.back().GetNotZ();
    LastHitTimeX[it] = Candidates.back().GetT();

    TrackLengthX[it] = TMS_TrackFinder::GetFinder().GetTrackLengthX()[it];
    TotalTrackEnergyX[it] = TMS_TrackFinder::GetFinder().GetTrackEnergyX()[it];
    OccupancyX[it] = double(HoughCandsX[it].size())/TotalHits;
    
    float earliest_hit_time = 1e32;
    float latest_hit_time = -1e32;
    // Get each hit in the track and save its energy
    for (unsigned int j = 0; j < Candidates.size(); ++j) {
      TrackHitEnergyX[it][j] = Candidates[j].GetE();
      TrackHitTimeX[it][j] = Candidates[j].GetT();
      TrackHitPosX[it][j][0] = Candidates[j].GetZ();
      TrackHitPosX[it][j][1] = Candidates[j].GetNotZ();
      
      float time = Candidates[j].GetT();
      
      if (time < earliest_hit_time) earliest_hit_time = time;
      if (time > latest_hit_time) latest_hit_time = time;
    }
    EarliestHitTimeX[it] = earliest_hit_time;
    LatestHitTimeX[it] = latest_hit_time;
    

    double maxenergy = 0;
    unsigned int nLastHits_temp = nLastHits;
    if ((unsigned int)nLastHits > Candidates.size()) nLastHits_temp = Candidates.size();
    for (unsigned int i = 0; i < nLastHits_temp; ++i) {
      double hitenergy = Candidates[nLastHits_temp-1-i].GetE();
      if (hitenergy > maxenergy) maxenergy = hitenergy;
    }
    if (maxenergy > EnergyCut) TrackStoppingX[it] = true;


    it++;
  }


/*  std::vector<TMS_Track> HoughCands3D = TMS_TrackFinder::GetFinder().GetHoughTracks3D(); //TODO this should be done by Liam with Reco_Tree
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
  std::vector<std::vector<TMS_Hit> > ClustersU = TMS_TrackFinder::GetFinder().GetClusterCandidatesU();
  std::vector<std::vector<TMS_Hit> > ClustersV = TMS_TrackFinder::GetFinder().GetClusterCandidatesV();
  std::vector<std::vector<TMS_Hit> > ClustersX = TMS_TrackFinder::GetFinder().GetClusterCandidatesX();
  nClustersU = ClustersU.size();
  nClustersV = ClustersV.size();
  nClustersX = ClustersX.size();
  if (nClustersU > __TMS_MAX_CLUSTERS__) {
    std::cerr << "Too many clusters in TMS_TreeWriter" << std::endl;
    std::cerr << "Hard-coded maximum: " << __TMS_MAX_CLUSTERS__ << std::endl;
    std::cerr << "nClustersU in event: " << nClustersU << std::endl;
    return;
  }
  if (nClustersV > __TMS_MAX_CLUSTERS__) {
    std::cerr << "Too many clusters in TMS_TreeWriter" << std::endl;
    std::cerr << "Hard-coded maximum: " << __TMS_MAX_CLUSTERS__ << std::endl;
    std::cerr << "nClustersV in event: " << nClustersV << std::endl;
    return;
  }
  if (nClustersX > __TMS_MAX_CLUSTERS__) {
    std::cerr << "Too many clusters in TMS_TreeWriter" << std::endl;
    std::cerr << "Hard-coded maximum: " << __TMS_MAX_CLUSTERS__ << std::endl;
    std::cerr << "nClustersX in event: " << nClustersX << std::endl;
    return;
  }

  // Sum up the cluster candidates
  int stdit = 0;
  // Calculate the cluster by cluster summaries
  // e.g. total energy in cluster, cluster position, and cluster standard deviation
  for (auto it = ClustersU.begin(); it != ClustersU.end(); ++it, ++stdit) {
    double total_energy = 0;
    // Mean of cluster in z and not z
    double mean_z = 0;
    double mean_notz = 0;
    // Mean of square of cluster in z and not z
    double mean2_z = 0;
    double mean2_notz = 0;
    double min_cluster_time = 1e10;
    int nhits = (*it).size();
//    std::cout << "Cluster One nHits: " << nhits << std::endl;
    for (int j = 0; j < nhits; ++j) {
      mean_z += (*it)[j].GetZ();
      mean_notz += (*it)[j].GetNotZ();
      mean2_z += (*it)[j].GetZ()*(*it)[j].GetZ();
      mean2_notz += (*it)[j].GetNotZ()*(*it)[j].GetNotZ();
      total_energy += (*it)[j].GetE();
      float time = (*it)[j].GetT();
      if (time < min_cluster_time) min_cluster_time = time;
//      std::cout << (*it)[j].GetZ() << " " << (*it)[j].GetNotZ() << std::endl;
      ClusterHitPosU[stdit][j][0] = (*it)[j].GetZ();
      ClusterHitPosU[stdit][j][1] = (*it)[j].GetNotZ();
      ClusterHitEnergyU[stdit][j] = (*it)[j].GetE();
      ClusterHitTimeU[stdit][j] = (*it)[j].GetT();
      ClusterHitSliceU[stdit][j] = (*it)[j].GetSlice();
    }
    mean_z /= nhits;
    mean_notz /= nhits;

    mean2_z /= nhits;
    mean2_notz /= nhits;

    ClusterEnergyU[stdit] = total_energy;
    ClusterTimeU[stdit] = min_cluster_time;
    nHitsInClusterU[stdit] = nhits;
    ClusterPosMeanU[stdit][0] = mean_z;
    ClusterPosMeanU[stdit][1] = mean_notz;
    // Calculate the standard deviation
    double std_dev_z = mean2_z-mean_z*mean_z;
    double std_dev_notz = mean2_z-mean_z*mean_z;
    if (std_dev_z > 0) std_dev_z = sqrt(std_dev_z);
    if (std_dev_notz > 0) std_dev_notz = sqrt(std_dev_notz);
    ClusterPosStdDevU[stdit][0] = std_dev_z;
    ClusterPosStdDevU[stdit][1] = std_dev_notz;
  }
  stdit = 0;
  for (auto it = ClustersV.begin(); it != ClustersV.end(); ++it, ++stdit) {
    double total_energy = 0;
    double mean_z = 0;
    double mean_notz = 0;
    double mean2_z = 0;
    double mean2_notz = 0;
    double min_cluster_time = 1e10;
    int nhits = (*it).size();
//    std::cout << "Cluster Other nHits: " << nhits << std::endl;
    for (int j = 0; j < nhits; ++j) {
      mean_z += (*it)[j].GetZ();
      mean_notz += (*it)[j].GetNotZ();
      mean2_z += (*it)[j].GetZ()*(*it)[j].GetZ();
      mean2_notz += (*it)[j].GetNotZ()*(*it)[j].GetNotZ();
      total_energy += (*it)[j].GetE();
      float time = (*it)[j].GetT();
      if (time < min_cluster_time) min_cluster_time = time;
//      std::cout << (*it)[j].GetZ() << " " << (*it)[j].GetNotZ() << std::endl;
      ClusterHitPosV[stdit][j][0] = (*it)[j].GetZ();
      ClusterHitPosV[stdit][j][1] = (*it)[j].GetNotZ();
      ClusterHitEnergyV[stdit][j] = (*it)[j].GetE();
      ClusterHitTimeV[stdit][j] = (*it)[j].GetT();
      ClusterHitSliceV[stdit][j] = (*it)[j].GetSlice();
    }
    mean_z /= nhits;
    mean_notz /= nhits;

    mean2_z /= nhits;
    mean2_notz /= nhits;

    ClusterEnergyV[stdit] = total_energy;
    ClusterTimeV[stdit] = min_cluster_time;
    nHitsInClusterV[stdit] = nhits;
    ClusterPosMeanV[stdit][0] = mean_z;
    ClusterPosMeanV[stdit][1] = mean_notz;

    double std_dev_z = mean2_z-mean_z*mean_z;
    double std_dev_notz = mean2_z-mean_z*mean_z;
    if (std_dev_z > 0) std_dev_z = sqrt(std_dev_z);
    if (std_dev_notz > 0) std_dev_notz = sqrt(std_dev_notz);
    ClusterPosStdDevV[stdit][0] = std_dev_z;
    ClusterPosStdDevV[stdit][1] = std_dev_notz;
  }
  stdit = 0;
  for (auto it = ClustersX.begin(); it != ClustersX.end(); ++it, ++stdit) {
    double total_energy = 0;
    // Mean of cluster in z and not z
    double mean_z = 0;
    double mean_notz = 0;
    // Mean of square of cluster in z and not z
    double mean2_z = 0;
    double mean2_notz = 0;
    double min_cluster_time = 1e10;
    int nhits = (*it).size();
//    std::cout << "Cluster One nHits: " << nhits << std::endl;
    for (int j = 0; j < nhits; ++j) {
      mean_z += (*it)[j].GetZ();
      mean_notz += (*it)[j].GetNotZ();
      mean2_z += (*it)[j].GetZ()*(*it)[j].GetZ();
      mean2_notz += (*it)[j].GetNotZ()*(*it)[j].GetNotZ();
      total_energy += (*it)[j].GetE();
      float time = (*it)[j].GetT();
      if (time < min_cluster_time) min_cluster_time = time;
//      std::cout << (*it)[j].GetZ() << " " << (*it)[j].GetNotZ() << std::endl;
      ClusterHitPosX[stdit][j][0] = (*it)[j].GetZ();
      ClusterHitPosX[stdit][j][1] = (*it)[j].GetNotZ();
      ClusterHitEnergyX[stdit][j] = (*it)[j].GetE();
      ClusterHitTimeX[stdit][j] = (*it)[j].GetT();
      ClusterHitSliceX[stdit][j] = (*it)[j].GetSlice();
    }
    mean_z /= nhits;
    mean_notz /= nhits;

    mean2_z /= nhits;
    mean2_notz /= nhits;

    ClusterEnergyX[stdit] = total_energy;
    ClusterTimeX[stdit] = min_cluster_time;
    nHitsInClusterX[stdit] = nhits;
    ClusterPosMeanX[stdit][0] = mean_z;
    ClusterPosMeanX[stdit][1] = mean_notz;
    // Calculate the standard deviation
    double std_dev_z = mean2_z-mean_z*mean_z;
    double std_dev_notz = mean2_z-mean_z*mean_z;
    if (std_dev_z > 0) std_dev_z = sqrt(std_dev_z);
    if (std_dev_notz > 0) std_dev_notz = sqrt(std_dev_notz);
    ClusterPosStdDevX[stdit][0] = std_dev_z;
    ClusterPosStdDevX[stdit][1] = std_dev_notz;
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
  RecoTrackN = Reco_Tracks.size();

  for (auto RecoTrack = Reco_Tracks.begin(); RecoTrack != Reco_Tracks.end(); ++RecoTrack, ++itTrack) {
    nHitsIn3DTrack[itTrack]         = (int) RecoTrack->Hits.size(); // Do we need to cast it? idk
//    std::cout << "TreeWriter number of hits: " << nHitsIn3DTrack[itTrack] << std::endl;
    RecoTrackEnergyRange[itTrack]   =       RecoTrack->EnergyRange;
    RecoTrackLength[itTrack]        =       RecoTrack->Length;
    RecoTrackEnergyDeposit[itTrack] =       RecoTrack->EnergyDeposit;
    for (int j = 0; j < 3; j++) {
      RecoTrackStartPos[itTrack][j]  = RecoTrack->Start[j];
      RecoTrackEndPos[itTrack][j]    = RecoTrack->End[j];
      RecoTrackDirection[itTrack][j] = RecoTrack->Direction[j];
    }
    for (unsigned int j = 0; j < RecoTrack->Hits.size(); ++j) {
      if (RecoTrack->Hits[j].GetBar().GetBarType() != TMS_Bar::kXBar) {
        RecoTrackHitPos[itTrack][j][0] = RecoTrack->Hits[j].GetRecoX();
        RecoTrackHitPos[itTrack][j][1] = RecoTrack->Hits[j].GetRecoY();
      } else if (RecoTrack->Hits[j].GetBar().GetBarType() == TMS_Bar::kXBar) {
        RecoTrackHitPos[itTrack][j][0] = RecoTrack->Hits[j].GetRecoX();
        RecoTrackHitPos[itTrack][j][1] = RecoTrack->Hits[j].GetNotZ();
      }
      RecoTrackHitPos[itTrack][j][2] = RecoTrack->Hits[j].GetZ();
//      std::cout << "TreeWriter hit position: " << RecoTrackHitPos[itTrack][j][0] << " " << RecoTrackHitPos[itTrack][j][1] << " " << RecoTrackHitPos[itTrack][j][2] << std::endl;
    }
    // Can manually compute direction if it hasn't been set
    if ( (RecoTrackDirection[itTrack][0] == 0) && (RecoTrackDirection[itTrack][1] == 0) && (RecoTrackDirection[itTrack][2] == 0) )
    { // If true^ it seems the direction hasn't been set
      for (int j = 0; j < 3; j++)
      { // Right now no need to make sure this is a unit vector
        RecoTrackDirection[itTrack][j] = RecoTrack->End[j] - RecoTrack->Start[j];
      }
    }
    
    // Now fill truth info
    double total_true_visible_energy = 0;
    double true_primary_visible_energy = -999;
    double true_secondary_visible_energy = -999;
    int true_primary_particle_index = -999;
    int true_secondary_particle_index = -999;
    auto particle_info = TMS_Utils::GetPrimaryIdsByEnergy(RecoTrack->Hits);
    total_true_visible_energy = particle_info.total_energy;
    if (particle_info.energies.size() > 0) {
      true_primary_visible_energy = particle_info.energies[0];
      true_primary_particle_index = particle_info.indices[0];
    }
    if (particle_info.energies.size() > 1) {
      true_secondary_visible_energy = particle_info.energies[1];
      true_secondary_particle_index = particle_info.indices[1];
    }
    // Now for the primary index, find the true starting and ending momentum and position
    // TODO add back
    if (true_primary_particle_index < 0) {
      // Do nothing, this means we didn't find a true particle associated with a reco track
      // This can't happen unless dark noise existed which is currently doesn't
      std::cout<<"Error: Found true_primary_particle_index < 0. There should be at least one particle creating energy (assuming no dark noise) but instead the index is: "<<true_primary_particle_index<<", with energy: "<<true_primary_visible_energy<<std::endl;
    }
    else if ((size_t)true_primary_particle_index >= TrueParticles.size()) {
      std::cout<<"Error: Found true_primary_particle_index >= TrueParticles.size(). This shouldn't happen since each index should point to a true particle"<<std::endl;
    }
    else {
      TMS_TrueParticle tp = TrueParticles[true_primary_particle_index];
      double start_z = RecoTrack->Start[2];
      double end_z = RecoTrack->End[2];
      const double max_z_distance = 1e9; // Want the closest possible starting and ending points, regardless of distance
      if (itTrack < __TMS_MAX_LINES__) {
        setMomentum(RecoTrackPrimaryParticleTrueMomentumTrackStart[itTrack], tp.GetMomentumAtZ(start_z, max_z_distance));
        setPosition(RecoTrackPrimaryParticleTruePositionTrackStart[itTrack], tp.GetPositionAtZ(start_z, max_z_distance));
        setMomentum(RecoTrackPrimaryParticleTrueMomentumTrackEnd[itTrack], tp.GetMomentumAtZ(end_z, max_z_distance));
        setPosition(RecoTrackPrimaryParticleTruePositionTrackEnd[itTrack], tp.GetPositionAtZ(end_z, max_z_distance));
      }
    }
    
    RecoTrackTrueVisibleEnergy[itTrack] = total_true_visible_energy;
    RecoTrackPrimaryParticleIndex[itTrack] = true_primary_particle_index;
    RecoTrackPrimaryParticleTrueVisibleEnergy[itTrack] = true_primary_visible_energy;
    RecoTrackSecondaryParticleIndex[itTrack] = true_secondary_particle_index;
    RecoTrackSecondaryParticleTrueVisibleEnergy[itTrack] = true_secondary_visible_energy;
    
  }

  Reco_Tree->Fill();
  Truth_Info->Fill();
}

// Reset the variables
void TMS_TreeWriter::Clear() {

  // Reset truth information
  EventNo = nParticles = NeutrinoPDG = LeptonPDG = Muon_TrueKE = Muon_TrueTrackLength = VertexIdOfMostEnergyInEvent = -999;
  VertexIdOfMostEnergyInEvent = VisibleEnergyFromUVertexInSlice = TotalVisibleEnergyFromVertex = VisibleEnergyFromVVerticesInSlice = -999;
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
  nLinesU = -999;
  nLinesV = -999;
  nLinesX = -999;
  for (int i = 0; i < __TMS_MAX_LINES__; ++i) {
    SlopeU[i] = -999;
    SlopeV[i] = -999;
    SlopeX[i] = -999;
    InterceptU[i] = -999;
    InterceptV[i] = -999;
    InterceptX[i] = -999;
    Slope_DownstreamU[i] = -999;
    Slope_DownstreamV[i] = -999;
    Slope_DownstreamX[i] = -999;
    Intercept_DownstreamU[i] = -999;
    Intercept_DownstreamV[i] = -999;
    Intercept_DownstreamX[i] = -999;
    Slope_UpstreamU[i] = -999;
    Slope_UpstreamV[i] = -999;
    Slope_UpstreamX[i] = -999;
    Intercept_UpstreamU[i] = -999;
    Intercept_UpstreamV[i] = -999;
    Intercept_UpstreamX[i] = -999;

    DirectionZU[i] = -999;
    DirectionXU[i] = -999;
    DirectionZU_Upstream[i] = -999;
    DirectionXU_Upstream[i] = -999;
    DirectionZU_Downstream[i] = -999;
    DirectionXU_Downstream[i] = -999;

    DirectionZV[i] = -999;
    DirectionXV[i] = -999;
    DirectionZV_Upstream[i] = -999;
    DirectionXV_Upstream[i] = -999;
    DirectionZV_Downstream[i] = -999;
    DirectionXV_Downstream[i] = -999;

    DirectionZX[i] = -999;
    DirectionYX[i] = -999;
    DirectionZX_Upstream[i] = -999;
    DirectionYX_Upstream[i] = -999;
    DirectionZX_Downstream[i] = -999;
    DirectionYX_Downstream[i] = -999;

    OccupancyU[i] = -999;
    OccupancyV[i] = -999;
    OccupancyX[i] = -999;
    TrackLengthU[i] = -999;
    TrackLengthV[i] = -999;
    TrackLengthX[i] = -999;
    TotalTrackEnergyU[i] = -999;
    TotalTrackEnergyV[i] = -999;
    TotalTrackEnergyX[i] = -999;
    FirstPlaneU[i] = -999;
    FirstPlaneV[i] = -999;
    FirstPlaneX[i] = -999;
    LastPlaneU[i] = -999;
    LastPlaneV[i] = -999;
    LastPlaneX[i] = -999;
    nHitsInTrackU[i] = -999;
    nHitsInTrackV[i] = -999;
    nHitsInTrackX[i] = -999;
    TrackStoppingU[i] = false;
    TrackStoppingV[i] = false;
    TrackStoppingX[i] = false;
  }
/*    Occupancy3D[i] = -999;
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
  }*/

  // Reset hit information
  nHits = -999;
  for (int i = 0; i < __TMS_MAX_HITS__; ++i) {
    for (int j = 0; j < 4; ++j) RecoHitPos[i][j] = -999;
    RecoHitEnergy[i] = -999;
  }

  // Reset Cluster info
  nClustersU = -999;
  nClustersV = -999;
  nClustersX = -999;
  for (int i = 0; i < __TMS_MAX_CLUSTERS__; ++i) {
    ClusterEnergyU[i] = -999;
    ClusterEnergyV[i] = -999;
    ClusterEnergyX[i] = -999;
    nHitsInClusterU[i] = -999;
    nHitsInClusterV[i] = -999;
    nHitsInClusterX[i] = -999;
    for (int j = 0; j < 2; ++j) {
      ClusterPosMeanU[i][j] = -999;
      ClusterPosStdDevU[i][j] = -999;
      
      ClusterPosMeanV[i][j] = -999;
      ClusterPosStdDevV[i][j] = -999;

      ClusterPosMeanX[i][j] = -999;
      ClusterPosStdDevX[i][j] = -999;
    }
    for (int j = 0; j < __TMS_MAX_LINE_HITS__; ++j) {
      ClusterHitPosU[i][j][0] = -999;
      ClusterHitPosU[i][j][1] = -999;
      ClusterHitEnergyU[i][j] = -999;

      ClusterHitPosV[i][j][0] = -999;
      ClusterHitPosV[i][j][1] = -999;
      ClusterHitEnergyV[i][j] = -999;

      ClusterHitPosX[i][j][0] = -999;
      ClusterHitPosX[i][j][1] = -999;
      ClusterHitEnergyX[i][j] = -999;
    }
  }

  // Reset track information
  nTracks = -999;
  for (int i = 0; i < __TMS_MAX_TRACKS__; ++i) {
    for (int j = 0; j < 3; ++j) {
      RecoTrackStartPos[i][j] = -999;
      RecoTrackDirection[i][j] = -999;
      RecoTrackEndPos[i][j] = -999;
    }
    for (int k = 0; k < __TMS_MAX_LINE_HITS__; ++k) {
      RecoTrackHitPos[i][k][0] = -999;
      RecoTrackHitPos[i][k][1] = -999;
      RecoTrackHitPos[i][k][2] = -999;
    }
    RecoTrackEnergyRange[i] = -999;
    RecoTrackEnergyDeposit[i] = -999;
    RecoTrackLength[i] = -999;
  }


}

