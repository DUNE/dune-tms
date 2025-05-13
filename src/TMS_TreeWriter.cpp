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

  Truth_Spill = new TTree("Truth_Spill", "Truth_Spill");
  Truth_Spill->SetDirectory(Output);
  Truth_Spill->SetAutoSave(__TMS_AUTOSAVE__);

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
  Branch_Lines->Branch("RunNo", &RunNo, "RunNo/I");
  Branch_Lines->Branch("nLinesU", &nLinesU, "nLinesU/I");
  Branch_Lines->Branch("nLinesV", &nLinesV, "nLinesV/I");
  Branch_Lines->Branch("nLinesX", &nLinesX, "nLinesX/I");
  Branch_Lines->Branch("nLinesY", &nLinesY, "nLinesY/I");
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

  Branch_Lines->Branch("SlopeY",                SlopeY,                 "SlopeY[nLinesY]/F");
  Branch_Lines->Branch("InterceptY",            InterceptY,             "InterceptY[nLinesY]/F");
  Branch_Lines->Branch("Slope_DownstreamY",     Slope_DownstreamY,      "Slope_DownstreamY[nLinesY]/F");
  Branch_Lines->Branch("Intercept_DownstreamY", Intercept_DownstreamY,  "Intercept_DownstreamY[nLinesY]/F");
  Branch_Lines->Branch("Slope_UpstreamY",       Slope_UpstreamY,        "Slope_UpstreamY[nLinesY]/F");
  Branch_Lines->Branch("Intercept_UpstreamY",   Intercept_UpstreamY,    "Intercept_UpstreamY[nLinesY]/F");

  Branch_Lines->Branch("DirectionZU", DirectionZU, "DirectionZU[nLinesU]/F");
  Branch_Lines->Branch("DirectionZV", DirectionZV, "DirectionZV[nLinesV]/F");
  Branch_Lines->Branch("DirectionZX", DirectionZX, "DirectionZX[nLinesX]/F");
  Branch_Lines->Branch("DirectionZY", DirectionZY, "DirectionZY[nLinesY]/F");
  Branch_Lines->Branch("DirectionXU", DirectionXU, "DirectionXU[nLinesU]/F");
  Branch_Lines->Branch("DirectionXV", DirectionXV, "DirectionXV[nLinesV]/F");
  Branch_Lines->Branch("DirectionYX", DirectionYX, "DirectionYX[nLinesX]/F");
  Branch_Lines->Branch("DirectionXY", DirectionXY, "DirectionXY[nLinesY]/F");
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
  Branch_Lines->Branch("DirectionZY_Upstream",    DirectionZY_Upstream,   "DirectionZY_Upstream[nLinesY]/F");
  Branch_Lines->Branch("DirectionZY_Downstream",  DirectionZY_Downstream, "DirectionZY_Downstream[nLinesY]/F");
  Branch_Lines->Branch("DirectionXY_Downstream",  DirectionXY_Downstream, "DirectionXY_Downstream[nLinesY]/F");
  Branch_Lines->Branch("DirectionZY_Upstream",  DirectionZY_Upstream, "DirectionZY_Upstream[nLinesY]/F");
  Branch_Lines->Branch("DirectionXY_Upstream",  DirectionXY_Upstream, "DirectionXY_Upstream[nLinesY]/F");

  Branch_Lines->Branch("FirstHoughHitU",     FirstHitU,     "FirstHoughHitU[nLinesU][2]/F");
  Branch_Lines->Branch("FirstHoughHitV",     FirstHitV,     "FirstHoughHitV[nLinesV][2]/F");
  Branch_Lines->Branch("FirstHoughHitX",     FirstHitX,     "FirstHoughHitX[nLinesX][2]/F");
  Branch_Lines->Branch("FirstHoughHitY",     FirstHitY,     "FirstHoughHitY[nLinesY][2]/F");
  Branch_Lines->Branch("LastHoughHitU",      LastHitU,      "LastHoughHitU[nLinesU][2]/F");
  Branch_Lines->Branch("LastHoughHitV",      LastHitV,      "LastHoughHitV[nLinesV][2]/F");
  Branch_Lines->Branch("LastHoughHitX",      LastHitX,      "LastHoughHitX[nLinesX][2]/F");
  Branch_Lines->Branch("LastHoughHitY",      LastHitY,      "LastHoughHitY[nLinesY][2]/F");
  Branch_Lines->Branch("FirstHoughHitTimeU", FirstHitTimeU, "FirstHoughHitTimeU[nLinesU]/F");
  Branch_Lines->Branch("FirstHoughHitTimeV", FirstHitTimeV, "FirstHoughHitTimeV[nLinesV]/F"); 
  Branch_Lines->Branch("FirstHoughHitTimeX", FirstHitTimeX, "FirstHoughHitTimeX[nLinesX]/F");
  Branch_Lines->Branch("FirstHoughHitTimeY", FirstHitTimeY, "FirstHoughHitTimeY[nLinesY]/F");
  Branch_Lines->Branch("LastHoughHitTimeU",  LastHitTimeU,  "LastHoughHitTimeU[nLinesU]/F");
  Branch_Lines->Branch("LastHoughHitTimeV",  LastHitTimeV,  "LastHoughHitTimeV[nLinesV]/F");
  Branch_Lines->Branch("LastHoughHitTImeX",  LastHitTimeX,  "LastHoughHitTimeX[nLinesX]/F");
  Branch_Lines->Branch("LastHoughHitTImeY",  LastHitTimeY,  "LastHoughHitTimeY[nLinesY]/F");
  Branch_Lines->Branch("HoughEarliestHitTimeU", EarliestHitTimeU, "HoughEarliestHitTimeU[nLinesU]/F");
  Branch_Lines->Branch("HoughEarliestHitTimeV", EarliestHitTimeV, "HoughEarliestHitTimeV[nLinesV]/F");
  Branch_Lines->Branch("HoughEarliestHitTimeX", EarliestHitTimeX, "HoughEarliestHitTimeX[nLinesX]/F");
  Branch_Lines->Branch("HoughEarliestHitTimeY", EarliestHitTimeY, "HoughEarliestHitTimeY[nLinesY]/F");
  Branch_Lines->Branch("HoughLatestHitTimeU",   LatestHitTimeU,   "HoughLatestHitTimeU[nLinesU]/F");
  Branch_Lines->Branch("HoughLatestHitTimeV",   LatestHitTimeV,   "HoughLatestHitTimeV[nLinesV]/F");
  Branch_Lines->Branch("HoughLatestHitTimeX",   LatestHitTimeX,   "HoughLatestHitTimeX[nLinesX]/F");
  Branch_Lines->Branch("HoughLatestHitTimeY",   LatestHitTimeY,   "HoughLatestHitTimeY[nLinesY]/F");
  Branch_Lines->Branch("FirstHoughPlaneU", FirstPlaneU,   "FirstHoughPlaneU[nLinesU]/I");
  Branch_Lines->Branch("FirstHoughPlaneV", FirstPlaneV,   "FirstHoughPlaneV[nLinesV]/I");
  Branch_Lines->Branch("FirstHoughPlaneX", FirstPlaneX,   "FirstHoughPlaneX[nLinesX]/I");
  Branch_Lines->Branch("FirstHoughPlaneY", FirstPlaneY,   "FirstHoughPlaneY[nLinesY]/I");
  Branch_Lines->Branch("LastHoughPlaneU",  LastPlaneU,    "LastHoughPlaneU[nLinesU]/I");
  Branch_Lines->Branch("LastHoughPlaneV",  LastPlaneV,    "LastHoughPlaneV[nLinesV]/I");
  Branch_Lines->Branch("LastHoughPlaneX",  LastPlaneX,    "LastHoughPlaneX[nLinesX]/I");
  Branch_Lines->Branch("LastHoughPlaneY",  LastPlaneY,    "LastHoughPlaneY[nLinesY]/I");
  Branch_Lines->Branch("TMSStart",         &TMSStart,     "TMSStart/O");
  Branch_Lines->Branch("TMSStartTime",     &TMSStartTime, "TMSStartTime/F");
  Branch_Lines->Branch("OccupancyU",       OccupancyU,    "OccupancyU[nLinesU]/F");
  Branch_Lines->Branch("OccupancyV",       OccupancyV,    "OccupancyV[nLinesV]/F");
  Branch_Lines->Branch("OccupancyX",       OccupancyX,    "OccupancyX[nLinesX]/F");
  Branch_Lines->Branch("OccupancyY",       OccupancyY,    "OccupancyY[nLinesY]/F");
  Branch_Lines->Branch("TrackLengthU",     TrackLengthU,  "TrackLengthU[nLinesU]/F");
  Branch_Lines->Branch("TrackLengthV",     TrackLengthV,  "TrackLengthV[nLinesV]/F");
  Branch_Lines->Branch("TrackLengthX",     TrackLengthX,  "TrackLengthX[nLinesX]/F");
  Branch_Lines->Branch("TrackLengthY",     TrackLengthY,  "TrackLengthY[nLinesY]/F");
  Branch_Lines->Branch("TotalTrackEnergyU", TotalTrackEnergyU, "TotalTrackEnergyU[nLinesU]/F");
  Branch_Lines->Branch("TotalTrackEnergyV", TotalTrackEnergyV, "TotalTrackEnergyV[nLinesV]/F");
  Branch_Lines->Branch("TotalTrackEnergyX", TotalTrackEnergyX, "TotalTrackEnergyX[nLinesX]/F");
  Branch_Lines->Branch("TotalTrackEnergyY", TotalTrackEnergyY, "TotalTrackEnergyY[nLinesY]/F");
  Branch_Lines->Branch("TrackStoppingU",    TrackStoppingU,    "TrackStoppingU[nLinesU]/O");
  Branch_Lines->Branch("TrackStoppingV",    TrackStoppingV,    "TrackStoppingV[nLinesV]/O");
  Branch_Lines->Branch("TrackStoppingX",    TrackStoppingX,    "TrackStoppingX[nLinesX]/O");
  Branch_Lines->Branch("TrackStoppingY",    TrackStoppingY,    "TrackStoppingY[nLinesY]/O");

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
  Branch_Lines->Branch("nHitsInTrackY",   &nHitsInTrackY,  "nHitsInTrackY[nLinesY]/I");
  Branch_Lines->Branch("TrackHitEnergyU", TrackHitEnergyU, "TrackHitEnergyU[10][200]/F");
  Branch_Lines->Branch("TrackHitEnergyV", TrackHitEnergyV, "TrackHitEnergyV[10][200]/F");
  Branch_Lines->Branch("TrackHitEnergyX", TrackHitEnergyX, "TrackHitEnergyX[10][200]/F");
  Branch_Lines->Branch("TrackHitEnergyY", TrackHitEnergyY, "TrackHitEnergyY[10][200]/F");
  Branch_Lines->Branch("TrackHitPosU",    TrackHitPosU,    "TrackHitPosU[10][200][2]/F");
  Branch_Lines->Branch("TrackHitPosV",    TrackHitPosV,    "TrackHitPosV[10][200][2]/F");
  Branch_Lines->Branch("TrackHitPosX",    TrackHitPosX,    "TrackHitPosX[10][200][2]/F");
  Branch_Lines->Branch("TrackHitPosY",    TrackHitPosY,    "TrackHitPosY[10][200][2]/F");
  Branch_Lines->Branch("TrackHitTimeU",   TrackHitTimeU,   "TrackHitTimeU[10][200]/F");
  Branch_Lines->Branch("TrackHitTimeV",   TrackHitTimeV,   "TrackHitTimeV[10][200]/F");
  Branch_Lines->Branch("TrackHitTimeX",   TrackHitTimeX,   "TrackHitTimeX[10][200]/F");
  Branch_Lines->Branch("TrackHitTimeY",   TrackHitTimeY,   "TrackHitTimeY[10][200]/F");

  // Cluster information
  Branch_Lines->Branch("nClustersU",        &nClustersU,       "nClustersU/I");
  Branch_Lines->Branch("nClustersV",        &nClustersV,       "nClustersV/I");
  Branch_Lines->Branch("nClusterX",         &nClustersX,       "nClustersX/I");
  Branch_Lines->Branch("nClusterY",         &nClustersY,       "nClustersY/I");
  Branch_Lines->Branch("ClusterEnergyU",    ClusterEnergyU,    "ClusterEnergyU[nClustersU]/F");
  Branch_Lines->Branch("ClusterEnergyV",    ClusterEnergyV,    "ClusterEnergyV[nClustersV]/F");
  Branch_Lines->Branch("ClusterEnergyX",    ClusterEnergyX,    "ClusterEnergyX[nClustersX]/F");
  Branch_Lines->Branch("ClusterEnergyY",    ClusterEnergyY,    "ClusterEnergyY[nClustersY]/F");
  Branch_Lines->Branch("ClusterTimeU",      ClusterTimeU,      "ClusterTimeU[nClustersU]/F");
  Branch_Lines->Branch("ClusterTimeV",      ClusterTimeV,      "ClusterTimeV[nClustersV]/F");
  Branch_Lines->Branch("ClusterTimeX",      ClusterTimeX,      "ClusterTimeX[nClustersX]/F");
  Branch_Lines->Branch("ClusterTimeY",      ClusterTimeY,      "ClusterTimeY[nClustersY]/F");
  Branch_Lines->Branch("ClusterPosMeanU",   ClusterPosMeanU,   "ClusterPosMeanU[25][2]/F");
  Branch_Lines->Branch("ClusterPosMeanV",   ClusterPosMeanV,   "ClusterPosMeanV[25][2]/F");
  Branch_Lines->Branch("ClusterPosMeanX",   ClusterPosMeanX,   "ClusterPosMeanX[25][2]/F");
  Branch_Lines->Branch("ClusterPosMeanY",   ClusterPosMeanY,   "ClusterPosMeanY[25][2]/F");
  Branch_Lines->Branch("ClusterPosStdDevU", ClusterPosStdDevU, "ClusterPosStdDevU[25][2]/F");
  Branch_Lines->Branch("ClusterPosStdDevV", ClusterPosStdDevV, "ClusterPosStdDevV[25][2]/F");
  Branch_Lines->Branch("ClusterPosStdDevX", ClusterPosStdDevX, "ClusterPosStdDevX[25][2]/F");
  Branch_Lines->Branch("ClusterPosStdDevY", ClusterPosStdDevY, "ClusterPosStdDevY[25][2]/F");
  Branch_Lines->Branch("nHitsInClusterU",   nHitsInClusterU,   "nHitsInClusterU[nClustersU]/I");
  Branch_Lines->Branch("nHitsInClusterV",   nHitsInClusterV,   "nHitsInClusterV[nClustersV]/I");
  Branch_Lines->Branch("nHitsInClusterX",   nHitsInClusterX,   "nHitsInClusterX[nClustersX]/I");
  Branch_Lines->Branch("nHitsInClusterY",   nHitsInClusterY,   "nHitsInClusterY[nClustersY]/I");
  Branch_Lines->Branch("ClusterHitPosU",    ClusterHitPosU,    "ClusterHitPosU[25][200][2]/F");
  Branch_Lines->Branch("ClusterHitPosV",    ClusterHitPosV,    "ClusterHitPosV[25][200][2]/F");
  Branch_Lines->Branch("ClusterHitPosX",    ClusterHitPosX,    "ClusterHitPosX[25][200][2]/F");
  Branch_Lines->Branch("ClusterHitPosY",    ClusterHitPosY,    "ClusterHitPosY[25][200][2]/F");
  Branch_Lines->Branch("ClusterHitEnergyU", ClusterHitEnergyU, "ClusterHitEnergyU[25][200]/F");
  Branch_Lines->Branch("ClusterHitEnergyV", ClusterHitEnergyV, "ClusterHitEnergyV[25][200]/F");
  Branch_Lines->Branch("ClusterHitEnergyX", ClusterHitEnergyX, "ClusterHitEnergyX[25][200]/F");
  Branch_Lines->Branch("ClusterHitEnergyY", ClusterHitEnergyY, "ClusterHitEnergyY[25][200]/F");
  Branch_Lines->Branch("ClusterHitTimeU",   ClusterHitTimeU,   "ClusterHitTimeU[25][200]/F");
  Branch_Lines->Branch("ClusterHitTimeV",   ClusterHitTimeV,   "ClusterHitTimeV[25][200]/F");
  Branch_Lines->Branch("ClusterHitTimeX",   ClusterHitTimeX,   "ClusterHitTimeX[25][200]/F");
  Branch_Lines->Branch("ClusterHitTimeY",   ClusterHitTimeY,   "ClusterHitTimeY[25][200]/F");
  Branch_Lines->Branch("ClusterHitSliceU",  ClusterHitSliceU,  "ClusterHitSliceU[25][200]/I");
  Branch_Lines->Branch("ClusterHitSliceV",  ClusterHitSliceV,  "ClusterHitSliceV[25][200]/I");
  Branch_Lines->Branch("ClusterHitSliceX",  ClusterHitSliceX,  "ClusterHitSliceX[25][200]/I");
  Branch_Lines->Branch("ClusterHitSliceY",  ClusterHitSliceY,  "ClusterHitSliceY[25][200]/I");

  // Hit information
  Branch_Lines->Branch("nHits",         &nHits,        "nHits/I");
  Branch_Lines->Branch("RecoHitPos",    RecoHitPos,    "RecoHitPos[nHits][4]/F");
  Branch_Lines->Branch("RecoHitEnergy", RecoHitEnergy, "RecoHitEnergy[nHits]/F");
  Branch_Lines->Branch("RecoHitPE", RecoHitPE, "RecoHitPE[nHits]/F");
  Branch_Lines->Branch("RecoHitBar",  RecoHitBar,  "RecoHitBar[nHits]/I");
  Branch_Lines->Branch("RecoHitBarType", RecoHitBarType, "RecoHitBarType[nHits]/I");
  Branch_Lines->Branch("RecoHitPlane",  RecoHitPlane,  "RecoHitPlane[nHits]/I");
  Branch_Lines->Branch("RecoHitSlice",  RecoHitSlice,  "RecoHitSlice[nHits]/I");

  // Track information
  // TODO: Fill these properly
  Reco_Tree->Branch("EventNo", &EventNo, "EventNo/I");
  Reco_Tree->Branch("SliceNo", &SliceNo, "SliceNo/I");
  Reco_Tree->Branch("SpillNo", &SpillNo, "SpillNo/I");
  Reco_Tree->Branch("RunNo", &RunNo, "RunNo/I");

  Reco_Tree->Branch("nTracks",        &nTracks,                 "nTracks/I");
  Reco_Tree->Branch("nHits",          nHitsIn3DTrack,           "nHits[nTracks]/I");
  Reco_Tree->Branch("TrackHitPos",    RecoTrackHitPos,          "TrackHitPos[nTracks][200][3]/F");
  Reco_Tree->Branch("nKalmanNodes",   nKalmanNodes,             "nKalmanNodes[nTracks]/I");

  Reco_Tree->Branch("KalmanErrorDetVol",KalmanErrorDetVol,      "KalmanErrorDetVol[nTracks]/I");
  Reco_Tree->Branch("nKalmanNodes_plus",   nKalmanNodes_plus,             "nKalmanNodes_plus[nTracks]/I");
  Reco_Tree->Branch("nKalmanNodes_minus",   nKalmanNodes_minus,             "nKalmanNodes_minus[nTracks]/I");

  Reco_Tree->Branch("KalmanPos",      RecoTrackKalmanPos,       "KalmanPos[nTracks][200][3]/F");
  Reco_Tree->Branch("RecoTrackKalmanFirstPlaneBarView", RecoTrackKalmanFirstPlaneBarView, "RecoTrackKalmanFirstPlaneBarView[nTracks][3]/I");
  Reco_Tree->Branch("RecoTrackKalmanLastPlaneBarView", RecoTrackKalmanLastPlaneBarView, "RecoTrackKalmanLastPlaneBarView[nTracks][3]/I");
  Reco_Tree->Branch("RecoTrackKalmanPlaneBarView", RecoTrackKalmanPlaneBarView, "RecoTrackKalmanPlaneBarView[nTracks][200][3]/I");
  Reco_Tree->Branch("KalmanTruePos",  RecoTrackKalmanTruePos,   "TrackHitTruePos[nTracks][200][3]/F");
  Reco_Tree->Branch("RecoTrackKalmanFirstPlaneBarViewTrue", RecoTrackKalmanFirstPlaneBarViewTrue, "RecoTrackKalmanFirstPlaneBarViewTrue[nTracks][3]/I");
  Reco_Tree->Branch("RecoTrackKalmanLastPlaneBarViewTrue", RecoTrackKalmanLastPlaneBarViewTrue, "RecoTrackKalmanLastPlaneBarViewTrue[nTracks][3]/I");
  Reco_Tree->Branch("RecoTrackKalmanPlaneBarViewTrue", RecoTrackKalmanPlaneBarViewTrue, "RecoTrackKalmanPlaneBarViewTrue[nTracks][200][3]/I");
  Reco_Tree->Branch("StartDirection", RecoTrackStartDirection,  "StartDirection[nTracks][3]/F");
  Reco_Tree->Branch("EndDirection",   RecoTrackEndDirection,    "EndDirection[nTracks][3]/F");
  Reco_Tree->Branch("StartPos",       RecoTrackStartPos,        "StartPos[nTracks][3]/F");
  Reco_Tree->Branch("EndPos",         RecoTrackEndPos,          "EndPos[nTracks][3]/F");
  Reco_Tree->Branch("EnergyRange",    RecoTrackEnergyRange,     "EnergyRange[nTracks]/F");
  Reco_Tree->Branch("EnergyDeposit",  RecoTrackEnergyDeposit,   "EnergyDeposit[nTracks]/F");
  Reco_Tree->Branch("Momentum",       RecoTrackMomentum,        "Momentum[nTracks]/F");
  Reco_Tree->Branch("Length",         RecoTrackLength,          "Length[nTracks]/F");
  Reco_Tree->Branch("Length_3D",      RecoTrackLength_3D,       "Length_3D[nTracks]/F");
  Reco_Tree->Branch("Charge",         RecoTrackCharge,          "Charge[nTracks]/I");

  Reco_Tree->Branch("Charge_Kalman",         RecoTrackCharge_Kalman,          "Charge_Kalman[nTracks]/I");
  
  Reco_Tree->Branch("Chi2",           RecoTrackChi2,             "Chi2[nTracks]/F");
  Reco_Tree->Branch("Chi2_plus",           RecoTrackChi2_plus,             "Chi2_plus[nTracks]/F");
  Reco_Tree->Branch("Chi2_minus",           RecoTrackChi2_minus,             "Chi2_minus[nTracks]/F");

  Reco_Tree->Branch("TrackHitEnergies", RecoTrackHitEnergies,   "TrackHitEnergies[nTracks][200]/F");
  Reco_Tree->Branch("TrackHitBarType", RecoTrackHitBarType,   "RecoTrackHitBarType[nTracks][200]/I");
  
  Reco_Tree->Branch("TimeSliceStartTime", &TimeSliceStartTime, "TimeSliceStartTime/F");
  Reco_Tree->Branch("TimeSliceEndTime",   &TimeSliceEndTime,   "TimeSliceEndTime/F");

  MakeTruthBranches(Truth_Info);
  MakeTruthBranches(Truth_Spill);
  
  // Truth information
  Truth_Info->Branch("nParticles", &nParticles, "nParticles/I");
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

  Truth_Info->Branch("LArOuterShellEnergy", &LArOuterShellEnergy, "LArOuterShellEnergy/F");
  Truth_Info->Branch("LArOuterShellEnergyFromVertex", &LArOuterShellEnergyFromVertex, "LArOuterShellEnergyFromVertex/F");
  Truth_Info->Branch("LArTotalEnergy", &LArTotalEnergy, "LArTotalEnergy/F");
  Truth_Info->Branch("LArTotalEnergyFromVertex", &LArTotalEnergyFromVertex, "LArTotalEnergyFromVertex/F");
  Truth_Info->Branch("TotalNonTMSEnergy", &TotalNonTMSEnergy, "TotalNonTMSEnergy/F");
  Truth_Info->Branch("TotalNonTMSEnergyFromVertex", &TotalNonTMSEnergyFromVertex, "TotalNonTMSEnergyFromVertex/F");
  
  Truth_Info->Branch("RecoTrackN", &RecoTrackN, "RecoTrackN/I");
  Truth_Info->Branch("RecoTrackTrueVisibleEnergy", RecoTrackTrueVisibleEnergy,
                     "RecoTrackTrueVisibleEnergy[RecoTrackN]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleIndex", RecoTrackPrimaryParticleIndex, "RecoTrackPrimaryParticleIndex[RecoTrackN]/I");
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueVisibleEnergy", RecoTrackPrimaryParticleTrueVisibleEnergy,
                     "RecoTrackPrimaryParticleTrueVisibleEnergy[RecoTrackN]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueNHits", RecoTrackPrimaryParticleTrueNHits,
                     "RecoTrackPrimaryParticleTrueNHits[RecoTrackN]/I");
  Truth_Info->Branch("RecoTrackSecondaryParticleIndex", RecoTrackSecondaryParticleIndex, "RecoTrackSecondaryParticleIndex[RecoTrackN]/I");
  Truth_Info->Branch("RecoTrackSecondaryParticleTrueVisibleEnergy", RecoTrackSecondaryParticleTrueVisibleEnergy,
                     "RecoTrackSecondaryParticleTrueVisibleEnergy[RecoTrackN]/F");
  Truth_Info->Branch("RecoTrackSecondaryParticleTrueNHits", RecoTrackSecondaryParticleTrueNHits,
                     "RecoTrackSecondaryParticleTrueNHits[RecoTrackN]/I");
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueMomentumTrackStart", RecoTrackPrimaryParticleTrueMomentumTrackStart,
                     "RecoTrackPrimaryParticleTrueMomentumTrackStart[RecoTrackN][4]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleTruePositionTrackStart", RecoTrackPrimaryParticleTruePositionTrackStart,
                     "RecoTrackPrimaryParticleTruePositionTrackStart[RecoTrackN][4]/F"); 
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueMomentumTrackEnd", RecoTrackPrimaryParticleTrueMomentumTrackEnd,
                     "RecoTrackPrimaryParticleTrueMomentumTrackEnd[RecoTrackN][4]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleTruePositionTrackEnd", RecoTrackPrimaryParticleTruePositionTrackEnd,
                     "RecoTrackPrimaryParticleTruePositionTrackEnd[RecoTrackN][4]/F"); 

  Truth_Info->Branch("RecoTrackNHits", RecoTrackNHits, "RecoTrackNHits[RecoTrackN]/I");
  Truth_Info->Branch("RecoTrackTrueHitPosition", RecoTrackTrueHitPosition,
                     "RecoTrackTrueHitPosition[RecoTrackN][200][4]/F");

  Truth_Info->Branch("RecoTrackPrimaryParticleTrueTrackLengthAsMeasured",
                      RecoTrackPrimaryParticleTrueTrackLengthAsMeasured,
                     "RecoTrackPrimaryParticleTrueTrackLengthAsMeasured[RecoTrackN]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueTrackLengthAsMeasuredIgnoreY",
                      RecoTrackPrimaryParticleTrueTrackLengthAsMeasuredIgnoreY,
                     "RecoTrackPrimaryParticleTrueTrackLengthAsMeasuredIgnoreY[RecoTrackN]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueTrackLengthRecoStart",
                      RecoTrackPrimaryParticleTrueTrackLengthRecoStart,
                     "RecoTrackPrimaryParticleTrueTrackLengthRecoStart[RecoTrackN]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueTrackLengthRecoStartIgnoreY",
                      RecoTrackPrimaryParticleTrueTrackLengthRecoStartIgnoreY,
                     "RecoTrackPrimaryParticleTrueTrackLengthRecoStartIgnoreY[RecoTrackN]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueTrackLengthInTMSIgnoreY",
                      RecoTrackPrimaryParticleTrueTrackLengthInTMSIgnoreY,
                     "RecoTrackPrimaryParticleTrueTrackLengthInTMSIgnoreY[RecoTrackN]/F");
                     
  Truth_Info->Branch("RecoTrackPrimaryParticlePDG", &RecoTrackPrimaryParticlePDG, "RecoTrackPrimaryParticlePDG[RecoTrackN]/I");
  Truth_Info->Branch("RecoTrackPrimaryParticleIsPrimary", &RecoTrackPrimaryParticleIsPrimary, "RecoTrackPrimaryParticleIsPrimary[RecoTrackN]/O");
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueMomentum", RecoTrackPrimaryParticleTrueMomentum,
                     "RecoTrackPrimaryParticleTrueMomentum[RecoTrackN][4]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleTruePositionStart", RecoTrackPrimaryParticleTruePositionStart,
                     "RecoTrackPrimaryParticleTruePositionStart[RecoTrackN][4]/F"); 
  Truth_Info->Branch("RecoTrackPrimaryParticleTruePositionEnd", RecoTrackPrimaryParticleTruePositionEnd,
                     "RecoTrackPrimaryParticleTruePositionEnd[RecoTrackN][4]/F"); 
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueTrackLength", RecoTrackPrimaryParticleTrueTrackLength,
                     "RecoTrackPrimaryParticleTrueTrackLength[RecoTrackN]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueTrackLengthIgnoreY", RecoTrackPrimaryParticleTrueTrackLengthIgnoreY,
                     "RecoTrackPrimaryParticleTrueTrackLengthIgnoreY[RecoTrackN]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueTrackLengthInTMS", RecoTrackPrimaryParticleTrueTrackLengthInTMS,
                     "RecoTrackPrimaryParticleTrueTrackLengthInTMS[RecoTrackN]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueMomentumEnteringTMS", RecoTrackPrimaryParticleTrueMomentumEnteringTMS,
                     "RecoTrackPrimaryParticleTrueMomentumEnteringTMS[RecoTrackN][4]/F"); 
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueMomentumLeavingTMS", RecoTrackPrimaryParticleTrueMomentumLeavingTMS,
                     "RecoTrackPrimaryParticleTrueMomentumLeavingTMS[RecoTrackN][4]/F"); 
  
  Truth_Info->Branch("RecoTrackPrimaryParticleTruePositionEnteringTMS", RecoTrackPrimaryParticleTruePositionEnteringTMS,
                     "RecoTrackPrimaryParticleTruePositionEnteringTMS[RecoTrackN][4]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleTruePositionLeavingTMS", RecoTrackPrimaryParticleTruePositionLeavingTMS,
                     "RecoTrackPrimaryParticleTruePositionLeavingTMS[RecoTrackN][4]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleTruePositionLeavingLAr", RecoTrackPrimaryParticleTruePositionLeavingLAr,
                     "RecoTrackPrimaryParticleTruePositionLeavingLAr[RecoTrackN][4]/F");
  Truth_Info->Branch("RecoTrackPrimaryParticleTrueMomentumLeavingLAr", RecoTrackPrimaryParticleTrueMomentumLeavingLAr,
                     "RecoTrackPrimaryParticleTrueMomentumLeavingLAr[RecoTrackN][4]/F");

  Truth_Info->Branch("RecoTrackPrimaryParticleTMSFiducialStart", RecoTrackPrimaryParticleTMSFiducialStart,
    "RecoTrackPrimaryParticleTMSFiducialStart[RecoTrackN]/O");
  Truth_Info->Branch("RecoTrackPrimaryParticleTMSFiducialTouch", RecoTrackPrimaryParticleTMSFiducialTouch,
    "RecoTrackPrimaryParticleTMSFiducialTouch[RecoTrackN]/O");
  Truth_Info->Branch("RecoTrackPrimaryParticleTMSFiducialEnd", RecoTrackPrimaryParticleTMSFiducialEnd,
    "RecoTrackPrimaryParticleTMSFiducialEnd[RecoTrackN]/O");
  Truth_Info->Branch("RecoTrackPrimaryParticleLArFiducialStart", RecoTrackPrimaryParticleLArFiducialStart,
    "RecoTrackPrimaryParticleLArFiducialStart[RecoTrackN]/O");
  Truth_Info->Branch("RecoTrackPrimaryParticleLArFiducialTouch", RecoTrackPrimaryParticleLArFiducialTouch,
    "RecoTrackPrimaryParticleLArFiducialTouch[RecoTrackN]/O");
  Truth_Info->Branch("RecoTrackPrimaryParticleLArFiducialEnd", RecoTrackPrimaryParticleLArFiducialEnd,
    "RecoTrackPrimaryParticleLArFiducialEnd[RecoTrackN]/O");

  Truth_Info->Branch("RecoTrackPrimaryParticleVtxId", RecoTrackPrimaryParticleVtxId,
                     "RecoTrackPrimaryParticleVtxId[RecoTrackN]/I");
  Truth_Info->Branch("RecoTrackPrimaryParticleVtxFiducialCut", RecoTrackPrimaryParticleVtxFiducialCut,
    "RecoTrackPrimaryParticleVtxFiducialCut[RecoTrackN]/O");
  Truth_Info->Branch("RecoTrackPrimaryParticleVtxShellEnergyCut", RecoTrackPrimaryParticleVtxShellEnergyCut,
    "RecoTrackPrimaryParticleVtxShellEnergyCut[RecoTrackN]/O");
  Truth_Info->Branch("RecoTrackPrimaryParticleVtxNDPhysicsCut", RecoTrackPrimaryParticleVtxNDPhysicsCut,
    "RecoTrackPrimaryParticleVtxNDPhysicsCut[RecoTrackN]/O");
                     
  Truth_Info->Branch("RecoTrackSecondaryParticlePDG", &RecoTrackSecondaryParticlePDG, "RecoTrackSecondaryParticlePDG[RecoTrackN]/I");
  Truth_Info->Branch("RecoTrackSecondaryParticleIsPrimary", &RecoTrackSecondaryParticleIsPrimary, "RecoTrackSecondaryParticleIsPrimary[RecoTrackN]/O");
  Truth_Info->Branch("RecoTrackSecondaryParticleTrueMomentum", RecoTrackSecondaryParticleTrueMomentum,
                     "RecoTrackSecondaryParticleTrueMomentum[RecoTrackN][4]/F");
  Truth_Info->Branch("RecoTrackSecondaryParticleTruePositionStart", RecoTrackSecondaryParticleTruePositionStart,
                     "RecoTrackSecondaryParticleTruePositionStart[RecoTrackN][4]/F"); 
  Truth_Info->Branch("RecoTrackSecondaryParticleTruePositionEnd", RecoTrackSecondaryParticleTruePositionEnd,
                     "RecoTrackSecondaryParticleTruePositionEnd[RecoTrackN][4]/F"); 
                     
  Truth_Info->Branch("TrueVisibleEnergyInSlice", TrueVisibleEnergyInSlice, "TrueVisibleEnergyInSlice[nTrueParticles]/F");
  Truth_Info->Branch("TrueNHitsInSlice", TrueNHitsInSlice, "TrueNHitsInSlice[nTrueParticles]/I");
  
  
  Truth_Info->Branch("NTrueHits", &NTrueHits);
  Truth_Info->Branch("TrueHitX", &TrueHitX, "TrueHitX[NTrueHits]/F");
  Truth_Info->Branch("TrueHitY", &TrueHitY, "TrueHitY[NTrueHits]/F");
  Truth_Info->Branch("TrueHitZ", &TrueHitZ, "TrueHitZ[NTrueHits]/F");
  Truth_Info->Branch("TrueHitT", &TrueHitT, "TrueHitT[NTrueHits]/F");
  Truth_Info->Branch("TrueHitE", &TrueHitE, "TrueHitE[NTrueHits]/F");
  Truth_Info->Branch("TrueHitPE", &TrueHitPE, "TrueHitPE[NTrueHits]/F");
  Truth_Info->Branch("TrueHitPEAfterFibers", &TrueHitPEAfterFibers, "TrueHitPEAfterFibers[NTrueHits]/F");
  Truth_Info->Branch("TrueHitPEAfterFibersLongPath", &TrueHitPEAfterFibersLongPath, "TrueHitPEAfterFibersLongPath[NTrueHits]/F");
  Truth_Info->Branch("TrueHitPEAfterFibersShortPath", &TrueHitPEAfterFibersShortPath, "TrueHitPEAfterFibersShortPath[NTrueHits]/F");
  Truth_Info->Branch("TrueNTrueParticles", &TrueNTrueParticles, "TrueNTrueParticles[NTrueHits]/I");
  Truth_Info->Branch("TrueLeptonicEnergy", &TrueLeptonicEnergy, "TrueLeptonicEnergy[NTrueHits]/F");
  Truth_Info->Branch("TrueHadronicEnergy", &TrueHadronicEnergy, "TrueHadronicEnergy[NTrueHits]/F");
  
  
  Truth_Info->Branch("TrueRecoHitX", &TrueRecoHitX, "TrueRecoHitX[NTrueHits]/F");
  Truth_Info->Branch("TrueRecoHitY", &TrueRecoHitY, "TrueRecoHitY[NTrueHits]/F");
  Truth_Info->Branch("TrueRecoHitZ", &TrueRecoHitZ, "TrueRecoHitZ[NTrueHits]/F");
  Truth_Info->Branch("TrueRecoHitTrackX", &TrueRecoHitTrackX, "TrueRecoHitTrackX[NTrueHits]/F");
  Truth_Info->Branch("TrueRecoHitTrackY", &TrueRecoHitTrackY, "TrueRecoHitTrackY[NTrueHits]/F");
  Truth_Info->Branch("TrueRecoHitTrackXUncertainty", &TrueRecoHitTrackXUncertainty, "TrueRecoHitTrackXUncertainty[NTrueHits]/F");
  Truth_Info->Branch("TrueRecoHitTrackYUncertainty", &TrueRecoHitTrackYUncertainty, "TrueRecoHitTrackYUncertainty[NTrueHits]/F");
  Truth_Info->Branch("TrueRecoHitNotZ", &TrueRecoHitNotZ, "TrueRecoHitNotZ[NTrueHits]/F");
  Truth_Info->Branch("TrueRecoHitT", &TrueRecoHitT, "TrueRecoHitT[NTrueHits]/F");
  Truth_Info->Branch("TrueRecoHitE", &TrueRecoHitE, "TrueRecoHitE[NTrueHits]/F");
  Truth_Info->Branch("TrueRecoHitPE", &TrueRecoHitPE, "TrueRecoHitPE[NTrueHits]/F");
  Truth_Info->Branch("TrueRecoHitEVis", &TrueRecoHitEVis, "TrueRecoHitEVis[NTrueHits]/F");
  Truth_Info->Branch("TrueRecoHitIsPedSupped", &TrueRecoHitIsPedSupped, "TrueRecoHitIsPedSupped[NTrueHits]/O");
  Truth_Info->Branch("TrueHitBar", &TrueHitBar, "TrueHitBar[NTrueHits]/I");
  Truth_Info->Branch("TrueHitView", &TrueHitView, "TrueHitView[NTrueHits]/I");
  Truth_Info->Branch("TrueHitPlane", &TrueHitPlane, "TrueHitPlane[NTrueHits]/I");
}

void TMS_TreeWriter::MakeTruthBranches(TTree* truth) {
  // Truth information saved per spill only
  truth->Branch("EventNo", &EventNo, "EventNo/I");
  truth->Branch("SpillNo", &SpillNo, "SpillNo/I");
  truth->Branch("RunNo", &RunNo, "RunNo/I");
  truth->Branch("IsCC", &IsCC, "IsCC/O");
  truth->Branch("Interaction", &Reaction);
  truth->Branch("TruthInfoIndex", &TruthInfoIndex, "TruthInfoIndex/I");
  truth->Branch("TruthInfoNSlices", &TruthInfoNSlices, "TruthInfoNSlices/I");
  truth->Branch("nPrimaryVertices", &nPrimaryVertices, "nPrimaryVertices/I");
  truth->Branch("HasPileup", &HasPileup, "HasPileup/O");
  truth->Branch("NeutrinoPDG", &NeutrinoPDG, "NeutrinoPDG/I");
  truth->Branch("NeutrinoP4", NeutrinoP4, "NeutrinoP4[4]/F");
  truth->Branch("NeutrinoX4", NeutrinoX4, "NeutrinoX4[4]/F");
  
  /*
  truth->Branch("TrueNonTMSNHits", &TrueNonTMSNHits, "TrueNonTMSNHits/I");
  truth->Branch("TrueNonTMSHitPos", TrueNonTMSHitPos, "TrueNonTMSHitPos[TrueNonTMSNHits][4]/F");
  truth->Branch("TrueNonTMSHitEnergy", TrueNonTMSHitEnergy, "TrueNonTMSHitEnergy[TrueNonTMSNHits]/F");
  truth->Branch("TrueNonTMSHitHadronicEnergy", TrueNonTMSHitHadronicEnergy, "TrueNonTMSHitHadronicEnergy[TrueNonTMSNHits]/F");
  truth->Branch("TrueNonTMSHitDx", TrueNonTMSHitDx, "TrueNonTMSHitDx[TrueNonTMSNHits]/F");
  truth->Branch("TrueNonTMSHitdEdx", TrueNonTMSHitdEdx, "TrueNonTMSHitdEdx[TrueNonTMSNHits]/F");
  truth->Branch("TrueNonTMSHitVertexID", TrueNonTMSHitVertexID, "TrueNonTMSHitVertexID[TrueNonTMSNHits]/I");
  */

  truth->Branch("nTrueParticles", &nTrueParticles, "nTrueParticles/I");
  truth->Branch("nTruePrimaryParticles", &nTruePrimaryParticles, "nTruePrimaryParticles/I");
  truth->Branch("nTrueForgottenParticles", &nTrueForgottenParticles, "nTrueForgottenParticles/I");
  truth->Branch("VertexID", VertexID, "VertexID[nTrueParticles]/I");
  truth->Branch("Parent", Parent, "Parent[nTrueParticles]/I");
  truth->Branch("TrackId", TrackId, "TrackId[nTrueParticles]/I");
  truth->Branch("PDG", PDG, "PDG[nTrueParticles]/I");
  truth->Branch("IsPrimary", IsPrimary, "IsPrimary[nTrueParticles]/O");
  truth->Branch("TrueVisibleEnergy", TrueVisibleEnergy, "TrueVisibleEnergy[nTrueParticles]/F");
  truth->Branch("TrueNHits", TrueNHits, "TrueNHits[nTrueParticles]/I");
  truth->Branch("TruePathLength", TruePathLength, "TruePathLength[nTrueParticles]/F");
  truth->Branch("TruePathLengthIgnoreY", TruePathLengthIgnoreY, "TruePathLengthIgnoreY[nTrueParticles]/F");
  truth->Branch("TruePathLengthInTMS", TruePathLengthInTMS, "TruePathLengthInTMS[nTrueParticles]/F");
  truth->Branch("TruePathLengthInTMSIgnoreY", TruePathLengthInTMSIgnoreY, "TruePathLengthInTMSIgnoreY[nTrueParticles]/F");
  
  truth->Branch("InteractionTMSFiducial", &InteractionTMSFiducial, "InteractionTMSFiducial/O");
  truth->Branch("InteractionTMSFirstTwoModules", &InteractionTMSFirstTwoModules, "InteractionTMSFirstTwoModules/O");
  truth->Branch("InteractionTMSThin", &InteractionTMSThin, "InteractionTMSThin/O");
  truth->Branch("InteractionLArFiducial", &InteractionLArFiducial, "InteractionLArFiducial/O");
  truth->Branch("TMSFiducialStart", TMSFiducialStart, "TMSFiducialStart[nTrueParticles]/O");
  truth->Branch("TMSFiducialTouch", TMSFiducialTouch, "TMSFiducialTouch[nTrueParticles]/O");
  truth->Branch("TMSFiducialEnd", TMSFiducialEnd, "TMSFiducialEnd[nTrueParticles]/O");
  truth->Branch("LArFiducialStart", LArFiducialStart, "LArFiducialStart[nTrueParticles]/O");
  truth->Branch("LArFiducialTouch", LArFiducialTouch, "LArFiducialTouch[nTrueParticles]/O");
  truth->Branch("LArFiducialEnd", LArFiducialEnd, "LArFiducialEnd[nTrueParticles]/O");
  
  truth->Branch("BirthMomentum", BirthMomentum, "BirthMomentum[nTrueParticles][4]/F");
  truth->Branch("BirthPosition", BirthPosition, "BirthPosition[nTrueParticles][4]/F");
  truth->Branch("DeathMomentum", DeathMomentum, "DeathMomentum[nTrueParticles][4]/F");
  truth->Branch("DeathPosition", DeathPosition, "DeathPosition[nTrueParticles][4]/F");
  
  // IsInside-based start/end
  truth->Branch("MomentumLArStart", MomentumLArStart, "MomentumLArStart[nTrueParticles][4]/F");
  truth->Branch("PositionLArStart", PositionLArStart, "PositionLArStart[nTrueParticles][4]/F");
  truth->Branch("MomentumLArEnd", MomentumLArEnd, "MomentumLArEnd[nTrueParticles][4]/F");
  truth->Branch("PositionLArEnd", PositionLArEnd, "PositionLArEnd[nTrueParticles][4]/F");
  truth->Branch("MomentumTMSStart", MomentumTMSStart, "MomentumTMSStart[nTrueParticles][4]/F");
  truth->Branch("PositionTMSStart", PositionTMSStart, "PositionTMSStart[nTrueParticles][4]/F");
  truth->Branch("MomentumTMSFirstTwoModulesEnd", MomentumTMSFirstTwoModulesEnd, "MomentumTMSFirstTwoModulesEnd[nTrueParticles][4]/F");
  truth->Branch("PositionTMSFirstTwoModulesEnd", PositionTMSFirstTwoModulesEnd, "PositionTMSFirstTwoModulesEnd[nTrueParticles][4]/F"); 
  truth->Branch("MomentumTMSThinEnd", MomentumTMSThinEnd, "MomentumTMSThinEnd[nTrueParticles][4]/F");
  truth->Branch("PositionTMSThinEnd", PositionTMSThinEnd, "PositionTMSThinEnd[nTrueParticles][4]/F"); 
  truth->Branch("MomentumTMSEnd", MomentumTMSEnd, "MomentumTMSEnd[nTrueParticles][4]/F");
  truth->Branch("PositionTMSEnd", PositionTMSEnd, "PositionTMSEnd[nTrueParticles][4]/F"); 
  
  // Z-based start/end
  truth->Branch("MomentumZIsLArEnd", MomentumZIsLArEnd, "MomentumZIsLArEnd[nTrueParticles][4]/F");
  truth->Branch("PositionZIsLArEnd", PositionZIsLArEnd, "PositionZIsLArEnd[nTrueParticles][4]/F");
  truth->Branch("MomentumZIsTMSStart", MomentumZIsTMSStart, "MomentumZIsTMSStart[nTrueParticles][4]/F");
  truth->Branch("PositionZIsTMSStart", PositionZIsTMSStart, "PositionZIsTMSStart[nTrueParticles][4]/F");
  truth->Branch("MomentumZIsTMSEnd", MomentumZIsTMSEnd, "MomentumZIsTMSEnd[nTrueParticles][4]/F");
  truth->Branch("PositionZIsTMSEnd", PositionZIsTMSEnd, "PositionZIsTMSEnd[nTrueParticles][4]/F");
  
  truth->Branch("TrueVtxN", &TrueVtxN, "TrueVtxN/I");
  // Position
  truth->Branch("TrueVtxX", TrueVtxX, "TrueVtxX[TrueVtxN]/F");
  truth->Branch("TrueVtxY", TrueVtxY, "TrueVtxY[TrueVtxN]/F");
  truth->Branch("TrueVtxZ", TrueVtxZ, "TrueVtxZ[TrueVtxN]/F");
  truth->Branch("TrueVtxT", TrueVtxT, "TrueVtxT[TrueVtxN]/F");
  // Vtx E
  truth->Branch("TrueVtxPx", TrueVtxPx, "TrueVtxPx[TrueVtxN]/F");
  truth->Branch("TrueVtxPy", TrueVtxPy, "TrueVtxPy[TrueVtxN]/F");
  truth->Branch("TrueVtxPz", TrueVtxPz, "TrueVtxPz[TrueVtxN]/F");
  truth->Branch("TrueVtxE", TrueVtxE, "TrueVtxE[TrueVtxN]/F");
  // Other info
  truth->Branch("TrueVtxPDG", TrueVtxPDG, "TrueVtxPDG[TrueVtxN]/I");
  truth->Branch("TrueVtxID", TrueVtxID, "TrueVtxID[TrueVtxN]/I");
  truth->Branch("TrueVtxReaction", &TrueVtxReaction);
  // Hadronic E
  truth->Branch("TrueVtxHadronicELarShell", TrueVtxHadronicELarShell, "TrueVtxHadronicELarShell[TrueVtxN]/F");
  truth->Branch("TrueVtxHadronicELAr", TrueVtxHadronicELAr, "TrueVtxHadronicELAr[TrueVtxN]/F");
  truth->Branch("TrueVtxHadronicETMS", TrueVtxHadronicETMS, "TrueVtxHadronicETMS[TrueVtxN]/F");
  truth->Branch("TrueVtxHadronicE", TrueVtxHadronicE, "TrueVtxHadronicE[TrueVtxN]/F");
  // Visible E
  truth->Branch("TrueVtxVisibleETMS", TrueVtxVisibleETMS, "TrueVtxVisibleETMS[TrueVtxN]/F");
  truth->Branch("TrueVtxVisibleELAr", TrueVtxVisibleELAr, "TrueVtxVisibleELAr[TrueVtxN]/F");
  truth->Branch("TrueVtxVisibleE", TrueVtxVisibleE, "TrueVtxVisibleE[TrueVtxN]/F");
  // Truth cuts
  truth->Branch("TrueVtxFiducialCut", TrueVtxFiducialCut, "TrueVtxFiducialCut[TrueVtxN]/O");
  truth->Branch("TrueVtxShellEnergyCut", TrueVtxShellEnergyCut, "TrueVtxShellEnergyCut[TrueVtxN]/O");
  truth->Branch("TrueVtxNDPhysicsCut", TrueVtxNDPhysicsCut, "TrueVtxNDPhysicsCut[TrueVtxN]/O");
}

static void setMomentum(float *branch, TVector3 momentum, double energy = -9999) {
    branch[0] = momentum.Px();
    branch[1] = momentum.Py();
    branch[2] = momentum.Pz();
    branch[3] = energy;
}

static void setMomentum(float *branch, TLorentzVector momentum) {
    branch[0] = momentum.Px();
    branch[1] = momentum.Py();
    branch[2] = momentum.Pz();
    branch[3] = momentum.E();
}

static void setPosition(float *branch, TLorentzVector position) {
    branch[0] = position.X();
    branch[1] = position.Y();
    branch[2] = position.Z();
    // We're saving floats, so remove giant offset or else we'll have trouble
    branch[3] = std::fmod(position.T(), TMS_Manager::GetInstance().Get_Nersc_Spill_Period());
}

void TMS_TreeWriter::Fill(TMS_Event &event) {
  // Clear old info
  Clear();

  // See if track is exiting or not
  int nLastHits = TMS_Manager::GetInstance().Get_Reco_STOPPING_nLastHits();
  double EnergyCut = TMS_Manager::GetInstance().Get_Reco_STOPPING_EnergyCut();
  
  FillTruthInfo(event);


  // Fill the truth info
  EventNo = event.GetEventNumber();
  SliceNo = event.GetSliceNumber();
  SpillNo = event.GetSpillNumber();
  RunNo = event.GetRunNumber();
  Reaction = event.GetReaction();

  NeutrinoPDG = event.GetNeutrinoPDG();
  NeutrinoP4[0] = event.GetNeutrinoP4().X();
  NeutrinoP4[1] = event.GetNeutrinoP4().Y();
  NeutrinoP4[2] = event.GetNeutrinoP4().Z();
  NeutrinoP4[3] = event.GetNeutrinoP4().T();
  NeutrinoX4[0] = event.GetNeutrinoVtx().X();
  NeutrinoX4[1] = event.GetNeutrinoVtx().Y();
  NeutrinoX4[2] = event.GetNeutrinoVtx().Z();
  NeutrinoX4[3] = event.GetNeutrinoVtx().T();
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

  Muon_TrueTrackLength= event.GetMuonTrueTrackLength();
  //Muon_TrueTrackLength = -999.99;
  //std::cout << Muon_TrueTrackLength << std::endl;
  Muon_TrueKE = event.GetMuonTrueKE();
  
  // Fill LAr hit outer shell energy info
  // Case 1: All energy in outer shell, useful only for single event interactions
  // Case 2: All energy from primary vertex, useful for pileup. We're assuming reco can distinguish
  double thickness = TMS_Manager::GetInstance().Get_LAR_OUTER_SHELL_THICKNESS(); // mm
  LArOuterShellEnergy = event.CalculateEnergyInLArOuterShell(thickness);
  LArOuterShellEnergyFromVertex = event.CalculateEnergyInLArOuterShell(thickness, VertexIdOfMostEnergyInEvent);
  LArTotalEnergy = event.CalculateEnergyInLAr();
  LArTotalEnergyFromVertex = event.CalculateEnergyInLAr(VertexIdOfMostEnergyInEvent);
  TotalNonTMSEnergy = event.CalculateTotalNonTMSEnergy();
  TotalNonTMSEnergyFromVertex = event.CalculateTotalNonTMSEnergy(VertexIdOfMostEnergyInEvent);

  // Fill the reco info
  std::vector<std::pair<bool, TF1*>> HoughLinesU = TMS_TrackFinder::GetFinder().GetHoughLinesU();
  std::vector<std::pair<bool, TF1*>> HoughLinesV = TMS_TrackFinder::GetFinder().GetHoughLinesV();
  std::vector<std::pair<bool, TF1*>> HoughLinesX = TMS_TrackFinder::GetFinder().GetHoughLinesX();
  std::vector<std::pair<bool, TF1*>> HoughLinesY = TMS_TrackFinder::GetFinder().GetHoughLinesY();
  // Also get the size of the hits to get a measure of relative goodness
  std::vector<std::vector<TMS_Hit> > HoughCandidatesU = TMS_TrackFinder::GetFinder().GetHoughCandidatesU();
  std::vector<std::vector<TMS_Hit> > HoughCandidatesV = TMS_TrackFinder::GetFinder().GetHoughCandidatesV();
  std::vector<std::vector<TMS_Hit> > HoughCandidatesX = TMS_TrackFinder::GetFinder().GetHoughCandidatesX();
  std::vector<std::vector<TMS_Hit> > HoughCandidatesY = TMS_TrackFinder::GetFinder().GetHoughCandidatesY();
  nLinesU = HoughCandidatesU.size();
  nLinesV = HoughCandidatesV.size();
  nLinesX = HoughCandidatesX.size();
  nLinesY = HoughCandidatesY.size();


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

  if (nLinesY > __TMS_MAX_LINES__) {
    std::cerr << "Exceeded max number of HoughLines to write to file" << std::endl;
    std::cerr << "Max lines: " << __TMS_MAX_LINES__ << std::endl;
    std::cerr << "Number of lines in event (x): " << nLinesY << std::endl;
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

  it = 0;
  for (auto &Lines: HoughLinesY) {
    // Get the slopes saved down
    InterceptY[it] = Lines.second->GetParameter(0);
    Intercept_UpstreamY[it] = TMS_TrackFinder::GetFinder().GetHoughLinesY_Upstream()[it].first;
    Intercept_DownstreamY[it] = TMS_TrackFinder::GetFinder().GetHoughLinesY_Downstream()[it].first;
    SlopeY[it] = Lines.second->GetParameter(1);
    Slope_UpstreamY[it] = TMS_TrackFinder::GetFinder().GetHoughLinesY_Upstream()[it].second;
    Slope_DownstreamY[it] = TMS_TrackFinder::GetFinder().GetHoughLinesY_Downstream()[it].second;

    // Calculate the z and x vectors by evaling the TF1 in thin and thick target
    double xlow = TMS_Const::TMS_Thin_Start;
    double xhi = TMS_Const::TMS_Thick_Start;
    double ylow = Lines.second->Eval(xlow);
    double yhi = Lines.second->Eval(xhi);
    // Do the same for up and downstream portions
    double ylow_upstream = Intercept_UpstreamY[it]+xlow*Slope_UpstreamY[it];
    double yhi_upstream = Intercept_UpstreamY[it]+xhi*Slope_UpstreamY[it];

    double ylow_downstream = Intercept_DownstreamY[it]+xlow*Slope_DownstreamY[it];
    double yhi_downstream = Intercept_DownstreamY[it]+xhi*Slope_DownstreamY[it];

    double xlen = xhi-xlow;
    double ylen = yhi-ylow;
    double len = sqrt(xlen*xlen+ylen*ylen);

    double ylen_upstream = yhi_upstream - ylow_upstream;
    double ylen_downstream = yhi_downstream - ylow_downstream;

    double len_upstream = sqrt(xlen*xlen+ylen_upstream*ylen_upstream);
    double len_downstream = sqrt(xlen*xlen+ylen_downstream*ylen_downstream);

    DirectionZY[it] = xlen/len;
    DirectionXY[it] = ylen/len;

    DirectionZY_Upstream[it] = xlen/len_upstream;
    DirectionXY_Upstream[it] = ylen_upstream/len_upstream;

    DirectionZY_Downstream[it] = xlen/len_downstream;
    DirectionXY_Downstream[it] = ylen_downstream/len_downstream;

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

  std::vector<std::vector<TMS_Hit> > HoughCandsY = TMS_TrackFinder::GetFinder().GetHoughCandidatesY();
  it = 0;
  for (auto &Candidates: HoughCandsY) {
    // Loop over hits
    for (auto &hit: Candidates) {
      if (FirstTrack == NULL) {
        FirstTrack = &hit;
      } else if (hit.GetZ() < FirstTrack->GetZ()) {
        FirstTrack = &hit;
      }
    }
    nHitsInTrackY[it] = Candidates.size();

    // Then save the hit info
    FirstPlaneY[it] = Candidates.front().GetPlaneNumber();
    FirstHitY[it][0] = Candidates.front().GetZ();
    FirstHitY[it][1] = Candidates.front().GetNotZ();
    FirstHitTimeY[it] = Candidates.front().GetT();

    LastPlaneY[it] = Candidates.back().GetPlaneNumber();
    LastHitY[it][0] = Candidates.back().GetZ();
    LastHitY[it][1] = Candidates.back().GetNotZ();
    LastHitTimeY[it] = Candidates.back().GetT();

    TrackLengthY[it] = TMS_TrackFinder::GetFinder().GetTrackLengthY()[it];
    TotalTrackEnergyY[it] = TMS_TrackFinder::GetFinder().GetTrackEnergyY()[it];
    OccupancyY[it] = double(HoughCandsY[it].size())/TotalHits;
    
    float earliest_hit_time = 1e32;
    float latest_hit_time = -1e32;
    // Get each hit in the track and save its energy
    for (unsigned int j = 0; j < Candidates.size(); ++j) {
      TrackHitEnergyY[it][j] = Candidates[j].GetE();
      TrackHitTimeY[it][j] = Candidates[j].GetT();
      TrackHitPosY[it][j][0] = Candidates[j].GetZ();
      TrackHitPosY[it][j][1] = Candidates[j].GetNotZ();
      
      float time = Candidates[j].GetT();
      
      if (time < earliest_hit_time) earliest_hit_time = time;
      if (time > latest_hit_time) latest_hit_time = time;
    }
    EarliestHitTimeY[it] = earliest_hit_time;
    LatestHitTimeY[it] = latest_hit_time;
    

    double maxenergy = 0;
    unsigned int nLastHits_temp = nLastHits;
    if ((unsigned int)nLastHits > Candidates.size()) nLastHits_temp = Candidates.size();
    for (unsigned int i = 0; i < nLastHits_temp; ++i) {
      double hitenergy = Candidates[nLastHits_temp-1-i].GetE();
      if (hitenergy > maxenergy) maxenergy = hitenergy;
    }
    if (maxenergy > EnergyCut) TrackStoppingY[it] = true;


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
  std::vector<std::vector<TMS_Hit> > ClustersY = TMS_TrackFinder::GetFinder().GetClusterCandidatesY();
  nClustersU = ClustersU.size();
  nClustersV = ClustersV.size();
  nClustersX = ClustersX.size();
  nClustersY = ClustersY.size();
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
  if (nClustersY > __TMS_MAX_CLUSTERS__) {
    std::cerr << "Too many clusters in TMS_TreeWriter" << std::endl;
    std::cerr << "Hard-coded maximum: " << __TMS_MAX_CLUSTERS__ << std::endl;
    std::cerr << "nClustersX in event: " << nClustersY << std::endl;
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
  stdit = 0;
  for (auto it = ClustersY.begin(); it != ClustersY.end(); ++it, ++stdit) {
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
      ClusterHitPosY[stdit][j][0] = (*it)[j].GetZ();
      ClusterHitPosY[stdit][j][1] = (*it)[j].GetNotZ();
      ClusterHitEnergyY[stdit][j] = (*it)[j].GetE();
      ClusterHitTimeY[stdit][j] = (*it)[j].GetT();
      ClusterHitSliceY[stdit][j] = (*it)[j].GetSlice();
    }
    mean_z /= nhits;
    mean_notz /= nhits;

    mean2_z /= nhits;
    mean2_notz /= nhits;

    ClusterEnergyY[stdit] = total_energy;
    ClusterTimeY[stdit] = min_cluster_time;
    nHitsInClusterY[stdit] = nhits;
    ClusterPosMeanY[stdit][0] = mean_z;
    ClusterPosMeanY[stdit][1] = mean_notz;
    // Calculate the standard deviation
    double std_dev_z = mean2_z-mean_z*mean_z;
    double std_dev_notz = mean2_z-mean_z*mean_z;
    if (std_dev_z > 0) std_dev_z = sqrt(std_dev_z);
    if (std_dev_notz > 0) std_dev_notz = sqrt(std_dev_notz);
    ClusterPosStdDevY[stdit][0] = std_dev_z;
    ClusterPosStdDevY[stdit][1] = std_dev_notz;
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
    RecoHitBarType[stdit] = (*it).GetBar().GetBarType();
    RecoHitPE[stdit] = (*it).GetPE();
    RecoHitBar[stdit] = (*it).GetBarNumber();
    RecoHitPlane[stdit] = (*it).GetPlaneNumber();
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
  
  TimeSliceStartTime = event.GetTimeSliceBounds().first;
  TimeSliceEndTime = event.GetTimeSliceBounds().second;

  for (auto RecoTrack = Reco_Tracks.begin(); RecoTrack != Reco_Tracks.end(); ++RecoTrack, ++itTrack) {
    nHitsIn3DTrack[itTrack]         =  RecoTrack->Hits.size(); // Do we need to cast it? idk
    nKalmanNodes[itTrack]           =  RecoTrack->KalmanNodes.size();
    nKalmanNodes_plus[itTrack]           =  RecoTrack->KalmanNodes_plus.size();
    nKalmanNodes_minus[itTrack]           =  RecoTrack->KalmanNodes_minus.size();
    KalmanErrorDetVol[itTrack]      =       RecoTrack->KalmanErrorDetVol;

//    std::cout << "TreeWriter number of hits: " << nHitsIn3DTrack[itTrack] << std::endl;
    RecoTrackEnergyRange[itTrack]   =       RecoTrack->EnergyRange;
    if (nLinesU<=0){
        RecoTrackLength[itTrack]        =       0.5 * (TrackLengthX[itTrack] + TrackLengthY[itTrack]); //RecoTrack->Length;// RecoTrack->Length;, 2d is better estimate than 3d because of y jumps
    }
    else {
        RecoTrackLength[itTrack]        =       0.5 * (TrackLengthU[itTrack] + TrackLengthV[itTrack]); //RecoTrack->Length;// RecoTrack->Length;, 2d is better estimate than 3d because of y jumps
    }

    RecoTrackLength_3D[itTrack]        =    RecoTrack->Length; 
    RecoTrackEnergyDeposit[itTrack] =       RecoTrack->EnergyDeposit;
    RecoTrackMomentum[itTrack]      =       RecoTrack->Momentum;
    RecoTrackCharge[itTrack]        =       RecoTrack->Charge;
    RecoTrackCharge_Kalman[itTrack]        =       RecoTrack->Charge_Kalman;
    RecoTrackChi2[itTrack]          =       RecoTrack->Chi2;
    RecoTrackChi2_plus[itTrack]          =       RecoTrack->Chi2_plus;
    RecoTrackChi2_minus[itTrack]          =       RecoTrack->Chi2_minus;
    

    for (int j = 0; j < 3; j++) {
      RecoTrackStartPos[itTrack][j]  = RecoTrack->Start[j];
      RecoTrackEndPos[itTrack][j]    = RecoTrack->End[j];
      RecoTrackStartDirection[itTrack][j] = RecoTrack->StartDirection[j];
      RecoTrackEndDirection[itTrack][j] = RecoTrack->EndDirection[j];
    }
    
    if (RecoTrack->KalmanNodes.size() > 0) {
      size_t last_index = RecoTrack->KalmanNodes.size() - 1;
      TMS_Bar first_bar(RecoTrack->KalmanNodes[0].RecoX, RecoTrack->KalmanNodes[0].RecoY,
                        RecoTrack->KalmanNodes[0].z);
      TMS_Bar last_bar(RecoTrack->KalmanNodes[last_index].RecoX, RecoTrack->KalmanNodes[last_index].RecoY,
                       RecoTrack->KalmanNodes[last_index].z);
      RecoTrackKalmanFirstPlaneBarView[itTrack][0] = first_bar.GetPlaneNumber();
      RecoTrackKalmanFirstPlaneBarView[itTrack][1] = first_bar.GetBarNumber();
      RecoTrackKalmanFirstPlaneBarView[itTrack][2] = first_bar.GetBarTypeNumber();

      RecoTrackKalmanLastPlaneBarView[itTrack][0] = last_bar.GetPlaneNumber();
      RecoTrackKalmanLastPlaneBarView[itTrack][1] = last_bar.GetBarNumber();
      RecoTrackKalmanLastPlaneBarView[itTrack][2] = last_bar.GetBarTypeNumber();

      TMS_Bar first_bar_true(RecoTrack->KalmanNodes[0].TrueX, RecoTrack->KalmanNodes[0].TrueY,
                             RecoTrack->KalmanNodes[0].z);
      TMS_Bar last_bar_true(RecoTrack->KalmanNodes[last_index].TrueX, RecoTrack->KalmanNodes[last_index].TrueY,
                            RecoTrack->KalmanNodes[last_index].z);
      RecoTrackKalmanFirstPlaneBarViewTrue[itTrack][0] = first_bar_true.GetPlaneNumber();
      RecoTrackKalmanFirstPlaneBarViewTrue[itTrack][1] = first_bar_true.GetBarNumber();
      RecoTrackKalmanFirstPlaneBarViewTrue[itTrack][2] = first_bar_true.GetBarTypeNumber();

      RecoTrackKalmanLastPlaneBarViewTrue[itTrack][0] = last_bar_true.GetPlaneNumber();
      RecoTrackKalmanLastPlaneBarViewTrue[itTrack][1] = last_bar_true.GetBarNumber();
      RecoTrackKalmanLastPlaneBarViewTrue[itTrack][2] = last_bar_true.GetBarTypeNumber();
    }

    for (unsigned int j = 0; j < RecoTrack->KalmanNodes.size(); ++j) {
      RecoTrackKalmanPos[itTrack][j][0] = RecoTrack->KalmanNodes[j].RecoX;
      RecoTrackKalmanPos[itTrack][j][1] = RecoTrack->KalmanNodes[j].RecoY;
      RecoTrackKalmanPos[itTrack][j][2] = RecoTrack->KalmanNodes[j].z;

      RecoTrackKalmanTruePos[itTrack][j][0] = RecoTrack->KalmanNodes[j].TrueX;
      RecoTrackKalmanTruePos[itTrack][j][1] = RecoTrack->KalmanNodes[j].TrueY;
      RecoTrackKalmanTruePos[itTrack][j][2] = RecoTrack->KalmanNodes[j].z;

      //RecoTrackKalmanPos[itTrack][j][0] = RecoTrack->KalmanNodes[j].CurrentState.x;
      //RecoTrackKalmanPos[itTrack][j][1] = RecoTrack->KalmanNodes[j].CurrentState.y;
      //RecoTrackKalmanPos[itTrack][j][2] = RecoTrack->KalmanNodes[j].CurrentState.z;

      TMS_Bar current_bar(RecoTrack->KalmanNodes[j].RecoX, RecoTrack->KalmanNodes[j].RecoY,
                          RecoTrack->KalmanNodes[j].z);
      RecoTrackKalmanPlaneBarView[itTrack][j][0] = current_bar.GetPlaneNumber();
      RecoTrackKalmanPlaneBarView[itTrack][j][1] = current_bar.GetBarNumber();
      RecoTrackKalmanPlaneBarView[itTrack][j][2] = current_bar.GetBarTypeNumber();

      TMS_Bar current_bar_true(RecoTrack->KalmanNodes[j].TrueX, RecoTrack->KalmanNodes[j].TrueY,
                               RecoTrack->KalmanNodes[j].z);
      RecoTrackKalmanPlaneBarViewTrue[itTrack][j][0] = current_bar_true.GetPlaneNumber();
      RecoTrackKalmanPlaneBarViewTrue[itTrack][j][1] = current_bar_true.GetBarNumber();
      RecoTrackKalmanPlaneBarViewTrue[itTrack][j][2] = current_bar_true.GetBarTypeNumber();
    }

    for (unsigned int j = 0; j < RecoTrack->Hits.size(); ++j) {
      RecoTrackHitEnergies[itTrack][j] = RecoTrack->Hits[j].GetE(); // Add the energy deposit from each hit
      RecoTrackHitBarType[itTrack][j] = RecoTrack->Hits[j].GetBar().GetBarType();

      // Here we check for bar orientation
      if (RecoTrack->Hits[j].GetBar().GetBarType() != TMS_Bar::kXBar) {
        RecoTrackHitPos[itTrack][j][0] = RecoTrack->Hits[j].GetRecoX(); // GetNotZ?
        RecoTrackHitPos[itTrack][j][1] = RecoTrack->Hits[j].GetRecoY();
      } else if (RecoTrack->Hits[j].GetBar().GetBarType() == TMS_Bar::kXBar) {
        RecoTrackHitPos[itTrack][j][0] = RecoTrack->Hits[j].GetRecoX();
        RecoTrackHitPos[itTrack][j][1] = RecoTrack->Hits[j].GetNotZ();
      }
      RecoTrackHitPos[itTrack][j][2] = RecoTrack->Hits[j].GetZ();
      
    }
    // Can manually compute direction if it hasn't been set
//    if ( (RecoTrackDirection[itTrack][0] == 0) && (RecoTrackDirection[itTrack][1] == 0) && (RecoTrackDirection[itTrack][2] == 0) )
//    { // If true it seems the direction hasn't been set
//      for (int j = 0; j < 3; j++)
//      { // Right now no need to make sure this is a unit vector
//        RecoTrackDirection[itTrack][j] = RecoTrack->End[j] - RecoTrack->Start[j];
//      }
//    }

    auto TrueParticles = event.GetTrueParticles();
    
    // Now fill truth info
    if (itTrack >= __TMS_MAX_LINES__) {
      std::cout<<"Warning: RecoTrackN < __TMS_MAX_LINES__. If this happens often, increase __TMS_MAX_LINES__"<<std::endl;
      RecoTrackN = __TMS_MAX_LINES__;
      continue;
    }
    double total_true_visible_energy = 0;
    double true_primary_visible_energy = -999;
    double true_secondary_visible_energy = -999;
    int true_primary_particle_index = -999;
    int true_secondary_particle_index = -999;
    auto particle_info = TMS_Utils::GetPrimaryIdsByEnergy(RecoTrack->Hits);
    total_true_visible_energy = particle_info.total_energy;
    if (particle_info.energies.size() > 0) {
      true_primary_visible_energy = particle_info.energies[0];
      true_primary_particle_index = event.GetTrueParticleIndex(particle_info.vertexids[0], particle_info.trackids[0]);
    }
    // Now for the primary index, find the true starting and ending momentum and position
    if (true_primary_particle_index < 0) {
      // Do nothing, this means we didn't find a true particle associated with a reco track
      // This can't happen unless dark noise existed which is currently doesn't
      std::cout<<"Warning: Found true_primary_particle_index < 0. There should be at least one true particle creating energy but instead the index is: "<<true_primary_particle_index<<", with energy: "<<true_primary_visible_energy<<std::endl;
      std::cout<<"This can happen with dark noise or if TMS_Event.OnlyPrimaryOrVisibleEnergy is false"<<std::endl;
    }
    else if ((size_t)true_primary_particle_index >= TrueParticles.size()) {
      // This can happen if TMS_Event.OnlyPrimaryOrVisibleEnergy is false
      std::cout<<"Warning: Found a true_primary_particle_index >= TrueParticles.size() case. "<<true_primary_particle_index<<" >= "<<TrueParticles.size()<<std::endl;
      true_primary_particle_index = -800000000 - true_primary_particle_index;
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
        
        // TODO needs fixing
        //if ( RecoTrack->Hits.size() !=  RecoTrack->nHits) std::cout<<"N hits mismatch: "<< RecoTrack->Hits.size() << " vs "<< RecoTrack->nHits<<std::endl;
        RecoTrackNHits[itTrack] = RecoTrack->Hits.size();
        for (size_t h = 0; h <  RecoTrack->Hits.size() && h < __TMS_MAX_LINE_HITS__; h++) {
          auto& hit = RecoTrack->Hits[h];
          RecoTrackTrueHitPosition[itTrack][h][0] = hit.GetTrueHit().GetX();
          RecoTrackTrueHitPosition[itTrack][h][1] = hit.GetTrueHit().GetY();
          RecoTrackTrueHitPosition[itTrack][h][2] = hit.GetTrueHit().GetZ();
          RecoTrackTrueHitPosition[itTrack][h][3] = hit.GetTrueHit().GetT();
        }

        // Now calulate the true track length from true start to true end
        RecoTrackPrimaryParticleTrueTrackLengthAsMeasured[itTrack] =
            TMS_Geom::GetInstance().GetTrackLength(tp.GetPositionPoints(start_z, end_z));
        RecoTrackPrimaryParticleTrueTrackLengthAsMeasuredIgnoreY[itTrack] =
            TMS_Geom::GetInstance().GetTrackLength(tp.GetPositionPoints(start_z, end_z), true);
        // Again from true start to particle's end
        // Picking 10x TMS_Const::TMS_Thick_End breaks geometry
        // Picking 2x TMS_Const::TMS_Thick_End causes infinite loop
        // Or maybe crazy slowdown in sand?
        // So let's go slightly beyond end of TMS
        const double LARGE_Z = TMS_Const::TMS_Thick_End + 1000;
        const double SMALL_Z = TMS_Const::LAr_Start_Exact[2] - 1000;
        RecoTrackPrimaryParticleTrueTrackLengthRecoStart[itTrack] =
            TMS_Geom::GetInstance().GetTrackLength(tp.GetPositionPoints(start_z, LARGE_Z));
        RecoTrackPrimaryParticleTrueTrackLengthRecoStartIgnoreY[itTrack] =
            TMS_Geom::GetInstance().GetTrackLength(tp.GetPositionPoints(start_z, LARGE_Z), true);
        // Again from true start to particle's end within tms
        RecoTrackPrimaryParticleTrueTrackLengthInTMS[itTrack] =
            TMS_Geom::GetInstance().GetTrackLength(tp.GetPositionPoints(start_z, LARGE_Z, true));
        RecoTrackPrimaryParticleTrueTrackLengthInTMSIgnoreY[itTrack] =
            TMS_Geom::GetInstance().GetTrackLength(tp.GetPositionPoints(start_z, LARGE_Z, true), true);
            
        RecoTrackPrimaryParticleTrueTrackLength[itTrack] =
            TMS_Geom::GetInstance().GetTrackLength(tp.GetPositionPoints(SMALL_Z, LARGE_Z));
        RecoTrackPrimaryParticleTrueTrackLengthIgnoreY[itTrack] =
            TMS_Geom::GetInstance().GetTrackLength(tp.GetPositionPoints(SMALL_Z, LARGE_Z), true);
            
        RecoTrackPrimaryParticlePDG[itTrack] = tp.GetPDG();
        RecoTrackPrimaryParticleIsPrimary[itTrack] = tp.IsPrimary();
        setMomentum(RecoTrackPrimaryParticleTrueMomentum[itTrack], tp.GetBirthMomentum());
        setPosition(RecoTrackPrimaryParticleTruePositionStart[itTrack], tp.GetBirthPosition());
        setMomentum(RecoTrackPrimaryParticleTruePositionEnd[itTrack], tp.GetDeathPosition());
        
        setMomentum(RecoTrackPrimaryParticleTrueMomentumEnteringTMS[itTrack], tp.GetMomentumEnteringTMS());
        setPosition(RecoTrackPrimaryParticleTruePositionEnteringTMS[itTrack], tp.GetPositionEnteringTMS());
        setMomentum(RecoTrackPrimaryParticleTrueMomentumLeavingTMS[itTrack], tp.GetMomentumLeavingTMS());
        setPosition(RecoTrackPrimaryParticleTruePositionLeavingTMS[itTrack], tp.GetPositionLeavingTMS());
        setMomentum(RecoTrackPrimaryParticleTrueMomentumLeavingLAr[itTrack], tp.GetMomentumLeavingLAr());
        setPosition(RecoTrackPrimaryParticleTruePositionLeavingLAr[itTrack], tp.GetPositionLeavingLAr());
        
        TVector3 location_birth = tp.GetBirthPosition().Vect();
        TVector3 location_death = tp.GetDeathPosition().Vect();
        RecoTrackPrimaryParticleTMSFiducialStart[itTrack] = TMS_Geom::GetInstance().IsInsideTMS(location_birth);
        RecoTrackPrimaryParticleTMSFiducialTouch[itTrack] = tp.EntersVolume(TMS_Geom::StaticIsInsideTMS);
        RecoTrackPrimaryParticleTMSFiducialEnd[itTrack] = TMS_Geom::GetInstance().IsInsideTMS(location_death);
        RecoTrackPrimaryParticleLArFiducialStart[itTrack] = TMS_Geom::GetInstance().IsInsideLAr(location_birth);
        RecoTrackPrimaryParticleLArFiducialTouch[itTrack] = tp.EntersVolume(TMS_Geom::StaticIsInsideLAr);
        RecoTrackPrimaryParticleLArFiducialEnd[itTrack] = TMS_Geom::GetInstance().IsInsideLAr(location_death);
        
        auto* vtx_info = event.GetVertexInfo(tp.GetVertexID());
        if (vtx_info != NULL) {
          RecoTrackPrimaryParticleVtxId[itTrack] = vtx_info->vtx_id;
          RecoTrackPrimaryParticleVtxFiducialCut[itTrack] = vtx_info->fiducial_cut;
          RecoTrackPrimaryParticleVtxShellEnergyCut[itTrack] = vtx_info->shell_energy_cut;
          RecoTrackPrimaryParticleVtxNDPhysicsCut[itTrack] = vtx_info->nd_physics_cut;
        }
        else {
          RecoTrackPrimaryParticleVtxId[itTrack] = -999999999.0;
          RecoTrackPrimaryParticleVtxFiducialCut[itTrack] = false;
          RecoTrackPrimaryParticleVtxShellEnergyCut[itTrack] = false;
          RecoTrackPrimaryParticleVtxNDPhysicsCut[itTrack] = false;
        }
      }
    }
    
    if (particle_info.energies.size() > 1) {
      true_secondary_visible_energy = particle_info.energies[1];
      true_secondary_particle_index = event.GetTrueParticleIndex(particle_info.vertexids[1], particle_info.trackids[1]);
    }
    if (true_secondary_particle_index < 0 || (size_t)true_secondary_particle_index  >= TrueParticles.size()) {
      true_secondary_particle_index = -999999999;
    }
    else {
      if (itTrack < __TMS_MAX_LINES__) {
        TMS_TrueParticle tp = TrueParticles[true_secondary_particle_index];
        RecoTrackSecondaryParticlePDG[itTrack] = tp.GetPDG();
        RecoTrackSecondaryParticleIsPrimary[itTrack] = tp.IsPrimary();
        setMomentum(RecoTrackSecondaryParticleTrueMomentum[itTrack], tp.GetBirthMomentum());
        setPosition(RecoTrackSecondaryParticleTruePositionStart[itTrack], tp.GetBirthPosition());
        setMomentum(RecoTrackSecondaryParticleTruePositionEnd[itTrack], tp.GetDeathPosition());
      }
    }
    
    RecoTrackTrueVisibleEnergy[itTrack] = total_true_visible_energy;
    RecoTrackPrimaryParticleIndex[itTrack] = true_primary_particle_index;
    RecoTrackPrimaryParticleTrueVisibleEnergy[itTrack] = true_primary_visible_energy;
    RecoTrackSecondaryParticleIndex[itTrack] = true_secondary_particle_index;
    RecoTrackSecondaryParticleTrueVisibleEnergy[itTrack] = true_secondary_visible_energy;
    
    
    RecoTrackPrimaryParticleTrueNHits[itTrack] = 0;
    RecoTrackSecondaryParticleTrueNHits[itTrack] = 0;
    if (particle_info.vertexids.size() > 0) {
      // Loop through all the hits, check if truth info matches primary (or secondary) particle, and then add to count
      for (size_t ih = 0; ih < RecoTrack->Hits.size(); ih++) {
        auto true_hit = RecoTrack->Hits[ih].GetTrueHit();
        for (size_t i = 0; i < true_hit.GetNTrueParticles(); i++) {
          if (true_hit.GetVertexIds(i) == particle_info.vertexids[0] && true_hit.GetPrimaryIds(i) == particle_info.trackids[0]) {
            RecoTrackPrimaryParticleTrueNHits[itTrack] += 1;
            // Only add 1 hit per true hit. True hits can have more than one instance of the same track id and vertex id after merging
            break; 
          }
          if (particle_info.vertexids.size() > 1) {
            if (true_hit.GetVertexIds(i) == particle_info.vertexids[1] && true_hit.GetPrimaryIds(i) == particle_info.trackids[1]) {
              RecoTrackSecondaryParticleTrueNHits[itTrack] += 1;
              // Only add 1 hit per true hit. True hits can have more than one instance of the same track id and vertex id after merging
              break; 
            }
          }
        }
      }
    }
    
  }
  
  
  // Clear branches
  NTrueHits = 0;
  int index = 0;
  for (auto& hit : event.GetHitsRaw()) {
    if (index >= __MAX_TRUE_TREE_ARRAY_LENGTH__) {
      std::cout<<"TMS_TreeWriter WARNING: Too many hits in event. Increase __MAX_TRUE_TREE_ARRAY_LENGTH__. Saving partial event"<<std::endl;
      break;
    }
    if (index < __MAX_TRUE_TREE_ARRAY_LENGTH__) {
      // In theory a reco hit should have many true hits based on how the merging worked
      // Plus true hits should have noise hits which don't have any parent info
      auto true_hit = hit.GetTrueHit();
      
      // Only save if more than 0.5 PE since reco hits are pedestal subtracted if < 3 PE currently
      if (true_hit.GetPE() > 0.5) {
      
        // True info
        NTrueHits += 1;
        TrueHitX[index] = true_hit.GetX();
        TrueHitY[index] = true_hit.GetY();
        TrueHitZ[index] = true_hit.GetZ();
        TrueHitT[index] = true_hit.GetT();
        TrueHitE[index] = true_hit.GetE();
        TrueHitPE[index] = true_hit.GetPE();
        TrueHitPEAfterFibers[index] = true_hit.GetPEAfterFibers();
        TrueHitPEAfterFibersLongPath[index] = true_hit.GetPEAfterFibersLongPath();
        TrueHitPEAfterFibersShortPath[index] = true_hit.GetPEAfterFibersShortPath();
        TrueHitBar[index] = hit.GetBarNumber();
        TrueHitView[index] = hit.GetBar().GetBarTypeNumber();
        TrueHitPlane[index] = hit.GetPlaneNumber();
        TrueNTrueParticles[index] = true_hit.GetNTrueParticles();
        TrueLeptonicEnergy[index] = true_hit.GetLeptonicEnergy();
        TrueHadronicEnergy[index] = true_hit.GetHadronicEnergy();
        
        // Reco info
        TrueRecoHitX[index] = hit.GetX();
        TrueRecoHitY[index] = hit.GetY();
        TrueRecoHitZ[index] = hit.GetZ();
        TrueRecoHitNotZ[index] = hit.GetNotZ();
        TrueRecoHitT[index] = hit.GetT();
        TrueRecoHitE[index] = hit.GetE();
        TrueRecoHitEVis[index] = hit.GetEVis();
        TrueRecoHitPE[index] = hit.GetPE();
        TrueRecoHitIsPedSupped[index] = hit.GetPedSup();
        TrueRecoHitTrackX[index] = hit.GetRecoX();
        TrueRecoHitTrackY[index] = hit.GetRecoY();
        TrueRecoHitTrackXUncertainty[index] = hit.GetRecoXUncertainty();
        TrueRecoHitTrackYUncertainty[index] = hit.GetRecoYUncertainty();
        
        index += 1;
      }
    }
    
  }

  Reco_Tree->Fill();
  Truth_Info->Fill();
}

void TMS_TreeWriter::FillTruthInfo(TMS_Event &event) {
  // Common code between Fill (which fills Truth_Info) and FillSpill (which fills Truth_Spill)
  
  // Fill the truth info
  EventNo = event.GetEventNumber();
  SliceNo = event.GetSliceNumber();
  SpillNo = event.GetSpillNumber();
  RunNo = event.GetRunNumber();
  Reaction = event.GetReaction();
  HasPileup = event.GetNVertices() != 1;
  nPrimaryVertices = event.GetNVertices();

  NeutrinoPDG = event.GetNeutrinoPDG();
  NeutrinoP4[0] = event.GetNeutrinoP4().X();
  NeutrinoP4[1] = event.GetNeutrinoP4().Y();
  NeutrinoP4[2] = event.GetNeutrinoP4().Z();
  NeutrinoP4[3] = event.GetNeutrinoP4().T();
  NeutrinoX4[0] = event.GetNeutrinoVtx().X();
  NeutrinoX4[1] = event.GetNeutrinoVtx().Y();
  NeutrinoX4[2] = event.GetNeutrinoVtx().Z();
  NeutrinoX4[3] = event.GetNeutrinoVtx().T();
  IsCC = (event.GetReaction().find("[CC]") != std::string::npos);

  TVector3 interaction_location = event.GetNeutrinoVtx().Vect(); 
  InteractionTMSFiducial = TMS_Geom::GetInstance().IsInsideTMS(interaction_location);
  InteractionTMSFirstTwoModules = TMS_Geom::GetInstance().IsInsideTMSFirstTwoModules(interaction_location);
  InteractionTMSThin = TMS_Geom::GetInstance().IsInsideTMSThin(interaction_location);
  InteractionLArFiducial = TMS_Geom::GetInstance().IsInsideLAr(interaction_location);
  
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
  nTruePrimaryParticles = 0;
  nTrueForgottenParticles = event.GetNTrueForgottenParticles();
  if (nTrueParticles > __TMS_MAX_TRUE_PARTICLES__) nTrueParticles = __TMS_MAX_TRUE_PARTICLES__;
  for (auto it = TrueParticles.begin(); it != TrueParticles.end(); ++it) {
    int index = it - TrueParticles.begin();
    
    if (index >= __TMS_MAX_TRUE_PARTICLES__) {
      std::cerr<<"WARNING: Found more particles than __TMS_MAX_TRUE_PARTICLES__. Stopping loop early. If this happens often, increase the max"<<std::endl;
      std::cerr<<"WARNING: In this case, the __TMS_MAX_TRUE_PARTICLES__ is "<<__TMS_MAX_TRUE_PARTICLES__<<" but need "<<TrueParticles.size()<<std::endl;
      break;
    }
    
  }
    
  nTrueParticles = TrueParticles.size();
  nTruePrimaryParticles = 0;
  nTrueForgottenParticles = event.GetNTrueForgottenParticles();
  if (nTrueParticles > __TMS_MAX_TRUE_PARTICLES__) nTrueParticles = __TMS_MAX_TRUE_PARTICLES__;
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
    IsPrimary[index] = (*it).IsPrimary();
    if ((*it).IsPrimary()) nTruePrimaryParticles += 1;
    TrueVisibleEnergy[index] = (*it).GetTrueVisibleEnergy(false);
    TrueNHits[index] = (*it).GetNTrueHits(false);
    TrueVisibleEnergyInSlice[index] = (*it).GetTrueVisibleEnergy(true);
    TrueNHitsInSlice[index] = (*it).GetNTrueHits(true);

    TVector3 location_birth = (*it).GetBirthPosition().Vect();
    TVector3 location_death = (*it).GetDeathPosition().Vect();
    TMSFiducialStart[index] = TMS_Geom::GetInstance().IsInsideTMS(location_birth);
    TMSFiducialTouch[index] = (*it).EntersVolume(TMS_Geom::StaticIsInsideTMS);
    TMSFiducialEnd[index] = TMS_Geom::GetInstance().IsInsideTMS(location_death);
    LArFiducialStart[index] = TMS_Geom::GetInstance().IsInsideLAr(location_birth);
    LArFiducialTouch[index] = (*it).EntersVolume(TMS_Geom::StaticIsInsideLAr);
    LArFiducialEnd[index] = TMS_Geom::GetInstance().IsInsideLAr(location_death);
    
    setMomentum(BirthMomentum[index], (*it).GetBirthMomentum(), (*it).GetBirthEnergy());
    setPosition(BirthPosition[index], (*it).GetBirthPosition());
    
    setMomentum(DeathMomentum[index], (*it).GetDeathMomentum(), (*it).GetDeathEnergy());
    setPosition(DeathPosition[index], (*it).GetDeathPosition());

    TruePathLength[index] = TMS_Geom::GetInstance().GetTrackLength((*it).GetPositionPoints(BirthPosition[index][2], DeathPosition[index][2]));
    TruePathLengthIgnoreY[index] =
        TMS_Geom::GetInstance().GetTrackLength((*it).GetPositionPoints(BirthPosition[index][2], DeathPosition[index][2]), true);
    TruePathLengthInTMS[index] =
        TMS_Geom::GetInstance().GetTrackLength((*it).GetPositionPoints(BirthPosition[index][2], DeathPosition[index][2], true));
    TruePathLengthInTMSIgnoreY[index] =
        TMS_Geom::GetInstance().GetTrackLength((*it).GetPositionPoints(BirthPosition[index][2], DeathPosition[index][2], true), true);

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
  
  TrueNonTMSNHits = event.GetNonTMSHits().size();
  if (TrueNonTMSNHits > __TMS_MAX_TRUE_NONTMS_HITS__) TrueNonTMSNHits = __TMS_MAX_TRUE_NONTMS_HITS__;
  int index = 0;
  int n_TrueNonTMSNHits_filled = 0;
  for (auto& hit : event.GetNonTMSHits()) {
    if (index >= __TMS_MAX_TRUE_NONTMS_HITS__) {
      std::cout<<"Warning: Found more nontms hits than __TMS_MAX_TRUE_NONTMS_HITS__. "
                 "If this happens often, increase limit from "<<__TMS_MAX_TRUE_NONTMS_HITS__<<std::endl;
      break;
    }
    //if (hit.GetE() < 0.5) continue; // Don't fill below energy threshold
    TrueNonTMSHitPos[index][0] = hit.GetX();
    TrueNonTMSHitPos[index][1] = hit.GetY();
    TrueNonTMSHitPos[index][2] = hit.GetZ();
    TrueNonTMSHitPos[index][3] = hit.GetT();
    TrueNonTMSHitEnergy[index] = hit.GetE();
    TrueNonTMSHitHadronicEnergy[index] = hit.GetHadronicEnergy();
    TrueNonTMSHitDx[index] = hit.GetdX();
    TrueNonTMSHitdEdx[index] = hit.GetdEdx();
    if (hit.GetNTrueParticles() > 1) {
      int target_vertex_id = hit.GetVertexIds(0);
      for (size_t v = 1; v < hit.GetNTrueParticles(); v++) {
        if (target_vertex_id != hit.GetVertexIds(v)) {
          std::cout<<"Fatal: found > 1 true hit GetNTrueParticles() with different vertex ids. Expecting exactly one"<<target_vertex_id<<" vs "<<hit.GetVertexIds(v)<<std::endl;
          throw std::runtime_error("Fatal: found > 1 true hit GetNTrueParticles()");
        }
      }
    }
    if (hit.GetNTrueParticles() == 1) {
      TrueNonTMSHitVertexID[index] = hit.GetVertexIds(0);
    }
    n_TrueNonTMSNHits_filled += 1;
    index += 1;
  }
  TrueNonTMSNHits = n_TrueNonTMSNHits_filled;

  auto vtx_info = event.GetVertexInfo();
  TrueVtxN = std::min((int)vtx_info.size(), __TMS_MAX_TRUE_VERTICES__);
  int true_vtx_index = 0;
  for (auto& itvtx : vtx_info) {
    auto& vtx = itvtx.second;
    if (true_vtx_index >= __TMS_MAX_TRUE_VERTICES__) {
      std::cout<<"Stopping loop. Vtx index hit max of __TMS_MAX_TRUE_VERTICES__ = "<<__TMS_MAX_TRUE_VERTICES__<<", increase limit"<<std::endl;
      break;
    }
    // Position
    TrueVtxX[true_vtx_index] = vtx.vtx.X();
    TrueVtxY[true_vtx_index] = vtx.vtx.Y();
    TrueVtxZ[true_vtx_index] = vtx.vtx.Z();
    TrueVtxT[true_vtx_index] = vtx.vtx.T();
    // Momentum
    TrueVtxPx[true_vtx_index] = vtx.p4.X();
    TrueVtxPy[true_vtx_index] = vtx.p4.Y();
    TrueVtxPz[true_vtx_index] = vtx.p4.Z();
    TrueVtxE[true_vtx_index] = vtx.p4.E();
    // Pdg
    TrueVtxPDG[true_vtx_index] = vtx.pdg;
    TrueVtxID[true_vtx_index] = vtx.vtx_id;
    TrueVtxReaction.push_back(vtx.reaction);
    // Hadronic E
    TrueVtxHadronicELarShell[true_vtx_index] = vtx.hadronic_energy_lar_shell;
    TrueVtxHadronicELAr[true_vtx_index] = vtx.hadronic_energy_lar;
    TrueVtxHadronicETMS[true_vtx_index] = vtx.hadronic_energy_tms;
    TrueVtxHadronicE[true_vtx_index] = vtx.hadronic_energy_total;
    // Visible E
    TrueVtxVisibleETMS[true_vtx_index] = vtx.true_visible_energy_tms;
    TrueVtxVisibleELAr[true_vtx_index] = vtx.true_visible_energy_lar;
    TrueVtxVisibleE[true_vtx_index] = vtx.true_visible_energy_total;
    // Truth cuts
    TrueVtxFiducialCut[true_vtx_index] = vtx.fiducial_cut;
    TrueVtxShellEnergyCut[true_vtx_index] = vtx.shell_energy_cut;
    TrueVtxNDPhysicsCut[true_vtx_index] = vtx.nd_physics_cut;
    // Finally update index
    true_vtx_index++;
  }
}

void TMS_TreeWriter::FillSpill(TMS_Event &event, int truth_info_entry_number, int truth_info_n_slices) {
  // Clear old info
  Clear();
  TruthInfoIndex = truth_info_entry_number;
  TruthInfoNSlices = truth_info_n_slices;
  
  FillTruthInfo(event);

  Truth_Spill->Fill();
}

// Reset the variables
void TMS_TreeWriter::Clear() {

  const float DEFAULT_CLEARING_FLOAT = -999999999;

  // Reset truth information

  EventNo = nParticles = NeutrinoPDG = LeptonPDG = Muon_TrueKE = Muon_TrueTrackLength = VertexIdOfMostEnergyInEvent = -999;
  VertexIdOfMostEnergyInEvent = VisibleEnergyFromVertexInSlice = TotalVisibleEnergyFromVertex = VisibleEnergyFromOtherVerticesInSlice = -999;
  LArOuterShellEnergy = LArTotalEnergy = TotalNonTMSEnergy = DEFAULT_CLEARING_FLOAT;
  LArOuterShellEnergyFromVertex = LArTotalEnergyFromVertex = TotalNonTMSEnergyFromVertex = DEFAULT_CLEARING_FLOAT;
  Reaction = "";
  IsCC = false;

  for (int i = 0; i < 4; ++i) {
    MuonP4[i]=DEFAULT_CLEARING_FLOAT;
    Muon_Vertex[i]=DEFAULT_CLEARING_FLOAT;
    Muon_Death[i]=DEFAULT_CLEARING_FLOAT;
    NeutrinoP4[i]=DEFAULT_CLEARING_FLOAT;
    NeutrinoX4[i]=DEFAULT_CLEARING_FLOAT;
    LeptonP4[i]=DEFAULT_CLEARING_FLOAT;
    LeptonX4[i]=DEFAULT_CLEARING_FLOAT;
  }

  // Reset line information
  TMSStart = false;
  nLinesU = DEFAULT_CLEARING_FLOAT;
  nLinesV = DEFAULT_CLEARING_FLOAT;
  nLinesX = DEFAULT_CLEARING_FLOAT;
  nLinesY = DEFAULT_CLEARING_FLOAT;
  for (int i = 0; i < __TMS_MAX_LINES__; ++i) {
    SlopeU[i] = DEFAULT_CLEARING_FLOAT;
    SlopeV[i] = DEFAULT_CLEARING_FLOAT;
    SlopeX[i] = DEFAULT_CLEARING_FLOAT;
    SlopeY[i] = DEFAULT_CLEARING_FLOAT;
    InterceptU[i] = DEFAULT_CLEARING_FLOAT;
    InterceptV[i] = DEFAULT_CLEARING_FLOAT;
    InterceptX[i] = DEFAULT_CLEARING_FLOAT;
    InterceptY[i] = DEFAULT_CLEARING_FLOAT;
    Slope_DownstreamU[i] = DEFAULT_CLEARING_FLOAT;
    Slope_DownstreamV[i] = DEFAULT_CLEARING_FLOAT;
    Slope_DownstreamX[i] = DEFAULT_CLEARING_FLOAT;
    Slope_DownstreamY[i] = DEFAULT_CLEARING_FLOAT;
    Intercept_DownstreamU[i] = DEFAULT_CLEARING_FLOAT;
    Intercept_DownstreamV[i] = DEFAULT_CLEARING_FLOAT;
    Intercept_DownstreamX[i] = DEFAULT_CLEARING_FLOAT;
    Intercept_DownstreamY[i] = DEFAULT_CLEARING_FLOAT;
    Slope_UpstreamU[i] = DEFAULT_CLEARING_FLOAT;
    Slope_UpstreamV[i] = DEFAULT_CLEARING_FLOAT;
    Slope_UpstreamX[i] = DEFAULT_CLEARING_FLOAT;
    Slope_UpstreamY[i] = DEFAULT_CLEARING_FLOAT;
    Intercept_UpstreamU[i] = DEFAULT_CLEARING_FLOAT;
    Intercept_UpstreamV[i] = DEFAULT_CLEARING_FLOAT;
    Intercept_UpstreamX[i] = DEFAULT_CLEARING_FLOAT;
    Intercept_UpstreamY[i] = DEFAULT_CLEARING_FLOAT;

    DirectionZU[i] = DEFAULT_CLEARING_FLOAT;
    DirectionXU[i] = DEFAULT_CLEARING_FLOAT;
    DirectionZU_Upstream[i] = DEFAULT_CLEARING_FLOAT;
    DirectionXU_Upstream[i] = DEFAULT_CLEARING_FLOAT;
    DirectionZU_Downstream[i] = DEFAULT_CLEARING_FLOAT;
    DirectionXU_Downstream[i] = DEFAULT_CLEARING_FLOAT;

    DirectionZV[i] = DEFAULT_CLEARING_FLOAT;
    DirectionXV[i] = DEFAULT_CLEARING_FLOAT;
    DirectionZV_Upstream[i] = DEFAULT_CLEARING_FLOAT;
    DirectionXV_Upstream[i] = DEFAULT_CLEARING_FLOAT;
    DirectionZV_Downstream[i] = DEFAULT_CLEARING_FLOAT;
    DirectionXV_Downstream[i] = DEFAULT_CLEARING_FLOAT;

    DirectionZX[i] = DEFAULT_CLEARING_FLOAT;
    DirectionYX[i] = DEFAULT_CLEARING_FLOAT;
    DirectionZX_Upstream[i] = DEFAULT_CLEARING_FLOAT;
    DirectionYX_Upstream[i] = DEFAULT_CLEARING_FLOAT;
    DirectionZX_Downstream[i] = DEFAULT_CLEARING_FLOAT;
    DirectionYX_Downstream[i] = DEFAULT_CLEARING_FLOAT;

    DirectionZY[i] = DEFAULT_CLEARING_FLOAT;
    DirectionXY[i] = DEFAULT_CLEARING_FLOAT;
    DirectionZY_Upstream[i] = DEFAULT_CLEARING_FLOAT;
    DirectionXY_Upstream[i] = DEFAULT_CLEARING_FLOAT;
    DirectionZY_Downstream[i] = DEFAULT_CLEARING_FLOAT;
    DirectionXY_Downstream[i] = DEFAULT_CLEARING_FLOAT;

    OccupancyU[i] = DEFAULT_CLEARING_FLOAT;
    OccupancyV[i] = DEFAULT_CLEARING_FLOAT;
    OccupancyX[i] = DEFAULT_CLEARING_FLOAT;
    OccupancyY[i] = DEFAULT_CLEARING_FLOAT;
    TrackLengthU[i] = DEFAULT_CLEARING_FLOAT;
    TrackLengthV[i] = DEFAULT_CLEARING_FLOAT;
    TrackLengthX[i] = DEFAULT_CLEARING_FLOAT;
    TrackLengthY[i] = DEFAULT_CLEARING_FLOAT;
    TotalTrackEnergyU[i] = DEFAULT_CLEARING_FLOAT;
    TotalTrackEnergyV[i] = DEFAULT_CLEARING_FLOAT;
    TotalTrackEnergyX[i] = DEFAULT_CLEARING_FLOAT;
    TotalTrackEnergyY[i] = DEFAULT_CLEARING_FLOAT;
    FirstPlaneU[i] = DEFAULT_CLEARING_FLOAT;
    FirstPlaneV[i] = DEFAULT_CLEARING_FLOAT;
    FirstPlaneX[i] = DEFAULT_CLEARING_FLOAT;
    FirstPlaneY[i] = DEFAULT_CLEARING_FLOAT;
    LastPlaneU[i] = DEFAULT_CLEARING_FLOAT;
    LastPlaneV[i] = DEFAULT_CLEARING_FLOAT;
    LastPlaneX[i] = DEFAULT_CLEARING_FLOAT;
    LastPlaneY[i] = DEFAULT_CLEARING_FLOAT;
    nHitsInTrackU[i] = DEFAULT_CLEARING_FLOAT;
    nHitsInTrackV[i] = DEFAULT_CLEARING_FLOAT;
    nHitsInTrackX[i] = DEFAULT_CLEARING_FLOAT;
    nHitsInTrackY[i] = DEFAULT_CLEARING_FLOAT;
    TrackStoppingU[i] = false;
    TrackStoppingV[i] = false;
    TrackStoppingX[i] = false;
    TrackStoppingY[i] = false;

    for (int k = 0; k < 3; k++) {
      RecoTrackKalmanFirstPlaneBarView[i][k] = DEFAULT_CLEARING_FLOAT;
      RecoTrackKalmanLastPlaneBarView[i][k] = DEFAULT_CLEARING_FLOAT;
      for (int j = 0; j < __TMS_MAX_LINE_HITS__; j++) {
        RecoTrackKalmanPlaneBarView[i][j][k] = DEFAULT_CLEARING_FLOAT;
        RecoTrackKalmanPlaneBarViewTrue[i][j][k] = DEFAULT_CLEARING_FLOAT;
      }
    }
  }
/*    Occupancy3D[i] = DEFAULT_CLEARING_FLOAT;
    TrackLength3D[i] = DEFAULT_CLEARING_FLOAT;
    TotalTrackEnergy3D[i] = DEFAULT_CLEARING_FLOAT;
    FirstPlane3D[i] = DEFAULT_CLEARING_FLOAT;
    LastPlane3D[i] = DEFAULT_CLEARING_FLOAT;
    nHitsInTrack3D[i] = DEFAULT_CLEARING_FLOAT;
    TrackStopping3D[i] = false;

    for (int j = 0; j < 2; ++j) {
      FirstHitOne[i][j] = DEFAULT_CLEARING_FLOAT;
      FirstHitOther[i][j] = DEFAULT_CLEARING_FLOAT;
      LastHitOne[i][j] = DEFAULT_CLEARING_FLOAT;
      LastHitOther[i][j] = DEFAULT_CLEARING_FLOAT;
    }
    for (int j = 0; j < 3; ++j) {
      FirstHit3D[i][j] = DEFAULT_CLEARING_FLOAT;
      LastHit3D[i][j] = DEFAULT_CLEARING_FLOAT;
    }
    for (int j = 0; j < __TMS_MAX_LINE_HITS__; ++j) {
      TrackHitEnergyOne[i][j]=DEFAULT_CLEARING_FLOAT;
      TrackHitEnergyOther[i][j]=DEFAULT_CLEARING_FLOAT;
      TrackHitPosOne[i][j][0]=DEFAULT_CLEARING_FLOAT;
      TrackHitPosOne[i][j][1]=DEFAULT_CLEARING_FLOAT;
      TrackHitPosOther[i][j][0]=DEFAULT_CLEARING_FLOAT;
      TrackHitPosOther[i][j][1]=DEFAULT_CLEARING_FLOAT;

      TrackHitEnergy3D[i][j] = DEFAULT_CLEARING_FLOAT;
      TrackHitPos3D[i][j][0] = DEFAULT_CLEARING_FLOAT;
      TrackHitPos3D[i][j][1] = DEFAULT_CLEARING_FLOAT;
      TrackHitPos3D[i][j][2] = DEFAULT_CLEARING_FLOAT;
    }
  }*/

  // Reset hit information
  nHits = DEFAULT_CLEARING_FLOAT;
  for (int i = 0; i < __TMS_MAX_HITS__; ++i) {
    for (int j = 0; j < 4; ++j) RecoHitPos[i][j] = DEFAULT_CLEARING_FLOAT;
    RecoHitEnergy[i] = DEFAULT_CLEARING_FLOAT;
    RecoHitBarType[i] = DEFAULT_CLEARING_FLOAT;
    RecoHitPE[i] = DEFAULT_CLEARING_FLOAT;
    RecoHitBar[i] = DEFAULT_CLEARING_FLOAT;
    RecoHitPlane[i] = DEFAULT_CLEARING_FLOAT;
    RecoHitSlice[i] = DEFAULT_CLEARING_FLOAT;
  }

  // Reset Cluster info
  nClustersU = DEFAULT_CLEARING_FLOAT;
  nClustersV = DEFAULT_CLEARING_FLOAT;
  nClustersX = DEFAULT_CLEARING_FLOAT;
  nClustersY = DEFAULT_CLEARING_FLOAT;
  for (int i = 0; i < __TMS_MAX_CLUSTERS__; ++i) {
    ClusterEnergyU[i] = DEFAULT_CLEARING_FLOAT;
    ClusterEnergyV[i] = DEFAULT_CLEARING_FLOAT;
    ClusterEnergyX[i] = DEFAULT_CLEARING_FLOAT;
    ClusterEnergyY[i] = DEFAULT_CLEARING_FLOAT;
    nHitsInClusterU[i] = DEFAULT_CLEARING_FLOAT;
    nHitsInClusterV[i] = DEFAULT_CLEARING_FLOAT;
    nHitsInClusterX[i] = DEFAULT_CLEARING_FLOAT;
    nHitsInClusterY[i] = DEFAULT_CLEARING_FLOAT;
    for (int j = 0; j < 2; ++j) {
      ClusterPosMeanU[i][j] = DEFAULT_CLEARING_FLOAT;
      ClusterPosStdDevU[i][j] = DEFAULT_CLEARING_FLOAT;
      
      ClusterPosMeanV[i][j] = DEFAULT_CLEARING_FLOAT;
      ClusterPosStdDevV[i][j] = DEFAULT_CLEARING_FLOAT;

      ClusterPosMeanX[i][j] = DEFAULT_CLEARING_FLOAT;
      ClusterPosStdDevX[i][j] = DEFAULT_CLEARING_FLOAT;

      ClusterPosMeanY[i][j] = DEFAULT_CLEARING_FLOAT;
      ClusterPosStdDevY[i][j] = DEFAULT_CLEARING_FLOAT;
    }
    for (int j = 0; j < __TMS_MAX_LINE_HITS__; ++j) {
      ClusterHitPosU[i][j][0] = DEFAULT_CLEARING_FLOAT;
      ClusterHitPosU[i][j][1] = DEFAULT_CLEARING_FLOAT;
      ClusterHitEnergyU[i][j] = DEFAULT_CLEARING_FLOAT;

      ClusterHitPosV[i][j][0] = DEFAULT_CLEARING_FLOAT;
      ClusterHitPosV[i][j][1] = DEFAULT_CLEARING_FLOAT;
      ClusterHitEnergyV[i][j] = DEFAULT_CLEARING_FLOAT;

      ClusterHitPosX[i][j][0] = DEFAULT_CLEARING_FLOAT;
      ClusterHitPosX[i][j][1] = DEFAULT_CLEARING_FLOAT;
      ClusterHitEnergyX[i][j] = DEFAULT_CLEARING_FLOAT;

      ClusterHitPosY[i][j][0] = DEFAULT_CLEARING_FLOAT;
      ClusterHitPosY[i][j][1] = DEFAULT_CLEARING_FLOAT;
      ClusterHitEnergyY[i][j] = DEFAULT_CLEARING_FLOAT;
    }
  }

  // Reset track information
  nTracks = DEFAULT_CLEARING_FLOAT;
  TimeSliceStartTime = DEFAULT_CLEARING_FLOAT;
  TimeSliceEndTime = DEFAULT_CLEARING_FLOAT;
  for (int i = 0; i < __TMS_MAX_TRACKS__; ++i) {
    for (int j = 0; j < 3; ++j) {
      RecoTrackStartPos[i][j] = DEFAULT_CLEARING_FLOAT;
      RecoTrackStartDirection[i][j] = DEFAULT_CLEARING_FLOAT;
      RecoTrackEndPos[i][j] = DEFAULT_CLEARING_FLOAT;
      RecoTrackEndDirection[i][j] = DEFAULT_CLEARING_FLOAT;
    }
    for (int k = 0; k < __TMS_MAX_LINE_HITS__; ++k) {
      RecoTrackHitPos[i][k][0] = DEFAULT_CLEARING_FLOAT;
      RecoTrackHitPos[i][k][1] = DEFAULT_CLEARING_FLOAT;
      RecoTrackHitPos[i][k][2] = DEFAULT_CLEARING_FLOAT;
     
      RecoTrackHitEnergies[i][k] = DEFAULT_CLEARING_FLOAT;
    }
    RecoTrackEnergyRange[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackEnergyDeposit[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackLength[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackLength_3D[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackCharge[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackCharge_Kalman[i] = DEFAULT_CLEARING_FLOAT;
  }
  
  RecoTrackN = 0;
  for (int i = 0; i < __TMS_MAX_LINES__; ++i) {
    RecoTrackTrueVisibleEnergy[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackPrimaryParticleIndex[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackPrimaryParticleTrueVisibleEnergy[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackPrimaryParticleTrueNHits[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackSecondaryParticleIndex[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackSecondaryParticleTrueVisibleEnergy[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackSecondaryParticleTrueNHits[i] = DEFAULT_CLEARING_FLOAT;

    RecoTrackPrimaryParticleTrueTrackLengthAsMeasured[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackPrimaryParticleTrueTrackLengthAsMeasuredIgnoreY[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackPrimaryParticleTrueTrackLengthRecoStart[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackPrimaryParticleTrueTrackLengthRecoStartIgnoreY[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackPrimaryParticleTrueTrackLengthInTMS[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackPrimaryParticleTrueTrackLengthInTMSIgnoreY[i] = DEFAULT_CLEARING_FLOAT;
    
    RecoTrackPrimaryParticlePDG[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackSecondaryParticlePDG[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackPrimaryParticleIsPrimary[i] = false;
    RecoTrackSecondaryParticleIsPrimary[i] = false;
    RecoTrackPrimaryParticleTrueTrackLength[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackPrimaryParticleTrueTrackLengthIgnoreY[i] = DEFAULT_CLEARING_FLOAT;
    
    RecoTrackPrimaryParticleTMSFiducialStart[i] = false;
    RecoTrackPrimaryParticleTMSFiducialTouch[i] = false;
    RecoTrackPrimaryParticleTMSFiducialEnd[i] = false;
    RecoTrackPrimaryParticleLArFiducialStart[i] = false;
    RecoTrackPrimaryParticleLArFiducialTouch[i] = false;
    RecoTrackPrimaryParticleLArFiducialEnd[i] = false;

    RecoTrackPrimaryParticleVtxId[i] = DEFAULT_CLEARING_FLOAT;
    RecoTrackPrimaryParticleVtxFiducialCut[i] = false;
    RecoTrackPrimaryParticleVtxShellEnergyCut[i] = false;
    RecoTrackPrimaryParticleVtxNDPhysicsCut[i] = false;

    for (int j = 0; j < 4; ++j) {
      RecoTrackPrimaryParticleTrueMomentumTrackStart[i][j] = DEFAULT_CLEARING_FLOAT;
      RecoTrackPrimaryParticleTruePositionTrackStart[i][j] = DEFAULT_CLEARING_FLOAT;
      RecoTrackPrimaryParticleTrueMomentumTrackEnd[i][j] = DEFAULT_CLEARING_FLOAT;
      RecoTrackPrimaryParticleTruePositionTrackEnd[i][j] = DEFAULT_CLEARING_FLOAT;
      
      RecoTrackPrimaryParticleTrueMomentum[i][j] = DEFAULT_CLEARING_FLOAT;
      RecoTrackPrimaryParticleTruePositionStart[i][j] = DEFAULT_CLEARING_FLOAT;
      RecoTrackPrimaryParticleTruePositionEnd[i][j] = DEFAULT_CLEARING_FLOAT;
      RecoTrackPrimaryParticleTrueMomentumEnteringTMS[i][j] = DEFAULT_CLEARING_FLOAT;
      RecoTrackPrimaryParticleTruePositionEnteringTMS[i][j] = DEFAULT_CLEARING_FLOAT;
      RecoTrackPrimaryParticleTrueMomentumLeavingTMS[i][j] = DEFAULT_CLEARING_FLOAT;
      RecoTrackPrimaryParticleTruePositionLeavingTMS[i][j] = DEFAULT_CLEARING_FLOAT;
      RecoTrackPrimaryParticleTrueMomentumLeavingLAr[i][j] = DEFAULT_CLEARING_FLOAT;
      RecoTrackPrimaryParticleTruePositionLeavingLAr[i][j] = DEFAULT_CLEARING_FLOAT;
      
      RecoTrackSecondaryParticleTrueMomentum[i][j] = DEFAULT_CLEARING_FLOAT;
      RecoTrackSecondaryParticleTruePositionStart[i][j] = DEFAULT_CLEARING_FLOAT;
      RecoTrackSecondaryParticleTruePositionEnd[i][j] = DEFAULT_CLEARING_FLOAT;
      
      for (int h = 0; h < __TMS_MAX_LINE_HITS__; h++) {
        RecoTrackTrueHitPosition[i][h][j] = DEFAULT_CLEARING_FLOAT;
      }
    }
  }
    
  nTrueParticles = 0;
  TrueNonTMSNHits = 0;
  for (int i = 0; i < __TMS_MAX_TRUE_PARTICLES__; ++i) {
    VertexID[i] = DEFAULT_CLEARING_FLOAT;
    Parent[i] = DEFAULT_CLEARING_FLOAT;
    TrackId[i] = DEFAULT_CLEARING_FLOAT;
    PDG[i] = DEFAULT_CLEARING_FLOAT;
    IsPrimary[i] = false;
    TrueVisibleEnergy[i] = DEFAULT_CLEARING_FLOAT;
    TrueNHits[i] = DEFAULT_CLEARING_FLOAT;
    TrueVisibleEnergyInSlice[i] = DEFAULT_CLEARING_FLOAT;
    TrueNHitsInSlice[i] = DEFAULT_CLEARING_FLOAT;
    TruePathLength[i] = DEFAULT_CLEARING_FLOAT;
    TruePathLengthIgnoreY[i] = DEFAULT_CLEARING_FLOAT;
    TruePathLengthInTMS[i] = DEFAULT_CLEARING_FLOAT;
    TruePathLengthInTMSIgnoreY[i] = DEFAULT_CLEARING_FLOAT;
    
    TMSFiducialStart[i] = false;
    TMSFiducialTouch[i] = false;
    TMSFiducialEnd[i] = false;
    LArFiducialStart[i] = false;
    LArFiducialTouch[i] = false;
    LArFiducialEnd[i] = false;
    
    TrueNonTMSHitEnergy[i] = DEFAULT_CLEARING_FLOAT;
    TrueNonTMSHitHadronicEnergy[i] = DEFAULT_CLEARING_FLOAT;
    TrueNonTMSHitDx[i] = DEFAULT_CLEARING_FLOAT;
    TrueNonTMSHitdEdx[i] = DEFAULT_CLEARING_FLOAT;
    TrueNonTMSHitVertexID[i] = -99999;

    for (int j = 0; j < 4; ++j) {
      TrueNonTMSHitPos[i][j] = DEFAULT_CLEARING_FLOAT;
      
      BirthMomentum[i][j] = DEFAULT_CLEARING_FLOAT;
      BirthPosition[i][j] = DEFAULT_CLEARING_FLOAT;
      DeathMomentum[i][j] = DEFAULT_CLEARING_FLOAT;
      DeathPosition[i][j] = DEFAULT_CLEARING_FLOAT;
      MomentumZIsLArEnd[i][j] = DEFAULT_CLEARING_FLOAT;
      PositionZIsLArEnd[i][j] = DEFAULT_CLEARING_FLOAT;
      MomentumLArStart[i][j] = DEFAULT_CLEARING_FLOAT;
      PositionLArStart[i][j] = DEFAULT_CLEARING_FLOAT;
      MomentumLArEnd[i][j] = DEFAULT_CLEARING_FLOAT;
      PositionLArEnd[i][j] = DEFAULT_CLEARING_FLOAT;
      MomentumZIsTMSStart[i][j] = DEFAULT_CLEARING_FLOAT;
      PositionZIsTMSStart[i][j] = DEFAULT_CLEARING_FLOAT;
      MomentumZIsTMSEnd[i][j] = DEFAULT_CLEARING_FLOAT;
      PositionZIsTMSEnd[i][j] = DEFAULT_CLEARING_FLOAT;
      MomentumTMSStart[i][j] = DEFAULT_CLEARING_FLOAT;
      PositionTMSStart[i][j] = DEFAULT_CLEARING_FLOAT;
      MomentumTMSFirstTwoModulesEnd[i][j] = DEFAULT_CLEARING_FLOAT;
      PositionTMSFirstTwoModulesEnd[i][j] = DEFAULT_CLEARING_FLOAT;
      MomentumTMSThinEnd[i][j] = DEFAULT_CLEARING_FLOAT;
      PositionTMSThinEnd[i][j] = DEFAULT_CLEARING_FLOAT;
      MomentumTMSEnd[i][j] = DEFAULT_CLEARING_FLOAT;
      PositionTMSEnd[i][j] = DEFAULT_CLEARING_FLOAT;
    }
  }


}

