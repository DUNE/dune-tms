#include "TMS_Manager.h"

TMS_Manager::TMS_Manager() {

  if (!std::getenv("TMS_DIR")) {
    std::cerr << "Need ${TMS_DIR} environment set for reconstruction, please export TMS_DIR" << std::endl;
    throw;
  }

  std::string filename;
  if (std::getenv("TMS_TOML") != NULL) filename = std::string(std::getenv("TMS_TOML"));
  else filename = std::string(std::getenv("TMS_DIR"))+"/config/TMS_Default_Config.toml";

  std::cout << "Creating TMS Manager instance using TOML: " << filename << std::endl;

  // Read the TOML file
  const auto data = toml::parse(filename);

  // The minimum hits needed to run reconstruction in a TMS event
  _RECO_MinHits = toml::find<int>(data, "Recon", "MinHits");
  
  _RECO_TIME_RunTimeSlicer = toml::find<bool>(data, "Recon", "Time", "RunTimeSlicer");
  _RECO_TIME_RunSimpleTimeSlicer = toml::find<bool>(data, "Recon", "Time", "RunSimpleTimeSlicer");
  _RECO_TIME_TimeSlicerThresholdStart = toml::find<double>(data, "Recon", "Time", "TimeSlicerThresholdStart");
  _RECO_TIME_TimeSlicerThresholdEnd = toml::find<double>(data, "Recon", "Time", "TimeSlicerThresholdEnd");
  _RECO_TIME_TimeSlicerSliceUnit = toml::find<double>(data, "Recon", "Time", "SliceUnit");
  _RECO_TIME_TimeSlicerEnergyWindowInUnits = toml::find<int>(data, "Recon", "Time", "TimeSlicerEnergyWindowInUnits");
  _RECO_TIME_TimeSlicerMinimumSliceWidthInUnits = toml::find<int>(data, "Recon", "Time", "TimeSlicerMinimumSliceWidthInUnits");
  _RECO_TIME_TimeSlicerMaxTime = toml::find<double>(data, "Recon", "Time", "TimeSlicerMaxTime");

  _RECO_DBSCAN_MinPoints = toml::find<int>(data, "Recon", "DBSCAN", "MinPoints");
  _RECO_DBSCAN_Epsilon = toml::find<double>(data, "Recon", "DBSCAN", "Epsilon");
  _RECO_DBSCAN_PreDBNeighbours = toml::find<int>(data, "Recon", "DBSCAN", "PreDBNeighbours");
  _RECO_DBSCAN_PreDBDistance = toml::find<double>(data, "Recon", "DBSCAN", "PreDBDistance");

  _RECO_HOUGH_MaxHough = toml::find<int>(data, "Recon", "Hough", "MaxTrans");
  _RECO_HOUGH_MinInterp = toml::find<double>(data, "Recon", "Hough", "MinInter");
  _RECO_HOUGH_MaxInterp = toml::find<double>(data, "Recon", "Hough", "MaxInter");
  _RECO_HOUGH_MinSlope = toml::find<double>(data, "Recon", "Hough", "MinSlope");
  _RECO_HOUGH_MaxSlope = toml::find<double>(data, "Recon", "Hough", "MaxSlope");
  _RECO_HOUGH_NSlope = toml::find<int>(data, "Recon", "Hough", "NSlope");
  _RECO_HOUGH_NInter = toml::find<int>(data, "Recon", "Hough", "NInter");
  _RECO_HOUGH_HitMult = toml::find<double>(data, "Recon", "Hough", "HitMult");
  _RECO_HOUGH_MergeTracks = toml::find<bool>(data, "Recon", "Hough", "MergeTracks");
  _RECO_HOUGH_RunAStar = toml::find<bool>(data, "Recon", "Hough", "RunAStarCleanup");
  _RECO_HOUGH_FirstCluster = toml::find<bool>(data, "Recon", "Hough", "FirstCluster");
  _RECO_HOUGH_MinDist = toml::find<double>(data, "Recon", "Hough", "MinDist");

  _RECO_EXTRAPOLATION_Extrapolation = toml::find<bool>(data, "Recon", "Extrapolation", "Extrapolation");
  _RECO_EXTRAPOLATION_ExtrapolateDist = toml::find<int>(data, "Recon", "Extrapolation", "ExtrapolateDist");
  _RECO_EXTRAPOLATION_ExtrapolateLimit = toml::find<int>(data, "Recon", "Extrapolation", "ExtrapolateLimit");
  _RECO_EXTRAPOLATION_NumBarsEnd = toml::find<int>(data, "Recon", "Extrapolation", "NumBarsEnd");
  _RECO_EXTRAPOLATION_NumBarsStart = toml::find<int>(data, "Recon", "Extrapolation", "NumBarsStart");

  _RECO_TRACKMATCH_PlaneLimit = toml::find<int>(data, "Recon", "TrackMatch3D", "PlaneLimit");
  _RECO_TRACKMATCH_BarLimit = toml::find<int>(data, "Recon", "TrackMatch3D", "BarLimit");
  _RECO_TRACKMATCH_TimeLimit = toml::find<int>(data, "Recon", "TrackMatch3D", "TimeLimit");
  _RECO_TRACKMATCH_XTimeLimit = toml::find<int>(data, "Recon", "TrackMatch3D", "XTimeLimit");
  _RECO_TRACKMATCH_YAnchor = toml::find<float>(data, "Recon", "TrackMatch3D", "YAnchor");
  _RECO_TRACKMATCH_TiltAngle = toml::find<double>(data, "Recon", "TrackMatch3D", "TiltAngle");
  _RECO_TRACKMATCH_YDifference = toml::find<float>(data, "Recon", "TrackMatch3D", "YDifference");
  _RECO_TRACKMATCH_DirectionDistance = toml::find<int>(data, "Recon", "TrackMatch3D", "DirectionDistance");

  _RECO_ASTAR_IsGreedy = toml::find<bool> (data, "Recon", "AStar", "IsGreedy");
  _RECO_ASTAR_CostMetric = toml::find<std::string> (data, "Recon", "AStar", "CostMetric");

  _RECO_STOPPING_nLastHits = toml::find<int>(data, "Recon", "Stopping", "nLastHits");
  _RECO_STOPPING_EnergyCut = toml::find<double>(data, "Recon", "Stopping", "EnergyCut");

  _RECO_TRACK_METHOD  = toml::find<std::string>(data, "Recon", "TrackMethod");
  _RECO_CLUSTERING    = toml::find<bool>  (data, "Recon", "Clustering");

  _RECO_KALMAN_RUN = toml::find<bool>(data, "Recon", "Kalman", "Run");
  
  _RECO_CALIBRATION_EnergyCalibration = toml::find<double>  (data, "Recon", "Calibration", "EnergyCalibration");
  
  _FIDUCIAL_TMS_START_X = toml::find<double>(data, "Fiducial", "TMS", "Start", "X");
  _FIDUCIAL_TMS_START_Y = toml::find<double>(data, "Fiducial", "TMS", "Start", "Y");
  _FIDUCIAL_TMS_START_Z = toml::find<double>(data, "Fiducial", "TMS", "Start", "Z");
  _FIDUCIAL_TMS_END_X = toml::find<double>(data, "Fiducial", "TMS", "End", "X");
  _FIDUCIAL_TMS_END_Y = toml::find<double>(data, "Fiducial", "TMS", "End", "Y");
  _FIDUCIAL_TMS_END_Z = toml::find<double>(data, "Fiducial", "TMS", "End", "Z");
  _ACTIVE_LAR_START_X = toml::find<double>(data, "Fiducial", "LAr", "Start", "X");
  _ACTIVE_LAR_START_Y = toml::find<double>(data, "Fiducial", "LAr", "Start", "Y");
  _ACTIVE_LAR_START_Z = toml::find<double>(data, "Fiducial", "LAr", "Start", "Z");
  _ACTIVE_LAR_END_X = toml::find<double>(data, "Fiducial", "LAr", "End", "X");
  _ACTIVE_LAR_END_Y = toml::find<double>(data, "Fiducial", "LAr", "End", "Y");
  _ACTIVE_LAR_END_Z = toml::find<double>(data, "Fiducial", "LAr", "End", "Z");

  _TRUTH_LIGHTWEIGHT = toml::find<bool> (data, "Truth", "LightWeight");

  _APPLICATIONS_DrawPDF =  toml::find<bool> (data, "Applications", "DrawPDF");
  
  _APPLICATIONS_MaximumNEvents = toml::find<int>(data, "Applications", "MaximumNEvents");

  _GEOMETRY_YMIDDLE = toml::find<double>(data, "Geometry", "YBarMiddle");
  
  _NERSC_SPILL_PERIOD = 1.2e9;
}
