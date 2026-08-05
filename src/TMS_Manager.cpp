#include "TMS_Manager.h"

#include <stdexcept>

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

  _RECO_TRACKSMOOTHING_UseTrackSmoothing = toml::find<bool>(data, "Recon", "TrackSmoothing", "UseTrackSmoothing");
  _RECO_TRACKSMOOTHING_TrackSmoothingStrategy = toml::find<std::string>(data, "Recon", "TrackSmoothing", "TrackSmoothingStrategy");
  _RECO_TRACKSMOOTHING_MaxYDistanceBetweenUVTransitionPoints =
        toml::find<double>(data, "Recon", "TrackSmoothing", "MaxYDistanceBetweenUVTransitionPoints");
  _RECO_TRACKSMOOTHING_UncertaintyGoodDirection = toml::find<double>(data, "Recon", "TrackSmoothing", "UncertaintyGoodDirection");
  _RECO_TRACKSMOOTHING_UncertaintyBadDirection = toml::find<double>(data, "Recon", "TrackSmoothing", "UncertaintyBadDirection");
  _RECO_TRACKSMOOTHING_UncertaintyForUVTransitionPoints = toml::find<double>(data, "Recon", "TrackSmoothing", "UncertaintyForUVTransitionPoints");

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
  _RECO_HOUGH_EndpointRescanFraction = toml::find<double>(data, "Recon", "Hough", "EndpointRescanFraction");
  _RECO_HOUGH_EndpointRescanMinimumHits = toml::find<int>(data, "Recon", "Hough", "EndpointRescanMinimumHits");
  _RECO_HOUGH_RefitInterceptTolerance = toml::find<double>(data, "Recon", "Hough", "RefitInterceptTolerance");
  _RECO_HOUGH_RefitSlopeTolerance = toml::find<double>(data, "Recon", "Hough", "RefitSlopeTolerance");
  _RECO_HOUGH_MergeEndpointPlaneGap = toml::find<int>(data, "Recon", "Hough", "MergeEndpointPlaneGap");
  _RECO_HOUGH_MergeEndpointBarGap = toml::find<int>(data, "Recon", "Hough", "MergeEndpointBarGap");
  _RECO_HOUGH_MergeInterceptDifference = toml::find<double>(data, "Recon", "Hough", "MergeInterceptDifference");
  _RECO_HOUGH_MergeSlopeDifference = toml::find<double>(data, "Recon", "Hough", "MergeSlopeDifference");
  _RECO_HOUGH_ContainmentHalfWidthX = toml::find<double>(data, "Recon", "Hough", "ContainmentHalfWidthX");
  _RECO_HOUGH_ContainmentHalfWidthY = toml::find<double>(data, "Recon", "Hough", "ContainmentHalfWidthY");

  _RECO_EXTRAPOLATION_Extrapolation = toml::find<bool>(data, "Recon", "Extrapolation", "Extrapolation");
  _RECO_EXTRAPOLATION_ExtrapolateDist = toml::find<int>(data, "Recon", "Extrapolation", "ExtrapolateDist");
  _RECO_EXTRAPOLATION_ExtrapolateLimit = toml::find<int>(data, "Recon", "Extrapolation", "ExtrapolateLimit");
  _RECO_EXTRAPOLATION_NumBarsEnd = toml::find<int>(data, "Recon", "Extrapolation", "NumBarsEnd");
  _RECO_EXTRAPOLATION_NumBarsStart = toml::find<int>(data, "Recon", "Extrapolation", "NumBarsStart");
  _RECO_EXTRAPOLATION_ContainmentWidthScaleX = toml::find<double>(data, "Recon", "Extrapolation", "ContainmentWidthScaleX");
  _RECO_EXTRAPOLATION_ContainmentWidthScaleY = toml::find<double>(data, "Recon", "Extrapolation", "ContainmentWidthScaleY");
  _RECO_EXTRAPOLATION_EndpointSeedHits = toml::find<int>(data, "Recon", "Extrapolation", "EndpointSeedHits");
  _RECO_EXTRAPOLATION_MultiCandidateThreshold = toml::find<int>(data, "Recon", "Extrapolation", "MultiCandidateThreshold");

  _RECO_TRACKMATCH_PlaneLimit = toml::find<int>(data, "Recon", "TrackMatch3D", "PlaneLimit");
  _RECO_TRACKMATCH_BarLimit = toml::find<int>(data, "Recon", "TrackMatch3D", "BarLimit");
  _RECO_TRACKMATCH_TimeLimit = toml::find<int>(data, "Recon", "TrackMatch3D", "TimeLimit");
  _RECO_TRACKMATCH_XTimeLimit = toml::find<int>(data, "Recon", "TrackMatch3D", "XTimeLimit");
  _RECO_TRACKMATCH_YAnchor = toml::find<float>(data, "Recon", "TrackMatch3D", "YAnchor");
  _RECO_TRACKMATCH_TiltAngle = toml::find<double>(data, "Recon", "TrackMatch3D", "TiltAngle");
  _RECO_TRACKMATCH_YDifference = toml::find<float>(data, "Recon", "TrackMatch3D", "YDifference");
  _RECO_TRACKMATCH_DirectionDistance = toml::find<int>(data, "Recon", "TrackMatch3D", "DirectionDistance");
  _RECO_TRACKMATCH_UVMaxSeparationBars = toml::find<double>(data, "Recon", "TrackMatch3D", "UVMaxSeparationBars");

  _RECO_ASTAR_IsGreedy = toml::find<bool> (data, "Recon", "AStar", "IsGreedy");
  _RECO_ASTAR_CostMetric = toml::find<std::string> (data, "Recon", "AStar", "CostMetric");
  _RECO_ASTAR_MergePlaneGap = toml::find<int>(data, "Recon", "AStar", "MergePlaneGap");
  _RECO_ASTAR_MergeBarGap = toml::find<int>(data, "Recon", "AStar", "MergeBarGap");
  _RECO_ASTAR_NeighbourPlaneWindow = toml::find<int>(data, "Recon", "AStar", "NeighbourPlaneWindow");
  _RECO_ASTAR_NeighbourBarWindow = toml::find<int>(data, "Recon", "AStar", "NeighbourBarWindow");
  _RECO_ASTAR_DownstreamGradientTolerance = toml::find<double>(data, "Recon", "AStar", "DownstreamGradientTolerance");
  _RECO_ASTAR_UpstreamGradientTolerance = toml::find<double>(data, "Recon", "AStar", "UpstreamGradientTolerance");

  _RECO_STOPPING_nLastHits = toml::find<int>(data, "Recon", "Stopping", "nLastHits");
  _RECO_STOPPING_EnergyCut = toml::find<double>(data, "Recon", "Stopping", "EnergyCut");

  _RECO_TRACK_METHOD  = toml::find<std::string>(data, "Recon", "TrackMethod");
  _RECO_CLUSTERING    = toml::find<bool>  (data, "Recon", "Clustering");

  _RECO_KALMAN_RUN = toml::find<bool>(data, "Recon", "Kalman", "Run");
  _RECO_KALMAN_ASSUMED_CHARGE = toml::find<double>(data, "Recon", "Kalman", "Assumed_Charge");
  _RECO_KALMAN_SLOPE_SEED_LAYERS = toml::find<int>(data, "Recon", "Kalman", "SlopeSeedLayers");
  _RECO_KALMAN_MAX_SEED_X_SLOPE = toml::find<double>(data, "Recon", "Kalman", "MaxSeedXSlope");
  _RECO_KALMAN_MAX_SEED_Y_SLOPE = toml::find<double>(data, "Recon", "Kalman", "MaxSeedYSlope");
  _RECO_KALMAN_SAME_Z_TOLERANCE_MM = toml::find<double>(data, "Recon", "Kalman", "SameZToleranceMm");
  _RECO_KALMAN_INITIAL_MOMENTUM_SEED = toml::find<double>(data, "Recon", "Kalman", "InitialMomentumSeed");
  _RECO_KALMAN_INITIAL_COVARIANCE_PATH_LENGTH_THRESHOLD_MM = toml::find<double>(data, "Recon", "Kalman", "InitialCovariancePathLengthThresholdMm");
  _RECO_KALMAN_INITIAL_COVARIANCE_X = toml::find<double>(data, "Recon", "Kalman", "InitialCovarianceX");
  _RECO_KALMAN_INITIAL_COVARIANCE_Y = toml::find<double>(data, "Recon", "Kalman", "InitialCovarianceY");
  _RECO_KALMAN_INITIAL_COVARIANCE_X_SLOPE = toml::find<double>(data, "Recon", "Kalman", "InitialCovarianceXSlope");
  _RECO_KALMAN_INITIAL_COVARIANCE_Y_SLOPE = toml::find<double>(data, "Recon", "Kalman", "InitialCovarianceYSlope");
  _RECO_KALMAN_INITIAL_COVARIANCE_Q_OVER_P = toml::find<double>(data, "Recon", "Kalman", "InitialCovarianceQoverP");
  _RECO_KALMAN_UNMEASURED_COORDINATE_SIGMA_MM = toml::find<double>(data, "Recon", "Kalman", "UnmeasuredCoordinateSigmaMm");

  if (_RECO_EXTRAPOLATION_EndpointSeedHits < 2) {
    throw std::runtime_error("Recon.Extrapolation.EndpointSeedHits must be at least 2");
  }
  if (_RECO_HOUGH_ContainmentHalfWidthX <= 0.0 || _RECO_HOUGH_ContainmentHalfWidthY <= 0.0 ||
      _RECO_EXTRAPOLATION_ContainmentWidthScaleX <= 0.0 || _RECO_EXTRAPOLATION_ContainmentWidthScaleY <= 0.0) {
    throw std::runtime_error("Reco containment-width settings must be positive");
  }
  if (_RECO_KALMAN_SLOPE_SEED_LAYERS < 1) {
    throw std::runtime_error("Recon.Kalman.SlopeSeedLayers must be at least 1");
  }
  if (_RECO_KALMAN_INITIAL_MOMENTUM_SEED <= 0.0) {
    throw std::runtime_error("Recon.Kalman.InitialMomentumSeed must be positive");
  }

  
  _RECO_CALIBRATION_EnergyCalibration = toml::find<double>  (data, "Recon", "Calibration", "EnergyCalibration");
  
  _GEOMETRY_YMIDDLE = toml::find<double>(data, "Geometry", "YBarMiddle");
  _FIDUCIAL_TMS_START_X = toml::find<double>(data, "Fiducial", "TMS", "Start", "X");
  _FIDUCIAL_TMS_START_Y = toml::find<double>(data, "Fiducial", "TMS", "Start", "Y");
  _FIDUCIAL_TMS_START_Z = toml::find<double>(data, "Fiducial", "TMS", "Start", "Z");
  _FIDUCIAL_TMS_END_X = toml::find<double>(data, "Fiducial", "TMS", "End", "X");
  _FIDUCIAL_TMS_END_Y = toml::find<double>(data, "Fiducial", "TMS", "End", "Y");
  _FIDUCIAL_TMS_END_Z = toml::find<double>(data, "Fiducial", "TMS", "End", "Z");
  _ACTIVE_LAR_START_X = toml::find<double>(data, "Active", "LAr", "Start", "X");
  _ACTIVE_LAR_START_Y = toml::find<double>(data, "Active", "LAr", "Start", "Y");
  _ACTIVE_LAR_START_Z = toml::find<double>(data, "Active", "LAr", "Start", "Z");
  _ACTIVE_LAR_END_X = toml::find<double>(data, "Active", "LAr", "End", "X");
  _ACTIVE_LAR_END_Y = toml::find<double>(data, "Active", "LAr", "End", "Y");
  _ACTIVE_LAR_END_Z = toml::find<double>(data, "Active", "LAr", "End", "Z");
  
  _LAR_FIDUCIAL_DOWNSTREAM_Z_CUT = toml::find<double>(data, "Fiducial", "LAr", "DownstreamZCut");
  _LAR_FIDUCIAL_XY_CUT = toml::find<double>(data, "Fiducial", "LAr", "XYCut");
  _LAR_OUTER_SHELL_THICKNESS = toml::find<double>(data, "Fiducial", "LAr", "OuterShellThickness");
  _ND_PHYSICS_MUON_LAR_SHELL_CUT_ENERGY = toml::find<double>(data, "Fiducial", "LAr", "CutEnergy");

  _TRUTH_LIGHTWEIGHT = toml::find<bool> (data, "Truth", "LightWeight");

  _APPLICATIONS_DrawPDF =  toml::find<bool> (data, "Applications", "DrawPDF");
  
  _APPLICATIONS_MaximumNEvents = toml::find<int>(data, "Applications", "MaximumNEvents");
  
  _NERSC_SPILL_PERIOD = 1.2e9;

  _GEOMETRY_git_tag = toml::find<std::string>(data, "Geometry", "GitTag");
  _GEOMETRY_git_branch = toml::find<std::string>(data, "Geometry", "GitBranch");
  _GEOMETRY_git_commit = toml::find<std::string>(data, "Geometry", "GitCommit");
  _GEOMETRY_NumberOfScintillatorPlanes = toml::find<int>(data, "Geometry", "NumberOfScintillatorPlanes");
  _GEOMETRY_NumberOfSteelPlatesThin = toml::find<int>(data, "Geometry", "NumberOfSteelPlatesThin");
  _GEOMETRY_NumberOfSteelPlatesThick = toml::find<int>(data, "Geometry", "NumberOfSteelPlatesThick");
  _GEOMETRY_NumberOfSteelPlatesDouble = toml::find<int>(data, "Geometry", "NumberOfSteelPlatesDouble");
  _META_version_major = toml::find<int>(data, "Meta", "Version", "Major");
  _META_version_minor = toml::find<int>(data, "Meta", "Version", "Minor");
  _META_version_patch = toml::find<int>(data, "Meta", "Version", "Patch");
}
