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

  // The minimum hits needed to run reconstruction
  _RECO_MinHits = toml::find<int>(data, "Recon", "MinHits");

  _RECO_DBSCAN_MinPoints = toml::find<int>(data, "Recon", "DBSCAN", "MinPoints");
  _RECO_DBSCAN_Epsilon = toml::find<double>(data, "Recon", "DBSCAN", "Epsilon");

  _RECO_HOUGH_MaxHough = toml::find<int>(data, "Recon", "Hough", "MaxTrans");
  _RECO_HOUGH_MinInterp = toml::find<double>(data, "Recon", "Hough", "MinInter");
  _RECO_HOUGH_MaxInterp = toml::find<double>(data, "Recon", "Hough", "MaxInter");
  _RECO_HOUGH_MinSlope = toml::find<double>(data, "Recon", "Hough", "MinSlope");
  _RECO_HOUGH_MaxSlope = toml::find<double>(data, "Recon", "Hough", "MaxSlope");
  _RECO_HOUGH_HitMult = toml::find<double>(data, "Recon", "Hough", "HitMult");
  _RECO_HOUGH_MergeTracks = toml::find<bool>(data, "Recon", "Hough", "MergeTracks");
  _RECO_HOUGH_RunAStar = toml::find<bool>(data, "Recon", "Hough", "RunAStarCleanup");

  _RECO_ASTAR_IsGreedy = toml::find<bool> (data, "Recon", "AStar", "IsGreedy");
  _RECO_ASTAR_CostMetric = toml::find<std::string> (data, "Recon", "AStar", "CostMetric");

  _RECO_TRACK_METHOD  = toml::find<std::string>(data, "Recon", "TrackMethod");
  _RECO_CLUSTERING    = toml::find<bool>  (data, "Recon", "Clustering");

  _TRUTH_LIGHTWEIGHT = toml::find<bool> (data, "Truth", "LightWeight");

  _APPLICATIONS_DrawPDF =  toml::find<bool> (data, "Applications", "DrawPDF");
}
