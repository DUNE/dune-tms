// Stage 2.5 diagnostic: runs the REAL TMS_SpacePointDBScan + TMS_SpacePointCluster
// (PCA) classes on the space points of one already-reconstructed real event,
// pulled directly out of an existing ConvertToTMSTree output file's Reco_Tree
// (SpacePointX/Y/Z/Time). This is not part of the eventual production pipeline
// (that's Stage 4) -- it's a standalone look at real clustering results before
// committing to Stage 3/4, dumping per-point cluster labels to CSV for external
// plotting.
//
// TMS_SpacePointDBScan/TMS_SpacePointCluster operate purely on
// TMS_SpacePoint::GetX()/GetY()/GetZ() (confirmed by reading their source) --
// hit indices are never dereferenced during clustering/PCA -- so it's valid to
// reconstruct TMS_SpacePoint objects here from just the tree's saved
// coordinates, with dummy hit indices, rather than re-deriving the full event
// from the original edep-sim file.
//
// TMS_SpacePointDBScan's anisotropic neighbor test needs each point's real
// z-ordered plane index and the real bar pitch (its quantization unit for the
// transverse tolerance), both of which require the actual surveyed geometry
// (see its header comment) -- so this tool now also loads the original
// edep-sim spill file's TGeoManager (same pattern as app/ShootRay.cpp), purely
// to derive these via TMS_Geom::GetPlaneIndexNearestZ()/GetMaxBarPitch().

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TGeoManager.h"
#include "TTree.h"

#include "TMS_Geom.h"
#include "TMS_SpacePoint.h"
#include "TMS_SpacePointCluster.h"
#include "TMS_SpacePointDBScan.h"

int main(int argc, char **argv) {
  if (argc != 9) {
    std::cerr << "Usage: " << argv[0]
              << " <edep_sim_geom_file> <input_reco_tree.root> <entry_number> <base_transverse_bars>"
                 " <transverse_bars_per_plane_gap> <max_plane_gap> <min_points> <output.csv>"
              << std::endl;
    return -1;
  }

  const std::string geom_filename = argv[1];
  const std::string input_filename = argv[2];
  const long long entry_number = std::stoll(argv[3]);
  const int base_transverse_bars = std::stoi(argv[4]);
  const int transverse_bars_per_plane_gap = std::stoi(argv[5]);
  const int max_plane_gap = std::stoi(argv[6]);
  const unsigned int min_points = std::stoul(argv[7]);
  const std::string output_csv = argv[8];

  TFile geom_input(geom_filename.c_str());
  if (geom_input.IsZombie()) {
    std::cerr << "Failed to open geometry source file: " << geom_filename << std::endl;
    return -1;
  }
  TGeoManager *geom = (TGeoManager *)geom_input.Get("EDepSimGeometry");
  if (!geom) {
    std::cerr << "Geometry source file is missing 'EDepSimGeometry': " << geom_filename << std::endl;
    return -1;
  }
  TMS_Geom::GetInstance().SetGeometry(geom);
  const double max_plane_pitch = TMS_Geom::GetInstance().GetMaxPlanePitch();
  const double bar_pitch = TMS_Geom::GetInstance().GetMaxBarPitch();
  if (max_plane_pitch <= 0 || bar_pitch <= 0) {
    std::cerr << "TMS_Geom found fewer than 2 surveyed planes or bars -- cannot derive a clustering tolerance."
              << std::endl;
    return -1;
  }
  const double worst_case_transverse = (base_transverse_bars + max_plane_gap * transverse_bars_per_plane_gap) * bar_pitch;
  const double broad_phase_radius =
      std::sqrt(worst_case_transverse * worst_case_transverse +
                std::pow(max_plane_pitch * (max_plane_gap + 1), 2));

  TFile input(input_filename.c_str());
  if (input.IsZombie()) {
    std::cerr << "Failed to open input file: " << input_filename << std::endl;
    return -1;
  }
  TTree *reco_tree = (TTree *)input.Get("Reco_Tree");
  if (!reco_tree) {
    std::cerr << "Input file is missing the required 'Reco_Tree' tree: " << input_filename << std::endl;
    return -1;
  }
  if (entry_number < 0 || entry_number >= reco_tree->GetEntries()) {
    std::cerr << "Entry " << entry_number << " out of range (tree has " << reco_tree->GetEntries()
              << " entries)" << std::endl;
    return -1;
  }

  const int kMaxSpacePoints = 10000;  // matches __TMS_MAX_SPACEPOINTS__
  int n_space_points = 0;
  std::vector<float> sp_x(kMaxSpacePoints), sp_y(kMaxSpacePoints), sp_z(kMaxSpacePoints),
      sp_time(kMaxSpacePoints);
  reco_tree->SetBranchAddress("nSpacePoints", &n_space_points);
  reco_tree->SetBranchAddress("SpacePointX", sp_x.data());
  reco_tree->SetBranchAddress("SpacePointY", sp_y.data());
  reco_tree->SetBranchAddress("SpacePointZ", sp_z.data());
  reco_tree->SetBranchAddress("SpacePointTime", sp_time.data());
  reco_tree->GetEntry(entry_number);

  std::cout << "Entry " << entry_number << ": " << n_space_points << " space points" << std::endl;

  std::vector<TMS_SpacePoint> space_points;
  space_points.reserve(n_space_points);
  for (int i = 0; i < n_space_points; ++i) {
    space_points.emplace_back(sp_x[i], sp_y[i], sp_z[i], /*x_idx=*/-1, /*y_idx=*/-1, sp_time[i]);
  }

  // A space point's Z is always exactly some plane's Z by BuildSpacePoints()'s
  // own convention (confirmed earlier this session), so the nearest-Z lookup
  // recovers the real plane index without needing the original hit objects.
  std::vector<int> plane_index;
  plane_index.reserve(n_space_points);
  for (int i = 0; i < n_space_points; ++i) {
    plane_index.push_back(TMS_Geom::GetInstance().GetPlaneIndexNearestZ(sp_z[i]));
  }

  TMS_SpacePointDBScan dbscan(space_points, plane_index, min_points, bar_pitch, base_transverse_bars,
                               transverse_bars_per_plane_gap, max_plane_gap, broad_phase_radius);
  std::vector<std::vector<int>> cluster_indices = dbscan.RunAndGetClusterIndices();

  std::vector<TMS_SpacePointCluster> clusters;
  clusters.reserve(cluster_indices.size());
  for (auto &indices : cluster_indices) {
    clusters.emplace_back(space_points, indices);
  }

  std::cout << "Found " << clusters.size() << " clusters (bar_pitch=" << bar_pitch << "mm, base_transverse_bars="
            << base_transverse_bars << ", transverse_bars_per_plane_gap=" << transverse_bars_per_plane_gap
            << ", max_plane_gap=" << max_plane_gap << ", min_points=" << min_points << ")" << std::endl;

  // Per-point cluster id: 0 = noise, 1..N = cluster index (1-based, matching
  // TMS_SpacePointDBScan's own convention).
  std::vector<int> point_cluster_id(n_space_points, 0);
  for (size_t c = 0; c < cluster_indices.size(); ++c) {
    for (int idx : cluster_indices[c]) point_cluster_id[idx] = static_cast<int>(c) + 1;
  }

  const double kLinearityThreshold = 0.8;
  const size_t kMinClusterSizeForTrack = 5;
  for (size_t c = 0; c < clusters.size(); ++c) {
    const auto &cl = clusters[c];
    std::cout << "  Cluster " << (c + 1) << ": n=" << cl.GetSize() << " linearity=" << cl.GetLinearity()
              << " IsTrackLike=" << (cl.IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack) ? "yes" : "no")
              << std::endl;
  }

  std::ofstream csv(output_csv);
  csv << "x,y,z,time,cluster_id,cluster_linearity,is_track_like\n";
  for (int i = 0; i < n_space_points; ++i) {
    int cid = point_cluster_id[i];
    double linearity = 0.0;
    bool is_track_like = false;
    if (cid > 0) {
      const auto &cl = clusters[cid - 1];
      linearity = cl.GetLinearity();
      is_track_like = cl.IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack);
    }
    csv << sp_x[i] << "," << sp_y[i] << "," << sp_z[i] << "," << sp_time[i] << "," << cid << ","
        << linearity << "," << (is_track_like ? 1 : 0) << "\n";
  }
  csv.close();
  std::cout << "Wrote " << output_csv << std::endl;

  // PCA axes for every track-like cluster, one row per cluster: centroid plus
  // all 3 eigenvalue/eigenvector pairs (descending), for plotting each axis
  // through the centroid with length ~3*sqrt(eigenvalue), same convention the
  // user draws these with elsewhere.
  const std::string pca_output_csv =
      (output_csv.size() > 4 && output_csv.compare(output_csv.size() - 4, 4, ".csv") == 0)
          ? output_csv.substr(0, output_csv.size() - 4) + "_pca.csv"
          : output_csv + "_pca.csv";
  std::ofstream pca_csv(pca_output_csv);
  pca_csv << "cluster_id,cx,cy,cz,eval0,ex0,ey0,ez0,eval1,ex1,ey1,ez1,eval2,ex2,ey2,ez2\n";
  for (size_t c = 0; c < clusters.size(); ++c) {
    const auto &cl = clusters[c];
    if (!cl.IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack)) continue;
    const auto &centroid = cl.GetCentroid();
    const auto &eigenvalues = cl.GetEigenvalues();
    const auto &eigenvectors = cl.GetEigenvectors();
    pca_csv << (c + 1) << "," << centroid[0] << "," << centroid[1] << "," << centroid[2];
    for (int rank = 0; rank < 3; ++rank) {
      pca_csv << "," << eigenvalues[rank] << "," << eigenvectors[rank][0] << "," << eigenvectors[rank][1] << ","
              << eigenvectors[rank][2];
    }
    pca_csv << "\n";
  }
  pca_csv.close();
  std::cout << "Wrote " << pca_output_csv << std::endl;

  return 0;
}
