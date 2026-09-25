// Whole-spill variant of RealEventClusterDisplay (Stage 2.5 diagnostic): runs the
// REAL TMS_SpacePointDBScan + TMS_SpacePointCluster (PCA) classes independently on
// EVERY already-reconstructed time-slice belonging to one spill, pulled directly out
// of an existing ConvertToTMSTree output file's Reco_Tree (SpacePointX/Y/Z/Time,
// SliceNo, SpillNo). Clustering is run separately per slice, exactly as the real
// per-slice pipeline does -- a single combined DBSCAN pass across the whole spill
// would risk spatially merging two unrelated tracks from different time slices that
// happen to sit at the same z/transverse position, since TMS_SpacePointDBScan has no
// time dimension. Results from every slice are then concatenated into one CSV, with
// a globally-unique cluster id (offset per slice) and the originating SliceNo/entry
// recorded per row, so all slices of a spill can be plotted together.
//
// See RealEventClusterDisplay.cpp for the geometry-loading rationale (loaded once
// here, not once per slice, since a spill can have 100+ slices).

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
#include "TMS_LayerGrouping.h"
#include "SpacePointLayerInput.h"

int main(int argc, char **argv) {
  if (argc != 9) {
    std::cerr << "Usage: " << argv[0]
              << " <edep_sim_geom_file> <input_reco_tree.root> <spill_number> <max_dz_mm>"
                 " <base_transverse_mm> <transverse_per_dz> <min_points> <output.csv>\n"
                 "  (DBSCAN tolerance, TMS_SpacePointDBScan::Tolerance; defaults 270 <bar pitch> 0.55)"
              << std::endl;
    return -1;
  }

  const std::string geom_filename = argv[1];
  const std::string input_filename = argv[2];
  const int target_spill = std::stoi(argv[3]);
  const double max_dz_mm = std::stod(argv[4]);
  const double base_transverse_mm = std::stod(argv[5]);
  const double transverse_per_dz = std::stod(argv[6]);
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
  const double bar_pitch = TMS_Geom::GetInstance().GetMaxBarPitch();
  if (bar_pitch <= 0) {
    std::cerr << "TMS_Geom found fewer than 2 surveyed bars -- cannot derive a clustering tolerance." << std::endl;
    return -1;
  }
  // DBSCAN tolerance from the command line (mm), as TMS_SpacePointDBScan::Tolerance.
  TMS_SpacePointDBScan::Tolerance dbscan_tolerance;
  dbscan_tolerance.MaxDzMM = max_dz_mm;
  dbscan_tolerance.BaseTransverseMM = base_transverse_mm;
  dbscan_tolerance.TransversePerDzMM = transverse_per_dz;

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

  const int kMaxSpacePoints = 10000;  // matches __TMS_MAX_SPACEPOINTS__
  int n_space_points = 0;
  int spill_no = 0, slice_no = 0;
  double slice_start = 0, slice_end = 0;
  std::vector<float> sp_x(kMaxSpacePoints), sp_y(kMaxSpacePoints), sp_z(kMaxSpacePoints),
      sp_time(kMaxSpacePoints);
  reco_tree->SetBranchAddress("nSpacePoints", &n_space_points);
  reco_tree->SetBranchAddress("SpacePointX", sp_x.data());
  reco_tree->SetBranchAddress("SpacePointY", sp_y.data());
  reco_tree->SetBranchAddress("SpacePointZ", sp_z.data());
  const SpacePointLayerInput sp_layer(reco_tree, kMaxSpacePoints);
  reco_tree->SetBranchAddress("SpacePointTime", sp_time.data());
  reco_tree->SetBranchAddress("SpillNo", &spill_no);
  reco_tree->SetBranchAddress("SliceNo", &slice_no);
  reco_tree->SetBranchAddress("TimeSliceStartTime", &slice_start);
  reco_tree->SetBranchAddress("TimeSliceEndTime", &slice_end);

  const double kLinearityThreshold = 0.8;
  const size_t kMinClusterSizeForTrack = 5;

  std::ofstream csv(output_csv);
  csv << "entry,slice,slice_start,slice_end,x,y,z,time,cluster_id,cluster_linearity,is_track_like\n";

  const std::string pca_output_csv =
      (output_csv.size() > 4 && output_csv.compare(output_csv.size() - 4, 4, ".csv") == 0)
          ? output_csv.substr(0, output_csv.size() - 4) + "_pca.csv"
          : output_csv + "_pca.csv";
  std::ofstream pca_csv(pca_output_csv);
  pca_csv << "entry,slice,slice_start,slice_end,cluster_id,n,cx,cy,cz,eval0,ex0,ey0,ez0,eval1,ex1,ey1,ez1,eval2,ex2,"
             "ey2,ez2\n";

  Long64_t n_entries = reco_tree->GetEntries();
  int global_cluster_offset = 0;
  int n_slices_processed = 0;
  int n_total_clusters = 0;
  int n_track_like_total = 0;
  long n_total_points = 0;

  for (Long64_t entry = 0; entry < n_entries; ++entry) {
    reco_tree->GetEntry(entry);
    if (spill_no != target_spill) continue;
    if (n_space_points <= 0) continue;

    std::vector<TMS_SpacePoint> space_points;
    space_points.reserve(n_space_points);
    for (int i = 0; i < n_space_points; ++i) {
      space_points.emplace_back(sp_x[i], sp_y[i], sp_z[i], /*x_idx=*/-1, /*y_idx=*/-1, sp_time[i],
                              sp_layer.Layer(i, sp_z[i]));
    }


    TMS_SpacePointDBScan dbscan(space_points, min_points, dbscan_tolerance);
    std::vector<std::vector<int>> cluster_indices = dbscan.RunAndGetClusterIndices();

    std::vector<TMS_SpacePointCluster> clusters;
    clusters.reserve(cluster_indices.size());
    for (auto &indices : cluster_indices) {
      clusters.emplace_back(space_points, indices);
    }

    std::vector<int> point_cluster_id(n_space_points, 0);
    for (size_t c = 0; c < cluster_indices.size(); ++c) {
      for (int idx : cluster_indices[c]) point_cluster_id[idx] = static_cast<int>(c) + 1 + global_cluster_offset;
    }

    for (int i = 0; i < n_space_points; ++i) {
      int cid = point_cluster_id[i];
      double linearity = 0.0;
      bool is_track_like = false;
      if (cid > 0) {
        const auto &cl = clusters[cid - 1 - global_cluster_offset];
        linearity = cl.GetLinearity();
        is_track_like = cl.IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack);
      }
      csv << entry << "," << slice_no << "," << slice_start << "," << slice_end << "," << sp_x[i] << "," << sp_y[i]
          << "," << sp_z[i] << "," << sp_time[i] << "," << cid << "," << linearity << ","
          << (is_track_like ? 1 : 0) << "\n";
    }

    for (size_t c = 0; c < clusters.size(); ++c) {
      const auto &cl = clusters[c];
      if (!cl.IsTrackLike(kLinearityThreshold, kMinClusterSizeForTrack)) continue;
      const auto &centroid = cl.GetCentroid();
      const auto &eigenvalues = cl.GetEigenvalues();
      const auto &eigenvectors = cl.GetEigenvectors();
      pca_csv << entry << "," << slice_no << "," << slice_start << "," << slice_end << ","
              << (c + 1 + global_cluster_offset) << "," << cl.GetSize() << "," << centroid[0] << "," << centroid[1]
              << "," << centroid[2];
      for (int rank = 0; rank < 3; ++rank) {
        pca_csv << "," << eigenvalues[rank] << "," << eigenvectors[rank][0] << "," << eigenvectors[rank][1] << ","
                << eigenvectors[rank][2];
      }
      pca_csv << "\n";
      ++n_track_like_total;
    }

    global_cluster_offset += static_cast<int>(clusters.size());
    n_total_clusters += static_cast<int>(clusters.size());
    n_total_points += n_space_points;
    ++n_slices_processed;
    std::cout << "  entry=" << entry << " slice=" << slice_no << " nSP=" << n_space_points << " -> "
              << clusters.size() << " clusters" << std::endl;
  }

  csv.close();
  pca_csv.close();

  std::cout << "Spill " << target_spill << ": " << n_slices_processed << " non-empty slices, " << n_total_points
            << " space points, " << n_total_clusters << " clusters (" << n_track_like_total << " track-like)"
            << std::endl;
  std::cout << "Wrote " << output_csv << " and " << pca_output_csv << std::endl;

  return 0;
}
