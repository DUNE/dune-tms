#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <filesystem>

// Root specific
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVector3.h>

#include "Truth_Info.h"
#include "Reco_Tree.h"
#include "Line_Candidates.h"

bool isTMSContained(TVector3 position, bool thin_only = false) {
  bool out = true;
    // Z positions of the first hits of the TMS
  const double TMS_Thin_Start = 11185;
  // Where do we transition to the thick region (first layer of scintillator before the change)
  const double TMS_Thick_Start = 14435;
  // Where do we transition to the double region
  const double TMS_Double_Start = 17495;
  // Where does the thick region end
  const double TMS_Double_End = 18535;
  const double TMS_Start_Bars_Only[] = {-3350, 240, TMS_Thin_Start};
  const double TMS_End_Bars_Only[] = {3350, -2950, TMS_Double_End};
  if (position.X() < TMS_Start_Bars_Only[0]) out = false;
  if (position.X() > TMS_End_Bars_Only[0]) out = false;
  if (position.Y() > TMS_Start_Bars_Only[1]) out = false;
  if (position.Y() < TMS_End_Bars_Only[1]) out = false;
  if (position.Z() < TMS_Start_Bars_Only[2]) out = false;
  if (thin_only) { if (position.Z() > TMS_Thick_Start) out = false; }
  else { if (position.Z() > TMS_End_Bars_Only[2]) out = false; }
  //std::cout<<position.X()<<", "<<position.Y()<<", "<<position.Z()<<": "<<out<<std::endl;
  return out;
}

Long64_t PrimaryLoop(Truth_Info& truth, Reco_Tree& reco, Line_Candidates& lc, int numEvents, TFile& outputFile) {
    // List all the hists here
    // Make sure to save them too
    int n_areal_density_bins = 50;
    double areal_density_max = 3000;
    int n_length_bins = 50;
    double length_max = 6000;
    TH1D hist_muon_true_track_length("muon_true_track_length", "True Areal Density;True Areal Density (g/cm2);N Muons", 
      n_areal_density_bins, 0, areal_density_max);
    TH1D hist_muon_true_track_length_ignore_y("muon_true_track_length_ignore_y", 
      "True Areal Density, Ignoring Y;True Areal Density (g/cm2);N Muons", n_areal_density_bins, 0, areal_density_max);
    TH1D hist_muon_true_track_length_ratio("muon_true_track_length_ratio", 
      "True Areal Density Ratio w/ and w/o Y;Ratio (ignore y/include y);N Muons", n_areal_density_bins, 0, 2);
      
    // Comparisons to reco tracks
    TH1D hist_reco_track_length("reco_track_length", "Reco Areal Density;Reco Areal Density (g/cm2);N Muons", 
      n_areal_density_bins, 0, areal_density_max);
    TH1D hist_true_track_length("true_track_length", "True Areal Density;True Areal Density (g/cm2);N Muons", 
      n_areal_density_bins, 0, areal_density_max);
    TH1D hist_true_track_length_original("true_track_length_original", "True Areal Density;True Areal Density (g/cm2);N Muons", 
      n_areal_density_bins, 0, areal_density_max);
    TH2D hist_reco_vs_true_track_length("reco_vs_true_track_length", 
      "Areal Density Smearing Matrix;True Areal Density (g/cm2);Reco Areal Density (g/cm2);N Muons", 
      n_areal_density_bins, 0, areal_density_max, n_areal_density_bins, 0, areal_density_max);
    TH2D hist_reco_u_vs_true_track_length("reco_u_vs_true_track_length", 
      "Areal Density Smearing Matrix;True Areal Density (g/cm2);Reco U Areal Density (g/cm2);N Muons", 
      n_areal_density_bins, 0, areal_density_max, n_areal_density_bins, 0, areal_density_max);
    TH2D hist_reco_v_vs_true_track_length("reco_v_vs_true_track_length", 
      "Areal Density Smearing Matrix;True Areal Density (g/cm2);Reco V Areal Density (g/cm2);N Muons", 
      n_areal_density_bins, 0, areal_density_max, n_areal_density_bins, 0, areal_density_max);
    TH2D hist_reco_v_vs_reco_track_length("reco_v_vs_reco_track_length", 
      "Areal Density Reco Comparison;Reco Areal Density (g/cm2);Reco V Areal Density (g/cm2);N Muons", 
      n_areal_density_bins, 0, areal_density_max, n_areal_density_bins, 0, areal_density_max);
    TH2D hist_reco_v_vs_reco_u_track_length("reco_v_vs_reco_u_track_length", 
      "Areal Density Reco Comparison;Reco U Areal Density (g/cm2);Reco V Areal Density (g/cm2);N Muons", 
      n_areal_density_bins, 0, areal_density_max, n_areal_density_bins, 0, areal_density_max);
    TH2D hist_reco_vs_true_track_length_original("reco_vs_true_track_length_original", 
      "Areal Density Smearing Matrix;True Areal Density (g/cm2);Reco Areal Density (g/cm2);N Muons", 
      n_areal_density_bins, 0, areal_density_max, n_areal_density_bins, 0, areal_density_max);
    TH2D hist_true_track_length_vs_density("true_track_length_vs_density", 
      "True Track Length Vs Areal Density;True Track Length (mm);True Areal Density (g/cm2);N Muons", 
      n_length_bins, 0, length_max, n_areal_density_bins, 0, areal_density_max);
    TH2D hist_true_track_length_vs_density_original("true_track_length_vs_density_original", 
      "True Track Length Vs Areal Density;True Track Length (mm);True Areal Density (g/cm2);N Muons", 
      n_length_bins, 0, length_max, n_areal_density_bins, 0, areal_density_max);
    TH2D hist_reco_track_length_vs_density("reco_track_length_vs_density", 
      "Reco Track Length Vs Areal Density;Reco Track Length (mm);Reco Areal Density (g/cm2);N Muons", 
      n_length_bins, 0, length_max, n_areal_density_bins, 0, areal_density_max);
      
    TH2D hist_reco_vs_true_track_actual_length("reco_vs_true_track_actual_length", 
      "Track Length Smearing Matrix;True Track Length (mm);Reco Track Length (mm);N Muons", 
      n_length_bins, 0, length_max, n_length_bins, 0, length_max);
      
    int n_points = 100;
    int sample_rate = 50;// 10;
    double y_min = -4000;
    double y_max = 1000;
    TH2D pseudo_event_yz_display_reco("pseudo_event_yz_display_reco", 
      "All Reco Track Hits;Reco Hit Z (mm);Reco Hit Y (mm);N Hits", 
      n_points, 11000, 19000, n_points, y_min, y_max);
    TH2D pseudo_event_yz_display_true("pseudo_event_yz_display_true", 
      "All True Track Hits;True Hit Z (mm);True Hit Y (mm);N Hits", 
      n_points, 11000, 19000, n_points, y_min, y_max);
      
    double x_max = 4000;
    TH2D pseudo_event_xz_display_reco("pseudo_event_xz_display_reco", 
      "All Reco Track Hits;Reco Hit Z (mm);Reco Hit Y (mm);N Hits", 
      n_points, 11000, 19000, n_points, -x_max, x_max);
    TH2D pseudo_event_xz_display_true("pseudo_event_xz_display_true", 
      "All True Track Hits;True Hit Z (mm);True Hit Y (mm);N Hits", 
      n_points, 11000, 19000, n_points, -x_max, x_max);
      
    int n_diff_bins = 101;
    TH2D hist_hit_reco_vs_true_x("hit_reco_vs_true_x", 
      "Track Hit Position Comparison X;True Hit Position X (mm);Reco Hit Position X (mm);N Hits", 
      n_points, -x_max, x_max, n_points, -x_max, x_max);
    TH2D hist_hit_reco_vs_true_y("hit_reco_vs_true_y", 
      "Track Hit Position Comparison Y;True Hit Position Y (mm);Reco Hit Position Y (mm);N Hits", 
      n_points, y_min, y_max, n_points, y_min, y_max);
    TH2D hist_hit_reco_vs_true_z("hit_reco_vs_true_z", 
      "Track Hit Position Comparison Z;True Hit Position Z (mm);Reco Hit Position Z (mm);N Hits", 
      n_points, 11000, 19000, n_points, 11000, 19000);
    TH1D hist_hit_x_diff("hit_x_diff", 
      "Track Hit Position Comparison X;Reco - True Hit Position X (mm);N Hits", 
      n_diff_bins, -1000, 1000);
    TH1D hist_hit_y_diff("hit_y_diff", 
      "Track Hit Position Comparison Y;Reco - True Hit Position Y (mm);N Hits", 
      n_diff_bins, -1000, 1000);
    TH1D hist_hit_z_diff("hit_z_diff", 
      "Track Hit Position Comparison Z;Reco - True Hit Position Z (mm);N Hits", 
      n_diff_bins, -1000, 1000);

    Long64_t entry_number = 0;
    // Now loop over the ttree
    for ( ; entry_number < truth.GetEntriesFast() && entry_number < reco.GetEntriesFast() \
      && (numEvents < 0 || entry_number < numEvents); entry_number++) {
      if (entry_number % 10000 == 0) std::cout<<"On entry: "<<entry_number<<std::endl;
      
      // Get the current entry
      // Currently reco and truth match
      truth.GetEntry(entry_number);
      reco.GetEntry(entry_number);
      lc.GetEntry(entry_number);
      
      for (int itrack = 0; itrack < truth.RecoTrackN; itrack++) {
        int particle_index = truth.RecoTrackPrimaryParticleIndex[itrack];
        int pdg = -999999999;
        if (particle_index >= 0 && particle_index < truth.nTrueParticles) pdg = truth.PDG[particle_index];
        //else std::cout<<"Found unusual RecoTrackPrimaryParticleIndex: "<<particle_index<<std::endl;
        auto end_array = truth.RecoTrackPrimaryParticleTruePositionTrackEnd[itrack];
        TVector3 endpoint(end_array[0], end_array[1], end_array[2]);
        auto start_array = truth.RecoTrackPrimaryParticleTruePositionTrackStart[itrack];
        TVector3 startpoint(start_array[0], start_array[1], start_array[2]);
        
        bool is_muon = std::abs(pdg) == 13;
        bool is_TMS_contained = isTMSContained(endpoint);
        bool is_TMS_contained_start = isTMSContained(startpoint);
        if (is_muon && is_TMS_contained && is_TMS_contained_start) {
          // Stuff we only want filled if the primary particle that makes up a reco track is a muon
          hist_muon_true_track_length.Fill(truth.RecoTrackPrimaryParticleTrueTrackLength[itrack]);
          hist_muon_true_track_length_ignore_y.Fill(truth.RecoTrackPrimaryParticleTrueTrackLengthIgnoreY[itrack]);
          double ratio = truth.RecoTrackPrimaryParticleTrueTrackLengthIgnoreY[itrack] / truth.RecoTrackPrimaryParticleTrueTrackLength[itrack];
          if (ratio > 1) std::cout<<"ratio: "<<ratio<<", "<<truth.RecoTrackPrimaryParticleTrueTrackLengthIgnoreY[itrack]<<\
            ", "<<truth.RecoTrackPrimaryParticleTrueTrackLength[itrack]<<std::endl;
          hist_muon_true_track_length_ratio.Fill(ratio);
          
          hist_reco_track_length.Fill(reco.Length[itrack]);
          hist_reco_vs_true_track_length.Fill(truth.RecoTrackPrimaryParticleTrueTrackLength[itrack], reco.Length[itrack]);
          hist_reco_u_vs_true_track_length.Fill(truth.RecoTrackPrimaryParticleTrueTrackLength[itrack], lc.TrackLengthU[itrack]);
          hist_reco_v_vs_true_track_length.Fill(truth.RecoTrackPrimaryParticleTrueTrackLength[itrack], lc.TrackLengthV[itrack]);
          hist_reco_v_vs_reco_track_length.Fill(reco.Length[itrack], lc.TrackLengthV[itrack]);
          hist_reco_v_vs_reco_u_track_length.Fill(lc.TrackLengthU[itrack], lc.TrackLengthV[itrack]);
          
          
          hist_reco_vs_true_track_length_original.Fill(truth.Muon_TrueTrackLength, reco.Length[itrack]);
          hist_true_track_length.Fill(truth.RecoTrackPrimaryParticleTrueTrackLength[itrack]);
          hist_true_track_length_original.Fill(truth.Muon_TrueTrackLength);
          
          TVector3 start(truth.Muon_Vertex[0], truth.Muon_Vertex[1], truth.Muon_Vertex[2]);
          TVector3 end(truth.Muon_Death[0], truth.Muon_Death[1], truth.Muon_Death[2]);
          double true_length = (end - start).Mag();
          hist_true_track_length_vs_density.Fill(true_length, truth.RecoTrackPrimaryParticleTrueTrackLength[itrack]);
          hist_true_track_length_vs_density_original.Fill(true_length, truth.Muon_TrueTrackLength);
          
          TVector3 reco_start(reco.StartPos[itrack][0], reco.StartPos[itrack][1], reco.StartPos[itrack][2]);
          TVector3 reco_end(reco.EndPos[itrack][0], reco.EndPos[itrack][1], reco.EndPos[itrack][2]);
          double reco_length = (reco_end - reco_start).Mag();
          hist_reco_track_length_vs_density.Fill(reco_length, reco.Length[itrack]);
          
          hist_reco_vs_true_track_actual_length.Fill(true_length, reco_length);
          
          if (entry_number % sample_rate == 0) {
            for (int h = 0; h < reco.nHits[itrack]; h++) {
              double x = reco.TrackHitPos[itrack][h][0];
              double y = reco.TrackHitPos[itrack][h][1];
              double z = reco.TrackHitPos[itrack][h][2];
              pseudo_event_yz_display_reco.Fill(z, y);
              pseudo_event_xz_display_reco.Fill(z, x);
            }
            for (int h = 0; h < truth.RecoTrackNHits[itrack]; h++) {
              double x = truth.RecoTrackTrueHitPosition[itrack][h][0];
              double y = truth.RecoTrackTrueHitPosition[itrack][h][1];
              double z = truth.RecoTrackTrueHitPosition[itrack][h][2];
              pseudo_event_yz_display_true.Fill(z, y);
              pseudo_event_xz_display_true.Fill(z, x);
              //std::cout<<"itrack: "<<itrack<<", h: "<<h<<", out of n hits: "<<truth.RecoTrackNHits[itrack]<<", xyz: "<<x<<","<<y<<","<<z<<std::endl;
            }
            for (int h = 0; h < truth.RecoTrackNHits[itrack]; h++) {
              double rx = reco.TrackHitPos[itrack][h][0];
              double ry = reco.TrackHitPos[itrack][h][1];
              double rz = reco.TrackHitPos[itrack][h][2];
              double tx = truth.RecoTrackTrueHitPosition[itrack][h][0];
              double ty = truth.RecoTrackTrueHitPosition[itrack][h][1];
              double tz = truth.RecoTrackTrueHitPosition[itrack][h][2];
              hist_hit_reco_vs_true_x.Fill(tx, rx);
              hist_hit_reco_vs_true_y.Fill(ty, ry);
              hist_hit_reco_vs_true_z.Fill(tz, rz);
              hist_hit_x_diff.Fill(rx - tx);
              hist_hit_y_diff.Fill(ry - ty);
              hist_hit_z_diff.Fill(rz - tz);
            }
          }
        }
      }
    }
    // Now save the hists
    outputFile.Write();
    
    auto entries_visited = entry_number;
    
    std::cout<<"Finished loop over "<<entries_visited<<" entries"<<std::endl;
    return entries_visited;
}

bool createDirectory(const std::string& path) {
    try {
        // Create directory and its parents if they don't exist
        std::filesystem::create_directories(path);
        return true;
    } catch(const std::exception& e) {
        std::cerr << "Error creating directory: " << e.what() << std::endl;
        return false;
    }
}

std::string getOutputFilename(const std::string& inputFilename) {
    // Find the position of the last occurrence of '/' in the input filename
    size_t pos = inputFilename.find_last_of('/');
    // Extract the filename without the directory structure
    std::string filename = (pos != std::string::npos) ? inputFilename.substr(pos + 1) : inputFilename;
    return filename;
}

int main(int argc, char* argv[]) {;
    // Check if the correct number of arguments is provided
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_filename> <num_events (optional)>" << std::endl;
        return 1;
    }

    // Extract input filename and number of events from command line arguments
    std::string inputFilename = argv[1];
    int numEvents = -1;
    if (argc > 2) numEvents = atoi(argv[2]);

    // Load the tree and make the Truth_Info object
    TFile TF(inputFilename.c_str());
    TTree* truth = (TTree*) TF.Get("Truth_Info");
    TTree* reco = (TTree*) TF.Get("Reco_Tree");
    TTree* line_candidates = (TTree*) TF.Get("Line_Candidates");
    bool missing_ttree = false;
    if (!truth) { 
      std::string message = inputFilename + " doesn't contain Truth_Info";
      std::cerr<<message<<std::endl;
      missing_ttree = true;
    }
    if (!reco) { 
      std::string message = inputFilename + " doesn't contain Reco_Tree";
      std::cerr<<message<<std::endl;
      missing_ttree = true;
    }
    if (!line_candidates) { 
      std::string message = inputFilename + " doesn't contain Line_Candidates";
      std::cerr<<message<<std::endl;
      missing_ttree = true;
    }
    if (missing_ttree) {
      throw std::runtime_error("Missing one or more ttree from file");
    }
    Truth_Info ti(truth);
    Reco_Tree ri(reco);
    Line_Candidates li(line_candidates);

    std::string directoryPath ="/exp/dune/data/users/"+ std::string(getenv("USER")) + "/dune-tms/Validation/Plot_True_Track_Length/";

    if (createDirectory(directoryPath)) {
        std::cout << "Directory created: " << directoryPath << std::endl;
    } else {
        std::cerr << "Failed to create directory" << std::endl;
    }
    
    // Create output filename
    std::string outputFilename = directoryPath + getOutputFilename(inputFilename);

    // Create TFile with the output filename
    TFile outputFile(outputFilename.c_str(), "RECREATE");
    
    PrimaryLoop(ti, ri, li, numEvents, outputFile);

    // Close the output file
    outputFile.Close();

    std::cout << "Output file created: " << outputFilename << std::endl;

    return 0;
}
