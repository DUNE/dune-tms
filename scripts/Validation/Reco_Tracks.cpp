#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <filesystem>

// Root specific
#include <TFile.h>
#include <TH1D.h>
#include <TVector3.h>

#include "Truth_Info.h"
#include "Reco_Tree.h"

bool isTMSContained(TVector3 position, bool thin_only = false) {
  bool out = true;
    // Z positions of the first hits of the TMS
  const double TMS_Thin_Start = 11362;
  // Where do we transition to the thick region (first layer of scintillator before the change)
  const double TMS_Thick_Start = 13500;
  // Where does the thick region end
  const double TMS_Thick_End = 18294;
  const double TMS_Start_Bars_Only[] = {-3350, 240, TMS_Thin_Start};
  const double TMS_End_Bars_Only[] = {3350, -2950, TMS_Thick_End};
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

Long64_t PrimaryLoop(Truth_Info& truth, Reco_Tree& reco, int numEvents, TFile& outputFile) {
    // List all the hists here
    // Make sure to save them too
    int n_zbins = 50;
    double z_start = 11000;
    double z_end = 18500;
    TH1D hist_reco_track_z("reco_track_z", "Reco Track Z_{Start};Z (mm);N Tracks", 
      n_zbins, z_start, z_end);
    TH1D hist_true_track_z("true_track_z", "True Track Z_{Start};Z (mm);N Tracks", 
      n_zbins, z_start, z_end);
    double dz = 400;
    TH1D hist_true_track_z_very_front("true_track_z_very_front", "True Track Z_{Start};Z (mm);N Tracks", 
      n_zbins, 11362 - dz, 11362 + dz);
    TH1D hist_true_track_z_very_back("true_track_z_very_back", "True Track Z_{Start};Z (mm);N Tracks", 
      n_zbins, 18294 - dz, 18400 + dz);
    TH1D hist_true_track_z_end("true_track_z_end", "True Track Z_{End};Z (mm);N Tracks", 
      n_zbins, z_start, z_end);
    TH1D hist_true_track_z_end_very_back("true_track_z_end_very_back", "True Track Z_{End};Z (mm);N Tracks", 
      n_zbins, 18294 - dz, 18400 + dz);
    TH1D hist_diff_track_z("diff_track_z", "Reco - True Track Z_{Start};dZ (mm);N Tracks", 
      n_zbins, -300, 300);
    int n_xbins = 50;
    double x_start = -4000;
    double x_end = 4000;
    TH1D hist_reco_track_x("reco_track_x", "Reco Track X_{Start};X (mm);N Tracks", 
      n_xbins, x_start, x_end);
    TH1D hist_true_track_x("true_track_x", "True Track X_{Start};X (mm);N Tracks", 
      n_xbins, x_start, x_end);
    TH1D hist_diff_track_x("diff_track_x", "Reco - True Track X_{Start};dX (mm);N Tracks", 
      n_xbins, -300, 300);
    int n_ybins = 50;
    double y_start = -3500;
    double y_end = 500;
    TH1D hist_reco_track_y("reco_track_y", "Reco Track Y_{Start};Y (mm);N Tracks", 
      n_ybins, y_start, y_end);
    TH1D hist_true_track_y("true_track_y", "True Track Y_{Start};Y (mm);N Tracks", 
      n_ybins, y_start, y_end);
    TH1D hist_diff_track_y("diff_track_y", "Reco - True Track Y_{Start};dY (mm);N Tracks", 
      n_ybins, -900, 900);

    int n_direction_bins = 50;
    double direction_start = -1;
    double direction_end = 1;
    TH1D hist_reco_track_direction("reco_track_direction", "Reco Track Direction dx/dz_{Start};Direction (mm);N Tracks", 
      n_direction_bins, direction_start, direction_end);
    TH1D hist_true_track_direction("true_track_direction", "True Track Direction dx/dz_{Start};Direction (mm);N Tracks", 
      n_direction_bins, direction_start, direction_end);
    TH1D hist_diff_track_direction("diff_track_direction", "Reco - True Track Direction dx/dz_{Start};Delta Direction (mm);N Tracks", 
      n_direction_bins, -0.5, 0.5);

    Long64_t entry_number = 0;
    // Now loop over the ttree
    for ( ; entry_number < truth.GetEntriesFast() && entry_number < reco.GetEntriesFast() \
      && (numEvents < 0 || entry_number < numEvents); entry_number++) {
      if (entry_number % 10000 == 0) std::cout<<"On entry: "<<entry_number<<std::endl;
      
      // Get the current entry
      // Currently reco and truth match
      truth.GetEntry(entry_number);
      reco.GetEntry(entry_number);
      
      for (int itrack = 0; itrack < truth.RecoTrackN; itrack++) {
        int particle_index = truth.RecoTrackPrimaryParticleIndex[itrack];
        int pdg = -999999999;
        if (particle_index >= 0 && particle_index < truth.nTrueParticles) pdg = truth.PDG[particle_index];
        //else std::cout<<"Found unusual RecoTrackPrimaryParticleIndex: "<<particle_index<<std::endl;
        auto end_array = truth.RecoTrackPrimaryParticleTruePositionTrackEnd[itrack];
        TVector3 endpoint(end_array[0], end_array[1], end_array[2]);
        auto start_array = truth.RecoTrackPrimaryParticleTruePositionTrackStart[itrack];
        TVector3 startpoint(start_array[0], start_array[1], start_array[2]);
        
        bool startLAr = truth.RecoTrackPrimaryParticleLArFiducialStart[itrack];
        bool endTMS = truth.RecoTrackPrimaryParticleTMSFiducialEnd[itrack];
        
        if (startpoint.Z() < 11362.0) std::cout<<"Found case outside bounds: "<<startpoint.Z()<<std::endl;
        
        bool is_muon = std::abs(pdg) == 13;
        if (is_muon && startLAr) {
          // Stuff we only want filled if the primary particle that makes up a reco track is a muon
          hist_reco_track_z.Fill(reco.StartPos[itrack][2]);
          hist_true_track_z.Fill(truth.RecoTrackPrimaryParticleTruePositionTrackStart[itrack][2]);
          hist_true_track_z_very_front.Fill(truth.RecoTrackPrimaryParticleTruePositionTrackStart[itrack][2]);
          hist_true_track_z_very_back.Fill(truth.RecoTrackPrimaryParticleTruePositionTrackStart[itrack][2]);
          hist_true_track_z_end.Fill(endpoint.Z());
          hist_true_track_z_end_very_back.Fill(endpoint.Z());
          
          hist_diff_track_z.Fill(reco.StartPos[itrack][2] - truth.RecoTrackPrimaryParticleTruePositionTrackStart[itrack][2]);
          hist_reco_track_x.Fill(reco.StartPos[itrack][0]);
          hist_true_track_x.Fill(truth.RecoTrackPrimaryParticleTruePositionTrackStart[itrack][0]);
          hist_diff_track_x.Fill(reco.StartPos[itrack][0] - truth.RecoTrackPrimaryParticleTruePositionTrackStart[itrack][0]);
          hist_reco_track_y.Fill(reco.StartPos[itrack][1]);
          hist_true_track_y.Fill(truth.RecoTrackPrimaryParticleTruePositionTrackStart[itrack][1]);
          hist_diff_track_y.Fill(reco.StartPos[itrack][1] - truth.RecoTrackPrimaryParticleTruePositionTrackStart[itrack][1]);
          
          double reco_direction = reco.Direction[itrack][0] / reco.Direction[itrack][2];
          hist_reco_track_direction.Fill(reco_direction);
          double true_direction = truth.RecoTrackPrimaryParticleTrueMomentumTrackStart[itrack][0] / truth.RecoTrackPrimaryParticleTrueMomentumTrackStart[itrack][2];
          hist_true_track_direction.Fill(true_direction);
          hist_diff_track_direction.Fill(reco_direction - true_direction);
        }
      }
    }
    // Now save the hists
    outputFile.Write();
    
    auto entries_visited = entry_number;
    
    std::cout<<"Finished loop over "<<entries_visited<<" entries"<<std::endl;
    return entries_visited;
}

std::string getExecutableName(const char* arg) {
    std::string executableName = arg; // Convert to std::string
    size_t lastSlash = executableName.find_last_of("/\\"); // Find last occurrence of '/' or '\'
    if (lastSlash != std::string::npos) {
        executableName = executableName.substr(lastSlash + 1); // Extract substring after the last slash
    }
    size_t extensionPos = executableName.find_last_of('.'); // Find last occurrence of '.'
    if (extensionPos != std::string::npos) {
        executableName = executableName.substr(0, extensionPos); // Remove extension
    }
    return executableName;
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
    if (missing_ttree) {
      throw std::runtime_error("Missing one or more ttree from file");
    }
    Truth_Info ti(truth);
    Reco_Tree ri(reco);

    std::string exeName = getExecutableName(argv[0]);
    std::string directoryPath ="/exp/dune/data/users/" + std::string(getenv("USER")) + "/dune-tms/Validation/" + exeName + "/";

    if (createDirectory(directoryPath)) {
        std::cout << "Directory created: " << directoryPath << std::endl;
    } else {
        std::cerr << "Failed to create directory" << std::endl;
        exit(1);
    }
    
    // Create output filename
    std::string outputFilename = directoryPath + getOutputFilename(inputFilename);

    // Create TFile with the output filename
    TFile outputFile(outputFilename.c_str(), "RECREATE");
    
    PrimaryLoop(ti, ri, numEvents, outputFile);
    
    // Close the output file
    outputFile.Close();

    std::cout << "Output file created: " << outputFilename << std::endl;

    return 0;
}
