#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <filesystem>
#include <math.h>       /* atan2 */
#define TAU 6.283185307179586

// Root specific
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVector3.h>

#include "Line_Candidates.h"
#include "Truth_Info.h"
#include "Reco_Tree.h"

#define IS_WITHIN(x, center, tolerance) (std::abs((x) - (center)) <= (tolerance))

int GetHitLocationCodeSingle(float x, bool isx) {
  bool zero = IS_WITHIN(x, 0, 1);
  bool is999 = IS_WITHIN(x, -999, 1) || IS_WITHIN(x, -9999, 1) || IS_WITHIN(x, -99999, 1) || IS_WITHIN(x, -999999, 1) || IS_WITHIN(x, -999999, 1);
  bool crazy_small = x < -4000;
  bool ltTMS = (isx) ? (x < -3300.0) : (x < 11362.0);
  bool gtTMS = (isx) ? (x > 3300.0) : (x > 18314.0);
  bool xlooksz = (isx) ? IS_WITHIN(x, 18314, 10) : false;
  int out = -99999;
  if (zero) out = 0;
  else if (is999) out = -2;
  else if (crazy_small) out = -3;
  else if (ltTMS) out = -1;
  else if (gtTMS) out = 1;
  else if (xlooksz) out = 3;
  else out = 2;
  return out;
}

int isTMSContained(TVector3 position, bool thin_only = false) {
  int out = 0;
    // Z positions of the first hits of the TMS
  const double TMS_Thin_Start = 11362;
  // Where do we transition to the thick region (first layer of scintillator before the change)
  const double TMS_Thick_Start = 13500;
  // Where does the thick region end
  const double TMS_Thick_End = 18294;
  const double TMS_Start_Bars_Only[] = {-3350, 240, TMS_Thin_Start};
  const double TMS_End_Bars_Only[] = {3350, -2950, TMS_Thick_End};
  if (position.X() < TMS_Start_Bars_Only[0]) out += 1 << 0;
  if (position.X() > TMS_End_Bars_Only[0]) out += 1 << 0;
  if (position.Y() > TMS_Start_Bars_Only[1]) out += 1 << 1;
  if (position.Y() < TMS_End_Bars_Only[1]) out += 1 << 1;
  if (position.Z() < TMS_Start_Bars_Only[2]) out += 1 << 2;
  if (thin_only) { if (position.Z() > TMS_Thick_Start) out += 1 << 2; }
  else { if (position.Z() > TMS_End_Bars_Only[2]) out  += 1 << 2; }
  return out;
}

Long64_t PrimaryLoop(Truth_Info& truth, Reco_Tree& reco, Line_Candidates& lc, int numEvents, TFile& outputFile) {
    // List all the hists here
    // Make sure to save them too
    
    bool has_kalman = reco.HasBranch("KalmanPos");
    std::cout<<"has_kalman status: "<<has_kalman<<std::endl;
    
    // Want to plot:
    // Reco vs RecoTrackTrueHitPosition
    // Reco track length vs:
    // RecoTrackPrimaryParticleTrueTrackLengthAsMeasured
    // RecoTrackPrimaryParticleTrueTrackLengthRecoStart
    // RecoTrackPrimaryParticleTrueTrackLengthInTMS
    // Using RecoTrackPrimaryParticleTMSFiducialStart/Touch and RecoTrackPrimaryParticleTMSFiducialEnd
    // Or RecoTrackPrimaryParticleLArFiducialStart, RecoTrackPrimaryParticleLArFiducialTouch, and RecoTrackPrimaryParticleLArFiducialEnd
    
    // For particles ending in TMS, TMSFiducialStart, TMSFiducialEnd
    // RecoTrackPrimaryParticleTrueMomentumEnteringTMS
    
    // Secondary particle info:
    // Flag if it has one or not
    // RecoTrackSecondaryParticleTrueVisibleEnergy
    // RecoTrackPrimaryParticleTrueVisibleEnergy
    
    // Reco tree vars:
    // StartPos, EndPos, nHits, TrackHitPos/KalmanPos
    // EnergyRange, EnergyDeposit, Length
    
    const double plane_pitch = 37;
    const int nz = (int)(19000 / (plane_pitch * 10));
    const int nzfine = (int)(4000 / (plane_pitch));
    
    // First draw some hists for validation
    TH1D hist_ntracks("n_tracks", "N Tracks;N Tracks;N Events", 
      10, -0.5, 9.5);
    TH1D hist_nhits("n_hits", "N Hits;N Hits;N Tracks", 
      25, 0, 100);
    TH1D hist_energy_range("energy_range", "Energy Range;Energy Range (MeV);N Tracks", 
      20, 0, 5000);
    TH1D hist_energy_deposit("energy_deposit", "Energy Deposit;Energy Deposit (MeV);N Tracks", 
      20, 0, 5000);
    TH1D hist_track_length("track_length", "Track Length;Track Areal Density (g/cm^2);N Tracks", 
      20, 0, 3000);
    TH2D hist_track_hits_xz("track_hits_xz", "Track Hit X vs Z;Z (mm);X (mm); N Hits",
      nz, 0, 19000, 100, -5000, 5000);
    TH2D hist_track_hits_yz("track_hits_yz", "Track Hit Y vs Z;Z (mm);Y (mm); N Hits",
      nz, 0, 19000, 100, -6000, 1500);
    TH2D hist_track_hits_xy("track_hits_xy", "Track Hit X vs Y;X (mm);Y (mm); N Hits",
      100, -5000, -5000, 100, -6000, 1500);
    TH2D hist_track_start_pos_xz("track_start_pos_xz", "Track Start Pos X vs Z;Z (mm);X (mm); N Hits",
      nz, 0, 19000, 100, -5000, 5000);
    TH2D hist_track_end_pos_xz("track_end_pos_xz", "Track Start Pos X vs Z;Z (mm);X (mm); N Hits",
      nz, 0, 19000, 100, -5000, 5000);
      
    // Now some reco-only validation
    TH2D hist_track_hits_xz_relative_to_start("track_hits_xz_relative_to_start", "Track Hit X vs Z, Relative to Start;Z (mm);X (mm); N Hits",
      nzfine, -2000, 2000, 100, -4000, 4000);
    TH2D hist_track_hits_xz_relative_to_third("track_hits_xz_relative_to_third", "Track Hit X vs Z, Relative to Third Hit;Z (mm);X (mm); N Hits",
      nzfine, -2000, 2000, 100, -4000, 4000);
    TH2D hist_track_hits_xz_relative_to_end("track_hits_xz_relative_to_end", "Track Hit X vs Z, Relative to End;Z (mm);X (mm); N Hits",
      nzfine, -2000, 2000, 100, -4000, 4000);
    TH2D hist_track_hits_xz_relative_to_third_last("track_hits_xz_relative_to_third_last", 
      "Track Hit X vs Z, Relative to Third Last Hit;Z (mm);X (mm); N Hits",
      nzfine, -2000, 2000, 100, -4000, 4000);
      
    TH1D hist_track_hits_z_relative_to_start("track_hits_z_relative_to_start", "Track Hit Z, Relative to Start;Z (mm); N Hits",
      nzfine, -2000, 2000);
    TH1D hist_track_hits_z_relative_to_third("track_hits_z_relative_to_third", "Track Hit Z, Relative to Third Hit;Z (mm); N Hits",
      nzfine, -2000, 2000);
      
    TH2D hist_comp_track_dz_vs_n_hits("comp_track_dz_vs_n_hits", 
      "Track dz vs n hits;Track End - Start Z (mm);N Hits; N",
      nzfine*2, 0, 8000, 120, 0, 120);
    TH2D hist_comp_track_dz_vs_length("comp_track_dz_vs_length", 
      "Track dz vs Areal Density;Track End - Start Z (mm);Track Areal Density (g/cm^2); N",
      nzfine*2, 0, 8000, 100, 0, 4000);
      
    TH1D hist_fiducial_hit_outside("fiducial_hit_outside", "Track Hit is Outside TMS Fiducial;Outside Fiducial (1 -> yes, 0 -> no); N Hits",
      2, -0.5, 1.5);
    TH1D hist_fiducial_coordinate_outside("fiducial_coordinate_outside", 
      "Track Hit Coordinate Outside Fiducial;Coordinate + Outside * 0.5 (x,y,z=1,2,3, 0.5 -> yes, 0 -> no); N Hits",
      6, 0.75,  3.75);
      
    TH1D hist_truthcomp_hit_res_x("truthcomp_hit_res_x", 
      "Track Hit X Resolution;Reco - True X (mm);N Hits",
      101, -500, 500);
    TH1D hist_truthcomp_hit_res_y("truthcomp_hit_res_y", 
      "Track Hit Y Resolution;Reco - True Y (mm);N Hits",
      101, -500, 500);
    TH1D hist_truthcomp_hit_res_z("truthcomp_hit_res_z", 
      "Track Hit Z Resolution;Reco - True Z (mm);N Hits",
      101, -500, 500);
      
    TH2D hist_truthcomp_hit_error_x("truthcomp_hit_error_x", 
      "Track Hit True vs Reco X;True X (mm);Reco X (mm);N Hits",
      100, -5000, 5000, 100, -5000, 5000);
    TH2D hist_truthcomp_hit_error_y("truthcomp_hit_error_y", 
      "Track Hit True vs Reco Y;True Y (mm);Reco Y (mm);N Hits",
      100, -5000, 3000, 100, -5000, 3000);
    TH2D hist_truthcomp_hit_error_z("truthcomp_hit_error_z", 
      "Track Hit True vs Reco Z;True Z (mm);Reco Z (mm);N Hits",
      100, 0, 19000, 100, 0, 19000);
      
    TH1D hist_endpoint_resolution_x("endpoint_resolution_x", "Track X Endpoint Resolution (TMS-ending muons only);Reco - True X (mm);N Tracks", 
      100, -1000, 1000);
    TH1D hist_endpoint_resolution_z("endpoint_resolution_z", "Track Z Endpoint Resolution (TMS-ending muons only);Reco - True Z (mm);N Tracks", 
      100, -1000, 1000);
    TH1D hist_endpoint_resolution_r_2d("endpoint_resolution_r_2d", "Track XZ Endpoint Resolution (TMS-ending muons only);dr (mm);N Tracks", 
      100, 0, 1000);
      
    TH2D hist_endpoint_error_x("endpoint_error_x", "Track X Endpoint Error Matrix (TMS-ending muons only);True X (mm);Reco X (mm);N Tracks", 
      100, -5000, 5000, 100, -5000, 5000);
    TH2D hist_endpoint_error_z("endpoint_error_z", "Track Z Endpoint Error Matrix (TMS-ending muons only);True Z (mm);Reco Z (mm);N Tracks", 
      100, 12000, 19000, 100, 12000, 19000);
    TH2D hist_endpoint_error_z_including_secondaries("endpoint_error_z_including_secondaries", "Track Z Endpoint Error Matrix (TMS-ending muons only, including those with secondary particles in reco track);True Z (mm);Reco Z (mm);N Tracks", 
      100, 12000, 19000, 100, 12000, 19000);
    TH2D hist_endpoint_error_z_all("endpoint_error_z_all", "Track Z Endpoint Error Matrix (all TMS-ending particles);True Z (mm);Reco Z (mm);N Tracks", 
      100, 12000, 19000, 100, 12000, 19000);
    TH2D hist_endpoint_error_z_all_with_secondary("endpoint_error_z_all_with_secondary", "Track Z Endpoint Error Matrix (all TMS-ending particles, reco track can have secondary true particle);True Z (mm);Reco Z (mm);N Tracks", 
      100, 12000, 19000, 100, 12000, 19000);
    TH2D hist_endpoint_error_z_all_using_second_z("endpoint_error_z_all_using_second_z", "Track Z Endpoint Error Matrix (all TMS-ending particles, using secondary's z if greater);True Z (mm);Reco Z (mm);N Tracks", 
      100, 12000, 19000, 100, 12000, 19000);
      
      // For matching to LAr, we care about the xz direction resolution, the x position resolution, y exiting information
      // and occupancy information
      // Mostly for LAr-starting, tms-ending muons
    TH1D hist_matching_angle_resolution("matching_angle_resolution",  
      "Matching Angle Resolution (LAr-start, TMS-ending muons only);XZ Angle Reco - True (deg);N Tracks", 
      101, -60, 60);
    TH1D hist_matching_angle_true("matching_angle_true",  
      "True Angle, TMS First Plane (LAr-start, TMS-ending muons only);XZ Angle True (deg);N Tracks", 
      101, -60, 60);
    TH1D hist_matching_angle_reco("matching_angle_reco",  
      "Reco Angle, TMS First Plane (LAr-start, TMS-ending muons only);XZ Angle Reco (deg);N Tracks", 
      101, -60, 60);
    TH1D hist_matching_direction_resolution("matching_direction_resolution",  
      "Matching Direction Resolution (LAr-start, TMS-ending muons only);XZ Direction Reco - True;N Tracks", 
      101, -1, 1);
    TH1D hist_matching_x_position_resolution("matching_x_position_resolution",  
      "Matching X Position Resolution (LAr-start, TMS-ending muons only);Track X Position Reco - True (mm);N Tracks", 
      101, -250, 250);
    TH1D hist_matching_y_position_resolution("matching_y_position_resolution",  
      "Matching Y Position Resolution (LAr-start, TMS-ending muons only);Track Y Position Reco - True (mm);N Tracks", 
      101, -500, 500);
    TH1D hist_matching_z_position_resolution("matching_z_position_resolution",  
      "Matching Z Position Resolution (LAr-start, TMS-ending muons only);Track Z Position Reco - True (mm);N Tracks", 
      101, -300, 300);
    TH1D hist_matching_z_position_reco("matching_z_position_reco",  
      "Matching Reco Z Position (LAr-start, TMS-ending muons only);Track Z Position Reco (mm);N Tracks", 
      101, 11000, 16000);
    TH1D hist_matching_z_position_true("matching_z_position_true",  
      "Matching True Z Position (LAr-start, TMS-ending muons only);Track Z Position True (mm);N Tracks", 
      101, 11000, 16000);
    TH2D hist_matching_xz_position_true("matching_xz_position_true",  
      "Matching True XZ Position (LAr-start, TMS-ending muons only);Track Z Position True (mm);Track X Position True (mm);N Tracks", 
      101, 11000, 16000, 101, -4000, 4000);
    // TODO check for track leaving using reco Y. Also check for reco in first two planes cut, and occupancy > some amount cut
      
    // TODO make reco eff for muons reco'd in TMS vs not

    Long64_t entry_number = 0;
    // Now loop over the ttree
    for ( ; entry_number < truth.GetEntriesFast() && entry_number < reco.GetEntriesFast()\
      && (numEvents < 0 || entry_number < numEvents); entry_number++) {
      if (entry_number % 10000 == 0) std::cout<<"On entry: "<<entry_number<<std::endl;
      
      // Get the current entry
      // Currently reco and truth match
      truth.GetEntry(entry_number);
      reco.GetEntry(entry_number);
      lc.GetEntry(entry_number);
      
      // Adjust to kalman if needed
      if (!has_kalman) {
        for (int itrack = 0; itrack < reco.nTracks; itrack++) {
          for (int ihit = 0; ihit < reco.nHits[itrack]; ihit++) {
            for (int i = 0; i < 3; i++) {
              reco.KalmanPos[itrack][ihit][i] = reco.TrackHitPos[itrack][ihit][i];
              reco.KalmanTruePos[itrack][ihit][i] = truth.RecoTrackTrueHitPosition[itrack][ihit][i];
            }
          }
        }
      }
      
      // TODO calculate plane number and then check per plane occupancy
      // Also related is total visible energy so compare that to true value
      
      hist_ntracks.Fill(reco.nTracks);
      for (int itrack = 0; itrack < reco.nTracks; itrack++) {
        hist_nhits.Fill(reco.nHits[itrack]);
        hist_energy_range.Fill(reco.EnergyRange[itrack]);
        hist_energy_deposit.Fill(reco.EnergyDeposit[itrack]);
        hist_track_length.Fill(reco.Length[itrack]);
        
        // Matching
        //std::cout<<"truth.RecoTrackPrimaryParticleLArFiducialStart[itrack]: "<<truth.RecoTrackPrimaryParticleLArFiducialStart[itrack]<<",\ttruth.RecoTrackPrimaryParticleTMSFiducialEnd[itrack]: "<<truth.RecoTrackPrimaryParticleTMSFiducialEnd[itrack]<<std::endl;
        if (truth.RecoTrackPrimaryParticleLArFiducialStart[itrack] && truth.RecoTrackPrimaryParticleTMSFiducialEnd[itrack]) {
          double direction_reco = reco.Direction[itrack][0] / reco.Direction[itrack][2];
          double direction_true = truth.RecoTrackPrimaryParticleTrueMomentumTrackStart[itrack][0] / truth.RecoTrackPrimaryParticleTrueMomentumTrackStart[itrack][2];
          double angle_reco = atan2(direction_reco, 1) * 360 / TAU;
          double angle_true = atan2(direction_true, 1) * 360 / TAU;
          hist_matching_angle_resolution.Fill(angle_reco - angle_true);
          hist_matching_angle_true.Fill(angle_true);
          hist_matching_angle_reco.Fill(angle_reco);
          hist_matching_direction_resolution.Fill(direction_reco - direction_true);
          double track_x_reco = reco.StartPos[itrack][0];
          double track_y_reco = reco.StartPos[itrack][1];
          double track_z_reco = reco.StartPos[itrack][2];
          /*
          double track_x_true = truth.RecoTrackPrimaryParticleTruePositionEnteringTMS[itrack][0];
          double track_y_true = truth.RecoTrackPrimaryParticleTruePositionEnteringTMS[itrack][1];
          double track_z_true = truth.RecoTrackPrimaryParticleTruePositionEnteringTMS[itrack][2]; */
          int particle_index = truth.RecoTrackPrimaryParticleIndex[itrack];
          double track_x_true = truth.PositionZIsTMSStart[particle_index][0];
          double track_y_true = truth.PositionZIsTMSStart[particle_index][1];
          double track_z_true = truth.PositionZIsTMSStart[particle_index][2];
          hist_matching_x_position_resolution.Fill(track_x_reco - track_x_true);
          hist_matching_y_position_resolution.Fill(track_y_reco - track_y_true);
          hist_matching_z_position_resolution.Fill(track_z_reco - track_z_true);
          hist_matching_z_position_reco.Fill(track_z_reco);
          hist_matching_z_position_true.Fill(track_z_true);
          hist_matching_xz_position_true.Fill(track_z_true, track_x_true);
        }
        
        //std::cout<<"start xyz "<<itrack<<": "<<reco.StartPos[itrack][0]<<","<<reco.StartPos[itrack][1]<<","<<reco.StartPos[itrack][2]<<std::endl;
        //std::cout<<"end xyz "<<itrack<<":   "<<reco.EndPos[itrack][0]<<","<<reco.EndPos[itrack][1]<<","<<reco.EndPos[itrack][2]<<std::endl;
        hist_track_start_pos_xz.Fill(reco.StartPos[itrack][2], reco.StartPos[itrack][0]);
        hist_track_end_pos_xz.Fill(reco.EndPos[itrack][2], reco.EndPos[itrack][0]);
        for (int ihit = 0; ihit < reco.nHits[itrack]; ihit++) {
          hist_track_hits_xz.Fill(reco.KalmanPos[itrack][ihit][2], reco.KalmanPos[itrack][ihit][0]);
          hist_track_hits_yz.Fill(reco.KalmanPos[itrack][ihit][2], reco.KalmanPos[itrack][ihit][1]);
          hist_track_hits_xy.Fill(reco.KalmanPos[itrack][ihit][0], reco.KalmanPos[itrack][ihit][1]);
          
          int last_index = reco.nHits[itrack] - 1;
          double z_relative_start = reco.KalmanPos[itrack][ihit][2] - reco.KalmanPos[itrack][0][2];
          double x_relative_start = reco.KalmanPos[itrack][ihit][0] - reco.KalmanPos[itrack][0][0];
          double z_relative_end = reco.KalmanPos[itrack][ihit][2] - reco.KalmanPos[itrack][last_index][2];
          double x_relative_end = reco.KalmanPos[itrack][ihit][0] - reco.KalmanPos[itrack][last_index][0];
          hist_track_hits_xz_relative_to_start.Fill(z_relative_start, x_relative_start);
          hist_track_hits_xz_relative_to_end.Fill(z_relative_end, x_relative_end);
          hist_track_hits_z_relative_to_start.Fill(z_relative_start);
          if (reco.nHits[itrack] >= 3) {
            double z_relative_start_m3 = reco.KalmanPos[itrack][ihit][2] - reco.KalmanPos[itrack][2][2];
            double x_relative_start_m3 = reco.KalmanPos[itrack][ihit][0] - reco.KalmanPos[itrack][2][0];
            double z_relative_end_m3 = reco.KalmanPos[itrack][ihit][2] - reco.KalmanPos[itrack][last_index - 2][2];
            double x_relative_end_m3 = reco.KalmanPos[itrack][ihit][0] - reco.KalmanPos[itrack][last_index - 2][0];
            hist_track_hits_xz_relative_to_third.Fill(z_relative_start_m3, x_relative_start_m3);
            hist_track_hits_xz_relative_to_third_last.Fill(z_relative_end_m3, x_relative_end_m3);
            hist_track_hits_z_relative_to_third.Fill(z_relative_start_m3);
          }
            
          TVector3 hit_position(reco.KalmanPos[itrack][ihit][0], reco.KalmanPos[itrack][ihit][1], reco.KalmanPos[itrack][ihit][2]);
          int hit_outside_fiducial = isTMSContained(hit_position);
          // Fill with 0 if isTMSContained returns 0 (no coordinate is outside fiducial), 1 otherwise
          hist_fiducial_hit_outside.Fill((hit_outside_fiducial == 0) ? 0 : 1);
          bool x_is_out = hit_outside_fiducial & (1 << 0);
          bool y_is_out = hit_outside_fiducial & (1 << 1);
          bool z_is_out = hit_outside_fiducial & (1 << 2);
          hist_fiducial_coordinate_outside.Fill(1 + ((x_is_out) ? 0.5 : 0));
          hist_fiducial_coordinate_outside.Fill(2 + ((y_is_out) ? 0.5 : 0));
          hist_fiducial_coordinate_outside.Fill(3 + ((z_is_out) ? 0.5 : 0));
          
          // These hists only make sense if the true hit position is within the TMS, so remove any that are z=0
          if (reco.KalmanTruePos[itrack][ihit][2] > 100) {
            hist_truthcomp_hit_res_x.Fill(reco.KalmanPos[itrack][ihit][0] - reco.KalmanTruePos[itrack][ihit][0]);
            hist_truthcomp_hit_res_y.Fill(reco.KalmanPos[itrack][ihit][1] - reco.KalmanTruePos[itrack][ihit][1]);
            hist_truthcomp_hit_res_z.Fill(reco.KalmanPos[itrack][ihit][2] - reco.KalmanTruePos[itrack][ihit][2]);
            hist_truthcomp_hit_error_x.Fill(reco.KalmanTruePos[itrack][ihit][0], reco.KalmanPos[itrack][ihit][0]);
            hist_truthcomp_hit_error_y.Fill(reco.KalmanTruePos[itrack][ihit][1], reco.KalmanPos[itrack][ihit][1]);
            hist_truthcomp_hit_error_z.Fill(reco.KalmanTruePos[itrack][ihit][2], reco.KalmanPos[itrack][ihit][2]);
          }
          
        } // end for loop over hits
        
        // dz is a decent proxy of the length of the track and/or the number of planes
        // so it should both be proportional to areal density
        // as well as n hits
        // So something that's off diagonal in n hits means that it may have skipped a hit
        double dz = reco.EndPos[itrack][2] - reco.StartPos[itrack][2];
        hist_comp_track_dz_vs_n_hits.Fill(dz, reco.nHits[itrack]);
        hist_comp_track_dz_vs_length.Fill(dz, reco.Length[itrack]);
        
        // https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
        bool ismuon = std::abs(truth.RecoTrackPrimaryParticlePDG[itrack]) == 13;
        bool no_secondary = truth.RecoTrackSecondaryParticlePDG[itrack] == -1000000000;
        bool tms_ending = truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][2] < 18315.7 && truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][2] > 13728.7 && std::abs(truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][0]) < 4446.96; //truth.RecoTrackPrimaryParticleTMSFiducialEnd[itrack];
        //std::cout<<"x reco: "<<reco.EndPos[itrack][0]<<", x truth: "<<truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][0]<<std::endl;
        //std::cout<<"z reco: "<<reco.EndPos[itrack][2]<<", z truth: "<<truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][2]<<std::endl;
        //std::cout<<"truth.RecoTrackPrimaryParticleTMSFiducialEnd[itrack]: "<<truth.RecoTrackPrimaryParticleTMSFiducialEnd[itrack]<<std::endl;
        //std::cout<<"truth.RecoTrackPrimaryParticlePDG[itrack]: "<<truth.RecoTrackPrimaryParticlePDG[itrack]<<std::endl;
        //std::cout<<"truth.RecoTrackSecondaryParticlePDG[itrack]: "<<truth.RecoTrackSecondaryParticlePDG[itrack]<<std::endl;
        if (tms_ending && ismuon && no_secondary) { // && no_secondary && tms_ending) {
          double x = reco.EndPos[itrack][0] - truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][0];
          double z = reco.EndPos[itrack][2] - truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][2];
          double r_2d = std::sqrt(x * x + z * z);
          //std::cout<<"x: "<<x<<", z: "<<z<<", r_2d: "<<r_2d<<std::endl;
          hist_endpoint_resolution_r_2d.Fill(r_2d);
          hist_endpoint_resolution_x.Fill(x);
          hist_endpoint_resolution_z.Fill(z);
          hist_endpoint_error_z.Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][2], reco.EndPos[itrack][2]);
          hist_endpoint_error_x.Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][0], reco.EndPos[itrack][0]);
        }
        if (tms_ending && ismuon) hist_endpoint_error_z_including_secondaries.Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][2], reco.EndPos[itrack][2]);
        if (tms_ending && no_secondary) hist_endpoint_error_z_all.Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][2], reco.EndPos[itrack][2]);
        if (tms_ending) hist_endpoint_error_z_all_with_secondary.Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][2], reco.EndPos[itrack][2]);
        if (tms_ending) { 
          double true_z = truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][2];
          if (!no_secondary && truth.RecoTrackSecondaryParticleTruePositionEnd[itrack][2] > true_z) true_z = truth.RecoTrackSecondaryParticleTruePositionEnd[itrack][2];
          //std::cout<<truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][2]<<", "<<truth.RecoTrackSecondaryParticleTruePositionEnd[itrack][2]<<std::endl;
          hist_endpoint_error_z_all_using_second_z.Fill(true_z, reco.EndPos[itrack][2]);
        }
      }
    } // End for loop over entries
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

int main(int argc, char* argv[]) {
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
    // All these have a large memory footprint, especially Line_Candidates
    // by declaring them static, it moves it from the stack to the heap, which has more memory allocated
    // Otherwise, we get a confusing seg fault before main starts
    // See https://stackoverflow.com/questions/20253267/segmentation-fault-before-main
    static Truth_Info ti(truth); 
    static Reco_Tree ri(reco);
    std::cout<<"About to load Line_Candidates"<<std::endl; 
    static Line_Candidates li(line_candidates);
    std::cout<<"Loaded Line_Candidates"<<std::endl;

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
    
    PrimaryLoop(ti, ri, li, numEvents, outputFile);
    
    // Close the output file
    outputFile.Close();

    std::cout << "Output file created: " << outputFilename << std::endl;

    return 0;
}
