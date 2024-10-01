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
#include <TCanvas.h>
#include <TStyle.h>
#include <TMarker.h>
#include <TString.h>

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

void DrawSlice(std::string outfilename, std::string reason, std::string message, Reco_Tree& reco, Line_Candidates& lc, Truth_Info& truth) { 
  std::string directoryPath ="/exp/dune/data/users/" + std::string(getenv("USER")) + "/dune-tms/Validation/EventDisplays/" + reason + "/";
  createDirectory(directoryPath);
  static TCanvas* canvas;
  if (canvas == NULL) {
    gStyle->SetOptStat(0);
    gROOT->SetBatch(kTRUE); 
    canvas = new TCanvas();
    canvas->SetLeftMargin(0.12);
    canvas->Divide(1, 2);
  }
  
  TH2D hist("", "Event Display;Z (mm);X (mm);N Hits", 100, 11000, 19000, 100, -4000, 4000);
  TH2D histy("", "Event Display;Z (mm);Y (mm);N Hits", 100, 11000, 19000, 100, -5000, 1000);
  
  
  std::vector<TMarker> markers;
  std::vector<TMarker> markersy;
  // Get approx time slice range
  float hit_time_earliest = 1e9;
  float hit_time_latest = -1e9;
  for (int ih = 0; ih < lc.nHits; ih++) {
    if (lc.RecoHitPos[ih][3] > hit_time_latest) hit_time_latest = lc.RecoHitPos[ih][3];
    if (lc.RecoHitPos[ih][3] < hit_time_earliest) hit_time_earliest = lc.RecoHitPos[ih][3];
    {
      float mx = lc.RecoHitPos[ih][2];
      float my = lc.RecoHitPos[ih][0];
      TMarker marker(mx, my, 21);
      markers.push_back(marker);
    }
    {
      float mx = lc.RecoHitPos[ih][2];
      float my = lc.RecoHitPos[ih][1];
      TMarker marker(mx, my, 21);
      markersy.push_back(marker);
    }
  }
  for (int it = 0; it < reco.nTracks; it++) {
    for (int ih = 0; ih < reco.nKalmanNodes[it]; ih++) {
      {
        float mx = truth.RecoTrackTrueHitPosition[it][ih][2];
        float my = truth.RecoTrackTrueHitPosition[it][ih][0];
        TMarker marker(mx, my, 21);
        marker.SetMarkerColor(kGreen);
        markers.push_back(marker);
      }
      {
        float mx = truth.RecoTrackTrueHitPosition[it][ih][2];
        float my = truth.RecoTrackTrueHitPosition[it][ih][1];
        TMarker marker(mx, my, 21);
        marker.SetMarkerColor(kGreen);
        markersy.push_back(marker);
      }
      {
        float mx = reco.TrackHitPos[it][ih][2];
        float my = reco.TrackHitPos[it][ih][0];
        TMarker marker(mx, my, 20);
        marker.SetMarkerColor(kRed);
        markers.push_back(marker);
      }
      {
        float mx = reco.TrackHitPos[it][ih][2];
        float my = reco.TrackHitPos[it][ih][1];
        TMarker marker(mx, my, 20);
        marker.SetMarkerColor(kRed);
        markersy.push_back(marker);
      }
    }
  }
  if (false) {
    int colors[] = {kGreen, kRed, kBlue, kOrange, kYellow, kMagenta, kGreen + 2};
    int ncolors = sizeof(colors)/sizeof(colors[0]);
    int n_offset = 0;
    for (int it = 0; it < truth.nParticles; it++) {
      // Only draw muons
      if (std::abs(truth.PDG[it]) == 13) {
        {
          float mx = truth.BirthPosition[it][2];
          float my = truth.BirthPosition[it][0];
          if (mx < 11000) {
            mx = truth.PositionZIsTMSStart[it][2];
            my = truth.PositionZIsTMSStart[it][0];
          }
          TMarker marker(mx, my, 22);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          markers.push_back(marker);
        }
        {
          float mx = truth.BirthPosition[it][2];
          float my = truth.BirthPosition[it][1];
          if (mx < 11000) {
            mx = truth.PositionZIsTMSStart[it][2];
            my = truth.PositionZIsTMSStart[it][1];
          }
          TMarker marker(mx, my, 22);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          markersy.push_back(marker);
        }
        {
          float mx = truth.DeathPosition[it][2];
          float my = truth.DeathPosition[it][0];
          if (mx > 19000) {
            mx = truth.PositionZIsTMSEnd[it][2];
            my = truth.PositionZIsTMSEnd[it][0];
          }
          TMarker marker(mx, my, 23);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          markers.push_back(marker);
        }
        {
          float mx = truth.DeathPosition[it][2];
          float my = truth.DeathPosition[it][1];
          if (mx > 19000) {
            mx = truth.PositionZIsTMSEnd[it][2];
            my = truth.PositionZIsTMSEnd[it][1];
          }
          TMarker marker(mx, my, 23);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          markersy.push_back(marker);
        }
        n_offset += 1;
      }
    }
  }
  
  // Now draw
  canvas->cd(1);
  hist.Draw("Axis");
  for (auto& marker : markers) {
    marker.Draw();
  }
  canvas->cd(2);
  histy.Draw("Axis");
  for (auto& marker : markersy) {
    marker.Draw();
  }
  
  canvas->Print((directoryPath + outfilename + ".png").c_str());
}

// If one of these flags is set, then it fails the cut
// set flag: flags |= FOO;
// clear: flags &= ~FOO;
// check that flag is set: if (flags & FOO)
// check that flag is not set: if ((flags & FOO) == 0)
enum TMSCutFlags {
  TMSCutAny =        1 << 1, // Fails any cut
  TMSCutStartFront = 1 << 2, 
  TMSCutAnyEndCuts = 1 << 3, // fails 4, 5, or 6
  TMSCutZCut =       1 << 4,
  TMSCutYCut =       1 << 5,
  TMSCutUVCut =      1 << 6
};

const double TMS_START_Z = 11362;
const double TMS_LAST_PLANE_Z = 17900; // mm // 18314.0 - 50; // End Z - 5cm
const double TMS_Y_TOP = 160.0; // Top of scintillator (I'm 99% sure starting from where U/V overlap)
const double TMS_Y_BOTTOM = -2850.0; // Bottom of scinillator (also measured from U/V overlap)
const double TMS_Y_MIDDLE = 0.5 * (TMS_Y_TOP + TMS_Y_BOTTOM);
const double TMS_X_EXTENT = 3300.0; // Assumes things are symmetrical right now
const double Y_BUFFER = 400; // mm, 40 cm
const double X_BUFFER = 30; // mm, 3cm bar thickness?
const double Z_BUFFER = 100; // mm, two thin planes
const double VIEW_ANGLE = 3.0 * TMath::DegToRad(); 

TVector3 GetRotatedVector(TVector3 v, bool u) {
  double angle = VIEW_ANGLE;
  if (u) angle *= -1; // Opposite direction
  // Step 1: Translate to middle of TMS
  v.SetY(v.Y() - TMS_Y_MIDDLE);
  // Step 2: Rotate the vector around the z-axis
  v.RotateZ(angle);
  // Step 3: Translate the vector back
  v.SetY(v.Y() + TMS_Y_MIDDLE);
  return v;
}

void FillHistWithTMSCutFlagInfo(TH1 &hist, int flags) {
  // Fill passing case
  if (flags == 0) hist.Fill(0);
  // Check and fill for each flag
  if (flags & TMSCutAny) hist.Fill(1);
  if (flags & TMSCutStartFront) hist.Fill(2);
  if (flags & TMSCutAnyEndCuts) hist.Fill(3);
  if (flags & TMSCutZCut) hist.Fill(4);
  if (flags & TMSCutYCut) hist.Fill(5);
  if (flags & TMSCutUVCut) hist.Fill(6);
}

int TMSStartPointThruFront(TVector3 startpoint) {
  // Check that it's in the first two planes or so
  int out = 0;
  if (startpoint.Z() < TMS_START_Z) 
    out |= TMSCutStartFront;
  if (startpoint.Z() > TMS_START_Z + Z_BUFFER) 
    out |= TMSCutStartFront;
    
  const int ANYCUTS = (TMSCutZCut | TMSCutYCut | TMSCutUVCut | TMSCutAnyEndCuts | TMSCutStartFront);
  if ((out & ANYCUTS)) 
    out |= TMSCutAny;
  
  return out;
}

int TMSEndpointUncontained(TVector3 endpoint) {
  // Returns flag representing which cut was missed
  // Zero means it passes all cuts
  int out = 0;
  
  // Check that muon is not in last scinillator plane because it may be exiting then
  if (endpoint.Z() > TMS_LAST_PLANE_Z) 
    out |= TMSCutZCut;
  // Also check that it's contained inside the start of the TMS when using truth info
  if (endpoint.Z() < TMS_START_Z) 
    out |= TMSCutZCut;
  
  // Check that muon is not within 40cm in y of start and bottom
  const double Y = endpoint.Y();
  if (TMS_Y_BOTTOM + Y_BUFFER > Y || Y > TMS_Y_TOP - Y_BUFFER)
    out |= TMSCutYCut;
    
  // Check that x is within hexagon made by U and V views
  // Do it by rotating to the views and checking each one separately
  // Y position is checked by y cut above
  auto endpoint_v = GetRotatedVector(endpoint, false);
  auto endpoint_u = GetRotatedVector(endpoint, true);
  if (std::abs(endpoint_v.X()) > TMS_X_EXTENT - X_BUFFER || std::abs(endpoint_u.X()) > TMS_X_EXTENT - X_BUFFER) 
    out |= TMSCutUVCut;
  
  // Set that all cuts passed if they did
  const int ANYENDCUTS = (TMSCutZCut | TMSCutYCut | TMSCutUVCut);
  if ((out & ANYENDCUTS)) 
    out |= TMSCutAnyEndCuts;
  
  // Set that all cuts passed if they did
  const int ANYCUTS = (TMSCutZCut | TMSCutYCut | TMSCutUVCut | TMSCutAnyEndCuts | TMSCutStartFront);
  if ((out & ANYCUTS)) 
    out |= TMSCutAny;
  return out;
}

int CheckTMSCuts(TVector3 startpoint, TVector3 endpoint) {
  return TMSStartPointThruFront(startpoint) | TMSEndpointUncontained(endpoint);
}

TVector3 calculatePositionAtZ(double z_new, TVector3 start, double xz_dir, double yz_dir = 1.0) {
  double z = z_new;
  double dz = z - start.Z();
  double x = start.X() + xz_dir * dz;
  double y = start.Y() + yz_dir * dz;
  TVector3 out(x, y, z);
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

bool LArHadronCut(TVector3 position, double distance = 300) {
  bool out = true;
  const double LAr_Start_Exact[] = {-3478.48, -2166.71, 4179.24};
  const double LAr_End_Exact[] = {3478.48, 829.282, 9135.88};
  out = out && (position.X() >= LAr_Start_Exact[0] + distance) && (position.X() <= LAr_End_Exact[0] - distance);
  out = out && (position.Y() >= LAr_Start_Exact[1] + distance) && (position.Y() <= LAr_End_Exact[1] - distance);
  out = out && (position.Z() >= LAr_Start_Exact[2] + distance) && (position.Z() <= LAr_End_Exact[2] - distance);
  return out;
}

int PDGtoIndex(int pdgCode) {
    // Unknown is -999999999
    if (pdgCode < -999999990) return -1;
    switch (pdgCode) {
        case 11:   return 0;   // e-
        case -11:  return 0;   // e+
        case 22:   return 0;  // gamma
        case 13:   return 1;   // mu-
        case -13:  return 2;   // mu+
        case 211:  return 3;   // pi+
        case -211: return 4;   // pi-
        case 321:  return 5;   // K+
        case -321: return 5;   // K-
        case 310:  return 5;   // K0
        case 130:  return 5;  // K0_L
        case 311:  return 5;  // K0_S
        case 2112: return 6;  // Neutron
        case 2212: return 7;  // Proton
        case -2212: return 7;  // anti-Proton
        default:   return 8;  // other
    }
}

int PDGtoIndexReduced(int pdgCode) {
    const char *particle_types[] = {"electron", "muon", "pion", "kaon", "neutron", "proton", "other", "unknown"};
    if (pdgCode < -999999990) return 7;
    switch (pdgCode) {
        case 11:   return 0;   // e-
        case -11:  return 0;   // e+
        case 22:   return 0;  // gamma
        case 13:   return 1;   // mu-
        case -13:  return 1;   // mu+
        case 211:  return 2;   // pi+
        case -211: return 2;   // pi-
        case 321:  return 3;   // K+
        case -321: return 3;   // K-
        case 310:  return 3;   // K0
        case 130:  return 3;  // K0_L
        case 311:  return 3;  // K0_S
        case 2112: return 4;  // Neutron
        case 2212: return 5;  // Proton
        case -2212: return 5;  // anti-Proton
        default:   return 6;  // other
    }
}

Long64_t PrimaryLoop(Truth_Info& truth, Reco_Tree& reco, Line_Candidates& lc, int numEvents, TFile& outputFile) {
    // List all the hists here
    // Make sure to save them too

    bool has_kalman = reco.HasBranch("KalmanPos");
    std::cout<<"has_kalman status: "<<has_kalman<<std::endl;
    has_kalman = true;

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
    TH1D hist_ntracks("basic_n_tracks", "N Tracks;N Tracks;N Events",
      10, -0.5, 9.5);
    TH1D hist_nhits("basic_n_hits", "N Hits per Track;N Hits;N Tracks",
      25, 0, 100);
    TH1D hist_energy_range("basic_energy_range", "Energy Range;Energy Range (MeV);N Tracks",
      20, 0, 5000);
    TH1D hist_energy_deposit("basic_energy_deposit", "Energy Deposit;Energy Deposit (MeV);N Tracks",
      20, 0, 5000);
    TH1D hist_track_length("basic_track_length", "Track Length;Track Areal Density (g/cm^2);N Tracks",
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
    TH2D hist_track_hits_xz_relative_to_start("track_relative__hits_xz_relative_to_start", "Track Hit X vs Z, Relative to Start;Z (mm);X (mm); N Hits",
      nzfine, -2000, 2000, 100, -4000, 4000);
    TH2D hist_track_hits_xz_relative_to_third("track_relative__hits_xz_relative_to_third", "Track Hit X vs Z, Relative to Third Hit;Z (mm);X (mm); N Hits",
      nzfine, -2000, 2000, 100, -4000, 4000);
    TH2D hist_track_hits_xz_relative_to_end("track_relative__hits_xz_relative_to_end", "Track Hit X vs Z, Relative to End;Z (mm);X (mm); N Hits",
      nzfine, -2000, 2000, 100, -4000, 4000);
    TH2D hist_track_hits_xz_relative_to_third_last("track_relative__hits_xz_relative_to_third_last",
      "Track Hit X vs Z, Relative to Third Last Hit;Z (mm);X (mm); N Hits",
      nzfine, -2000, 2000, 100, -4000, 4000);

    TH1D hist_track_hits_z_relative_to_start("track_relative__hits_z_relative_to_start", "Track Hit Z, Relative to Start;Z (mm); N Hits",
      nzfine, -2000, 2000);
    TH1D hist_track_hits_z_relative_to_third("track_relative__hits_z_relative_to_third", "Track Hit Z, Relative to Third Hit;Z (mm); N Hits",
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
    TH1D hist_endpoint_resolution_z_within_10cm("endpoint_resolution_z_within_10cm", "Track Z Endpoint Resolution (TMS-ending muons only);Within 10cm, 0: No, 1: Yes;N Tracks",
      2, -0.5, 1.5);

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
    TH1D hist_matching_direction_yz_resolution("matching_direction_yz_resolution",
      "Matching Direction YZ Resolution (LAr-start, TMS-ending muons only);YZ Direction Reco - True;N Tracks",
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
    TH2D hist_matching_yz_position_true("matching_yz_position_true",
      "Matching True YZ Position (LAr-start, TMS-ending muons only);Track Z Position True (mm);Track Y Position True (mm);N Tracks",
      101, 11000, 16000, 101, -4000, 1000);
    TH2D hist_matching_xy_position_true("matching_xy_position_true",
      "Matching True XY Position (LAr-start, TMS-ending muons only);Track X Position True (mm);Track Y Position True (mm);N Tracks",
      101, -4000, 4000, 101, -4000, 1000);
    TH2D hist_matching_xz_position_reco("matching_xz_position_reco",
      "Matching Reco XZ Position (LAr-start, TMS-ending muons only);Track Z Position Reco (mm);Track X Position Reco (mm);N Tracks",
      101, 11000, 16000, 101, -4000, 4000);
    TH2D hist_matching_yz_position_reco("matching_yz_position_reco",
      "Matching Reco YZ Position (LAr-start, TMS-ending muons only);Track Z Position Reco (mm);Track Y Position Reco (mm);N Tracks",
      101, 11000, 16000, 101, -4000, 1000);
    TH2D hist_matching_xy_position_reco("matching_xy_position_reco",
      "Matching Reco XY Position (LAr-start, TMS-ending muons only);Track X Position Reco (mm);Track Y Position Reco (mm);N Tracks",
      101, -4000, 4000, 101, -4000, 1000);


    TH1D hist_matching_reco_extrap_x_vs_true_lar("matching_reco_extrap_x_vs_true_lar",
      "Matching X Position Reco TMS Track Extrap - True LAr (LAr-start, TMS-ending muons only);Reco Track X Extrap - True LAr X (mm);N Tracks",
      101, -500, 500);
    TH1D hist_matching_reco_extrap_x_vs_true_lar_no_direction("matching_reco_extrap_x_vs_true_lar_no_direction",
      "Matching X Position Reco TMS Track - True LAr (no direction extrapolation);Reco Track Start X - True LAr X (mm);N Tracks",
      101, -500, 500);
    TH1D hist_matching_reco_extrap_y_vs_true_lar("matching_reco_extrap_y_vs_true_lar",
      "Matching Y Position Reco TMS Track Extrap - True LAr (LAr-start, TMS-ending muons only);Reco Track Y Extrap - True LAr Y (mm);N Tracks",
      101, -3000, 3000);
    TH1D hist_matching_reco_extrap_y_vs_true_lar_no_direction("matching_reco_extrap_y_vs_true_lar_no_direction",
      "Matching Y Position Reco TMS Track - True LAr (no direction extrapolation);Reco Track Start Y - True LAr Y (mm);N Tracks",
      101, -3000, 3000);

    TH1D hist_matching_true_extrap_x_vs_true_lar("matching_true_extrap_x_vs_true_lar",
      "Matching X Position True TMS Track Extrap - True LAr (LAr-start, TMS-ending muons only);True Track X Extrap - True LAr X (mm);N Tracks",
      101, -500, 500);
    TH1D hist_matching_true_extrap_y_vs_true_lar("matching_true_extrap_y_vs_true_lar",
      "Matching Y Position True TMS Track Extrap - True LAr (LAr-start, TMS-ending muons only);True Track Y Extrap - True LAr Y (mm);N Tracks",
      101, -500, 500);
      
    TH1D hist_containment_TMS_reco_endpoint_cut("containment_TMS_reco_endpoint_cut",
      "Raw containment cut flags using reco muon endpoint;;N Tracks",
      7, 0, 7);
    TH1D hist_containment_TMS_reco_endpoint_cut_eff_numerator("containment_TMS_reco_endpoint_cut_eff_numerator",
      "Reco Tracks which are removed by TMS containment cuts;;N Tracks",
      7, 0, 7);
    TH1D hist_containment_TMS_reco_endpoint_cut_eff_denominator("containment_TMS_reco_endpoint_cut_eff_denominator",
      "Raw containment cut flags using reco muon endpoint;;N Tracks",
      7, 0, 7);
    TH1D hist_containment_TMS_true_endpoint_cut("containment_TMS_true_endpoint_cut",
      "Raw containment cut flags using true muon endpoint for reco tracks;;N Tracks",
      7, 0, 7);
    TH1D hist_containment_TMS_true_endpoint_cut_eff_numerator("containment_TMS_true_endpoint_cut_eff_numerator",
      "Raw containment cut flags using true muon endpoint for reco tracks;;N Tracks",
      7, 0, 7);
    TH1D hist_containment_TMS_true_endpoint_cut_eff_denominator("containment_TMS_true_endpoint_cut_eff_denominator",
      "Raw containment cut flags using true muon endpoint for reco tracks;;N Tracks",
      7, 0, 7);
    const char *flags[] = {"Passes all cuts", "Failed Any Cut", "Fail First Two Plane", "Fail Any Endpoint Cut", "Fail Last Plane Cut", "Fail Y Cut", "Fail UV Cut", " "};
    const int nflags = sizeof(flags) / sizeof(flags[0]);
    hist_containment_TMS_reco_endpoint_cut.SetNdivisions(nflags);
    hist_containment_TMS_reco_endpoint_cut_eff_numerator.SetNdivisions(nflags);
    hist_containment_TMS_reco_endpoint_cut_eff_denominator.SetNdivisions(nflags);
    hist_containment_TMS_true_endpoint_cut.SetNdivisions(nflags);
    hist_containment_TMS_true_endpoint_cut_eff_numerator.SetNdivisions(nflags);
    hist_containment_TMS_true_endpoint_cut_eff_denominator.SetNdivisions(nflags);
    for (int i = 0; i < nflags; i++) {
      hist_containment_TMS_reco_endpoint_cut.GetXaxis()->ChangeLabel(i+1, 350, -1, 11, -1, -1, flags[i]);
      hist_containment_TMS_reco_endpoint_cut_eff_numerator.GetXaxis()->ChangeLabel(i+1, 350, -1, 11, -1, -1, flags[i]);
      hist_containment_TMS_reco_endpoint_cut_eff_denominator.GetXaxis()->ChangeLabel(i+1, 350, -1, 11, -1, -1, flags[i]);
      hist_containment_TMS_true_endpoint_cut.GetXaxis()->ChangeLabel(i+1, 350, -1, 11, -1, -1, flags[i]);
      hist_containment_TMS_true_endpoint_cut_eff_numerator.GetXaxis()->ChangeLabel(i+1, 350, -1, 11, -1, -1, flags[i]);
      hist_containment_TMS_true_endpoint_cut_eff_denominator.GetXaxis()->ChangeLabel(i+1, 350, -1, 11, -1, -1, flags[i]);
    }
    
    TH1D hist_containment_reco_endpoint_all_tracks("containment_reco_endpoint_all_tracks",
      "All Track Endpoint Reco Z Position;Z (mm);N Tracks",
      101, 11000, 19000);
    TH1D hist_containment_reco_endpoint_passes_cuts("containment_reco_endpoint_passes_cuts",
      "Track Endpoint Reco Z Position, passes containment;Z (mm);N Tracks",
      101, 11000, 19000);
    TH2D hist_containment_check_tms_cut_all_tracks("containment_check_tms_cut_all_tracks",
      "Track XY Position;Track X Position (mm);Track Y Position (mm);N Tracks",
      101, -4000, 4000, 101, -4000, 1000);
    TH2D hist_containment_check_tms_cut_pass_uv_cut("containment_check_tms_cut_pass_uv_cut",
      "XY for Tracks That Pass UV Cut;Track X Position (mm);Track Y Position (mm);N Tracks",
      101, -4000, 4000, 101, -4000, 1000);
    TH2D hist_containment_check_tms_cut_pass_all_cuts("containment_check_tms_cut_pass_all_cuts",
      "XY for Tracks That Pass All Cuts;Track X Position (mm);Track Y Position (mm);N Tracks",
      101, -4000, 4000, 101, -4000, 1000);
    // TODO check reco eff vs muon entering energy given this cut
    // But also only calculate primary reco eff for nd physics sample

    // TODO Make direction extrap as well
    // Make plots about resolution assuming we don't take into account direction, just to see.
    // Like for y the direction might be so bad that it's making the extrapolation worse

    // TODO check for track leaving using reco Y. Also check for reco in first two planes cut, and occupancy > some amount cut

    // TODO make reco eff for muons reco'd in TMS vs not
    double muon_energy_bins[] = {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 5.1};
    int muon_energy_bins_n = sizeof(muon_energy_bins) / sizeof(double) - 1;
    TH1D hist_recoeff_muon_ke_tms_enter_numerator("recoeff_muon_ke_tms_enter_numerator",
      "Reco Eff Muon KE at Start of TMS Numerator;Muon KE at TMS Start (GeV);N Tracks",
      muon_energy_bins_n, muon_energy_bins);// 50, 0, 5); //
    TH1D hist_recoeff_muon_ke_tms_enter_denominator("recoeff_muon_ke_tms_enter_denominator",
      "Reco Eff Muon KE at Start of TMS Denominator;Muon KE at TMS Start (GeV);N Tracks",
      muon_energy_bins_n, muon_energy_bins); //50, 0, 5); //

    const char *pdg[] = {"e^{+/-}, #gamma", "#mu^{-}", "#mu^{+}", "#pi^{+}", "#pi^{-}", "K", "n", "p", "other", "unknown"};
    const int npdg = sizeof(pdg) / sizeof(pdg[0]);
    const int npdg_primary = npdg - 1;

    TH1D hist_recotrack_primary_particle_pdg("recotrack_primary_particle_pdg",
      "Reco Track Primary Particle;Particle;N Tracks", npdg_primary, -0.5, npdg_primary-0.5);
    TH1D hist_recotrack_secondary_particle_pdg("recotrack_secondary_particle_pdg",
      "Reco Track Secondary Particle;Particle;N Tracks", npdg, -0.5, npdg-0.5);
    TH1D hist_recotrack_has_secondary_particle("recotrack_has_secondary_particle",
      "Reco Track Has A Secondary Particle;0: No, 1: Yes;N Tracks", 2, -0.5, 1.5);
    // RecoTrackPrimaryParticleTrueVisibleEnergy / RecoTrackTrueVisibleEnergy
    TH1D hist_recotrack_primary_particle_energy_purity("recotrack_primary_particle_energy_purity",
      "Reco Track Primary Particle Visible Energy Purity;Purity (Primary Part. True Vis E in Track / All True Vis E in Track);N Tracks",
      50, 0.0, 1.0);
    // RecoTrackPrimaryParticleTrueVisibleEnergy / TrueVisibleEnergy
    TH1D hist_recotrack_primary_particle_energy_efficiency("recotrack_primary_particle_energy_efficiency",
      "Reco Track Primary Particle Visible Energy Efficiency;Efficiency (Primary Part. True Vis E in Track / All True Vis E for Part);N Tracks",
      50, 0.0, 1.0);

    hist_recotrack_primary_particle_pdg.SetNdivisions(npdg_primary);
    hist_recotrack_secondary_particle_pdg.SetNdivisions(npdg);
    for (int i = 0; i < npdg; i++) {
      hist_recotrack_primary_particle_pdg.GetXaxis()->ChangeLabel(i+1, -1, -1, -1, -1, -1, pdg[i]);
      hist_recotrack_secondary_particle_pdg.GetXaxis()->ChangeLabel(i+1, -1, -1, -1, -1, -1, pdg[i]);
    }

    //hist_recotrack_primary_particle_pdg.LabelsDeflate("X");
    //hist_recotrack_primary_particle_pdg.LabelsDeflate("Y");
    //hist_recotrack_primary_particle_pdg.LabelsOption("v");
    // TODO add hist to count n times a particle was in a reco track. Found a case of 3 times

    const char *particle_types[] = {"Electron", "Muon", "Pion", "Kaon", "Neutron", "Proton", "Other", "Unknown"};
    const int n_particle_types = sizeof(particle_types) / sizeof(particle_types[0]);
    std::vector<TH1D> hist_particledepth_track_endpoints_stack;
    std::vector<TH1D> hist_particledepth_true_endpoints_stack;
    std::vector<TH2D> hist_truthmatching_visible_energy_comparison;
    std::vector<TH2D> hist_truthmatching_x_comparison;
    std::vector<TH2D> hist_truthmatching_y_comparison;
    std::vector<TH2D> hist_truthmatching_z_comparison;
    std::vector<TH2D> hist_truthmatching_e_comparison;
    for (auto particle_type : particle_types) {
      {
        std::string name = "particledepth_track_endpoints_stack_" + std::string(particle_type);
        std::string title = "Reco Track Endpoint for " + std::string(particle_type) + ";Z (mm);N Tracks";
        TH1D hist(name.c_str(), title.c_str(), 51, 11000, 19000);
        hist_particledepth_track_endpoints_stack.push_back(hist);

        std::string name2 = "particledepth_true_endpoints_stack_" + std::string(particle_type);
        std::string title2 = "True Particle Endpoint for Primary Particle of Reco Track for " + std::string(particle_type) + ";Z (mm);N Tracks";
        TH1D hist2(name2.c_str(), title2.c_str(), 101, 0, 25000);
        hist_particledepth_true_endpoints_stack.push_back(hist2);

        std::string name3 = "truthmatching_visible_energy_comparison_" + std::string(particle_type);
        std::string title3 = "True vs Reco Visible Energy " + std::string(particle_type) + ";True Vis E (MeV);Reco Vis E (MeV);N Tracks";
        TH2D hist3(name3.c_str(), title3.c_str(), 100, 0, 1000, 100, 0, 300);
        hist_truthmatching_visible_energy_comparison.push_back(hist3);
      }

      {
        std::string name = "truthmatching_x_comparison_" + std::string(particle_type);
        std::string title = "True vs Reco Track True X " + std::string(particle_type) + ";True X (mm);Reco Track True X (mmMeV);N Tracks";
        TH2D hist(name.c_str(), title.c_str(), 100, -4000, 4000, 100, -4000, 4000);
        hist_truthmatching_x_comparison.push_back(hist);
      }

      {
        std::string name = "truthmatching_y_comparison_" + std::string(particle_type);
        std::string title = "True vs Reco Track True Y " + std::string(particle_type) + ";True Y (mm);Reco Track True Y (mm);N Tracks";
        TH2D hist(name.c_str(), title.c_str(), 100, -4000, 1000, 100, -4000, 1000);
        hist_truthmatching_y_comparison.push_back(hist);
      }

      {
        std::string name = "truthmatching_z_comparison_" + std::string(particle_type);
        std::string title = "True vs Reco Track True Z " + std::string(particle_type) + ";True Z (mm);Reco Track True X (mm);N Tracks";
        TH2D hist(name.c_str(), title.c_str(), 100, 0, 25000, 100, 0, 25000);
        hist_truthmatching_z_comparison.push_back(hist);
      }

      {
        std::string name = "truthmatching_e_comparison_" + std::string(particle_type);
        std::string title = "True vs Reco Track True E " + std::string(particle_type) + ";True E (MeV);Reco Track True E (MeV);N Tracks";
        TH2D hist(name.c_str(), title.c_str(), 100, 0, 25000, 100, 0, 25000);
        hist_truthmatching_e_comparison.push_back(hist);
      }
    }

    int last_spill_seen = -1;

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
          reco.nKalmanNodes[itrack] = reco.nHits[itrack];
        }
      }
      // Fix n kalman nodes which is always zero.
      for (int itrack = 0; itrack < reco.nTracks; itrack++) {
        int n_nonzero = 0;
        for (int ihit = 0; ihit < reco.nHits[itrack]; ihit++) {
          // Since kalman pos need to be above z = 11000, we can use this to find zero points safely
          if (reco.KalmanPos[itrack][ihit][2] > 1000) n_nonzero++;
        }
        reco.nKalmanNodes[itrack] = n_nonzero;
        //std::cout<<"track num: "<<itrack<<"\tn hits: "<<reco.nHits[itrack]<<"\tn kalman nodes: "<<reco.nKalmanNodes[itrack]<<std::endl;
      }

      // TODO calculate plane number and then check per plane occupancy
      // Also related is total visible energy so compare that to true value
      
       
      //std::cout<<entry_number<<": "<<reco.SpillNo<<", "<<reco.SliceNo<<", "<<reco.EventNo<<std::endl;
      if (lc.nHits > 100) {
        //DrawSlice(TString::Format("test_%d", entry_number).Data(), "test_event_displays", "test job", reco, lc, truth);
      }
     // if (entry_number > 700) exit(1); // TODO delete

      // Calculate reco eff, denominators
      // Only fill once per spill since all particles per spill are saved, not just this time slice's
      if (last_spill_seen != reco.SpillNo) {
        last_spill_seen = reco.SpillNo;
        for (int ipart = 0; ipart < truth.nTrueParticles; ipart++) {
          bool lar_starting = truth.LArFiducialStart[ipart];
          TVector3 start_pos(truth.BirthPosition[ipart][0], truth.BirthPosition[ipart][1], truth.BirthPosition[ipart][2]);
          bool lar_hadron_cut = LArHadronCut(start_pos);
          bool tms_ending = truth.TMSFiducialEnd[ipart];
          bool ismuon = abs(truth.PDG[ipart]) == 13;
          //if (ismuon && tms_ending && lar_starting) std::cout<<"True start pos: ("<<start_pos.X()<<","<<start_pos.Y()<<","<<start_pos.Y()<<"), hadron cut: "<<lar_hadron_cut<<", spill number: "<<reco.SpillNo<<std::endl;
          if (lar_starting && lar_hadron_cut && tms_ending && ismuon) {
            double muon_starting_ke = truth.MomentumTMSStart[ipart][3] * 1e-3;
            if (muon_starting_ke > 5) muon_starting_ke = 5.05; // Overflow bin
            hist_recoeff_muon_ke_tms_enter_denominator.Fill(muon_starting_ke);
          }
        }
      }

      // Calculate the reco eff, numerators
      for (int itrack = 0; itrack < reco.nTracks; itrack++) {
        bool lar_starting = truth.RecoTrackPrimaryParticleLArFiducialStart[itrack];
        TVector3 start_pos(truth.RecoTrackPrimaryParticleTruePositionStart[itrack][0],
                 truth.RecoTrackPrimaryParticleTruePositionStart[itrack][1],
                 truth.RecoTrackPrimaryParticleTruePositionStart[itrack][2]);
        bool lar_hadron_cut = LArHadronCut(start_pos);
        bool tms_ending = truth.RecoTrackPrimaryParticleTMSFiducialEnd[itrack];
        bool ismuon = abs(truth.RecoTrackPrimaryParticlePDG[itrack]) == 13;
        //if (ismuon && tms_ending && lar_starting) std::cout<<"Reco track true start pos: ("<<start_pos.X()<<","<<start_pos.Y()<<","<<start_pos.Y()<<"), hadron cut: "<<lar_hadron_cut<<", spill number: "<<reco.SpillNo<<", slice number: "<<reco.SliceNo<<", track num: "<<itrack<<std::endl;
        if (lar_starting && lar_hadron_cut && tms_ending && ismuon) {
          double muon_starting_ke = truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[itrack][3] * 1e-3;
          if (muon_starting_ke > 5) muon_starting_ke = 5.05; // Overflow bin
          hist_recoeff_muon_ke_tms_enter_numerator.Fill(muon_starting_ke);
        }
      }
      
      // Check the TMS cuts using reco information
      for (int itrack = 0; itrack < reco.nTracks; itrack++) {
        // Fill using reco track endpoint
        TVector3 track_startpoint(reco.StartPos[itrack][0], reco.StartPos[itrack][1], reco.StartPos[itrack][2]);
        TVector3 track_endpoint(reco.EndPos[itrack][0], reco.EndPos[itrack][1], reco.EndPos[itrack][2]);
        auto flags = CheckTMSCuts(track_startpoint, track_endpoint);
        FillHistWithTMSCutFlagInfo(hist_containment_TMS_reco_endpoint_cut, flags);
        FillHistWithTMSCutFlagInfo(hist_containment_TMS_reco_endpoint_cut_eff_numerator, flags);
        for (int i = 0; i < hist_containment_TMS_reco_endpoint_cut_eff_denominator.GetXaxis()->GetNbins(); i++)
          hist_containment_TMS_reco_endpoint_cut_eff_denominator.Fill(i);
          
        hist_containment_reco_endpoint_all_tracks.Fill(track_endpoint.Z());
        if (flags == 0) hist_containment_reco_endpoint_passes_cuts.Fill(track_endpoint.Z());
        
        // Fill same but using true endpoint of the reco track
        int particle_index = truth.RecoTrackPrimaryParticleIndex[itrack];
        if (particle_index >= 0 && particle_index < truth.nTrueParticles) {
          TVector3 true_track_startpoint(truth.PositionZIsTMSStart[particle_index][0], 
                                         truth.PositionZIsTMSStart[particle_index][1], 
                                         truth.PositionZIsTMSStart[particle_index][2]);
          TVector3 true_track_endpoint(truth.DeathPosition[itrack][0], 
                                        truth.DeathPosition[itrack][1], 
                                        truth.DeathPosition[itrack][2]);
          auto true_flags = CheckTMSCuts(true_track_startpoint, true_track_endpoint);
          FillHistWithTMSCutFlagInfo(hist_containment_TMS_true_endpoint_cut, true_flags);
          FillHistWithTMSCutFlagInfo(hist_containment_TMS_true_endpoint_cut_eff_numerator, true_flags);
          for (int i = 0; i < hist_containment_TMS_true_endpoint_cut_eff_denominator.GetXaxis()->GetNbins(); i++)
            hist_containment_TMS_true_endpoint_cut_eff_denominator.Fill(i);
          
          
          
          
          // Check which tracks pass the UV cut
          hist_containment_check_tms_cut_all_tracks.Fill(true_track_endpoint.X(), true_track_endpoint.Y());
          if ((true_flags & TMSCutUVCut) == 0) 
            hist_containment_check_tms_cut_pass_uv_cut.Fill(true_track_endpoint.X(), true_track_endpoint.Y());
          if (true_flags == 0) 
            hist_containment_check_tms_cut_pass_all_cuts.Fill(true_track_endpoint.X(), true_track_endpoint.Y());
        }
      }

      // Loop over tracks and record information about particles that were reconstructed
      for (int itrack = 0; itrack < reco.nTracks; itrack++) {
        int primary_pdg = PDGtoIndex(truth.RecoTrackPrimaryParticlePDG[itrack]);
        if (primary_pdg == 8) std::cout<<"pdg missing "<<truth.RecoTrackPrimaryParticlePDG[itrack]<<std::endl;
        int secondary_pdg = PDGtoIndex(truth.RecoTrackSecondaryParticlePDG[itrack]);
        // Not all particles are tracked, so check if there's a nonzero secondary particle visible energy and tag as unknown
        if (secondary_pdg == -1 && truth.RecoTrackSecondaryParticleTrueVisibleEnergy[itrack] > 0) secondary_pdg = 9;

        // Now fill hists:
        if (truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][2] > 11000) 
          hist_recotrack_primary_particle_pdg.Fill(primary_pdg, 1);
        if (truth.RecoTrackSecondaryParticleTruePositionEnd[itrack][2] > 11000) {
          hist_recotrack_secondary_particle_pdg.Fill(secondary_pdg, 1);
          hist_recotrack_has_secondary_particle.Fill((secondary_pdg < 0) ? 0 : 1);
        }
        // TODO are eff and purity calculations right? They seem low
        double pur = truth.RecoTrackPrimaryParticleTrueVisibleEnergy[itrack] / truth.RecoTrackTrueVisibleEnergy[itrack];
        if (pur < 1 && truth.RecoTrackSecondaryParticleTrueVisibleEnergy[itrack] == 0)
          std::cout<<"Warning: Found purity < 1 where there was no secondary energy"<<std::endl;
        hist_recotrack_primary_particle_energy_purity.Fill(pur);
        int particle_index = truth.RecoTrackPrimaryParticleIndex[itrack];
        if (particle_index >= 0 && particle_index < truth.nTrueParticles) {
          double eff = truth.RecoTrackPrimaryParticleTrueVisibleEnergy[itrack] / truth.TrueVisibleEnergy[particle_index];
          //std::cout<<"eff: "<<eff<<",\tprimary: "<<truth.RecoTrackPrimaryParticleTrueVisibleEnergy[itrack]<<" MeV,\tTotal: "<<truth.TrueVisibleEnergy[particle_index]<<" MeV"<<std::endl;
          hist_recotrack_primary_particle_energy_efficiency.Fill(eff);
          //std::cout<<"pur: "<<pur<<", eff: "<<eff<<", secondary_pdg: "<<secondary_pdg<<", "<<truth.RecoTrackSecondaryParticlePDG[itrack]<<std::endl;
        }
        //else std::cout<<"pur: "<<pur<<", secondary_pdg: "<<secondary_pdg<<", "<<truth.RecoTrackSecondaryParticlePDG[itrack]<<std::endl;
        
        bool is_good_truth_matching = true;
        if (truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][2] < 11000) is_good_truth_matching = false;

        int index = PDGtoIndexReduced(truth.RecoTrackPrimaryParticlePDG[itrack]);
        if (is_good_truth_matching) hist_particledepth_track_endpoints_stack[index].Fill(reco.EndPos[itrack][2]);
        //if (truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][2] > truth.RecoTrackPrimaryParticleTruePositionStart[itrack][2])
        if (is_good_truth_matching) hist_particledepth_true_endpoints_stack[index].Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][2]);
        

        if (particle_index >= 0 && particle_index < truth.nTrueParticles && is_good_truth_matching) { //  && !truth.IsPrimary[particle_index]
          double true_visible_energy = truth.TrueVisibleEnergy[particle_index];
          // RecoTrackPrimaryParticleTrueVisibleEnergy is true vis energy on track only
          double reco_visible_energy = truth.RecoTrackPrimaryParticleTrueVisibleEnergy[itrack];
          hist_truthmatching_visible_energy_comparison[index].Fill(true_visible_energy, reco_visible_energy);

          double true_x = truth.BirthPosition[particle_index][0];
          double reco_x = truth.RecoTrackPrimaryParticleTruePositionTrackStart[itrack][0];
          hist_truthmatching_x_comparison[index].Fill(true_x, reco_x);
          hist_truthmatching_y_comparison[index].Fill(truth.BirthPosition[particle_index][1],
            truth.RecoTrackPrimaryParticleTruePositionTrackStart[itrack][1]);
          hist_truthmatching_z_comparison[index].Fill(truth.BirthPosition[particle_index][2],
            truth.RecoTrackPrimaryParticleTruePositionTrackStart[itrack][2]);
          hist_truthmatching_e_comparison[index].Fill(truth.BirthMomentum[particle_index][3],
            truth.RecoTrackPrimaryParticleTrueMomentumTrackStart[itrack][3]);
          //std::cout<<"Birth energy: "<<truth.BirthMomentum[particle_index][3]<<"MeV,\tTrueMomentum: "<<truth.RecoTrackPrimaryParticleTrueMomentum[itrack][3]<<"MeV,\tTrackStart: "<<truth.RecoTrackPrimaryParticleTrueMomentumTrackStart[itrack][3]<<" MeV"<<std::endl;

          // Print info if more than 10% different, bad matching
          double truth_e = truth.BirthMomentum[particle_index][3];
          double reco_truth_e = truth.RecoTrackPrimaryParticleTrueMomentumTrackStart[itrack][3];
          bool should_print_info = false;
          if (truth_e > 3000) {
            if (reco_truth_e < 2000) should_print_info = true;
          }
          if (should_print_info && false) {
            std::cout<<"entry: "<<entry_number;
            std::cout<<"\tparticle_index: "<<particle_index;
            std::cout<<"\titrack: "<<itrack;
            std::cout<<"\tPDG: "<<truth.PDG[particle_index];
            std::cout<<"\ttrack pdg: "<<truth.RecoTrackPrimaryParticlePDG[itrack];
            std::cout<<"\tBirth energy: "<<truth.BirthMomentum[particle_index][3];
            std::cout<<"MeV\tRecoTrackPrimaryParticleTrueMomentum: "<<truth.RecoTrackPrimaryParticleTrueMomentum[itrack][3];
            std::cout<<" MeV,\tTrackStart: "<<truth.RecoTrackPrimaryParticleTrueMomentumTrackStart[itrack][3]<<" MeV\t";
            std::cout<<"\tbirth position z: "<<truth.BirthPosition[particle_index][2];
            std::cout<<"\tdeath position z: "<<truth.DeathPosition[particle_index][2];
            std::cout<<"\ttrack birth position z: "<<truth.RecoTrackPrimaryParticleTruePositionStart[itrack][2];
            std::cout<<"\ttrack death position z: "<<truth.RecoTrackPrimaryParticleTruePositionEnd[itrack][2];
            std::cout<<"\tparticle primary: "<<truth.IsPrimary[particle_index];
            std::cout<<"\ttrack is primary: "<<truth.RecoTrackPrimaryParticleIsPrimary[itrack];
            std::cout<<"\tsecondary index: "<<truth.RecoTrackSecondaryParticleIndex[itrack];
            std::cout<<std::endl;
          }

        }
      }

      hist_ntracks.Fill(reco.nTracks);
      for (int itrack = 0; itrack < reco.nTracks; itrack++) {
        hist_nhits.Fill(reco.nKalmanNodes[itrack]);
        hist_energy_range.Fill(reco.EnergyRange[itrack]);
        hist_energy_deposit.Fill(reco.EnergyDeposit[itrack]);
        hist_track_length.Fill(reco.Length[itrack]);

        // Matching
        //std::cout<<"truth.RecoTrackPrimaryParticleLArFiducialStart[itrack]: "<<truth.RecoTrackPrimaryParticleLArFiducialStart[itrack]<<",\ttruth.RecoTrackPrimaryParticleTMSFiducialEnd[itrack]: "<<truth.RecoTrackPrimaryParticleTMSFiducialEnd[itrack]<<std::endl;
        if (truth.RecoTrackPrimaryParticleLArFiducialStart[itrack] && truth.RecoTrackPrimaryParticleTMSFiducialEnd[itrack]) {
          //std::cout<<"Direction ("<<reco.Direction[itrack][0]<<","<<reco.Direction[itrack][1]<<","<<reco.Direction[itrack][2]<<")"<<std::endl;
          double direction_reco = reco.Direction[itrack][0] / reco.Direction[itrack][2];
          double direction_reco_y = reco.Direction[itrack][1] / reco.Direction[itrack][2];
          double direction_true = truth.RecoTrackPrimaryParticleTrueMomentumTrackStart[itrack][0] / truth.RecoTrackPrimaryParticleTrueMomentumTrackStart[itrack][2];
          double direction_true_y = truth.RecoTrackPrimaryParticleTrueMomentumTrackStart[itrack][1] / truth.RecoTrackPrimaryParticleTrueMomentumTrackStart[itrack][2];
          double angle_reco = atan2(direction_reco, 1) * 360 / TAU;
          double angle_true = atan2(direction_true, 1) * 360 / TAU;
          hist_matching_angle_resolution.Fill(angle_reco - angle_true);
          hist_matching_angle_true.Fill(angle_true);
          hist_matching_angle_reco.Fill(angle_reco);
          hist_matching_direction_resolution.Fill(direction_reco - direction_true);
          hist_matching_direction_yz_resolution.Fill(direction_reco_y - direction_true_y);
          double track_x_reco = reco.StartPos[itrack][0];
          double track_y_reco = reco.StartPos[itrack][1];
          double track_z_reco = reco.StartPos[itrack][2];

          hist_matching_xz_position_reco.Fill(track_z_reco, track_x_reco);
          hist_matching_yz_position_reco.Fill(track_z_reco, track_y_reco);
          hist_matching_xy_position_reco.Fill(track_x_reco, track_y_reco);
          /*
          double track_x_true = truth.RecoTrackPrimaryParticleTruePositionEnteringTMS[itrack][0];
          double track_y_true = truth.RecoTrackPrimaryParticleTruePositionEnteringTMS[itrack][1];
          double track_z_true = truth.RecoTrackPrimaryParticleTruePositionEnteringTMS[itrack][2]; */
          int particle_index = truth.RecoTrackPrimaryParticleIndex[itrack];
          if (particle_index >= 0 && particle_index < truth.nTrueParticles) {
            double track_x_true = truth.PositionZIsTMSStart[particle_index][0];
            double track_y_true = truth.PositionZIsTMSStart[particle_index][1];
            double track_z_true = truth.PositionZIsTMSStart[particle_index][2];
            hist_matching_x_position_resolution.Fill(track_x_reco - track_x_true);
            hist_matching_y_position_resolution.Fill(track_y_reco - track_y_true);
            hist_matching_z_position_resolution.Fill(track_z_reco - track_z_true);
            hist_matching_z_position_reco.Fill(track_z_reco);
            hist_matching_z_position_true.Fill(track_z_true);
            hist_matching_xz_position_true.Fill(track_z_true, track_x_true);
            hist_matching_yz_position_true.Fill(track_z_true, track_y_true);
            hist_matching_xy_position_true.Fill(track_x_true, track_y_true);


            const double lar_end_face = 9135.88;
            TVector3 track_start_position(track_x_reco, track_y_reco, track_z_reco);
            TVector3 reco_tms_position = calculatePositionAtZ(lar_end_face, track_start_position, direction_reco, direction_reco_y);

            hist_matching_reco_extrap_x_vs_true_lar.Fill(reco_tms_position.X() - truth.PositionZIsLArEnd[particle_index][0]);
            hist_matching_reco_extrap_x_vs_true_lar_no_direction.Fill(track_start_position.X() - truth.PositionZIsLArEnd[particle_index][0]);
            hist_matching_reco_extrap_y_vs_true_lar.Fill(reco_tms_position.Y() - truth.PositionZIsLArEnd[particle_index][1]);
            hist_matching_reco_extrap_y_vs_true_lar_no_direction.Fill(track_start_position.Y() - truth.PositionZIsLArEnd[particle_index][1]);

            TVector3 true_track_start_position(track_x_true, track_y_true, track_z_true);
            TVector3 true_tms_position = calculatePositionAtZ(lar_end_face, true_track_start_position, direction_true, direction_true_y);
            hist_matching_true_extrap_x_vs_true_lar.Fill(true_tms_position.X() - truth.PositionZIsLArEnd[particle_index][0]);
            hist_matching_true_extrap_y_vs_true_lar.Fill(true_tms_position.Y() - truth.PositionZIsLArEnd[particle_index][1]);

          }
          else {
            std::cout<<"Found particle_index "<<particle_index<<" which is outside 0-"<<truth.nTrueParticles<<std::endl;
          }

        }

        //std::cout<<"start xyz "<<itrack<<": "<<reco.StartPos[itrack][0]<<","<<reco.StartPos[itrack][1]<<","<<reco.StartPos[itrack][2]<<std::endl;
        //std::cout<<"end xyz "<<itrack<<":   "<<reco.EndPos[itrack][0]<<","<<reco.EndPos[itrack][1]<<","<<reco.EndPos[itrack][2]<<std::endl;
        hist_track_start_pos_xz.Fill(reco.StartPos[itrack][2], reco.StartPos[itrack][0]);
        hist_track_end_pos_xz.Fill(reco.EndPos[itrack][2], reco.EndPos[itrack][0]);
        for (int ihit = 0; ihit < reco.nKalmanNodes[itrack]; ihit++) {
          hist_track_hits_xz.Fill(reco.KalmanPos[itrack][ihit][2], reco.KalmanPos[itrack][ihit][0]);
          hist_track_hits_yz.Fill(reco.KalmanPos[itrack][ihit][2], reco.KalmanPos[itrack][ihit][1]);
          hist_track_hits_xy.Fill(reco.KalmanPos[itrack][ihit][0], reco.KalmanPos[itrack][ihit][1]);

          int last_index = reco.nKalmanNodes[itrack] - 1;
          double z_relative_start = reco.KalmanPos[itrack][ihit][2] - reco.KalmanPos[itrack][0][2];
          double x_relative_start = reco.KalmanPos[itrack][ihit][0] - reco.KalmanPos[itrack][0][0];
          double z_relative_end = reco.KalmanPos[itrack][ihit][2] - reco.KalmanPos[itrack][last_index][2];
          double x_relative_end = reco.KalmanPos[itrack][ihit][0] - reco.KalmanPos[itrack][last_index][0];
          hist_track_hits_xz_relative_to_start.Fill(z_relative_start, x_relative_start);
          hist_track_hits_xz_relative_to_end.Fill(z_relative_end, x_relative_end);
          hist_track_hits_z_relative_to_start.Fill(z_relative_start);
          if (reco.nKalmanNodes[itrack] >= 3) {
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
        hist_comp_track_dz_vs_n_hits.Fill(dz, reco.nKalmanNodes[itrack]);
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
          hist_endpoint_resolution_z_within_10cm.Fill((abs(z) < 100) ? 1 : 0);
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

    // Create output filename, GIT_BRANCH_NAME + "_" +
    std::string outputFilename = directoryPath + getOutputFilename(inputFilename);

    // Create TFile with the output filename
    TFile outputFile(outputFilename.c_str(), "RECREATE");

    PrimaryLoop(ti, ri, li, numEvents, outputFile);

    // Close the output file
    outputFile.Close();

    std::cout << "Output file created: " << outputFilename << std::endl;

    return 0;
}
