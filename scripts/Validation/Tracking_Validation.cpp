#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <filesystem>
#include <tuple>
#include <chrono> // for high res clock
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

std::string save_location = "";

const double CM = 0.1; // cm per mm
const double DEG = 360 / TAU;

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
  if (save_location == "") {
    throw std::runtime_error("Do not have a save directory for the event displays in DrawSlice");
  }
  std::string directoryPath = save_location + "EventDisplays/" + reason + "/";
  createDirectory(directoryPath);
  static TCanvas* canvas;
  if (canvas == NULL) {
    gStyle->SetOptStat(0);
    gROOT->SetBatch(kTRUE); 
    canvas = new TCanvas();
    canvas->SetLeftMargin(0.12);
    canvas->Divide(1, 2);
  }
  
  float buffer = 0;
  TH2D hist("", "X vs Z View;Z (mm);X (mm);N Hits", 100, 11000 - buffer, 19000 + buffer, 100, -4000 - buffer, 4000 + buffer);
  TH2D histy("", "Y vs Z View;Z (mm);Y (mm);N Hits", 100, 11000 - buffer, 19000 + buffer, 100, -5000 - buffer, 1000 + buffer);
  
  
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
      marker.SetMarkerStyle(21);
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
        marker.SetMarkerStyle(33);
        marker.SetMarkerSize(0.75);
        markers.push_back(marker);
      }
      {
        float mx = truth.RecoTrackTrueHitPosition[it][ih][2];
        float my = truth.RecoTrackTrueHitPosition[it][ih][1];
        TMarker marker(mx, my, 21);
        marker.SetMarkerColor(kGreen);
        marker.SetMarkerStyle(33);
        marker.SetMarkerSize(0.75);
        markersy.push_back(marker);
      }
      {
        float mx = reco.TrackHitPos[it][ih][2];
        float my = reco.TrackHitPos[it][ih][0];
        TMarker marker(mx, my, 20);
        marker.SetMarkerColor(kRed);
        marker.SetMarkerSize(0.5);
        markers.push_back(marker);
      }
      {
        float mx = reco.TrackHitPos[it][ih][2];
        float my = reco.TrackHitPos[it][ih][1];
        TMarker marker(mx, my, 20);
        marker.SetMarkerColor(kRed);
        marker.SetMarkerSize(0.5);
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

bool LArFiducialCut(TVector3 position, double distance = 500, double downstream = 1500) {
  // Returns true if inside fiducial, false otherwise
  bool out = true;
  const double LAr_Start_Exact[] = {-3478.48, -2166.71, 4179.24};
  const double LAr_End_Exact[] = {3478.48, 829.282, 9135.88};
  out = out && (position.X() >= LAr_Start_Exact[0] + distance) && (position.X() <= LAr_End_Exact[0] - distance);
  out = out && (position.Y() >= LAr_Start_Exact[1] + distance) && (position.Y() <= LAr_End_Exact[1] - distance);
  out = out && (position.Z() >= LAr_Start_Exact[2] + distance) && (position.Z() <= LAr_End_Exact[2] - downstream);
  return out;
}

int PDGtoIndex(int pdgCode) {
    // const char *pdg[] = {"e^{+/-}, #gamma", "#mu^{-}", "#mu^{+}", "#pi^{+}", "#pi^{-}", "K", "n", "p", "other", "unknown"};
    // Unknown is -999999999
    if (pdgCode < -999999990) return 9;
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
    // const char *particle_types[] = {"electron", "muon", "pion", "kaon", "neutron", "proton", "other", "unknown"};
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

std::unordered_map<std::string, TH1*> mapForGetHist;

std::tuple<std::string, int, double, double> GetBinning(std::string axis_name) {
  if (axis_name == "ntracks") return std::make_tuple("N Tracks", 10, -0.5, 9.5);
  if (axis_name == "n0-120") return std::make_tuple("N", 24, 0, 120);
  if (axis_name == "n0-500") return std::make_tuple("N", 25, 0, 500);
  if (axis_name == "EventNo") return std::make_tuple("Event Number", 100, 0, 5000);
  if (axis_name == "SliceNo") return std::make_tuple("Slice Number", 61, -0.5, 60.5);
  if (axis_name == "SpillNo") return std::make_tuple("Spill Number", 100, 0, 300);
  if (axis_name == "X") return std::make_tuple("X (cm)", 100, -400, 400);
  if (axis_name == "Y") return std::make_tuple("Y (cm)", 100, -500, 100);
  if (axis_name == "Z") return std::make_tuple("Z (cm)", 100, 1100, 1900);
  if (axis_name == "direction_xz") return std::make_tuple("XZ Direction", 31, -2, 2);
  if (axis_name == "dx") return std::make_tuple("dX (cm)", 100, -100, 100);
  if (axis_name == "dy") return std::make_tuple("dY (cm)", 100, -100, 100);
  if (axis_name == "dz") return std::make_tuple("dZ (cm)", 100, -100, 100);
  if (axis_name == "pdg") return std::make_tuple("Particle", 10, -0.5, 9.5);
  if (axis_name == "angle_tms_enter") return std::make_tuple("Angle (deg)", 30, -60, 60);
  std::cerr<<"Fatal: Add axis to GetBinning. Did not understand axis name "<<axis_name<<std::endl;
  throw std::runtime_error("Unable to understand axis name");
}

double muon_ke_bins[] = {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 5.1};
int n_muon_ke_bins = sizeof(muon_ke_bins) / sizeof(double) - 1;

std::tuple<bool, std::string, int, double*> GetComplexBinning(std::string axis_name) {
  if (axis_name == "ke_tms_enter") return std::make_tuple(true, "Muon KE Entering TMS (GeV)", n_muon_ke_bins, muon_ke_bins); 
  return std::make_tuple(false, "", 0, (double*)NULL);
}

void AdjustAxis(TH1* hist, std::string xaxis, std::string yaxis = "", std::string zaxis = "") {
  if (xaxis == "pdg") {
    const char *pdg[] = {"e^{+/-}, #gamma", "#mu^{-}", "#mu^{+}", "#pi^{+}", "#pi^{-}", "K", "n", "p", "other", "unknown"};
    const int npdg = sizeof(pdg) / sizeof(pdg[0]);
    hist->SetNdivisions(npdg);
    for (int ib = 0; ib < npdg; ib++) {
      hist->GetXaxis()->ChangeLabel(ib+1, -1, -1, -1, -1, -1, pdg[ib]);
    }
  }
}

TH1* MakeHist(std::string directory_and_name, std::string title, std::string xaxis, std::string yaxis = "", std::string zaxis = "") {
  if (zaxis != "") {
    // 3d hist case
    throw std::runtime_error("3d hists are not implmented yet");
  }
  else if (yaxis != "") {
    // 2d hist case
    throw std::runtime_error("2d hists are not implmented yet");
  }
  else {
    // 1d hist case
    TH1D* out;
    std::string xaxis_title;
    int xaxis_nbins;
    double* xaxis_bins;
    bool has_complex_binning;
    std::tie(has_complex_binning, xaxis_title,  xaxis_nbins, xaxis_bins) = GetComplexBinning(xaxis);
    if (has_complex_binning) {
    
      auto complete_title = title + ";" + xaxis_title;
      out = new TH1D(directory_and_name.c_str(), complete_title.c_str(), xaxis_nbins, xaxis_bins);
    } else {
      double xaxis_start;
      double xaxis_end;
      std::tie(xaxis_title,  xaxis_nbins, xaxis_start, xaxis_end) = GetBinning(xaxis);
      auto complete_title = title + ";" + xaxis_title;
      out = new TH1D(directory_and_name.c_str(), complete_title.c_str(), xaxis_nbins, xaxis_start, xaxis_end);
    }
    // Add special naming here
    AdjustAxis(out, xaxis);
    return out;
  }
}

double total_lookup_time = 0;
double total_no_make_time = 0;
double total_make_time = 0;
int n_lookups = 0;
int n_makes = 0;
int n_no_makes = 0;

TH1* GetHist(std::string directory_and_name, std::string title, std::string xaxis, std::string yaxis = "", std::string zaxis = "") {
    auto time_start = std::chrono::high_resolution_clock::now();
    bool did_make = false;
    if (mapForGetHist.find(directory_and_name) == mapForGetHist.end()) {
        // Object doesn't exist, create it and store it
        mapForGetHist[directory_and_name] = MakeHist(directory_and_name, title, xaxis, yaxis, zaxis);
      auto time_stop_make_only = std::chrono::high_resolution_clock::now();
      auto duration_make_only = std::chrono::duration_cast<std::chrono::microseconds>(time_stop_make_only - time_start).count();
      total_make_time += duration_make_only;
      n_makes += 1;
      did_make = true;
    }
    auto out = mapForGetHist[directory_and_name];
    auto time_stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(time_stop - time_start).count();
    total_lookup_time += duration;
    n_lookups += 1;
    if (!did_make) {
      total_no_make_time += duration;
      n_no_makes += 1;
    }
    if (n_lookups % 1000000 == 0) {
      double avg_lookup_time = total_lookup_time / n_lookups;
      double avg_make_time = total_make_time / n_makes;
      double avg_no_make_time = total_no_make_time / n_no_makes;
      
      std::cout<<"Avg lookup time: "<<avg_lookup_time<<" ("<<total_lookup_time<<"/"<<n_lookups<<")"<<std::endl;
      std::cout<<"Avg make time: "<<avg_make_time<<" ("<<total_make_time<<"/"<<n_makes<<")"<<std::endl;
      std::cout<<"Avg no make time: "<<avg_no_make_time<<" ("<<total_no_make_time<<"/"<<n_no_makes<<")"<<std::endl;
    }
    return out;
}

Long64_t PrimaryLoop(Truth_Info& truth, Reco_Tree& reco, Line_Candidates& lc, int numEvents, TFile& outputFile) {
    // List all the hists here
    // Make sure to save them too

    bool has_kalman = reco.HasBranch("KalmanPos");
    std::cout<<"has_kalman status: "<<has_kalman<<std::endl;
    has_kalman = false; // TODO this needs to be set based on the input file

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

    int last_spill_seen = -1;
    
    auto time_start = std::chrono::high_resolution_clock::now();

    Long64_t entry_number = 0;
    // Now loop over the ttree
    for ( ; entry_number < truth.GetEntriesFast() && entry_number < reco.GetEntriesFast()\
      && (numEvents < 0 || entry_number < numEvents); entry_number++) {
      if (entry_number % 10000 == 0) std::cout<<"On entry: "<<entry_number<<std::endl;
      
      if (entry_number > 10000) break;
      
      // Get the current entry
      // Currently reco and truth match
      truth.GetEntry(entry_number);
      reco.GetEntry(entry_number);
      lc.GetEntry(entry_number);
      
      // Check if we're on a new spill
      // Some things should only be filled per spill
      bool on_new_spill = false;
      if (last_spill_seen != reco.SpillNo) {
        on_new_spill = true;
        last_spill_seen = reco.SpillNo;
      }
      
      // Fill some basic "raw" variables from the reco tree
      GetHist("basic__raw__EventNo", "EventNo", "EventNo")->Fill(reco.EventNo);
      GetHist("basic__raw__SliceNo", "SliceNo", "SliceNo")->Fill(reco.SliceNo);
      GetHist("basic__raw__SpillNo", "SpillNo", "SpillNo")->Fill(reco.SpillNo);
   
      GetHist("basic__raw__ntracks", "N Reco Tracks", "ntracks")->Fill(reco.nTracks);
      for (int it = 0; it < reco.nTracks; it++) {
        GetHist("basic__raw__nHits", "nHits", "n0-120")->Fill(reco.nHits[it]);
        GetHist("basic__raw__nKalmanNodes", "nKalmanNodes", "n0-120")->Fill(reco.nKalmanNodes[it]);
        for (int ih = 0; ih < reco.nHits[it]; ih++) {
          GetHist("basic__raw__TrackHitPos_X", "TrackHitPos X", "X")->Fill(reco.TrackHitPos[it][ih][0] * CM);
          GetHist("basic__raw__TrackHitPos_Y", "TrackHitPos Y", "Y")->Fill(reco.TrackHitPos[it][ih][1] * CM);
          GetHist("basic__raw__TrackHitPos_Z", "TrackHitPos Z", "Z")->Fill(reco.TrackHitPos[it][ih][2] * CM);
        }
        for (int ih = 0; ih < reco.nKalmanNodes[it]; ih++) {
          GetHist("basic__raw__KalmanPos_X", "KalmanPos X", "X")->Fill(reco.KalmanPos[it][ih][0] * CM);
          GetHist("basic__raw__KalmanPos_Y", "KalmanPos Y", "Y")->Fill(reco.KalmanPos[it][ih][1] * CM);
          GetHist("basic__raw__KalmanPos_Z", "KalmanPos Z", "Z")->Fill(reco.KalmanPos[it][ih][2] * CM);
          GetHist("basic__raw__KalmanTruePos_X", "KalmanTruePos X", "X")->Fill(reco.KalmanTruePos[it][ih][0] * CM);
          GetHist("basic__raw__KalmanTruePos_Y", "KalmanTruePos Y", "Y")->Fill(reco.KalmanTruePos[it][ih][1] * CM);
          GetHist("basic__raw__KalmanTruePos_Z", "KalmanTruePos Z", "Z")->Fill(reco.KalmanTruePos[it][ih][2] * CM);
        }
        GetHist("basic__raw__StartDirection_X", "StartDirection X", "dx")->Fill(reco.StartDirection[it][0] * CM);
        GetHist("basic__raw__StartDirection_Y", "StartDirection Y", "dy")->Fill(reco.StartDirection[it][1] * CM);
        GetHist("basic__raw__StartDirection_Z", "StartDirection Z", "dz")->Fill(reco.StartDirection[it][2] * CM);
        GetHist("basic__raw__StartDirection_XZ", "StartDirection", "direction_xz")->Fill(reco.StartDirection[it][0] / reco.StartDirection[it][2]);
        GetHist("basic__raw__EndDirection_X", "EndDirection X", "dx")->Fill(reco.EndDirection[it][0] * CM);
        GetHist("basic__raw__EndDirection_Y", "EndDirection Y", "dy")->Fill(reco.EndDirection[it][1] * CM);
        GetHist("basic__raw__EndDirection_Z", "EndDirection Z", "dz")->Fill(reco.EndDirection[it][2] * CM);
        GetHist("basic__raw__EndDirection_XZ", "EndDirections", "direction_xz")->Fill(reco.EndDirection[it][0] / reco.EndDirection[it][2]);
      }
      // TODO finish adding these vars
      // TODO add "fixes" to direction, etc and see if that fixes things
      
      /*
   // Used to be Direction, now is StartDirection, check for both options depending on the file
   //if (HasBranch("Direction")) fChain->SetBranchAddress("Direction", Direction, &b_Direction);
   //else fChain->SetBranchAddress("StartDirection", Direction, &b_Direction);
   fChain->SetBranchAddress("StartDirection", StartDirection, &b_StartDirection);
   fChain->SetBranchAddress("EndDirection", EndDirection, &b_EndDirection);
   fChain->SetBranchAddress("StartPos", StartPos, &b_StartPos);
   fChain->SetBranchAddress("EndPos", EndPos, &b_EndPos);
   fChain->SetBranchAddress("EnergyRange", EnergyRange, &b_EnergyRange);
   fChain->SetBranchAddress("EnergyDeposit", EnergyDeposit, &b_EnergyDeposit);
   fChain->SetBranchAddress("Momentum", Momentum, &b_Momentum);
   fChain->SetBranchAddress("Length", Length, &b_Length);
   fChain->SetBranchAddress("Charge", Charge, &b_Charge);
   fChain->SetBranchAddress("TrackHitEnergies", TrackHitEnergies, &b_RecoTrackHitEnergies); */
      

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
        
        if (reco.nHits[itrack] > 0) {
          const int LOOKBACK_WINDOW = 10;
          int n_nodes = n_nonzero > LOOKBACK_WINDOW ? LOOKBACK_WINDOW : n_nonzero;
          n_nodes -= 1;
          reco.StartDirection[itrack][0] = reco.KalmanPos[itrack][0][0] - reco.KalmanPos[itrack][n_nodes][0];
          reco.StartDirection[itrack][1] = reco.KalmanPos[itrack][0][1] - reco.KalmanPos[itrack][n_nodes][1];
          reco.StartDirection[itrack][2] = reco.KalmanPos[itrack][0][2] - reco.KalmanPos[itrack][n_nodes][2];

          int last_index = reco.nKalmanNodes[itrack] - 1;
          reco.EndDirection[itrack][0] = reco.KalmanPos[itrack][last_index][0] - reco.KalmanPos[itrack][last_index - n_nodes][0];
          reco.EndDirection[itrack][1] = reco.KalmanPos[itrack][last_index][1] - reco.KalmanPos[itrack][last_index - n_nodes][1];
          reco.EndDirection[itrack][2] = reco.KalmanPos[itrack][last_index][2] - reco.KalmanPos[itrack][last_index - n_nodes][2];
        }
      }
      
      for (int it = 0; it < reco.nTracks; it++) {
        GetHist("basic__fixed__nKalmanNodes", "nKalmanNodes", "n0-120")->Fill(reco.nKalmanNodes[it]);
        for (int ih = 0; ih < reco.nKalmanNodes[it]; ih++) {
          GetHist("basic__fixed__KalmanPos_X", "KalmanPos X", "X")->Fill(reco.KalmanPos[it][ih][0] * CM);
          GetHist("basic__fixed__KalmanPos_Y", "KalmanPos Y", "Y")->Fill(reco.KalmanPos[it][ih][1] * CM);
          GetHist("basic__fixed__KalmanPos_Z", "KalmanPos Z", "Z")->Fill(reco.KalmanPos[it][ih][2] * CM);
          GetHist("basic__fixed__KalmanTruePos_X", "KalmanTruePos X", "X")->Fill(reco.KalmanTruePos[it][ih][0] * CM);
          GetHist("basic__fixed__KalmanTruePos_Y", "KalmanTruePos Y", "Y")->Fill(reco.KalmanTruePos[it][ih][1] * CM);
          GetHist("basic__fixed__KalmanTruePos_Z", "KalmanTruePos Z", "Z")->Fill(reco.KalmanTruePos[it][ih][2] * CM);
        }
        GetHist("basic__fixed__StartDirection_X", "StartDirection X", "dx")->Fill(reco.StartDirection[it][0] * CM);
        GetHist("basic__fixed__StartDirection_Y", "StartDirection Y", "dy")->Fill(reco.StartDirection[it][1] * CM);
        GetHist("basic__fixed__StartDirection_Z", "StartDirection Z", "dz")->Fill(reco.StartDirection[it][2] * CM);
        GetHist("basic__fixed__StartDirection_XZ", "StartDirection", "direction_xz")->Fill(reco.StartDirection[it][0] / reco.StartDirection[it][2]);
        GetHist("basic__fixed__EndDirection_X", "EndDirection X", "dx")->Fill(reco.EndDirection[it][0] * CM);
        GetHist("basic__fixed__EndDirection_Y", "EndDirection Y", "dy")->Fill(reco.EndDirection[it][1] * CM);
        GetHist("basic__fixed__EndDirection_Z", "EndDirection Z", "dz")->Fill(reco.EndDirection[it][2] * CM);
        GetHist("basic__fixed__EndDirection_XZ", "EndDirections", "direction_xz")->Fill(reco.EndDirection[it][0] / reco.EndDirection[it][2]);
        
        
        
      }
      
      if (on_new_spill) {
        GetHist("basic__truth__nTrueParticles", "nTrueParticles", "n0-500")->Fill(truth.nTrueParticles);
        GetHist("basic__truth__nTruePrimaryParticles", "nTruePrimaryParticles", "n0-500")->Fill(truth.nTruePrimaryParticles);
        GetHist("basic__truth__nTrueForgottenParticles", "nTrueForgottenParticles", "n0-120")->Fill(truth.nTrueForgottenParticles);
        for (int ip = 0; ip < truth.nTrueParticles; ip++) {
          GetHist("basic__truth__PDG", "PDG", "pdg")->Fill(PDGtoIndex(truth.PDG[ip]));
          if (truth.IsPrimary[ip]) GetHist("basic__truth__PDG_Primary", "PDG Primary Particles", "pdg")->Fill(PDGtoIndex(truth.PDG[ip]));
          if (!truth.IsPrimary[ip]) GetHist("basic__truth__PDG_Secondary", "PDG Secondary Particles", "pdg")->Fill(PDGtoIndex(truth.PDG[ip]));
          
          bool ismuon = abs(truth.PDG[ip]) == 13;
          // TODO add truth cut that we started in LAr fiducial and ended in TMS
          if (ismuon) {
            double muon_starting_ke = truth.MomentumTMSStart[ip][3] * 1e-3;
            if (muon_starting_ke > 5.1) muon_starting_ke = 5.05;
            GetHist("reco_eff__muon_ke_tms_enter_denominator", "Reco Eff Muon KE Entering TMS: Denominator", 
              "ke_tms_enter")->Fill(muon_starting_ke);
            double muon_starting_angle = std::atan2(truth.MomentumTMSStart[ip][0], truth.MomentumTMSStart[ip][2]) * DEG;
            GetHist("reco_eff__angle_tms_enter_denominator", "Reco Eff Muon Angle Entering TMS: Denominator", 
              "angle_tms_enter")->Fill(muon_starting_angle);
          }
        }
      }
      
      // Fill numerators of reco_eff plots here
      for (int it = 0; it < reco.nTracks; it++) {
        // TODO add cut that it starts in LAr and ends in TMS
        bool ismuon = abs(truth.RecoTrackPrimaryParticlePDG[it]) == 13;
        if (ismuon) {
          double muon_starting_ke = truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][3] * 1e-3;
          if (muon_starting_ke > 5.1) muon_starting_ke = 5.05;
          GetHist("reco_eff__muon_ke_tms_enter_numerator", "Reco Eff Muon KE Entering TMS: Numerator", 
            "ke_tms_enter")->Fill(muon_starting_ke);
          double muon_starting_angle = std::atan2(truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][0],
                                                  truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][2]) * DEG;
          GetHist("reco_eff__angle_tms_enter_numerator", "Reco Eff Muon Angle Entering TMS: Numerator", 
            "angle_tms_enter")->Fill(muon_starting_angle);
        }
      }
      
      // Example of drawing for a reason
      if (reco.nTracks > 2 && false) {
        DrawSlice(TString::Format("entry_%lld", entry_number).Data(), "high_reco_track_multiplicity", 
                  TString::Format("n tracks = %d", reco.nTracks).Data(), reco, lc, truth);
      }

      // TODO calculate plane number and then check per plane occupancy
      // Also related is total visible energy so compare that to true value
      
       
      //std::cout<<entry_number<<": "<<reco.SpillNo<<", "<<reco.SliceNo<<", "<<reco.EventNo<<std::endl;
      if (lc.nHits > 100) {
        //DrawSlice(TString::Format("final_%d", entry_number).Data(), "final", "test job", reco, lc, truth);
      }
      //if (entry_number > 700) exit(1); // TODO delete

    } // End for loop over entries
    
    
    auto time_stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(time_stop - time_start).count();
    
    // Now save the hists
    outputFile.Write();

    auto entries_visited = entry_number;
    double avg_time = duration / ((double) entries_visited);

    std::cout<<"Loop took "<<avg_time<<"us per event. "<<duration<<"us total for "<<entries_visited<<" entries"<<std::endl;
    return entries_visited;
}

std::string getOutputFilename(const std::string& inputFilename) {
    // Find the position of the last occurrence of '/' in the input filename
    size_t pos = inputFilename.find_last_of('/');
    // Extract the filename without the directory structure
    std::string filename = (pos != std::string::npos) ? inputFilename.substr(pos + 1) : inputFilename;
    return filename;
}

std::string getOutputDirname(const std::string& outputFilename) {
    // Find the position of the last occurrence of '.' in the output filename
    size_t pos = outputFilename.find_last_of('.');
    std::string filename = (pos != std::string::npos) ? outputFilename.substr(0, pos) : outputFilename;
    return filename + "_images/";
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
    
    save_location = getOutputDirname(outputFilename);

    // Create TFile with the output filename
    TFile outputFile(outputFilename.c_str(), "RECREATE");

    PrimaryLoop(ti, ri, li, numEvents, outputFile);

    // Close the output file
    outputFile.Close();

    std::cout << "Output file created: " << outputFilename << std::endl;

    return 0;
}
