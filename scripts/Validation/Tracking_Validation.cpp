#include <chrono> // for high res clock
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h> /* atan2 */
#include <stdexcept>
#include <tuple>
#include <vector>
#define TAU 6.283185307179586

// Root specific
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPaveText.h>
#include <TString.h>
#include <TStyle.h>
#include <TChain.h>
#include <TVector3.h>

// dune specific
#include <DUNEStyle.h>

#include "Line_Candidates.h"
#include "Reco_Tree.h"
#include "Truth_Info.h"
#include "Truth_Spill.h"

#define IS_WITHIN(x, center, tolerance)                                        \
  (std::abs((x) - (center)) <= (tolerance))
#define print(var) std::cout << #var << " = " << var << std::endl;

std::string save_location = "";

const bool only_nd_physics = false;

const double CM = 0.1; // cm per mm
const double DEG = 360 / TAU;
const double GEV = 1e-3; // GeV per MEV

#define length_to_energy_clarence(l) (l * 1.75 + 82) * 1e-3
#define default_length_to_energy(l) (l * 1.75) * 1e-3
// #define length_to_energy(l) default_length_to_energy(l) // No fit
// #define GEOM_V3 // for old geom
// #define length_to_energy(l) (l * 1.75 * 0.951 + 76.8) * 1e-3 // Old geom
// #define length_to_energy(l) (l * 1.75 * 0.9428 + 18.73) * 1e-3 // New geom
// #define length_to_energy(l) (1.75*l*0.9931 + 101.2445)*1e-3 // 2025-04-10
// processing, hybrid, et al. #define length_to_energy(l) (l*1.75*0.9460
// + 91.80)*1e-3 // New geom, no
// kalman
// #define length_to_energy(l) default_length_to_energy(l) // default for
// fitting #define length_to_energy(l) (1.75*l*0.9528 + 56.4912)*1e-3 //
// 2025-04-10 stereo, simple length3d
#define lar_length_to_energy(l) (l) * 1e-3

//#define GEOM_V3

namespace energy_function {
enum energy_function {
  none,
  old_geom,
  new_geom,
  length3d_20250410_stereo,
  f7_PDR_3_2b,
  f7_5v6h_4_2d,
  f7_5v6h_4_2c,
  f32b,
  f42c,
  f42d,
  test_reco_truth_override,
  test_kalman_refit,
};
}
energy_function::energy_function energy_function_to_use = energy_function::none;
double length_to_energy(double l) {
  if (energy_function_to_use == energy_function::none)
    return default_length_to_energy(l);
  if (energy_function_to_use == energy_function::old_geom)
    return (l * 1.75 * 0.951 + 76.8) * 1e-3;
  if (energy_function_to_use == energy_function::new_geom)
    return (l * 1.75 * 0.9428 + 18.73) * 1e-3;
  if (energy_function_to_use == energy_function::length3d_20250410_stereo)
    return (1.75 * l * 0.9528 + 56.4912) * 1e-3;
  if (energy_function_to_use == energy_function::f7_PDR_3_2b)
    return (1.75 * l * 1.0736 + 62.9978) * 1e-3;
  if (energy_function_to_use == energy_function::f7_5v6h_4_2d)
    return (1.75 * l * 0.9461 + 94.7684) * 1e-3;
  if (energy_function_to_use == energy_function::f7_5v6h_4_2c)
    return (1.75 * l * 1.0501 + 47.5804) * 1e-3;
  if (energy_function_to_use == energy_function::f32b)
    // return (1.75 * l * 1.0583 + 78.4743) * 1e-3;
    return (1.75 * l * 0.9258 + 46.6447) * 1e-3;
  if (energy_function_to_use == energy_function::f42c)
    // return (1.75 * l * 1.0501 + 47.5804) * 1e-3;
    return (1.75 * l * 0.9233 + 18.6727) * 1e-3;
  if (energy_function_to_use == energy_function::f42d)
    // return (1.75 * l * 0.9708 + 71.8783) * 1e-3;
    return (1.75 * l * 0.9208 + 22.7190) * 1e-3;
  if (energy_function_to_use == energy_function::test_reco_truth_override)
    // return (1.75*l*0.9110 + 119.4421)*1e-3;
    return (1.75 * l * 0.9132 + 117.5892) * 1e-3;
  if (energy_function_to_use == energy_function::test_kalman_refit)
    return (1.75*l*0.9577 + 63.4506)*1e-3;
  return default_length_to_energy(l);
}

void SetEnergyFunctionBasedOnInputFile(const std::string &inputFilename) {
  auto starting_function = energy_function_to_use;
  if (inputFilename.find("MicroProdN3p1_2024-12-18_candidate") !=
      std::string::npos)
    energy_function_to_use = energy_function::old_geom;
  if (inputFilename.find("MicroProdN1p2_2024-12-18_candidate") !=
      std::string::npos)
    energy_function_to_use = energy_function::old_geom;
  if (inputFilename.find("7_PDR_3_2b") != std::string::npos)
    energy_function_to_use = energy_function::f7_PDR_3_2b;
  if (inputFilename.find("7_5v6h_4_2d") != std::string::npos)
    energy_function_to_use = energy_function::f7_5v6h_4_2d;
  if (inputFilename.find("7_5v6h_4_2c") != std::string::npos)
    energy_function_to_use = energy_function::f7_5v6h_4_2c;
  if (inputFilename.find("32b") != std::string::npos)
    energy_function_to_use = energy_function::f32b;
  if (inputFilename.find("42c") != std::string::npos)
    energy_function_to_use = energy_function::f42c;
  if (inputFilename.find("42d") != std::string::npos)
    energy_function_to_use = energy_function::f42d;
  if (inputFilename.find("2025-10-18_finalize_kalman_v2") != std::string::npos)
    energy_function_to_use = energy_function::test_kalman_refit;
  if (inputFilename.find("2025-07-23_test_jdkio_reco_truth_override") !=
      std::string::npos)
    energy_function_to_use = energy_function::test_reco_truth_override;
  if (energy_function_to_use != starting_function)
    std::cout << "Set energy_function_to_use to " << energy_function_to_use
              << " based on inputFilename: " << inputFilename << std::endl;
  else
    std::cout << "Using default energy function" << std::endl;
}

const double MINIMUM_VISIBLE_ENERGY = 5; // MeV

int GetHitLocationCodeSingle(float x, bool isx) {
  bool zero = IS_WITHIN(x, 0, 1);
  bool is999 = IS_WITHIN(x, -999, 1) || IS_WITHIN(x, -9999, 1) ||
               IS_WITHIN(x, -99999, 1) || IS_WITHIN(x, -999999, 1) ||
               IS_WITHIN(x, -999999, 1);
  bool crazy_small = x < -4000;
  bool ltTMS = (isx) ? (x < -3300.0) : (x < 11362.0);
  bool gtTMS = (isx) ? (x > 3300.0) : (x > 18314.0);
  bool xlooksz = (isx) ? IS_WITHIN(x, 18314, 10) : false;
  int out = -99999;
  if (zero)
    out = 0;
  else if (is999)
    out = -2;
  else if (crazy_small)
    out = -3;
  else if (ltTMS)
    out = -1;
  else if (gtTMS)
    out = 1;
  else if (xlooksz)
    out = 3;
  else
    out = 2;
  return out;
}

enum LArSubRegion { LArFull, LArFiducial, LArUsingShellBuffer };

bool IsInLAr(TVector3 position, LArSubRegion subregion) {
  // LAr_detector_dims = {'x': (-350 , 350 ),'y': (-46.979-171.21,
  // -46.979+171.21),'z': (667-251.24,667+251.24)} # Active
  // Lar_fiducial_dimension = {'x': (-350+50 , 350-50 ),'y': (-46.979-171.21+50,
  // -46.979+171.21-50),'z': (667-251.24+50,667+251.24-150)} # Fiducial
  bool out = true;
  double buffer_xy = 0;
  double buffer_z = 0;
  if (subregion == LArFiducial) {
    buffer_xy = 50;
    buffer_z = 150;
  }
  if (subregion == LArUsingShellBuffer) {
    buffer_xy = 300;
    buffer_z = 30;
  }
  const double lar_y_start = -46.979 - 171.21;
  const double lar_y_end = -46.979 + 171.21;
  const double lar_z_start = 667 - 251.24;
  const double lar_z_end = 667 + 251.24;
  if (!(-350 + buffer_xy < position.X() && position.X() < 350 - buffer_xy))
    out = false;
  if (!(lar_y_start + buffer_xy < position.Y() &&
        position.Y() < lar_y_end - buffer_xy))
    out = false;
  if (!(lar_z_start + buffer_z < position.Z() &&
        position.Z() < lar_z_end - buffer_z))
    out = false;
  return out;
}

bool IsInTMSFull(TVector3 position) {
  bool out = true;
  if (-4500 < position.X() && position.X() < 4500)
    out = false;
  if (-4000 < position.Y() && position.Y() < 1000)
    out = false;
  if (11000 < position.Z() && position.Z() < 19000)
    out = false;
  return out;
}

std::string getExecutableName(const char *arg) {
  std::string executableName = arg; // Convert to std::string
  size_t lastSlash =
      executableName.find_last_of("/\\"); // Find last occurrence of '/' or '\'
  if (lastSlash != std::string::npos) {
    executableName = executableName.substr(
        lastSlash + 1); // Extract substring after the last slash
  }
  size_t extensionPos =
      executableName.find_last_of('.'); // Find last occurrence of '.'
  if (extensionPos != std::string::npos) {
    executableName = executableName.substr(0, extensionPos); // Remove extension
  }
  return executableName;
}

bool createDirectory(const std::string &path) {
  try {
    // Create directory and its parents if they don't exist
    std::filesystem::create_directories(path);
    return true;
  } catch (const std::exception &e) {
    std::cerr << "Error creating directory: " << e.what() << std::endl;
    return false;
  }
}

#include "DrawSlice.cxx"
#include "HistManager.cxx"

// If one of these flags is set, then it fails the cut
// set flag: flags |= FOO;
// clear: flags &= ~FOO;
// check that flag is set: if (flags & FOO)
// check that flag is not set: if ((flags & FOO) == 0)
enum TMSCutFlags {
  TMSCutAny = 1 << 1, // Fails any cut
  TMSCutStartFront = 1 << 2,
  TMSCutAnyEndCuts = 1 << 3, // fails 4, 5, or 6
  TMSCutZCut = 1 << 4,
  TMSCutYCut = 1 << 5,
  TMSCutUVCut = 1 << 6
};

const double TMS_START_Z = 11362;
#ifdef GEOM_V3
const double TMS_LAST_PLANE_Z = 17900; // mm // 18314.0 - 50; // End Z - 5cm
const double TMS_THICK_START = 13500;  // mm
#else
const double TMS_LAST_PLANE_Z = 18400;
const double TMS_THICK_START = 14480;        // mm
const double TMS_DOUBLE_THICK_START = 17560; // mm
#endif
const double TMS_Y_TOP =
    160.0; // Top of scintillator (I'm 99% sure starting from where U/V overlap)
const double TMS_Y_BOTTOM =
    -2850.0; // Bottom of scinillator (also measured from U/V overlap)
const double TMS_Y_MIDDLE = 0.5 * (TMS_Y_TOP + TMS_Y_BOTTOM);
const double TMS_X_EXTENT = 3300.0; // Assumes things are symmetrical right now
const double TMS_X_START = -TMS_X_EXTENT;
const double Y_BUFFER = 400; // mm, 40 cm
const double X_BUFFER = 30;  // mm, 3cm bar thickness?
const double Z_BUFFER = 100; // mm, two thin planes
const double VIEW_ANGLE = 3.0 * TMath::DegToRad();

TVector3 GetRotatedVector(TVector3 v, bool u) {
  double angle = VIEW_ANGLE;
  if (u)
    angle *= -1; // Opposite direction
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
  if (flags == 0)
    hist.Fill(0);
  // Check and fill for each flag
  if (flags & TMSCutAny)
    hist.Fill(1);
  if (flags & TMSCutStartFront)
    hist.Fill(2);
  if (flags & TMSCutAnyEndCuts)
    hist.Fill(3);
  if (flags & TMSCutZCut)
    hist.Fill(4);
  if (flags & TMSCutYCut)
    hist.Fill(5);
  if (flags & TMSCutUVCut)
    hist.Fill(6);
}

int TMSStartPointThruFront(TVector3 startpoint) {
  // Check that it's in the first two planes or so
  int out = 0;
  if (startpoint.Z() < TMS_START_Z)
    out |= TMSCutStartFront;
  if (startpoint.Z() > TMS_START_Z + Z_BUFFER)
    out |= TMSCutStartFront;

  const int ANYCUTS = (TMSCutZCut | TMSCutYCut | TMSCutUVCut |
                       TMSCutAnyEndCuts | TMSCutStartFront);
  if ((out & ANYCUTS))
    out |= TMSCutAny;

  return out;
}

int TMSEndpointUncontained(TVector3 endpoint) {
  // Returns flag representing which cut was missed
  // Zero means it passes all cuts
  int out = 0;

  // Check that muon is not in last scinillator plane because it may be exiting
  // then
  if (endpoint.Z() > TMS_LAST_PLANE_Z)
    out |= TMSCutZCut;
  // Also check that it's contained inside the start of the TMS when using truth
  // info
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
  if (std::abs(endpoint_v.X()) > TMS_X_EXTENT - X_BUFFER ||
      std::abs(endpoint_u.X()) > TMS_X_EXTENT - X_BUFFER)
    out |= TMSCutUVCut;

  // Set that all cuts passed if they did
  const int ANYENDCUTS = (TMSCutZCut | TMSCutYCut | TMSCutUVCut);
  if ((out & ANYENDCUTS))
    out |= TMSCutAnyEndCuts;

  // Set that all cuts passed if they did
  const int ANYCUTS = (TMSCutZCut | TMSCutYCut | TMSCutUVCut |
                       TMSCutAnyEndCuts | TMSCutStartFront);
  if ((out & ANYCUTS))
    out |= TMSCutAny;
  return out;
}

int CheckTMSCuts(TVector3 startpoint, TVector3 endpoint) {
  return TMSStartPointThruFront(startpoint) | TMSEndpointUncontained(endpoint);
}

TVector3 calculatePositionAtZ(double z_new, TVector3 start, double xz_dir,
                              double yz_dir = 1.0) {
  double z = z_new;
  double dz = z - start.Z();
  double x = start.X() + xz_dir * dz;
  double y = start.Y() + yz_dir * dz;
  TVector3 out(x, y, z);
  return out;
}

int isTMSContained(TVector3 position, bool thin_only = false) {
  // Returns 0 if contained. Otherwise returns code that says which cut was
  // missed
  int out = 0;
  // Z positions of the first hits of the TMS
  const double TMS_Thin_Start = 11362;
  // Where do we transition to the thick region (first layer of scintillator
  // before the change)
  const double TMS_Thick_Start = 13500;
  // Where does the thick region end
  const double TMS_Thick_End = 18000 - 100;
  const double TMS_Start_Bars_Only[] = {-3350, 240, TMS_Thin_Start};
  const double TMS_End_Bars_Only[] = {3350, -2950, TMS_Thick_End};
  if (position.X() < TMS_Start_Bars_Only[0])
    out += 1 << 0;
  if (position.X() > TMS_End_Bars_Only[0])
    out += 1 << 0;
  if (position.Y() > TMS_Start_Bars_Only[1])
    out += 1 << 1;
  if (position.Y() < TMS_End_Bars_Only[1])
    out += 1 << 1;
  if (position.Z() < TMS_Start_Bars_Only[2])
    out += 1 << 2;
  if (thin_only) {
    if (position.Z() > TMS_Thick_Start)
      out += 1 << 2;
  } else {
    if (position.Z() > TMS_End_Bars_Only[2])
      out += 1 << 2;
  }
  return out;
}

bool LArHadronCut(TVector3 position, double distance = 300) {
  bool out = true;
  const double LAr_Start_Exact[] = {-3478.48, -2166.71, 4179.24};
  const double LAr_End_Exact[] = {3478.48, 829.282, 9135.88};
  out = out && (position.X() >= LAr_Start_Exact[0] + distance) &&
        (position.X() <= LAr_End_Exact[0] - distance);
  out = out && (position.Y() >= LAr_Start_Exact[1] + distance) &&
        (position.Y() <= LAr_End_Exact[1] - distance);
  out = out && (position.Z() >= LAr_Start_Exact[2] + distance) &&
        (position.Z() <= LAr_End_Exact[2] - distance);
  return out;
}

bool LArFiducialCut(TVector3 position, double distance = 500,
                    double downstream = 1500) {
  // Returns true if inside fiducial, false otherwise
  bool out = true;
  const double LAr_Start_Exact[] = {-3478.48, -2166.71, 4179.24};
  const double LAr_End_Exact[] = {3478.48, 829.282, 9135.88};
  out = out && (position.X() >= LAr_Start_Exact[0] + distance) &&
        (position.X() <= LAr_End_Exact[0] - distance);
  out = out && (position.Y() >= LAr_Start_Exact[1] + distance) &&
        (position.Y() <= LAr_End_Exact[1] - distance);
  out = out && (position.Z() >= LAr_Start_Exact[2] + distance) &&
        (position.Z() <= LAr_End_Exact[2] - downstream);
  return out;
}

bool NDPhysicsMuon(Truth_Info &truth, Reco_Tree &reco, int track) {
  bool out = true;
  // Must start in the LAr fiducial volume
  if (!truth.RecoTrackPrimaryParticleLArFiducialStart[track])
    out = false;
  // Must end in the TMS
  if (!truth.RecoTrackPrimaryParticleTMSFiducialEnd[track])
    out = false;
  // Must be muon
  if (std::abs(truth.RecoTrackPrimaryParticlePDG[track]) != 13)
    out = false;
  // Hadonic E in LAr shell > 30 MeV
  if (truth.LArOuterShellEnergyFromVertex > 30)
    out = false;
  (void)reco;
  return out;
}

bool NDPhysicsSlice(Truth_Info &truth, Reco_Tree &reco) {
  // Checks that the entire slice qualifies as an ND physics slice
  // Ie. all reco tracks are muons from LAr to TMS
  bool out = true;
  for (int track = 0; track < reco.nTracks; track++) {
    out = out && NDPhysicsMuon(truth, reco, track);
  }
  return out;
}

struct ValidationContext {
  Truth_Info &truth;
  Truth_Spill &truth_spill;
  Reco_Tree &reco;
  Line_Candidates &lc;
  bool on_new_spill;
  int current_tree;
  int current_spill;
  Long64_t entry_number;
  std::map<std::string, int> &particle_fingerprints_reconstructed;
};

#include "Reco_Eff.cxx"

Long64_t PrimaryLoop(Truth_Info &truth, Truth_Spill &truth_spill,
                     Reco_Tree &reco, Line_Candidates &lc,
                     int numEvents, TFile &outputFile) {
  // List all the hists here
  // Make sure to save them too

  bool has_kalman = reco.HasBranch("KalmanPos");
  std::cout << "has_kalman status: " << has_kalman << std::endl;
  has_kalman = false; // TODO this needs to be set based on the input file

  // Want to plot:
  // Reco vs RecoTrackTrueHitPosition
  // Reco track length vs:
  // RecoTrackPrimaryParticleTrueTrackLengthAsMeasured
  // RecoTrackPrimaryParticleTrueTrackLengthRecoStart
  // RecoTrackPrimaryParticleTrueTrackLengthInTMS
  // Using RecoTrackPrimaryParticleTMSFiducialStart/Touch and
  // RecoTrackPrimaryParticleTMSFiducialEnd Or
  // RecoTrackPrimaryParticleLArFiducialStart,
  // RecoTrackPrimaryParticleLArFiducialTouch, and
  // RecoTrackPrimaryParticleLArFiducialEnd

  // For particles ending in TMS, TMSFiducialStart, TMSFiducialEnd
  // RecoTrackPrimaryParticleTrueMomentumEnteringTMS

  // Secondary particle info:
  // Flag if it has one or not
  // RecoTrackSecondaryParticleTrueVisibleEnergy
  // RecoTrackPrimaryParticleTrueVisibleEnergy

  // Reco tree vars:
  // StartPos, EndPos, nHits, TrackHitPos/KalmanPos
  // EnergyRange, EnergyDeposit, Length

  int last_tree_seen = -1;
  int last_spill_seen = -1;

  auto time_start = std::chrono::high_resolution_clock::now();
  std::map<int, double> manually_calculated_shell_hadronic_energy_per_vertex;

  std::map<std::string, int> particle_fingerprints_reconstructed;
  std::map<std::pair<int, int>, Long64_t> truth_spill_entry_by_tree_and_spill;
  for (Long64_t spill_entry = 0; spill_entry < truth_spill.GetEntriesFast();
       ++spill_entry) {
    truth_spill.GetEntry(spill_entry);
    const int spill_tree =
        truth_spill.fChain ? truth_spill.fChain->GetTreeNumber() : -1;
    truth_spill_entry_by_tree_and_spill.emplace(
        std::make_pair(spill_tree, truth_spill.SpillNo), spill_entry);
  }
  std::cout << "Indexed " << truth_spill_entry_by_tree_and_spill.size()
            << " Truth_Spill entries" << std::endl;

  Long64_t entry_number = 0;
  // Now loop over the ttree
  for (; entry_number < truth.GetEntriesFast() &&
         entry_number < reco.GetEntriesFast() &&
         (numEvents < 0 || entry_number < numEvents);
       entry_number++) {
    if (entry_number % 10000 == 0)
      std::cout << "On entry: " << entry_number << std::endl;

    // Get the current entry
    // Currently reco and truth match
    truth.GetEntry(entry_number);
    reco.GetEntry(entry_number);
    lc.GetEntry(entry_number);

    if (only_nd_physics && !NDPhysicsSlice(truth, reco))
      continue;

    // Check if we're on a new spill
    // Some things should only be filled per spill
    bool on_new_spill = false;
    int current_tree = truth.fChain ? truth.fChain->GetTreeNumber() : -1;
    int current_spill = truth.SpillNo;
    bool has_truth_spill = false;
    const auto spill_entry = truth_spill_entry_by_tree_and_spill.find(
        std::make_pair(current_tree, current_spill));
    if (spill_entry != truth_spill_entry_by_tree_and_spill.end()) {
      truth_spill.GetEntry(spill_entry->second);
      has_truth_spill = true;
    }
    if (last_tree_seen != current_tree || last_spill_seen != current_spill) {
      on_new_spill = true;
      last_tree_seen = current_tree;
      last_spill_seen = current_spill;
      // Empty on each new spill
      particle_fingerprints_reconstructed.clear();
      if (!has_truth_spill) {
        std::cout << "[Tracking_Validation] Missing Truth_Spill for tree "
                  << current_tree << ", SpillNo " << current_spill
                  << std::endl;
      }
    }
    if (!has_truth_spill)
      continue;

#include "Basic.cxx"

    // Truth matching information
    REGISTER_AXIS(completeness,
                  std::make_tuple("Primary on Track / in Slice", 20, 0, 1.01));
    REGISTER_AXIS(cleanliness,
                  std::make_tuple("Track Primary / Total", 20, 0, 1.01));
    REGISTER_AXIS(nhits_in_track,
                  std::make_tuple("N Hits in Track", 200, 0, 200));
    REGISTER_AXIS(energy_in_track,
                  std::make_tuple("Energy in Track (MeV)", 100, 0, 1000));
    REGISTER_AXIS(n_hits_per_plane,
                  std::make_tuple("N Hits per Plane", 5, 0.5, 5.5));
    for (int it = 0; it < reco.nTracks; it++) {
      GetHist("reco_track__primary_pdg", "Reco Track Primary Particle PDG",
              "pdg", "#N Tracks")
          ->Fill(PDGtoIndex(truth.RecoTrackPrimaryParticlePDG[it]));
      GetHist("reco_track__secondary_pdg", "Reco Track Secondary Particle PDG",
              "pdg", "#N Tracks")
          ->Fill(PDGtoIndex(truth.RecoTrackSecondaryParticlePDG[it]));
      if (NDPhysicsMuon(truth, reco, it) &&
          truth.RecoTrackTrueVisibleEnergy[it] > 0 && truth.RecoTrackNHits[it] > 0) {
        double cleanliness_energy =
            truth.RecoTrackPrimaryParticleTrueVisibleEnergy[it] /
            truth.RecoTrackTrueVisibleEnergy[it];
        double cleanliness_nhits = truth.RecoTrackPrimaryParticleTrueNHits[it] /
                                   (double)truth.RecoTrackNHits[it];
        GetHist("reco_track__cleanliness_energy",
                "Reco Track Cleanliness, Visible Energy", "cleanliness",
                "#N Tracks")
            ->Fill(cleanliness_energy);
        GetHist("reco_track__cleanliness_nhits",
                "Reco Track Cleanliness, N Hits", "cleanliness", "#N Tracks")
            ->Fill(cleanliness_nhits);

        bool ismuon = std::abs(truth.RecoTrackPrimaryParticlePDG[it]) == 13;
        if (ismuon && reco.nTracks == 1) {
          std::map<double, int> count;
          for (int ih = 0; ih < lc.nHits; ih++) {
            count[lc.RecoHitPos[ih][2]] += 1;
            if (count[lc.RecoHitPos[ih][2]] == 2)
              GetHist("reco_track__debugging__z_pos_with_two_or_more_hit_"
                      "planes",
                      "Z Pos with Two or More Hits in Plane", "Z", "#Count")
                  ->Fill(lc.RecoHitPos[ih][2] * CM);
          }
          for (auto &plane : count) {
            GetHist("reco_track__debugging__n_hits_in_plane",
                    "Hits per plane: All Reco Hits", "n_hits_per_plane",
                    "#Count")
                ->Fill(plane.second);
            GetHist("reco_track__debugging__n_hits_in_plane_both_"
                    "nostack_1_all_hits",
                    "Hits per plane: All Reco Hits", "n_hits_per_plane",
                    "#Count")
                ->Fill(plane.second);
          }

          count.clear();
          for (int ih = 0; ih < reco.nHits[0]; ih++) {
            count[reco.TrackHitPos[0][ih][2]] += 1;
            if (count[reco.TrackHitPos[0][ih][2]] == 2)
              GetHist("reco_track__debugging__z_pos_with_two_or_more_hit_"
                      "planes_in_track",
                      "Z Pos with Two or More Hits in Plane in Track", "Z",
                      "#Count")
                  ->Fill(reco.TrackHitPos[0][ih][2] * CM);
          }
          for (auto &plane : count) {
            GetHist("reco_track__debugging__n_hits_in_plane_track",
                    "Hits per plane: In Track", "n_hits_per_plane",
                    "#Count")
                ->Fill(plane.second);
            GetHist("reco_track__debugging__n_hits_in_plane_both_"
                    "nostack_2_track",
                    "Hits per plane: In Tracks", "n_hits_per_plane", "#Count")
                ->Fill(plane.second);
          }

          if (cleanliness_energy < 0.5)
            DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                      "completeness_cleanliness/cleanliness/under_50_percent",
                      TString::Format(
                          "Cleanliness by E;is under 50%%;Cleanliness = %.2f",
                          cleanliness_energy)
                          .Data(),
                      reco, lc, truth, DrawSliceN::handfull);
          if (cleanliness_energy >= 0.67 && cleanliness_energy <= 0.73)
            DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                      "completeness_cleanliness/cleanliness/around_70_percent",
                      TString::Format(
                          "Cleanliness by E;is around 70%%;Cleanliness = %.2f",
                          cleanliness_energy)
                          .Data(),
                      reco, lc, truth, DrawSliceN::handfull);
          if (cleanliness_energy >= 0.98)
            DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                      "completeness_cleanliness/cleanliness/above_98_percent",
                      TString::Format(
                          "Cleanliness by E;is above 98%%;Cleanliness = %.2f",
                          cleanliness_energy)
                          .Data(),
                      reco, lc, truth, DrawSliceN::handfull);
        }
      }
    }

    GetHist("basic__truth__LeptonX4_X", "Lepton X4 X", "X")
        ->Fill(truth.LeptonX4[0] * CM);
    GetHist("basic__truth__LeptonX4_Y", "Lepton X4 Y", "Y_full")
        ->Fill(truth.LeptonX4[1] * CM);
    GetHist("basic__truth__LeptonX4_Z", "Lepton X4 Z", "Z_full")
        ->Fill(truth.LeptonX4[2] * CM);

    if (truth.HasBranch("LArOuterShellEnergy")) {

      REGISTER_AXIS(
          LArOuterShellEnergy,
          std::make_tuple("LAr Visible Hadronic Energy in 30cm Shell (MeV)", 51,
                          0, 10000));
      if (on_new_spill)
        GetHist("basic__truth__LArOuterShellEnergy",
                "LAr Visible Hadronic Energy in 30cm Shell",
                "LArOuterShellEnergy")
            ->Fill(truth.LArOuterShellEnergy);
      REGISTER_AXIS(
          LArOuterShellEnergyFromVertex,
          std::make_tuple("LAr Visible Hadronic Energy in 30cm Shell (MeV)", 20,
                          0, 200));
      GetHist("basic__truth__LArOuterShellEnergyFromVertex",
              "LAr Visible Hadronic Energy from Primary Vertex in 30cm Shell",
              "LArOuterShellEnergyFromVertex")
          ->Fill(truth.LArOuterShellEnergyFromVertex);

      REGISTER_AXIS(LArTotalEnergy,
                    std::make_tuple("LAr Visible Energy (MeV)", 51, 0, 10000));
      if (on_new_spill)
        GetHist("basic__truth__LArTotalEnergy", "LAr Visible Energy",
                "LArTotalEnergy")
            ->Fill(truth.LArTotalEnergy);
      REGISTER_AXIS(LArTotalEnergyFromVertex,
                    std::make_tuple("LAr Visible Energy (MeV)", 51, 0, 10000));
      GetHist("basic__truth__LArTotalEnergyFromVertex",
              "LAr Visible Energy from Primary Vertex",
              "LArTotalEnergyFromVertex")
          ->Fill(truth.LArTotalEnergyFromVertex);

      REGISTER_AXIS(
          TotalNonTMSEnergy,
          std::make_tuple("Total E Deposited Outside TMS (MeV)", 51, 0, 10000));
      if (on_new_spill)
        GetHist("basic__truth__TotalNonTMSEnergy",
                "Total E Deposited Outside TMS", "TotalNonTMSEnergy")
            ->Fill(truth.TotalNonTMSEnergy);
      REGISTER_AXIS(TotalNonTMSEnergyFromVertex,
                    std::make_tuple("Total E Deposited from Outside TMS (MeV)",
                                    51, 0, 10000));
      GetHist("basic__truth__TotalNonTMSEnergyFromVertex",
              "Total E Deposited from Primary Vertex Outside TMS",
              "TotalNonTMSEnergyFromVertex")
          ->Fill(truth.TotalNonTMSEnergyFromVertex);

      if (on_new_spill)
        GetHist("nd_physics_cut__lar_outer_shell__energy_in_shell",
                "LAr Visible Energy in 30cm Shell", "LArOuterShellEnergy")
            ->Fill(truth.LArOuterShellEnergy);
      GetHist("nd_physics_cut__lar_outer_shell__energy_in_shell_from_vertex",
              "LAr Visible Energy from Primary Vertex in 30cm Shell",
              "LArOuterShellEnergyFromVertex")
          ->Fill(truth.LArOuterShellEnergyFromVertex);
      if (on_new_spill)
        GetHist("nd_physics_cut__lar_outer_shell__passes_cut_e_total",
                "Passes Total Hadronic E in ND-LAr Outer Shell < 30 MeV Cut",
                "yesno")
            ->Fill((truth.LArOuterShellEnergy < 30) ? 0 : 1);
      GetHist(
          "nd_physics_cut__lar_outer_shell__passes_cut_e_from_vertex",
          "Passes Hadronic E from Vertex in ND-LAr Outer Shell < 30 MeV Cut",
          "yesno")
          ->Fill((truth.LArOuterShellEnergyFromVertex < 30) ? 0 : 1);

      bool passes_ndlar_outer_shell_cut = truth.LArTotalEnergyFromVertex < 30;
      bool passes_ndlar_fiducial_cut = false;

      REGISTER_AXIS(lar_x, std::make_tuple("X (cm)", 51, -400, 400));
      REGISTER_AXIS(lar_y, std::make_tuple("Y (cm)", 51, -250, 150));
      REGISTER_AXIS(lar_z, std::make_tuple("Z (cm)", 51, 300, 1000));
      {
        TVector3 birth_pos(truth.LeptonX4[0], truth.LeptonX4[1],
                           truth.LeptonX4[2]);
        passes_ndlar_fiducial_cut = LArFiducialCut(birth_pos);

        if (passes_ndlar_fiducial_cut)
          GetHist("nd_physics_cut__lar_outer_shell__fiducial_passes_cut_e_from_"
                  "vertex",
                  "Passes Hadronic E from Vertex in ND-LAr Outer Shell < 30 "
                  "MeV Cut for Fiducial LAr Interactions",
                  "yesno")
              ->Fill((truth.LArOuterShellEnergyFromVertex < 30) ? 0 : 1);
        if (passes_ndlar_fiducial_cut)
          GetHist("nd_physics_cut__lar_outer_shell__fiducial_energy_in_shell_"
                  "from_vertex",
                  "LAr Visible Energy from Primary Vertex in 30cm Shell for "
                  "Fiducial LAr Interactions",
                  "LArOuterShellEnergyFromVertex")
              ->Fill(truth.LArOuterShellEnergyFromVertex);

        GetHist("nd_physics_cut__lar_fiducial__X_nostack_1_nocut", "X: All",
                "lar_x")
            ->Fill(truth.LeptonX4[0] * CM);
        GetHist("nd_physics_cut__lar_fiducial__Y_nostack_1_nocut", "Y: All",
                "lar_y")
            ->Fill(truth.LeptonX4[1] * CM);
        GetHist("nd_physics_cut__lar_fiducial__Z_nostack_1_nocut", "Z: All",
                "lar_z")
            ->Fill(truth.LeptonX4[2] * CM);
        if (passes_ndlar_fiducial_cut) {
          GetHist("nd_physics_cut__lar_fiducial__X_nostack_2_withcut",
                  "X: Passes", "lar_x")
              ->Fill(truth.LeptonX4[0] * CM);
          GetHist("nd_physics_cut__lar_fiducial__Y_nostack_2_withcut",
                  "Y: Passes", "lar_y")
              ->Fill(truth.LeptonX4[1] * CM);
          GetHist("nd_physics_cut__lar_fiducial__Z_nostack_2_withcut",
                  "Z: Passes", "lar_z")
              ->Fill(truth.LeptonX4[2] * CM);
        }

        GetHist("nd_physics_cut__lar_outer_shell__shell_X_nostack_1_all",
                "X: All", "lar_x")
            ->Fill(truth.LeptonX4[0] * CM);
        GetHist("nd_physics_cut__lar_outer_shell__shell_Y_nostack_1_all",
                "Y: All", "lar_y")
            ->Fill(truth.LeptonX4[1] * CM);
        GetHist("nd_physics_cut__lar_outer_shell__shell_Z_nostack_1_all",
                "Z: All", "lar_z")
            ->Fill(truth.LeptonX4[2] * CM);
        if (passes_ndlar_outer_shell_cut) {
          GetHist("nd_physics_cut__lar_outer_shell__shell_X_nostack_2_passes",
                  "X: Passes", "lar_x")
              ->Fill(truth.LeptonX4[0] * CM);
          GetHist("nd_physics_cut__lar_outer_shell__shell_Y_nostack_2_passes",
                  "Y: Passes", "lar_y")
              ->Fill(truth.LeptonX4[1] * CM);
          GetHist("nd_physics_cut__lar_outer_shell__shell_Z_nostack_2_passes",
                  "Z: Passes", "lar_z")
              ->Fill(truth.LeptonX4[2] * CM);
        } else {
          GetHist("nd_physics_cut__lar_outer_shell__shell_X_nostack_3_fails",
                  "X: Fails", "lar_x")
              ->Fill(truth.LeptonX4[0] * CM);
          GetHist("nd_physics_cut__lar_outer_shell__shell_Y_nostack_3_fails",
                  "Y: Fails", "lar_y")
              ->Fill(truth.LeptonX4[1] * CM);
          GetHist("nd_physics_cut__lar_outer_shell__shell_Z_nostack_3_fails",
                  "Z: Fails", "lar_z")
              ->Fill(truth.LeptonX4[2] * CM);
        }

        if (passes_ndlar_fiducial_cut) {
          GetHist("nd_physics_cut__lar_outer_shell__fid_shell_X_nostack_1_all",
                  "X: All", "lar_x")
              ->Fill(truth.LeptonX4[0] * CM);
          GetHist("nd_physics_cut__lar_outer_shell__fid_shell_Y_nostack_1_all",
                  "Y: All", "lar_y")
              ->Fill(truth.LeptonX4[1] * CM);
          GetHist("nd_physics_cut__lar_outer_shell__fid_shell_Z_nostack_1_all",
                  "Z: All", "lar_z")
              ->Fill(truth.LeptonX4[2] * CM);
          if (passes_ndlar_outer_shell_cut) {
            GetHist(
                "nd_physics_cut__lar_outer_shell__fid_shell_X_nostack_2_passes",
                "X: Passes", "lar_x")
                ->Fill(truth.LeptonX4[0] * CM);
            GetHist(
                "nd_physics_cut__lar_outer_shell__fid_shell_Y_nostack_2_passes",
                "Y: Passes", "lar_y")
                ->Fill(truth.LeptonX4[1] * CM);
            GetHist(
                "nd_physics_cut__lar_outer_shell__fid_shell_Z_nostack_2_passes",
                "Z: Passes", "lar_z")
                ->Fill(truth.LeptonX4[2] * CM);
          } else {
            GetHist(
                "nd_physics_cut__lar_outer_shell__fid_shell_X_nostack_3_fails",
                "X: Fails", "lar_x")
                ->Fill(truth.LeptonX4[0] * CM);
            GetHist(
                "nd_physics_cut__lar_outer_shell__fid_shell_Y_nostack_3_fails",
                "Y: Fails", "lar_y")
                ->Fill(truth.LeptonX4[1] * CM);
            GetHist(
                "nd_physics_cut__lar_outer_shell__fid_shell_Z_nostack_3_fails",
                "Z: Fails", "lar_z")
                ->Fill(truth.LeptonX4[2] * CM);
          }
        }
      }

      // First find the max true reco'd muon starting ke
      int highest_muon_track = -1;
      double max_muon_starting_ke = -999999;
      for (int it = 0; it < reco.nTracks; it++) {
        if (std::abs(truth.RecoTrackPrimaryParticlePDG[it]) == 13) {
          double muon_starting_ke =
              truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][3] *
              1e-3;
          if (muon_starting_ke > max_muon_starting_ke) {
            max_muon_starting_ke = muon_starting_ke;
            highest_muon_track = it;
          }
        }
      }
      if (highest_muon_track >= 0) {
        passes_ndlar_fiducial_cut =
            truth.RecoTrackPrimaryParticleLArFiducialStart[highest_muon_track];
      }
      // Now fill the plots if we found one
      if (max_muon_starting_ke >= 0) {
        GetHist(
            "nd_physics_cut__lar_outer_shell__muon_ke_nostack_1_no_cut_width",
            "Muon KE: No Cut", "ke_tms_enter")
            ->Fill(max_muon_starting_ke);
        if (passes_ndlar_fiducial_cut)
          GetHist("nd_physics_cut__lar_outer_shell__muon_ke_nostack_2_fiducial_"
                  "cut_width",
                  "Muon KE: ND LAr Fiducial Volume Cut", "ke_tms_enter")
              ->Fill(max_muon_starting_ke);
        if (passes_ndlar_outer_shell_cut)
          GetHist("nd_physics_cut__lar_outer_shell__muon_ke_nostack_3_shell_"
                  "cut_width",
                  "Muon KE: Outer Shell < 30 MeV Cut", "ke_tms_enter")
              ->Fill(max_muon_starting_ke);
        if (passes_ndlar_outer_shell_cut && passes_ndlar_fiducial_cut)
          GetHist("nd_physics_cut__lar_outer_shell__muon_ke_nostack_4_shell_"
                  "fiducial_cuts_width",
                  "Muon KE: Shell and Fiducial Cuts", "ke_tms_enter")
              ->Fill(max_muon_starting_ke);
      }

    } // end if (truth.HasBranch("LArOuterShellEnergy")) {

    REGISTER_AXIS(
        hit_position_resolution_x,
        std::make_tuple("Hit Position Resolution X (Reco - True) (cm)", 21, -4,
                        4));
    REGISTER_AXIS(
        hit_position_resolution_x_wide,
        std::make_tuple("Hit Position Resolution X (Reco - True) (cm)", 31, -30,
                        30));
    REGISTER_AXIS(
        hit_position_resolution_y,
        std::make_tuple("Hit Position Resolution Y (Reco - True) (cm)", 81, -40,
                        40));
    REGISTER_AXIS(
        slope_yz_resolution_track_start,
        std::make_tuple("Starting YZ Slope Resolution (Reco - True) (deg)", 21,
                        -40, 40));
    REGISTER_AXIS(
        slope_yz_resolution_track_end,
        std::make_tuple("Ending YZ Slope Resolution (Reco - True) (deg)", 21,
                        -40, 40));
    if (reco.nTracks == 1) {
      if (lc.nLinesX == 1) {
        // mm max distance in z
        const double z_eps = 50;
        double post_line_z = -1e9;
        double pre_line_z = 1e9;
        std::map<int, int> mapOfXLineIndices;
        for (int iH = 0; iH < reco.nHits[0]; iH++) {
          double true_z = truth.RecoTrackTrueHitPosition[0][iH][2];
          mapOfXLineIndices[iH] = -1;
          for (int ih = 0; ih < lc.nHitsInTrackX[0]; ih++) {
            double reco_z = lc.TrackHitPosX[0][ih][0];
            if (reco_z < pre_line_z)
              pre_line_z = reco_z;
            if (reco_z > post_line_z)
              post_line_z = reco_z;
            double dz = std::abs(reco_z - true_z);
            if (dz < z_eps) {
              mapOfXLineIndices[iH] = ih;
              break;
            }
          }
        }

        for (int iH = 0; iH < reco.nHits[0]; iH++) {
          int ih = mapOfXLineIndices[iH];
          double true_y = truth.RecoTrackTrueHitPosition[0][iH][1];
          if (ih >= 0) {
            double reco_y = lc.TrackHitPosX[0][ih][1];
            double dy = reco_y - true_y;
            GetHist("resolution__reco_track__lines_2d__line_x_hit_resolution_y",
                    "Reco Track Hit Resolution Y", "hit_position_resolution_y",
                    "#N Hits")
                ->Fill(dy * CM);
            GetHist("resolution__reco_track__lines_2d__line_x_hit_resolution_y_"
                    "v2_nostack_1_has_x",
                    "Reco Track Hit Resolution Y: Has X info",
                    "hit_position_resolution_y", "#N Hits")
                ->Fill(dy * CM);
          } else {
            double reco_z = reco.TrackHitPos[0][iH][2];
            double dy = reco.TrackHitPos[0][iH][1] - true_y;

            if (reco_z < pre_line_z)
              GetHist("resolution__reco_track__lines_2d__line_x_hit_resolution_"
                      "y_v2_nostack_3_pre_line",
                      "Reco Track Hit Resolution Y: Before X Line",
                      "hit_position_resolution_y", "#N Hits")
                  ->Fill(dy * CM);
            else if (reco_z > post_line_z)
              GetHist("resolution__reco_track__lines_2d__line_x_hit_resolution_"
                      "y_v2_nostack_4_post_line",
                      "Reco Track Hit Resolution Y: After X Line",
                      "hit_position_resolution_y", "#N Hits")
                  ->Fill(dy * CM);
            else
              GetHist("resolution__reco_track__lines_2d__line_x_hit_resolution_"
                      "y_v2_nostack_2_missing_x",
                      "Reco Track Hit Resolution Y: Missing X info",
                      "hit_position_resolution_y", "#N Hits")
                  ->Fill(dy * CM);
          }
        }
      }

      for (int ih = 0; ih < reco.nKalmanNodes[0]; ih++) {
        double dx =
            reco.KalmanPos[0][ih][0] - reco.KalmanTruePos[0][ih][0];
        GetSpecialHist("special__kalman_node_resolution_x",
                       "resolution__reco_track__kalman_node_resolution_x",
                       "Reco Kalman Node Resolution X",
                       "hit_position_resolution_x_wide", "#N Hits")
            ->Fill(dx * CM);
        double dy =
            reco.KalmanPos[0][ih][1] - reco.KalmanTruePos[0][ih][1];
        GetSpecialHist("special__kalman_node_resolution_y",
                       "resolution__reco_track__kalman_node_resolution_y",
                       "Reco Kalman Node Resolution Y",
                       "hit_position_resolution_y", "#N Hits")
            ->Fill(dy * CM);
      }

      double avg_offset = 0;
      double largest_dx_deviation = 0;
      for (int ih = 0; ih < reco.nHits[0]; ih++) {
        double dx = reco.TrackHitPos[0][ih][0] -
                    truth.RecoTrackTrueHitPosition[0][ih][0];
        GetSpecialHist("special__track_hit_resolution_x",
                       "resolution__reco_track__hit_resolution_x",
                       "Reco Track Hit Resolution X",
                       "hit_position_resolution_x", "#N Hits")
            ->Fill(dx * CM);
        if (abs(dx * CM) > abs(largest_dx_deviation))
          largest_dx_deviation = dx * CM;
        double dy = reco.TrackHitPos[0][ih][1] -
                    truth.RecoTrackTrueHitPosition[0][ih][1];
        GetSpecialHist("special__track_hit_resolution_y",
                       "resolution__reco_track__hit_resolution_y",
                       "Reco Track Hit Resolution Y",
                       "hit_position_resolution_y", "#N Hits")
            ->Fill(dy * CM);
        GetHist("resolution__reco_track__hit_resolution_y_comparison_nostack_"
                "with_offset",
                "Reco Track Hit Resolution Y: With Arb Y Offset",
                "hit_position_resolution_y")
            ->Fill(dy * CM);
        avg_offset += dy;

        if (ih == 0) {
          GetHist("resolution__reco_track__hit_resolution_x_first_hit",
                  "Reco Track Hit Resolution X, First Hit Only",
                  "hit_position_resolution_x")
              ->Fill(dx * CM);
          GetHist("resolution__reco_track__hit_resolution_y_first_hit",
                  "Reco Track Hit Resolution Y, First Hit Only",
                  "hit_position_resolution_y")
              ->Fill(dy * CM);
        }
        if (ih == reco.nHits[0] - 1) {
          GetHist("resolution__reco_track__hit_resolution_x_last_hit",
                  "Reco Track Hit Resolution X, Last Hit Only",
                  "hit_position_resolution_x")
              ->Fill(dx * CM);
          GetHist("resolution__reco_track__hit_resolution_y_last_hit",
                  "Reco Track Hit Resolution Y, Last Hit Only",
                  "hit_position_resolution_y")
              ->Fill(dy * CM);
        }
        if (ih >= reco.nHits[0] * 0.33 && ih <= reco.nHits[0] * 0.66) {
          GetHist("resolution__reco_track__hit_resolution_x_middle_third",
                  "Reco Track Hit Resolution X, Middle Third Hits Only",
                  "hit_position_resolution_x")
              ->Fill(dx * CM);
          GetHist("resolution__reco_track__hit_resolution_y_middle_third",
                  "Reco Track Hit Resolution Y, Middle Third Hits Only",
                  "hit_position_resolution_y")
              ->Fill(dy * CM);
        }
      }

      if (largest_dx_deviation > 4)
        DrawSlice(
            TString::Format("entry_%lld", entry_number).Data(),
            "resolution/poor_hit_x_resolution",
            TString::Format("Found poor dx;dx = %.1f cm", largest_dx_deviation)
                .Data(),
            reco, lc, truth, DrawSliceN::many);
      if (largest_dx_deviation > 10)
        DrawSlice(
            TString::Format("entry_%lld", entry_number).Data(),
            "resolution/poor_hit_x_resolution_worse",
            TString::Format("Found poor dx;dx = %.1f cm", largest_dx_deviation)
                .Data(),
            reco, lc, truth, DrawSliceN::many);
      if (largest_dx_deviation > 20)
        DrawSlice(
            TString::Format("entry_%lld", entry_number).Data(),
            "resolution/poor_hit_x_resolution_worst",
            TString::Format("Found poor dx;dx = %.1f cm", largest_dx_deviation)
                .Data(),
            reco, lc, truth, DrawSliceN::many);

      avg_offset /= reco.nHits[0];
      bool draw_slice = false;
      bool draw_slice_large = false;
      for (int ih = 0; ih < reco.nHits[0]; ih++) {
        double dy = (reco.TrackHitPos[0][ih][1] -
                     truth.RecoTrackTrueHitPosition[0][ih][1]) -
                    avg_offset;
        GetHist("resolution__reco_track__hit_resolution_y_offset_removed",
                "Reco Track Hit Resolution Y with Offset Removed",
                "hit_position_resolution_y")
            ->Fill(dy * CM);
        GetHist("resolution__reco_track__hit_resolution_y_comparison_nostack_"
                "without_offset",
                "Reco Track Hit Resolution Y: Without Arb Y Offset",
                "hit_position_resolution_y")
            ->Fill(dy * CM);
        if (dy * CM > 30 && !draw_slice) {
          DrawSlice(
              TString::Format("entry_%lld", entry_number).Data(),
              "resolution/poor_hit_y_resolution",
              TString::Format("Found poor dy;dy = %.1f cm", dy * CM).Data(),
              reco, lc, truth, DrawSliceN::many);
          draw_slice = true;
        }
        if (dy * CM > 60 && !draw_slice_large) {
          DrawSlice(
              TString::Format("entry_%lld", entry_number).Data(),
              "resolution/poor_hit_y_resolution_large_deviation",
              TString::Format("Found poor dy;dy = %.1f cm", dy * CM).Data(),
              reco, lc, truth, DrawSliceN::many);
          draw_slice_large = true;
        }
      }
      double slope_start_reco = DEG * std::atan2(reco.StartDirection[0][1],
                                                 reco.StartDirection[0][2]);
      double slope_start_true =
          DEG * std::atan2(
                    truth.RecoTrackPrimaryParticleTrueMomentumTrackStart[0][1],
                    truth.RecoTrackPrimaryParticleTrueMomentumTrackStart[0][2]);
      double slope_end_reco =
          DEG * std::atan2(reco.EndDirection[0][1], reco.EndDirection[0][2]) -
          180;
      double slope_end_true =
          DEG *
          std::atan2(truth.RecoTrackPrimaryParticleTrueMomentumTrackEnd[0][1],
                     truth.RecoTrackPrimaryParticleTrueMomentumTrackEnd[0][2]);

      GetHist("resolution__reco_track__slope_yz_start",
              "Reco Track Starting YZ Slope Resolution",
              "slope_yz_resolution_track_start")
          ->Fill(slope_start_reco - slope_start_true);
      GetHist("resolution__reco_track__slope_yz_end",
              "Reco Track Ending YZ Slope Resolution",
              "slope_yz_resolution_track_end")
          ->Fill(slope_end_reco - slope_end_true);

      if (std::abs(slope_start_reco - slope_start_true) > 30)
        DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                  "resolution/slope/start_resolution_above_30",
                  TString::Format("resolution = %f",
                                  slope_start_reco - slope_start_true)
                      .Data(),
                  reco, lc, truth, DrawSliceN::handfull);
      if (std::abs(slope_end_reco - slope_end_true) > 30)
        DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                  "resolution/slope/end_resolution_above_30",
                  TString::Format("resolution = %f",
                                  slope_start_reco - slope_start_true)
                      .Data(),
                  reco, lc, truth, DrawSliceN::handfull);
    }
    // Draw the first n for simple comparisons
    if (entry_number < 1000)
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                "first_n_events",
                TString::Format("n tracks = %d", reco.nTracks).Data(), reco, lc,
                truth, DrawSliceN::tons);
    if (entry_number == 128 || entry_number == 292)
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                "selected_events",
                TString::Format("n tracks = %d", reco.nTracks).Data(), reco, lc,
                truth, DrawSliceN::tons);

    if (lc.nLinesX > 0)
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                "track_multiplicity/has_nonzero_nLinesX",
                TString::Format("n tracks = %d", reco.nTracks).Data(), reco, lc,
                truth, DrawSliceN::many);

    // Example of drawing for a reason
    if (reco.nTracks > 2) {
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                "track_multiplicity/high_reco_track_multiplicity",
                TString::Format("n tracks = %d", reco.nTracks).Data(), reco, lc,
                truth, DrawSliceN::few);
    }

    if (reco.nTracks == 1) {
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                "track_multiplicity/single_track/all",
                TString::Format("n tracks = %d", reco.nTracks).Data(), reco, lc,
                truth, DrawSliceN::tons);

      bool ismuon = abs(truth.RecoTrackPrimaryParticlePDG[0]) == 13;
      if (ismuon) {
        DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                  "single_muon/all",
                  TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                  lc, truth, DrawSliceN::tons);
        if (truth.RecoTrackPrimaryParticleTMSFiducialEnd[0])
          DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                    "single_muon/tms_end",
                    TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                    lc, truth, DrawSliceN::tons);

        double muon_starting_angle =
            std::atan2(
                truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[0][1],
                truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[0][2]) *
            DEG;
        if (std::abs(muon_starting_angle) > 20) {
          DrawSlice(
              TString::Format("entry_%lld", entry_number).Data(),
              "track_multiplicity/single_track/high_angle",
              TString::Format("Track angle = %.1f deg", muon_starting_angle)
                  .Data(),
              reco, lc, truth, DrawSliceN::few);
        }
      }
    }
    if (reco.nTracks == 2) {
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                "track_multiplicity/two_track",
                TString::Format("n tracks = %d", reco.nTracks).Data(), reco, lc,
                truth, DrawSliceN::many);
    }
    if (reco.nTracks > 2) {
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                "track_multiplicity/3_or_more_track",
                TString::Format("n tracks = %d", reco.nTracks).Data(), reco, lc,
                truth, DrawSliceN::many);
    }

    // TODO calculate plane number and then check per plane occupancy
    // Also related is total visible energy so compare that to true value

    // std::cout<<entry_number<<": "<<reco.SpillNo<<", "<<reco.SliceNo<<",
    // "<<reco.EventNo<<std::endl;
    if (lc.nHits > 100) {
      // DrawSlice(TString::Format("final_%d", entry_number).Data(), "final",
      // "test job", reco, lc, truth);
    }
    // if (entry_number > 700) exit(1); // TODO delete

#include "Charge.cxx"
#include "EnergyResolution.cxx"
#include "EventRates.cxx"
#include "HitInformation.cxx"
    ValidationContext reco_eff_context{
        truth,         truth_spill,        reco,
        lc,            on_new_spill,       current_tree,
        current_spill, entry_number,       particle_fingerprints_reconstructed};
    FillRecoEff(reco_eff_context);
#include "TimeSlicing.cxx"
#include "Track_Resolution.cxx"
#include "TruthVtx.cxx"

  } // End for loop over entries

  auto time_stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                      time_stop - time_start)
                      .count();

  FinalizeHists();

  // Now save the hists
  outputFile.Write();

  auto entries_visited = entry_number;
  double avg_time = duration / ((double)entries_visited);

  std::cout << "Loop took " << avg_time << "us per event. " << duration
            << "us total for " << entries_visited << " entries" << std::endl;
  return entries_visited;
}

std::string getOutputFilename(const std::string &inputFilename) {
  // Find the position of the last occurrence of '/' in the input filename
  size_t pos = inputFilename.find_last_of('/');
  // Extract the filename without the directory structure
  std::string filename = (pos != std::string::npos)
                             ? inputFilename.substr(pos + 1)
                             : inputFilename;
  // Check for special case where user specified dir/*root
  if (filename.find("*.root") != std::string::npos) {
    // Make output filename out of remaining directory if possible
    if (pos != std::string::npos)
      filename = getOutputFilename(inputFilename.substr(0, pos));
    else
      filename = "validation.root"; // No choice but a default file
  }
  if (filename.size() >= 4 &&
      filename.substr(filename.size() - 4) == ".txt")
    filename = filename.substr(0, filename.size() - 4);
  // Now add .root if needed
  if (filename.find(".root") == std::string::npos)
    filename += ".root";
  if (only_nd_physics) {
    filename = filename.substr(0, filename.size() - 5);
    filename += "_only_nd_physics.root";
    std::cout << "Saving in filename: " << filename << std::endl;
  }
  return filename;
}

std::string getOutputDirname(const std::string &outputFilename) {
  // Find the position of the last occurrence of '.' in the output filename
  size_t pos = outputFilename.find_last_of('.');
  std::string filename = (pos != std::string::npos)
                             ? outputFilename.substr(0, pos)
                             : outputFilename;
  return filename + "_images/";
}

bool HasWildcard(const std::string &path) {
  return path.find('*') != std::string::npos ||
         path.find('?') != std::string::npos ||
         path.find('[') != std::string::npos;
}

std::string Trim(const std::string &text) {
  const auto first = text.find_first_not_of(" \t\r\n");
  if (first == std::string::npos)
    return "";
  const auto last = text.find_last_not_of(" \t\r\n");
  return text.substr(first, last - first + 1);
}

bool HasValidationTrees(const std::filesystem::path &path) {
  TFile file(path.string().c_str(), "READ");
  if (file.IsZombie())
    return false;

  return file.Get("Truth_Info") && file.Get("Truth_Spill") &&
         file.Get("Reco_Tree") &&
         file.Get("Line_Candidates");
}

struct InputFiles {
  std::vector<std::string> files;
  std::vector<std::string> patterns;

  std::string Description() const {
    if (!patterns.empty())
      return patterns.front();
    if (files.size() == 1)
      return files.front();
    return std::to_string(files.size()) + " validation ROOT files";
  }
};

std::vector<std::string> ReadInputFileList(const std::filesystem::path &path) {
  std::ifstream list(path);
  if (!list) {
    throw std::runtime_error("Could not open input file list: " +
                             path.string());
  }

  std::vector<std::string> files;
  std::string line;
  while (std::getline(list, line)) {
    const std::string file = Trim(line);
    if (file.empty() || file[0] == '#')
      continue;
    files.push_back(file);
  }

  if (files.empty()) {
    throw std::runtime_error("Input file list is empty: " + path.string());
  }

  return files;
}

InputFiles ResolveInputFiles(const std::string &inputFilename) {
  InputFiles inputs;

  if (HasWildcard(inputFilename)) {
    inputs.patterns.push_back(inputFilename);
    std::cout << "Using input pattern: " << inputFilename << std::endl;
    return inputs;
  }

  std::error_code ec;
  const std::filesystem::path inputPath(inputFilename);
  if (!std::filesystem::exists(inputPath, ec)) {
    throw std::runtime_error("Input path does not exist: " + inputFilename);
  }

  if (std::filesystem::is_regular_file(inputPath, ec)) {
    if (inputPath.extension() == ".txt") {
      inputs.files = ReadInputFileList(inputPath);
      std::cout << "Using " << inputs.files.size()
                << " input files from file list: " << inputPath << std::endl;
      return inputs;
    }

    inputs.files.push_back(inputPath.string());
    std::cout << "Using single input file: " << inputPath << std::endl;
    return inputs;
  }

  if (!std::filesystem::is_directory(inputPath, ec)) {
    throw std::runtime_error("Input path is neither a file nor a directory: " +
                             inputFilename);
  }

  std::cout << "Scanning input directory recursively: " << inputPath << std::endl;
  for (const auto &entry : std::filesystem::recursive_directory_iterator(
           inputPath, std::filesystem::directory_options::skip_permission_denied,
           ec)) {
    if (ec) {
      std::cerr << "Warning while scanning " << inputPath << ": " << ec.message()
                << std::endl;
      ec.clear();
      continue;
    }

    if (!entry.is_regular_file(ec))
      continue;

    const auto &path = entry.path();
    if (path.extension() != ".root")
      continue;

    if (HasValidationTrees(path))
      inputs.files.push_back(path.string());
  }

  std::sort(inputs.files.begin(), inputs.files.end());

  if (inputs.files.empty()) {
    throw std::runtime_error("No ROOT files with Truth_Info, Truth_Spill, "
                             "Reco_Tree, and Line_Candidates were found under: " +
                             inputFilename);
  }

  std::cout << "Found " << inputs.files.size()
            << " validation input files under " << inputPath << std::endl;
  return inputs;
}

void AddInputsToChain(TChain &chain, const InputFiles &inputs) {
  if (!inputs.patterns.empty()) {
    for (const auto &pattern : inputs.patterns)
      chain.Add(pattern.c_str());
    return;
  }

  for (const auto &file : inputs.files)
    chain.Add(file.c_str());
}

int main(int argc, char *argv[]) {
  // Check if the correct number of arguments is provided
  std::string exeName = getExecutableName(argv[0]);
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " <input file|glob|directory> <num_events (optional)> "
                 "<max num slices to draw (optional)> \\\n"
                 "<output filename (defaults to "
                 "/exp/dune/data/users/$USER/dune-tms/Validation/" +
                     exeName + "/<inputdirname>.root)"
              << std::endl;
    return 1;
  }

  // Extract input filename and number of events from command line arguments
  std::string inputFilename = argv[1];
  InputFiles inputFiles = ResolveInputFiles(inputFilename);
  int numEvents = -1;
  if (argc > 2)
    numEvents = atoi(argv[2]);
  // Can cut down on the number of slices to draw this way
  if (argc > 3)
    DrawSliceN::max_slices = atoi(argv[3]);

  // Load the trees and make the MakeClass objects
  TChain *truth = new TChain("Truth_Info");
  TChain *truth_spill = new TChain("Truth_Spill");
  TChain *reco = new TChain("Reco_Tree");
  TChain *line_candidates = new TChain("Line_Candidates");
  AddInputsToChain(*truth, inputFiles);
  AddInputsToChain(*truth_spill, inputFiles);
  AddInputsToChain(*reco, inputFiles);
  AddInputsToChain(*line_candidates, inputFiles);
  bool missing_ttree = false;
  if (!truth || truth->GetEntriesFast() == 0) {
    std::string message = inputFiles.Description() + " doesn't contain Truth_Info";
    std::cerr << message << std::endl;
    missing_ttree = true;
  }
  if (!truth_spill || truth_spill->GetEntriesFast() == 0) {
    std::string message =
        inputFiles.Description() + " doesn't contain Truth_Spill";
    std::cerr << message << std::endl;
    missing_ttree = true;
  }
  if (!reco || reco->GetEntriesFast() == 0) {
    std::string message = inputFiles.Description() + " doesn't contain Reco_Tree";
    std::cerr << message << std::endl;
    missing_ttree = true;
  }
  if (!line_candidates || line_candidates->GetEntriesFast() == 0) {
    std::string message =
        inputFiles.Description() + " doesn't contain Line_Candidates";
    std::cerr << message << std::endl;
    missing_ttree = true;
  }
  if (missing_ttree) {
    throw std::runtime_error("Missing one or more ttree from file");
  }
  // All these have a large memory footprint, especially Line_Candidates
  // by declaring them static, it moves it from the stack to the heap, which has
  // more memory allocated Otherwise, we get a confusing seg fault before main
  // starts See
  // https://stackoverflow.com/questions/20253267/segmentation-fault-before-main
  static Truth_Info ti(truth);
  static Truth_Spill tsi(truth_spill);
  static Reco_Tree ri(reco);
  std::cout << "About to load Line_Candidates" << std::endl;
  static Line_Candidates li(line_candidates);
  std::cout << "Loaded Line_Candidates" << std::endl;

  std::string directoryPath = "/exp/dune/data/users/" +
                              std::string(getenv("USER")) +
                              "/dune-tms/Validation/" + exeName + "/";

  if (createDirectory(directoryPath)) {
    std::cout << "Directory created: " << directoryPath << std::endl;
  } else {
    std::cerr << "Failed to create directory" << std::endl;
    exit(1);
  }

  SetEnergyFunctionBasedOnInputFile(inputFilename);

  // Create output filename, GIT_BRANCH_NAME + "_" +
  std::string outputFilename;
  if (argc > 4)
    outputFilename = argv[4];
  else
    outputFilename = directoryPath + getOutputFilename(inputFilename);

  save_location = getOutputDirname(outputFilename);

  // Create TFile with the output filename
  TFile outputFile(outputFilename.c_str(), "RECREATE");

  PrimaryLoop(ti, tsi, ri, li, numEvents, outputFile);

  // Close the output file
  outputFile.Close();

  std::cout << "Output file created: " << outputFilename << std::endl;

  return 0;
}
