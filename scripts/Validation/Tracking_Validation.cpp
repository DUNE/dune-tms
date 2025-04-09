#include <chrono> // for high res clock
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <math.h> /* atan2 */
#include <tuple>
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
#include <TVector3.h>

// dune specific
#include <DUNEStyle.h>

#include "Line_Candidates.h"
#include "Reco_Tree.h"
#include "Truth_Info.h"

#define IS_WITHIN(x, center, tolerance)                                        \
  (std::abs((x) - (center)) <= (tolerance))
#define print(var) std::cout << #var << " = " << var << std::endl;

std::string save_location = "";

const double CM = 0.1; // cm per mm
const double DEG = 360 / TAU;
const double GEV = 1e-3; // GeV per MEV

#define length_to_energy_clarence(l) (l * 1.75 + 82) * 1e-3
#define default_length_to_energy(l) (l * 1.75) * 1e-3
// #define length_to_energy(l) default_length_to_energy(l) // No fit
#define GEOM_V3 // for old geom
#ifdef GEOM_V3
#define length_to_energy(l) (l * 1.75 * 0.951 + 76.8) * 1e-3 // Old geom
#else
#define length_to_energy(l) (l * 1.75 * 0.9428 + 18.73) * 1e-3 // New geom
// #define length_to_energy(l) (l*1.75*0.9460 + 91.80)*1e-3 // New geom, no
// kalman
#endif
#define lar_length_to_energy(l) (l) * 1e-3

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
#else
const double TMS_LAST_PLANE_Z = 18400;
#endif
const double TMS_Y_TOP =
    160.0; // Top of scintillator (I'm 99% sure starting from where U/V overlap)
const double TMS_Y_BOTTOM =
    -2850.0; // Bottom of scinillator (also measured from U/V overlap)
const double TMS_Y_MIDDLE = 0.5 * (TMS_Y_TOP + TMS_Y_BOTTOM);
const double TMS_X_EXTENT = 3300.0; // Assumes things are symmetrical right now
const double Y_BUFFER = 400;        // mm, 40 cm
const double X_BUFFER = 30;         // mm, 3cm bar thickness?
const double Z_BUFFER = 100;        // mm, two thin planes
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

Long64_t PrimaryLoop(Truth_Info &truth, Reco_Tree &reco, Line_Candidates &lc,
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

  int last_spill_seen = -1;

  auto time_start = std::chrono::high_resolution_clock::now();
  std::map<int, double> manually_calculated_shell_hadronic_energy_per_vertex;

  std::map<int, int> particle_indices_reconstructed;

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

    // Check if we're on a new spill
    // Some things should only be filled per spill
    bool on_new_spill = false;
    if (last_spill_seen != reco.SpillNo) {
      on_new_spill = true;
      last_spill_seen = reco.SpillNo;
      // Empty on each new spill
      particle_indices_reconstructed.clear();
    }

    // Fill some basic "raw" variables from the reco tree
    GetHist("basic__raw__EventNo", "EventNo", "EventNo")->Fill(reco.EventNo);
    GetHist("basic__raw__SliceNo", "SliceNo", "SliceNo")->Fill(reco.SliceNo);
    GetHist("basic__raw__SpillNo", "SpillNo", "SpillNo")->Fill(reco.SpillNo);

    GetHist("basic__raw__ntracks", "N Reco Tracks", "ntracks")
        ->Fill(reco.nTracks);
    for (int it = 0; it < reco.nTracks; it++) {
      GetHist("basic__raw__nHits", "nHits", "n0-120")->Fill(reco.nHits[it]);
      GetHist("basic__raw__nKalmanNodes", "nKalmanNodes", "n0-120")
          ->Fill(reco.nKalmanNodes[it]);
      for (int ih = 0; ih < reco.nHits[it]; ih++) {
        GetHist("basic__raw__TrackHitPos_X", "TrackHitPos X", "X")
            ->Fill(reco.TrackHitPos[it][ih][0] * CM);
        GetHist("basic__raw__TrackHitPos_Y", "TrackHitPos Y", "Y")
            ->Fill(reco.TrackHitPos[it][ih][1] * CM);
        GetHist("basic__raw__TrackHitPos_Z", "TrackHitPos Z", "Z")
            ->Fill(reco.TrackHitPos[it][ih][2] * CM);
      }
      for (int ih = 0; ih < reco.nKalmanNodes[it]; ih++) {
        GetHist("basic__raw__KalmanPos_X", "KalmanPos X", "X")
            ->Fill(reco.KalmanPos[it][ih][0] * CM);
        GetHist("basic__raw__KalmanPos_Y", "KalmanPos Y", "Y")
            ->Fill(reco.KalmanPos[it][ih][1] * CM);
        GetHist("basic__raw__KalmanPos_Z", "KalmanPos Z", "Z")
            ->Fill(reco.KalmanPos[it][ih][2] * CM);
        GetHist("basic__raw__KalmanTruePos_X", "KalmanTruePos X", "X")
            ->Fill(reco.KalmanTruePos[it][ih][0] * CM);
        GetHist("basic__raw__KalmanTruePos_Y", "KalmanTruePos Y", "Y")
            ->Fill(reco.KalmanTruePos[it][ih][1] * CM);
        GetHist("basic__raw__KalmanTruePos_Z", "KalmanTruePos Z", "Z")
            ->Fill(reco.KalmanTruePos[it][ih][2] * CM);
      }
      GetHist("basic__raw__StartDirection_X", "StartDirection X", "dx")
          ->Fill(reco.StartDirection[it][0]);
      GetHist("basic__raw__StartDirection_Y", "StartDirection Y", "dy")
          ->Fill(reco.StartDirection[it][1]);
      GetHist("basic__raw__StartDirection_Z", "StartDirection Z", "dz")
          ->Fill(reco.StartDirection[it][2]);
      GetHist("basic__raw__EndDirection_X", "EndDirection X", "dx")
          ->Fill(reco.EndDirection[it][0]);
      GetHist("basic__raw__EndDirection_Y", "EndDirection Y", "dy")
          ->Fill(reco.EndDirection[it][1]);
      GetHist("basic__raw__EndDirection_Z", "EndDirection Z", "dz")
          ->Fill(reco.EndDirection[it][2]);

      REGISTER_AXIS(DirectionSanityCheck,
                    std::make_tuple("Direction Mag", 51, 0, 10));
      double start_direction_mag =
          std::sqrt(reco.StartDirection[it][0] * reco.StartDirection[it][0] +
                    reco.StartDirection[it][1] * reco.StartDirection[it][1] +
                    reco.StartDirection[it][2] * reco.StartDirection[it][2]);
      double end_direction_mag =
          std::sqrt(reco.EndDirection[it][0] * reco.EndDirection[it][0] +
                    reco.EndDirection[it][1] * reco.EndDirection[it][1] +
                    reco.EndDirection[it][2] * reco.EndDirection[it][2]);
      //const double epsilon = 1e-6;
      // if (std::abs(start_direction_mag - 1) > epsilon)
      //   std::cout<<"Warning: Found start_direction_mag outside of 1.
      //   "<<start_direction_mag<<std::endl;
      // if (std::abs(end_direction_mag - 1) > epsilon)
      //   std::cout<<"Warning: Found end_direction_mag outside of 1.
      //   "<<end_direction_mag<<std::endl;
      GetHist("basic__sanity__StartDirectionMag", "StartDirection Mag",
              "DirectionSanityCheck")
          ->Fill(start_direction_mag);
      GetHist("basic__sanity__EndDirectionMag", "EndDirection Mag",
              "DirectionSanityCheck")
          ->Fill(end_direction_mag);
      if (std::abs(reco.StartDirection[it][1]) > 1)
        std::cout << "big y dir: " << reco.StartDirection[it][0] << ","
                  << reco.StartDirection[it][1] << ","
                  << reco.StartDirection[it][2] << "\t" << start_direction_mag
                  << std::endl;

      const double big_epsilon = 1e-1;
      bool has_y_start_direction_zero =
          (std::abs(reco.StartDirection[it][1]) < big_epsilon);
      if (has_y_start_direction_zero)
        DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                  "y_start_direction_zero",
                  TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                  lc, truth, DrawSliceN::many);
      bool has_y_end_direction_zero =
          (std::abs(reco.EndDirection[it][1]) < big_epsilon);
      if (has_y_end_direction_zero)
        DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                  "y_end_direction_zero",
                  TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                  lc, truth, DrawSliceN::many);

      GetHist("basic__raw__StartPos_X", "StartPos X", "X")
          ->Fill(reco.StartPos[it][0] * CM);
      GetHist("basic__raw__StartPos_Y", "StartPos Y", "Y")
          ->Fill(reco.StartPos[it][1] * CM);
      GetHist("basic__raw__StartPos_Z", "StartPos Z", "Z")
          ->Fill(reco.StartPos[it][2] * CM);
      GetHist("basic__raw__EndPos_X", "EndPos X", "X")
          ->Fill(reco.EndPos[it][0] * CM);
      GetHist("basic__raw__EndPos_Y", "EndPos Y", "Y")
          ->Fill(reco.EndPos[it][1] * CM);
      GetHist("basic__raw__EndPos_Z", "EndPos Z", "Z")
          ->Fill(reco.EndPos[it][2] * CM);

      GetHist("basic__raw__Charge", "Charge", "charge")->Fill(reco.Charge[it]);
      GetHist("basic__raw__Length", "Areal Density, AKA Length",
              "areal_density")
          ->Fill(reco.Length[it]);
      GetHist("basic__raw__Momentum", "Momentum", "momentum")
          ->Fill(reco.Momentum[it]);
      GetHist("basic__raw__EnergyRange", "EnergyRange", "energy_range")
          ->Fill(reco.EnergyRange[it]);
      GetHist("basic__raw__EnergyDeposit", "EnergyDeposit", "energy_deposit")
          ->Fill(reco.EnergyDeposit[it]);
    }
    // TODO finish adding these vars
    // TODO add "fixes" to direction, etc and see if that fixes things

    /*
 // Used to be Direction, now is StartDirection, check for both options
 depending on the file
 //if (HasBranch("Direction")) fChain->SetBranchAddress("Direction", Direction,
 &b_Direction);
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
 fChain->SetBranchAddress("TrackHitEnergies", TrackHitEnergies,
 &b_RecoTrackHitEnergies); */

    // Adjust to kalman if needed
    if (!has_kalman) {
      for (int itrack = 0; itrack < reco.nTracks; itrack++) {
        for (int ihit = 0; ihit < reco.nHits[itrack]; ihit++) {
          for (int i = 0; i < 3; i++) {
            reco.KalmanPos[itrack][ihit][i] = reco.TrackHitPos[itrack][ihit][i];
            reco.KalmanTruePos[itrack][ihit][i] =
                truth.RecoTrackTrueHitPosition[itrack][ihit][i];
          }
        }
        reco.nKalmanNodes[itrack] = reco.nHits[itrack];
      }
    }
    // Fix n kalman nodes which is always zero.
    for (int itrack = 0; itrack < reco.nTracks; itrack++) {
      int n_nonzero = 0;
      for (int ihit = 0; ihit < reco.nHits[itrack]; ihit++) {
        // Since kalman pos need to be above z = 11000, we can use this to find
        // zero points safely
        if (reco.KalmanPos[itrack][ihit][2] > 1000)
          n_nonzero++;
      }
      reco.nKalmanNodes[itrack] = n_nonzero;

      /*if (reco.nHits[itrack] > 0) {
        const int LOOKBACK_WINDOW = 10;
        int n_nodes = n_nonzero > LOOKBACK_WINDOW ? LOOKBACK_WINDOW : n_nonzero;
        n_nodes -= 1;
        reco.StartDirection[itrack][0] = reco.KalmanPos[itrack][0][0] -
      reco.KalmanPos[itrack][n_nodes][0]; reco.StartDirection[itrack][1] =
      reco.KalmanPos[itrack][0][1] - reco.KalmanPos[itrack][n_nodes][1];
        reco.StartDirection[itrack][2] = reco.KalmanPos[itrack][0][2] -
      reco.KalmanPos[itrack][n_nodes][2];

        int last_index = reco.nKalmanNodes[itrack] - 1;
        reco.EndDirection[itrack][0] = reco.KalmanPos[itrack][last_index][0] -
      reco.KalmanPos[itrack][last_index - n_nodes][0];
        reco.EndDirection[itrack][1] = reco.KalmanPos[itrack][last_index][1] -
      reco.KalmanPos[itrack][last_index - n_nodes][1];
        reco.EndDirection[itrack][2] = reco.KalmanPos[itrack][last_index][2] -
      reco.KalmanPos[itrack][last_index - n_nodes][2];
      }*/
    }

    for (int it = 0; it < reco.nTracks; it++) {
      GetHist("basic__fixed__nKalmanNodes", "nKalmanNodes", "n0-120")
          ->Fill(reco.nKalmanNodes[it]);
      for (int ih = 0; ih < reco.nKalmanNodes[it]; ih++) {
        GetHist("basic__fixed__KalmanPos_X", "KalmanPos X", "X")
            ->Fill(reco.KalmanPos[it][ih][0] * CM);
        GetHist("basic__fixed__KalmanPos_Y", "KalmanPos Y", "Y")
            ->Fill(reco.KalmanPos[it][ih][1] * CM);
        GetHist("basic__fixed__KalmanPos_Z", "KalmanPos Z", "Z")
            ->Fill(reco.KalmanPos[it][ih][2] * CM);
        GetHist("basic__fixed__KalmanTruePos_X", "KalmanTruePos X", "X")
            ->Fill(reco.KalmanTruePos[it][ih][0] * CM);
        GetHist("basic__fixed__KalmanTruePos_Y", "KalmanTruePos Y", "Y")
            ->Fill(reco.KalmanTruePos[it][ih][1] * CM);
        GetHist("basic__fixed__KalmanTruePos_Z", "KalmanTruePos Z", "Z")
            ->Fill(reco.KalmanTruePos[it][ih][2] * CM);
      }
      GetHist("basic__fixed__StartDirection_X", "StartDirection X", "dx")
          ->Fill(reco.StartDirection[it][0]);
      GetHist("basic__fixed__StartDirection_Y", "StartDirection Y", "dy")
          ->Fill(reco.StartDirection[it][1]);
      GetHist("basic__fixed__StartDirection_Z", "StartDirection Z", "dz")
          ->Fill(reco.StartDirection[it][2]);
      GetHist("basic__fixed__EndDirection_X", "EndDirection X", "dx")
          ->Fill(reco.EndDirection[it][0]);
      GetHist("basic__fixed__EndDirection_Y", "EndDirection Y", "dy")
          ->Fill(reco.EndDirection[it][1]);
      GetHist("basic__fixed__EndDirection_Z", "EndDirection Z", "dz")
          ->Fill(reco.EndDirection[it][2]);

      // Check for reco hits outside the TMS scint volume
      for (int ih = 0; ih < reco.nKalmanNodes[it]; ih++) {
        int hit_flag_x = -1;
        int hit_flag_y = -1;
        int hit_flag_z = -1;
        if (std::abs(reco.KalmanPos[it][ih][2]) < 0.01)
          hit_flag_z = 0;
        if (reco.KalmanPos[it][ih][2] < -990)
          hit_flag_z = 1;
        if (reco.KalmanPos[it][ih][2] < 11362)
          hit_flag_z = 2;
        if (reco.KalmanPos[it][ih][2] > 18313)
          hit_flag_z = 3;
        // if (hit_flag_z == 3) std::cout<<"reco.KalmanPos[it][ih][2]:
        // "<<reco.KalmanPos[it][ih][2]<<std::endl;
        if (reco.KalmanPos[it][ih][1] < -2950)
          hit_flag_y = 2;
        if (reco.KalmanPos[it][ih][1] > 240)
          hit_flag_y = 3;
        if (reco.KalmanPos[it][ih][0] < -3350)
          hit_flag_x = 2;
        if (reco.KalmanPos[it][ih][0] > 3350)
          hit_flag_x = 3;
        // Only check x and y == 0 if z is strange
        if (hit_flag_z != -1) {
          if (std::abs(reco.KalmanPos[it][ih][0]) < 0.01)
            hit_flag_x = 0;
          if (std::abs(reco.KalmanPos[it][ih][1]) < 0.01)
            hit_flag_y = 0;
        }
        if (hit_flag_z != -1)
          GetHist("basic__sanity__Unusual_Reco_Hit_Locations_Z",
                  "Unusual Reco Hit Locations Z", "unusual_hit_locations")
              ->Fill(hit_flag_z);
        if (hit_flag_x != -1)
          GetHist("basic__sanity__Unusual_Reco_Hit_Locations_X",
                  "Unusual Reco Hit Locations X", "unusual_hit_locations")
              ->Fill(hit_flag_x);
        if (hit_flag_y != -1)
          GetHist("basic__sanity__Unusual_Reco_Hit_Locations_Y",
                  "Unusual Reco Hit Locations Y", "unusual_hit_locations")
              ->Fill(hit_flag_y);
      }
    }

    REGISTER_AXIS(HitEnergy, std::make_tuple("Hit Energy (MeV)", 51, 0, 100));
    REGISTER_AXIS(
        TrueVisibleEnergy,
        std::make_tuple("Particle True Visible Energy (MeV)", 51, 0, 1000));
    REGISTER_AXIS(HitHadronicEnergy,
                  std::make_tuple("Hit Hadronic Energy (MeV)", 51, 0, 100));
    REGISTER_AXIS(HitEnergy_zoom,
                  std::make_tuple("Hit Energy (MeV)", 20, 0, 1));
    REGISTER_AXIS(Hitdedx, std::make_tuple("Hit dEdx (MeV / cm)", 51, 0, 100));
    REGISTER_AXIS(TotalEnergy,
                  std::make_tuple("Total Hit Energy (GeV)", 51, 0, 100));
    REGISTER_AXIS(
        TotalHadEnergy,
        std::make_tuple("Total Hit Hadronic Energy (GeV)", 51, 0, 100));
    REGISTER_AXIS(TotalRatio, std::make_tuple("Ratio E_{Had}/E", 51, 0, 1.0));
    REGISTER_AXIS(
        TotalEnergyPerVertex,
        std::make_tuple("Total Hit Energy Per Vertex (GeV)", 51, 0, 20));
    REGISTER_AXIS(TotalHadronicEnergyPerVertex,
                  std::make_tuple("Total Hit Hadronic Energy Per Vertex (GeV)",
                                  51, 0, 20));
    REGISTER_AXIS(
        ShellHadronicEnergyPerVertex,
        std::make_tuple("Hit Hadronic Energy Per Vertex In Shell (GeV)", 51, 0,
                        20));
    REGISTER_AXIS(
        ShellHadronicEnergyPerVertexZoom,
        std::make_tuple("Hit Hadronic Energy Per Vertex In Shell (MeV)", 20, 0,
                        200));
    if (on_new_spill) {
      GetHist("basic__truth__nTrueParticles", "nTrueParticles", "n0-500")
          ->Fill(truth.nTrueParticles);
      GetHist("basic__truth__nTruePrimaryParticles", "nTruePrimaryParticles",
              "n0-500")
          ->Fill(truth.nTruePrimaryParticles);
      GetHist("basic__truth__nTrueForgottenParticles",
              "nTrueForgottenParticles", "n0-120")
          ->Fill(truth.nTrueForgottenParticles);

      std::map<int, bool> vertex_id_to_fiducial;
      for (int ip = 0; ip < truth.nTrueParticles; ip++) {
        int pdg = truth.PDG[ip];
        if (std::abs(pdg) != 13) {
          int vid = truth.VertexID[ip];
          TVector3 position(truth.BirthPosition[ip][0] * CM,
                            truth.BirthPosition[ip][1] * CM,
                            truth.BirthPosition[ip][2] * CM);
          bool result = IsInLAr(position, LArFiducial);
          vertex_id_to_fiducial[vid] = result;
        }
      }

      float total_energy = 0;
      float total_hadronic_energy = 0;
      std::map<int, float> total_hadronic_energy_per_vertex;
      std::map<int, float> total_energy_per_vertex;
      std::map<int, float> total_hadronic_energy_per_vertex_in_shell;

      if (truth.HasBranch("TrueNonTMSNHits")) {
        for (int ih = 0; ih < truth.TrueNonTMSNHits; ih++) {
          // Only care about above 0.5 MeV
          bool fiducial = false;
          int vid = truth.TrueNonTMSHitVertexID[ih];
          if (vertex_id_to_fiducial.find(vid) != vertex_id_to_fiducial.end())
            fiducial = vertex_id_to_fiducial[vid];
          if (truth.TrueNonTMSHitEnergy[ih] >= 0.5 && fiducial) {
            GetHist("basic__truth__nontms_hits__hit_energy",
                    "TrueNonTMSHitEnergy", "HitEnergy")
                ->Fill(truth.TrueNonTMSHitEnergy[ih]);
            GetHist("basic__truth__nontms_hits__hit_energy_zoom",
                    "TrueNonTMSHitEnergy", "HitEnergy_zoom")
                ->Fill(truth.TrueNonTMSHitEnergy[ih]);
            GetHist("basic__truth__nontms_hits__hit_hadronic_energy",
                    "TrueNonTMSHitHadronicEnergy", "HitHadronicEnergy")
                ->Fill(truth.TrueNonTMSHitHadronicEnergy[ih]);
            GetHist("basic__truth__nontms_hits__hit_dedx", "TrueNonTMSHitdEdx",
                    "Hitdedx")
                ->Fill(truth.TrueNonTMSHitdEdx[ih] /
                       CM); // in MeV / mm natively
            GetHist("basic__truth__nontms_hits__hit_x", "TrueNonTMSHitPosX",
                    "X")
                ->Fill(truth.TrueNonTMSHitPos[ih][0] * CM);
            GetHist("basic__truth__nontms_hits__hit_y", "TrueNonTMSHitPosY",
                    "Y_full")
                ->Fill(truth.TrueNonTMSHitPos[ih][1] * CM);
            GetHist("basic__truth__nontms_hits__hit_z", "TrueNonTMSHitPosZ",
                    "Z_full")
                ->Fill(truth.TrueNonTMSHitPos[ih][2] * CM);

            total_energy += truth.TrueNonTMSHitEnergy[ih];
            total_hadronic_energy += truth.TrueNonTMSHitHadronicEnergy[ih];
            total_energy_per_vertex[truth.TrueNonTMSHitVertexID[ih]] +=
                truth.TrueNonTMSHitEnergy[ih];
            total_hadronic_energy_per_vertex[truth.TrueNonTMSHitVertexID[ih]] +=
                truth.TrueNonTMSHitHadronicEnergy[ih];
            // Get position in CM
            TVector3 position(truth.TrueNonTMSHitPos[ih][0] * CM,
                              truth.TrueNonTMSHitPos[ih][1] * CM,
                              truth.TrueNonTMSHitPos[ih][2] * CM);
            // Want to check if in shell, so check if inside lar but not inside
            // the lar using the shell buffer size of 30 cm
            bool hit_in_shell = IsInLAr(position, LArFull) &&
                                !IsInLAr(position, LArUsingShellBuffer);
            // std::cout<<"IsInLAr(position, LArFull): "<<IsInLAr(position,
            // LArFull)<<", IsInLAr(position, LArUsingShellBuffer):
            // "<<IsInLAr(position, LArUsingShellBuffer)<<", IsInLAr(position,
            // LArFull) && !IsInLAr(position, LArUsingShellBuffer):
            // "<<(IsInLAr(position, LArFull) && !IsInLAr(position,
            // LArUsingShellBuffer))<<std::endl;
            if (hit_in_shell)
              total_hadronic_energy_per_vertex_in_shell
                  [truth.TrueNonTMSHitVertexID[ih]] +=
                  truth.TrueNonTMSHitHadronicEnergy[ih];

            if (truth.TrueNonTMSHitdEdx[ih] / CM >= 5)
              GetHist("basic__truth__nontms_hits__high_dedx_hit_energy",
                      "TrueNonTMSHitEnergy", "HitEnergy")
                  ->Fill(truth.TrueNonTMSHitEnergy[ih]);
            if (truth.TrueNonTMSHitEnergy[ih] >= 0.5)
              GetHist("basic__truth__nontms_hits__high_e_hit_dedx",
                      "TrueNonTMSHitdEdx", "Hitdedx")
                  ->Fill(truth.TrueNonTMSHitdEdx[ih] /
                         CM); // in MeV / mm natively

            double weight = truth.TrueNonTMSHitEnergy[ih];
            GetHist("basic__truth__nontms_hits__energy_weighted_hit_x",
                    "TrueNonTMSHitPosX", "X")
                ->Fill(truth.TrueNonTMSHitPos[ih][0] * CM, weight);
            GetHist("basic__truth__nontms_hits__energy_weighted_hit_y",
                    "TrueNonTMSHitPosY", "Y_full")
                ->Fill(truth.TrueNonTMSHitPos[ih][1] * CM, weight);
            GetHist("basic__truth__nontms_hits__energy_weighted_hit_z",
                    "TrueNonTMSHitPosZ", "Z_full")
                ->Fill(truth.TrueNonTMSHitPos[ih][2] * CM, weight);

            if (truth.TrueNonTMSHitEnergy[ih] >= 3) {
              GetHist("basic__truth__nontms_hits__above_3MeV_hit_x",
                      "TrueNonTMSHitPosX", "X")
                  ->Fill(truth.TrueNonTMSHitPos[ih][0] * CM);
              GetHist("basic__truth__nontms_hits__above_3MeV_hit_y",
                      "TrueNonTMSHitPosY", "Y_full")
                  ->Fill(truth.TrueNonTMSHitPos[ih][1] * CM);
              GetHist("basic__truth__nontms_hits__above_3MeV_hit_z",
                      "TrueNonTMSHitPosZ", "Z_full")
                  ->Fill(truth.TrueNonTMSHitPos[ih][2] * CM);
            }

            if (truth.TrueNonTMSHitdEdx[ih] / CM >= 5) {
              GetHist("basic__truth__nontms_hits__dedx_above_5MeVpcm_hit_x",
                      "TrueNonTMSHitPosX", "X")
                  ->Fill(truth.TrueNonTMSHitPos[ih][0] * CM);
              GetHist("basic__truth__nontms_hits__dedx_above_5MeVpcm_hit_y",
                      "TrueNonTMSHitPosY", "Y_full")
                  ->Fill(truth.TrueNonTMSHitPos[ih][1] * CM);
              GetHist("basic__truth__nontms_hits__dedx_above_5MeVpcm_hit_z",
                      "TrueNonTMSHitPosZ", "Z_full")
                  ->Fill(truth.TrueNonTMSHitPos[ih][2] * CM);
            }
          }
        }
      }

      GetHist("basic__truth__nontms_hits__total_energy",
              "Total LAr Hit Energy per Spill", "TotalEnergy")
          ->Fill(total_energy * GEV);
      GetHist("basic__truth__nontms_hits__total_hadronic_energy",
              "Total LAr Hit Hadronic Energy per Spill", "TotalHadEnergy")
          ->Fill(total_hadronic_energy * GEV);
      GetHist("basic__truth__nontms_hits__total_ratio",
              "Ratio of Hadronic to Total LAr Hit Energy per Spill",
              "TotalRatio")
          ->Fill(total_hadronic_energy / total_energy);
      for (auto it : total_energy_per_vertex) {
        GetHist("basic__truth__nontms_hits__per_vertex_total_energy",
                "Total LAr Hit Energy per Vertex", "TotalEnergyPerVertex")
            ->Fill(it.second * GEV);
      }
      for (auto it : total_hadronic_energy_per_vertex) {
        GetHist("basic__truth__nontms_hits__per_vertex_total_hadronic_energy",
                "Total LAr Hit Hadronic Energy per Vertex",
                "TotalHadronicEnergyPerVertex")
            ->Fill(it.second * GEV);
      }
      for (auto it : total_hadronic_energy_per_vertex_in_shell) {
        GetHist("basic__truth__nontms_hits__per_vertex_shell_hadronic_energy",
                "Shell LAr Hit Hadronic Energy per Vertex",
                "ShellHadronicEnergyPerVertex")
            ->Fill(it.second * GEV);
        GetHist(
            "basic__truth__nontms_hits__per_vertex_shell_hadronic_energy_zoom",
            "Shell LAr Hit Hadronic Energy per Vertex Zoomed",
            "ShellHadronicEnergyPerVertexZoom")
            ->Fill(it.second);
        manually_calculated_shell_hadronic_energy_per_vertex[it.first] =
            it.second;
        GetHist(
            "basic__truth__nontms_hits__per_vertex_shell_hadronic_energy_cut",
            "Shell LAr Hit Hadronic Energy per Vertex Cut Passes", "yesno")
            ->Fill((it.second < 30) ? 0 : 1);
      }

      for (int ip = 0; ip < truth.nTrueParticles; ip++) {
        GetHist("basic__truth__PDG", "PDG", "pdg")
            ->Fill(PDGtoIndex(truth.PDG[ip]));
        if (truth.IsPrimary[ip])
          GetHist("basic__truth__PDG_Primary", "PDG Primary Particles", "pdg")
              ->Fill(PDGtoIndex(truth.PDG[ip]));
        if (!truth.IsPrimary[ip])
          GetHist("basic__truth__PDG_Secondary", "PDG Secondary Particles",
                  "pdg")
              ->Fill(PDGtoIndex(truth.PDG[ip]));
        GetHist("basic__truth__TrueVisibleEnergy", "TrueVisibleEnergy",
                "TrueVisibleEnergy")
            ->Fill(truth.TrueVisibleEnergy[ip]);
        GetHist("basic__truth__RecoTrackPrimaryParticleTruePositionStart_X",
                "RecoTrackPrimaryParticleTruePositionStart X", "X")
            ->Fill(truth.RecoTrackPrimaryParticleTruePositionStart[ip][0] * CM);
        GetHist("basic__truth__RecoTrackPrimaryParticleTruePositionStart_Y",
                "RecoTrackPrimaryParticleTruePositionStart Y", "Y_full")
            ->Fill(truth.RecoTrackPrimaryParticleTruePositionStart[ip][1] * CM);
        GetHist("basic__truth__RecoTrackPrimaryParticleTruePositionStart_Z",
                "RecoTrackPrimaryParticleTruePositionStart Z", "Z_full")
            ->Fill(truth.RecoTrackPrimaryParticleTruePositionStart[ip][2] * CM);
        GetHist("basic__truth__RecoTrackPrimaryParticleTruePositionEnd_X",
                "RecoTrackPrimaryParticleTruePositionEnd X", "X")
            ->Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[ip][0] * CM);
        GetHist("basic__truth__RecoTrackPrimaryParticleTruePositionEnd_Y",
                "RecoTrackPrimaryParticleTruePositionEnd Y", "Y_full")
            ->Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[ip][1] * CM);
        GetHist("basic__truth__RecoTrackPrimaryParticleTruePositionEnd_Z",
                "RecoTrackPrimaryParticleTruePositionEnd Z", "Z_full")
            ->Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[ip][2] * CM);

        GetHist("basic__truth__RecoTrackPrimaryParticleTruePositionEnd_Z",
                "RecoTrackPrimaryParticleTruePositionEnd Z", "yesno")
            ->Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[ip][2] * CM);

        GetHist("basic__truth__RecoTrackPrimaryParticleLArFiducialStart",
                "RecoTrackPrimaryParticleLArFiducialStart", "yesno")
            ->Fill(truth.RecoTrackPrimaryParticleLArFiducialStart[ip] ? 0 : 1);
        GetHist("basic__truth__RecoTrackPrimaryParticleTMSFiducialEnd",
                "RecoTrackPrimaryParticleTMSFiducialEnd", "yesno")
            ->Fill(truth.RecoTrackPrimaryParticleTMSFiducialEnd[ip] ? 0 : 1);
      }
    }

    {
      int vid = truth.VertexIdOfMostEnergyInEvent;
      if (manually_calculated_shell_hadronic_energy_per_vertex.find(vid) !=
          manually_calculated_shell_hadronic_energy_per_vertex.end()) {
        double shell = truth.LArOuterShellEnergyFromVertex;
        double shell_manual =
            manually_calculated_shell_hadronic_energy_per_vertex[vid];
        GetHist("basic__truth__nontms_hits__comparison",
                "Comparison of Shell Energies", "ShellHadronicEnergyPerVertex",
                "ShellHadronicEnergyPerVertex")
            ->Fill(shell, shell_manual);
        GetHist("basic__truth__nontms_hits__comparison_zoom",
                "Comparison of Shell Energies",
                "ShellHadronicEnergyPerVertexZoom",
                "ShellHadronicEnergyPerVertexZoom")
            ->Fill(shell, shell_manual);
      }
    }

    // Truth matching information
    REGISTER_AXIS(
        completeness,
        std::make_tuple("Track Completeness (primary on track / primary)", 20,
                        0, 1.01));
    REGISTER_AXIS(
        cleanliness,
        std::make_tuple("Track Cleanliness (primary on track / total on track)",
                        20, 0, 1.01));
    REGISTER_AXIS(nhits_in_track,
                  std::make_tuple("N Hits in Track", 200, 0, 200));
    REGISTER_AXIS(energy_in_track,
                  std::make_tuple("Energy in Track (MeV)", 100, 0, 1000));
    REGISTER_AXIS(n_hits_per_plane,
                  std::make_tuple("N Hits per Plane", 5, 0.5, 5.5));
    for (int it = 0; it < reco.nTracks; it++) {
      GetHist("basic__reco_track__primary_pdg",
              "Reco Track Primary Particle PDG", "pdg")
          ->Fill(PDGtoIndex(truth.RecoTrackPrimaryParticlePDG[it]));
      GetHist("basic__reco_track__secondary_pdg",
              "Reco Track Secondary Particle PDG", "pdg")
          ->Fill(PDGtoIndex(truth.RecoTrackSecondaryParticlePDG[it]));
      int particle_index = truth.RecoTrackPrimaryParticleIndex[it];
      if (particle_index < 0 || particle_index >= truth.nTrueParticles) {
        std::cout << "Found a particle index outside the range: "
                  << particle_index
                  << ", nTrueParticles=" << truth.nTrueParticles << std::endl;
      } else {
        // Completeness = true primary on track / true primary, lower -> more
        // missed hits Cleanliness = true primary on track / true all on track,
        // lower -> more contamination
        double completeness_energy =
            truth.RecoTrackPrimaryParticleTrueVisibleEnergy[it] /
            truth.TrueVisibleEnergy[particle_index];
        double cleanliness_energy =
            truth.RecoTrackPrimaryParticleTrueVisibleEnergy[it] /
            truth.RecoTrackTrueVisibleEnergy[it];
        double completeness_nhits =
            truth.RecoTrackPrimaryParticleTrueNHits[it] /
            (double)
                truth.TrueNHits[particle_index]; // Also can do TrueNHitsInSlice
        double cleanliness_nhits = truth.RecoTrackPrimaryParticleTrueNHits[it] /
                                   (double)truth.RecoTrackNHits[it];
        GetHist("basic__reco_track__completeness_energy",
                "Reco Track Completeness, Visible Energy", "completeness")
            ->Fill(completeness_energy);
        GetHist("basic__reco_track__cleanliness_energy",
                "Reco Track Cleanliness, Visible Energy", "cleanliness")
            ->Fill(cleanliness_energy);
        GetHist("basic__reco_track__completeness_nhits",
                "Reco Track Completeness, N Hits", "completeness")
            ->Fill(completeness_nhits);
        GetHist("basic__reco_track__cleanliness_nhits",
                "Reco Track Cleanliness, N Hits", "cleanliness")
            ->Fill(cleanliness_nhits);

        GetHist("basic__reco_track__debugging__completeness_nostack_1_n_hits_"
                "in_track",
                "N Hits in Track: Reco Track", "nhits_in_track")
            ->Fill(truth.RecoTrackPrimaryParticleTrueNHits[it]);
        GetHist("basic__reco_track__debugging__completeness_nostack_2_n_hits_"
                "in_slice",
                "N Hits in Track: Slice", "nhits_in_track")
            ->Fill(truth.TrueNHitsInSlice[particle_index]);
        GetHist("basic__reco_track__debugging__completeness_nostack_3_n_hits_"
                "in_truth",
                "N Hits in Track: Truth", "nhits_in_track")
            ->Fill(truth.TrueNHits[particle_index]);

        GetHist("basic__reco_track__debugging__completeness_energy_nostack_1_"
                "energy_in_track",
                "Energy in Track: Reco Track", "energy_in_track")
            ->Fill(truth.RecoTrackPrimaryParticleTrueVisibleEnergy[it]);
        GetHist("basic__reco_track__debugging__completeness_energy_nostack_2_"
                "energy_in_truth",
                "Energy in Track: Truth", "energy_in_track")
            ->Fill(truth.TrueVisibleEnergy[particle_index]);

        bool ismuon = std::abs(truth.RecoTrackPrimaryParticlePDG[it]) == 13;
        if (ismuon) {
          GetHist("basic__reco_track__debugging__muon_only_completeness_"
                  "nostack_1_n_hits_in_track",
                  "N Hits in Track: Reco Track", "nhits_in_track")
              ->Fill(truth.RecoTrackPrimaryParticleTrueNHits[it]);
          GetHist("basic__reco_track__debugging__muon_only_completeness_"
                  "nostack_2_n_hits_in_slice",
                  "N Hits in Track: Slice", "nhits_in_track")
              ->Fill(truth.TrueNHitsInSlice[particle_index]);
          GetHist("basic__reco_track__debugging__muon_only_completeness_"
                  "nostack_3_n_hits_in_truth",
                  "N Hits in Track: Truth", "nhits_in_track")
              ->Fill(truth.TrueNHits[particle_index]);

          GetHist("basic__reco_track__debugging__muon_only_completeness_energy_"
                  "nostack_1_energy_in_track",
                  "Energy in Track: Reco Track", "energy_in_track")
              ->Fill(truth.RecoTrackPrimaryParticleTrueVisibleEnergy[it]);
          GetHist("basic__reco_track__debugging__muon_only_completeness_energy_"
                  "nostack_2_energy_in_truth",
                  "Energy in Track: Truth", "energy_in_track")
              ->Fill(truth.TrueVisibleEnergy[particle_index]);
        }

        if (ismuon && reco.nTracks == 1) {
          GetHist("basic__reco_track__debugging__single_muon_completeness_"
                  "nostack_1_n_hits_in_track",
                  "N Hits in Track: Reco Track", "nhits_in_track")
              ->Fill(truth.RecoTrackPrimaryParticleTrueNHits[it]);
          GetHist("basic__reco_track__debugging__single_muon_completeness_"
                  "nostack_2_n_hits_in_slice",
                  "N Hits in Track: Slice", "nhits_in_track")
              ->Fill(truth.TrueNHitsInSlice[particle_index]);
          GetHist("basic__reco_track__debugging__single_muon_completeness_"
                  "nostack_3_n_hits_in_truth",
                  "N Hits in Track: Truth", "nhits_in_track")
              ->Fill(truth.TrueNHits[particle_index]);
          {
            std::map<double, int> count;
            for (int ih = 0; ih < lc.nHits; ih++)
              count[lc.RecoHitPos[ih][2]] += 1;
            for (auto &plane : count) {
              GetHist("basic__reco_track__debugging__n_hits_in_plane",
                      "Hits per plane: All Reco Hits", "n_hits_per_plane")
                  ->Fill(plane.second);
              GetHist("basic__reco_track__debugging__n_hits_in_plane_both_"
                      "nostack_1_all_hits",
                      "Hits per plane: All Reco Hits", "n_hits_per_plane")
                  ->Fill(plane.second);
            }
          }
          {
            std::map<double, int> count;
            for (int ih = 0; ih < reco.nHits[0]; ih++)
              count[reco.TrackHitPos[0][ih][2]] += 1;
            for (auto &plane : count) {
              GetHist("basic__reco_track__debugging__n_hits_in_plane_track",
                      "Hits per plane: In Track", "n_hits_per_plane")
                  ->Fill(plane.second);
              GetHist("basic__reco_track__debugging__n_hits_in_plane_both_"
                      "nostack_2_track",
                      "Hits per plane: In Tracks", "n_hits_per_plane")
                  ->Fill(plane.second);
            }
          }
        }

        bool has_more_reco_energy =
            truth.RecoTrackPrimaryParticleTrueVisibleEnergy[it] >
            truth.TrueVisibleEnergy[particle_index];
        bool has_more_reco_nhits = truth.RecoTrackPrimaryParticleTrueNHits[it] >
                                   truth.TrueNHits[particle_index];
        bool has_zero_true_nhits = truth.TrueNHits[particle_index] == 0;

        if (has_zero_true_nhits)
          DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                    "completeness/zero_true_nhits",
                    TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                    lc, truth, DrawSliceN::handfull);
        if (has_more_reco_nhits)
          DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                    "completeness/has_more_reco_nhits",
                    TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                    lc, truth, DrawSliceN::handfull);
        if (has_more_reco_energy)
          DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                    "completeness/has_more_reco_energy",
                    TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                    lc, truth, DrawSliceN::handfull);

        if (completeness_energy < 0.5)
          DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                    "completeness/under_50_percent",
                    TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                    lc, truth, DrawSliceN::handfull);
        if (completeness_energy >= 0.67 && completeness_energy <= 0.73)
          DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                    "completeness/around_70_percent",
                    TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                    lc, truth, DrawSliceN::handfull);
        if (completeness_energy >= 0.98)
          DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                    "completeness/above_98_percent",
                    TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                    lc, truth, DrawSliceN::handfull);
        if (cleanliness_energy < 0.5)
          DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                    "cleanliness/under_50_percent",
                    TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                    lc, truth, DrawSliceN::handfull);
        if (cleanliness_energy >= 0.67 && cleanliness_energy <= 0.73)
          DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                    "cleanliness/around_70_percent",
                    TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                    lc, truth, DrawSliceN::handfull);
        if (cleanliness_energy >= 0.98)
          DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                    "cleanliness/above_98_percent",
                    TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                    lc, truth, DrawSliceN::handfull);
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
      int index = -99999;
      double max_muon_starting_ke = -999999;
      for (int it = 0; it < reco.nTracks; it++) {
        if (std::abs(truth.RecoTrackPrimaryParticlePDG[it]) == 13) {
          double muon_starting_ke =
              truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][3] *
              1e-3;
          if (muon_starting_ke > max_muon_starting_ke) {
            max_muon_starting_ke = muon_starting_ke;
            index = truth.RecoTrackPrimaryParticleIndex[it];
          }
        }
      }
      if (index >= 0) {
        TVector3 birth_pos(truth.BirthPosition[index][0],
                           truth.BirthPosition[index][1],
                           truth.BirthPosition[index][2]);
        passes_ndlar_fiducial_cut = LArFiducialCut(birth_pos);
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

    REGISTER_AXIS(energy_resolution,
                  std::make_tuple("Energy Resolution (Reco - True) / True", 21,
                                  -0.4, 0.4));
    for (int it = 0; it < reco.nTracks; it++) {
      bool ismuon = abs(truth.RecoTrackPrimaryParticlePDG[it]) == 13;
      if (ismuon) {
        double muon_starting_ke =
            truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][3] * 1e-3;
        double estimated_reco_ke = (82 + 1.75 * reco.Length[it]) * 1e-3;
        double resolution = estimated_reco_ke - muon_starting_ke;
        double fractional_resolution = resolution / muon_starting_ke;

        GetHist("resolution__muon_ke", "Muon KE Resolution",
                "energy_resolution")
            ->Fill(fractional_resolution);
        if (muon_starting_ke <= 0.5)
          GetHist("resolution__muon_ke_binned_area_norm_nostack_1_le05",
                  "Muon KE Resolution: True P_{#mu}/GeV<0.5",
                  "energy_resolution")
              ->Fill(fractional_resolution);
        if (0.5 < muon_starting_ke && muon_starting_ke <= 1.5)
          GetHist("resolution__muon_ke_binned_area_norm_nostack_2_gt05_le15",
                  "Muon KE Resolution: 0.5<True P_{#mu}/GeV<1.5",
                  "energy_resolution")
              ->Fill(fractional_resolution);
        if (1.5 < muon_starting_ke && muon_starting_ke <= 2.5)
          GetHist("resolution__muon_ke_binned_area_norm_nostack_3_gt15_le25",
                  "Muon KE Resolution: 1.5<True P_{#mu}/GeV<2.5",
                  "energy_resolution")
              ->Fill(fractional_resolution);
        if (2.5 < muon_starting_ke && muon_starting_ke <= 3.5)
          GetHist("resolution__muon_ke_binned_area_norm_nostack_4_gt25_le35",
                  "Muon KE Resolution: 2.5<True P_{#mu}/GeV<3.5",
                  "energy_resolution")
              ->Fill(fractional_resolution);
        if (3.5 < muon_starting_ke && muon_starting_ke <= 4.5)
          GetHist("resolution__muon_ke_binned_area_norm_nostack_5_gt35_le45",
                  "Muon KE Resolutionn: 3.5<True P_{#mu}/GeV<4.5",
                  "energy_resolution")
              ->Fill(fractional_resolution);
        if (muon_starting_ke > 4.5)
          GetHist("resolution__muon_ke_binned_area_norm_nostack_6_gt45",
                  "Muon KE Resolution: 4.5<True P_{#mu}/GeV",
                  "energy_resolution")
              ->Fill(fractional_resolution);
      }
    }

    REGISTER_AXIS(
        hit_position_resolution_x,
        std::make_tuple("Hit Position Resolution X (Reco - True) (cm)", 21, -4,
                        4));
    REGISTER_AXIS(
        hit_position_resolution_y,
        std::make_tuple("Hit Position Resolution Y (Reco - True) (cm)", 21, -40,
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
      double avg_offset = 0;
      for (int ih = 0; ih < reco.nHits[0]; ih++) {
        double dx = reco.TrackHitPos[0][ih][0] -
                    truth.RecoTrackTrueHitPosition[0][ih][0];
        GetHist("resolution__reco_track__hit_resolution_x",
                "Reco Track Hit Resolution X", "hit_position_resolution_x")
            ->Fill(dx * CM);
        double dy = reco.TrackHitPos[0][ih][1] -
                    truth.RecoTrackTrueHitPosition[0][ih][1];
        GetHist("resolution__reco_track__hit_resolution_y",
                "Reco Track Hit Resolution Y", "hit_position_resolution_y")
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
          DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                    "poor_hit_y_resolution",
                    TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                    lc, truth, DrawSliceN::many);
          draw_slice = true;
        }
        if (dy * CM > 60 && !draw_slice_large) {
          DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                    "poor_hit_y_resolution_large_deviation",
                    TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                    lc, truth, DrawSliceN::many);
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
                  "slope/start_resolution_above_30",
                  TString::Format("resolution = %f",
                                  slope_start_reco - slope_start_true)
                      .Data(),
                  reco, lc, truth, DrawSliceN::handfull);
      if (std::abs(slope_end_reco - slope_end_true) > 30)
        DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                  "slope/end_resolution_above_30",
                  TString::Format("resolution = %f",
                                  slope_start_reco - slope_start_true)
                      .Data(),
                  reco, lc, truth, DrawSliceN::handfull);
    }

    // Want to compare
    if (entry_number == 341 || entry_number == 118)
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(), "special",
                TString::Format("n tracks = %d", reco.nTracks).Data(), reco, lc,
                truth, DrawSliceN::many);
    // Draw the first n for simple comparisons
    if (entry_number < 1000)
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                "first_n_events",
                TString::Format("n tracks = %d", reco.nTracks).Data(), reco, lc,
                truth, DrawSliceN::tons);

    // Example of drawing for a reason
    if (reco.nTracks > 2) {
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                "high_reco_track_multiplicity",
                TString::Format("n tracks = %d", reco.nTracks).Data(), reco, lc,
                truth, DrawSliceN::few);
    }

    /*if (reco.nTracks == 1) {
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
    "lar_starting_single_track", TString::Format("n tracks = %d",
    reco.nTracks).Data(), reco, lc, truth, DrawSliceN::tons);
    }*/ // TODO finish nd physics sample only tracks here

    if (reco.nTracks == 1) {
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                "single_track",
                TString::Format("n tracks = %d", reco.nTracks).Data(), reco, lc,
                truth, DrawSliceN::tons);

      bool ismuon = abs(truth.RecoTrackPrimaryParticlePDG[0]) == 13;
      if (ismuon) {
        DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                  "single_muon",
                  TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                  lc, truth, DrawSliceN::tons);
        if (truth.RecoTrackPrimaryParticleTMSFiducialEnd[0])
          DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                    "single_muon_tms_end",
                    TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                    lc, truth, DrawSliceN::tons);

        double muon_starting_angle =
            std::atan2(
                truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[0][1],
                truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[0][2]) *
            DEG;
        if (std::abs(muon_starting_angle) > 20) {
          DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                    "single_track_high_angle",
                    TString::Format("n tracks = %d", reco.nTracks).Data(), reco,
                    lc, truth, DrawSliceN::few);
        }
      }
    }
    if (reco.nTracks == 2) {
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(), "two_track",
                TString::Format("n tracks = %d", reco.nTracks).Data(), reco, lc,
                truth, DrawSliceN::many);
    }
    if (reco.nTracks > 2) {
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                "3_or_more_track",
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

#include "EnergyResolution.cxx"
#include "EventRates.cxx"
#include "HitInformation.cxx"
#include "Reco_Eff.cxx"
#include "TimeSlicing.cxx"
#include "Track_Resolution.cxx"
#include "TruthVtx.cxx"

  } // End for loop over entries

  auto time_stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                      time_stop - time_start)
                      .count();

  NormalizeHists();

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

int main(int argc, char *argv[]) {
  // Check if the correct number of arguments is provided
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " <input_filename> <num_events (optional)> <num slices to "
                 "draw (optional)>"
              << std::endl;
    return 1;
  }

  // Extract input filename and number of events from command line arguments
  std::string inputFilename = argv[1];
  int numEvents = -1;
  if (argc > 2)
    numEvents = atoi(argv[2]);
  // Can cut down on the number of slices to draw this way
  if (argc > 3)
    DrawSliceN::max_slices = atoi(argv[3]);

  // Load the tree and make the Truth_Info object
  TFile TF(inputFilename.c_str());
  TTree *truth = (TTree *)TF.Get("Truth_Info");
  TTree *reco = (TTree *)TF.Get("Reco_Tree");
  TTree *line_candidates = (TTree *)TF.Get("Line_Candidates");
  bool missing_ttree = false;
  if (!truth) {
    std::string message = inputFilename + " doesn't contain Truth_Info";
    std::cerr << message << std::endl;
    missing_ttree = true;
  }
  if (!reco) {
    std::string message = inputFilename + " doesn't contain Reco_Tree";
    std::cerr << message << std::endl;
    missing_ttree = true;
  }
  if (!line_candidates) {
    std::string message = inputFilename + " doesn't contain Line_Candidates";
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
  static Reco_Tree ri(reco);
  std::cout << "About to load Line_Candidates" << std::endl;
  static Line_Candidates li(line_candidates);
  std::cout << "Loaded Line_Candidates" << std::endl;

  std::string exeName = getExecutableName(argv[0]);
  std::string directoryPath = "/exp/dune/data/users/" +
                              std::string(getenv("USER")) +
                              "/dune-tms/Validation/" + exeName + "/";

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
