#ifndef _TMSCONSTANTS_H_SEEN_
#define _TMSCONSTANTS_H_SEEN_

#include <string>
#include "TDatabasePDG.h"


// Hidden from us inside edep-sim EDepSimRooTrackerKinematicsGenerator.hh
// Number of maximum particles in an edep-sim event
#define __EDEP_SIM_MAX_PART__ 4000

#define __TMS_BAD_NUMBER__ -999.99

// Constants
namespace TMS_KinConst {
  static double GetMass(int pdg) {
    TDatabasePDG *database = TDatabasePDG::Instance();
    double out = 0;
    // The only particles not in the database are nuclei, and zero is a fine guess for their mass
    auto particle = database->GetParticle(pdg);
    if (particle) {
      out = particle->Mass() * 1000; // Convert from GeV to MeV
    }
    return out;
  }

  const double mass_mu = GetMass(13); // Used in Kalman filter

}

// General constants for the TMS
// Lots of these are hard-coded geometry constants that *NEED* to be updated for each production IF detectors move
namespace TMS_Const {

  // Number of planes, check against geometry
  const int TMS_nThinPlanes = 51;
  const int TMS_nThickPlanes = 34;
  const int TMS_nDoublePlanes = 9;

  // Dead region (area between LAr and TMS) track length contribution, in g/cm2
  const double Dead_Region_Track_Length = 24.35;

  // Material densities, in g/cm3
  // Check these against the density in the geometry file
  const double LAr_density = 1.3954;
  const double TMS_Steel_density = 7.85;
  const double TMS_Scint_density = 1.05;

  // Z positions of the first hits of the TMS
  const double TMS_Thin_Start = 11134;
  // Where do we transition to the thick region (first layer of scintillator before the change)
  const double TMS_Thick_Start = 14435;
  // Where does the thick region end
  const double TMS_Double_Start = 17495;
  // Where does the double region end
  const double TMS_Double_End = 18535;


  // Approximate starting and end positions of TMS detector in geometry for plotting hits, in {x,y,z}
  // in mm by default!
  const double TMS_Start[] = {-4000, -3500, 11000};
  const double TMS_End[] = {4000, 500, 19000};

  // From scanning all hits x,y,z position in the LAr active volume
  const double LAr_Start_Exact[] = {-3478.48, -2166.71, 4179.24};
  const double LAr_End_Exact[] = {3478.48, 829.282, 9135.88};

  // More exact locations of bars
  // This seems to contain the steel as well
  const double TMS_Start_Exact[] = {-3730, -3702.23, TMS_Thin_Start};
  const double TMS_End_Exact[] = {3730, 997.77, TMS_Double_End};

  // Plot TrueHitX,Y,Z and zoom in to see where the last hits are
  const double TMS_Start_Bars_Only[] = {-3488, -3076.23, TMS_Thin_Start};  // this is the bottom with horizontal bars, stereo: -3002.23
  const double TMS_End_Bars_Only[] = {3488, 371.77, TMS_Double_End};       // this is the top with horizontal bars, stereo: 297.77

  // Gap for TMS region that is thin iron layer (mm)
  const double TMS_Thin_gap = 65;
  // Gap for TMS region that is thick iron layer (mm)
  const double TMS_Thick_gap = 90;
  // Gap for TMS region that is double thick iron layer (mm)
  const double TMS_Double_gap = 130;

  // z
  // TMS scintillator width (10 mm)
  const double TMS_Scint_Width = 17;
  // TMS aluminium enclosure width (1mm)  // This is not implemented currently
  //const double TMS_Enclosure_Width = 1;
  // TMS steel width in thin region (15 mm);
  const double TMS_Thin_Steel_Width = 15;
  // TMS steel width in thick region (40 mm);
  const double TMS_Thick_Steel_Width = 40;
  // TMS steel width in double thick region (80 mm);
  const double TMS_Double_Steel_Width = 80;

  // Offsets to put the TMS in the middle
  const double TMS_Det_Offset[] = { 0., 0., 0. };

  // Start and end of top dead region in x
  const double TMS_Dead_Top[] = {1749, 1769};
  // Start and end of central dead region in x; THIS DOESNT SEEM TO BE PRESENT?
  const double TMS_Dead_Center[] = {0, 0};
  // Start and end of bottom dead region in x
  const double TMS_Dead_Bottom[] = {-1769, -1749};

  // Some distance after the end of the TMS
  const double TMS_End_z = TMS_Double_End+200;

  // Volume name of TMS related hits
  const std::string TMS_VolumeName = "TMS";
  // Volume name for edep-sim SegmentDetectors
  const std::string TMS_EDepSim_VolumeName = "volTMS";
  // To find in z
  const std::string TMS_ModuleLayerName = "modulelayervol";
  const std::string TMS_ModuleLayerName1 = "modulelayervol1"; // u orientation
  const std::string TMS_ModuleLayerName2 = "modulelayervol2"; // v orientation
  const std::string TMS_ModuleLayerName3 = "modulelayervol3"; // x orientation
  const std::string TMS_ModuleLayerName4 = "modulelayervol4"; // y orientation
  // To find scintillator "box"
  const std::string TMS_ModuleName = "ModuleBoxvol";
  // To find scintillator "box"
  const std::string TMS_ScintLayerName = "scinBoxlvTMS";
  const std::string TMS_ScintLayerOrthoName = "scinBoxlv_orthoTMS";
  const std::string TMS_ScintLayerParallelName = "scinBoxlv_parallelTMS";
  // The top layer name
  const std::string TMS_TopLayerName = "volWorld";
  // The detector enclosure
  const std::string TMS_DetEnclosure = "volDetEnclosure";
  // The rock volume
  const std::string TMS_Rock = "rockBox_lv";
  // LAr active region
  const std::string LAr_ActiveName = "volTPCActive";

  const double TMS_Small_Num = 1.E-5;

  const int nModulesPerSubModule = 32;
  const int nModules = 12;
  const int nPlanes = 94;

}

#endif
