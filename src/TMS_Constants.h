#ifndef _TMSCONSTANTS_H_SEEN_
#define _TMSCONSTANTS_H_SEEN_

#include <string>

// Constants
namespace TMS_KinConst {
  const double mass_mu = 105.6583755; // Muon mass in MeV/c2
}

// General constants for the TMS
// Lots of these are hard-coded geometry constants that *NEED* to be updated for each production IF detectors move
namespace TMS_Const {

  const int TMS_nThinPlanes = 40;
  const int TMS_nThickPlanes = 60;

  // Dead region (area between LAr and TMS) track length contribution, in g/cm2
  const double Dead_Region_Track_Length = 24.35;

  // Material densities, in g/cm3
  // Check these against the density in the geometry file
  const double LAr_density = 1.3954;
  const double TMS_Steel_density = 7.85;
  const double TMS_Scint_density = 1.05;

  // Z positions of the first hits of the TMS
  const double TMS_Thin_Start = 11362;
  // Where do we transition to the thick region (first layer of scintillator before the change)
  const double TMS_Thick_Start = 13500;
  // Where does the thick region end
  const double TMS_Thick_End = 18294;

  // Approximate starting and end positions of TMS detector in geometry for plotting hits, in {x,y,z}
  // in mm by default!
  const double TMS_Start[] = {-4000, -3500, 11000};
  const double TMS_End[] = {4000, 500, 19000};

  // More exact locations of bars
  const double TMS_Start_Exact[] = {-3520, -3864, TMS_Thin_Start};
  const double TMS_End_Exact[] = {3520, 1159, TMS_Thick_End};

  // Gap for TMS region that is thin iron layer (mm)
  const double TMS_Thin_gap = 55;
  // Gap for TMS region that is thick iron layer (mm)
  const double TMS_Thick_gap = 80;

  // z
  // TMS scintillator width (10 mm)
  const double TMS_Scint_Width = 10;
  //TMS steel width in thin region (15 mm);
  const double TMS_Thin_Steel_Width = 15;
  //TMS steel width in thick region (40 mm);
  const double TMS_Thick_Steel_Width = 40;

  // Offsets to put the TMS in the middle
  const double TMS_Det_Offset[] = { 0., 0., 0. };

  // Start and end of top dead region in x
  const double TMS_Dead_Top[] = {1749, 1769};
  // Start and end of central dead region in x; THIS DOESNT SEEM TO BE PRESENT?
  const double TMS_Dead_Center[] = {0, 0};
  // Start and end of bottom dead region in x
  const double TMS_Dead_Bottom[] = {-1769, -1749};

  // Some distance after the end of the TMS
  const double TMS_End_z = TMS_Thick_End+200;

  // Volume name of TMS related hits
  const std::string TMS_VolumeName = "TMS";
  // Volume name for edep-sim SegmentDetectors
  const std::string TMS_EDepSim_VolumeName = "volTMS";
  // To find in z
  const std::string TMS_ModuleLayerName = "modulelayervol_PV";
  // To find scintillator "box"
  const std::string TMS_ModuleName = "ModuleBoxvol_PV";
  // To find scintillator "box"
  const std::string TMS_ScintLayerName = "scinBoxlvTMS_PV";
  // The top layer name
  const std::string TMS_TopLayerName = "volWorld";

  const std::string LAr_ActiveName = "volTPCActive";
}

#endif
