#include "TMS_Bar.h"

#include <stdexcept>

// Construct a bar from a hit
TMS_Bar::TMS_Bar(TG4HitSegment &edep_seg) {

  BarNumber = -1;
  GlobalBarNumber = -1;
  PlaneNumber = -1;

  xw = -1;
  yw = -1;
  zw = -1;
  AxisReadoutCenter = -99999000;

  // Get the average distance
  double avgx = ( edep_seg.GetStart().X() + edep_seg.GetStop().X() ) / 2;
  double avgy = ( edep_seg.GetStart().Y() + edep_seg.GetStop().Y() ) / 2;
  double avgz = ( edep_seg.GetStart().Z() + edep_seg.GetStop().Z() ) / 2;

  x = avgx;
  y = avgy;
  z = avgz;

  // Find the bar in the geometry
  // Updates the x,y,z,xw,yw,zw
  FindModules(x, y, z);
}

// Construct a bar from x, y, z
TMS_Bar::TMS_Bar(double X, double Y, double Z) : PlaneNumber(-1), BarNumber(-1), GlobalBarNumber(-1),
    x(X), y(Y), z(Z), xw(-1), yw(-1), zw(-1), AxisReadoutCenter(-99999000)
{
  // Find the bar in the geometry
  // Updates the x,y,z,xw,yw,zw
  FindModules(x, y, z);
}

int TMS_Bar::GetBarTypeNumber() const {
  if (GetPlaneNumber() < 0) return -999999999;
  return (int) BarOrient;
}

// Find which bar a given x,y,z position corresponds to
// Maybe this function should be moved to the singleton instead
bool TMS_Bar::FindModules(double xval, double yval, double zval) {

  // Use the ROOT geometry to figure it out if available
  TGeoManager *geom = TMS_Geom::GetInstance().GetGeometry();

  // Find which node this position is equivalent too
  auto node = TMS_Geom::GetInstance().FindNode(xval,yval,zval);
  if (node == NULL) return false;
  std::string NodeName = std::string(node->GetName());

  // Volume-name substrings that identify each part of the tree; config-driven
  // (config/TMS_Default_Config.toml [Geometry.VolumeNames]) rather than compiled,
  // so a GDML renaming doesn't need a recompile.
  TMS_Manager &manager = TMS_Manager::GetInstance();
  const std::string &ModuleLayerName = manager.Get_GEOMETRY_VOLUME_ModuleLayer();
  const std::string &ModuleLayerNameU = manager.Get_GEOMETRY_VOLUME_ModuleLayerU();
  const std::string &ModuleLayerNameV = manager.Get_GEOMETRY_VOLUME_ModuleLayerV();
  const std::string &ModuleLayerNameX = manager.Get_GEOMETRY_VOLUME_ModuleLayerX();
  const std::string &ModuleLayerNameY = manager.Get_GEOMETRY_VOLUME_ModuleLayerY();
  const std::string &ModuleName = manager.Get_GEOMETRY_VOLUME_Module();
  const std::string &ScintLayerName = manager.Get_GEOMETRY_VOLUME_ScintLayer();
  const std::string &ScintLayerOrthoName = manager.Get_GEOMETRY_VOLUME_ScintLayerOrtho();
  const std::string &ScintLayerParallelName = manager.Get_GEOMETRY_VOLUME_ScintLayerParallel();
  const std::string &TopLayerName = manager.Get_GEOMETRY_VOLUME_TopLayer();
  const std::string &DetEnclosureName = manager.Get_GEOMETRY_VOLUME_DetEnclosure();

  // cd up in the geometry to find the right name
  while (NodeName.find(DetEnclosureName) == std::string::npos && NodeName.find(TopLayerName) == std::string::npos) {
    // We've found the plane number
    if (NodeName.find(ModuleLayerName) != std::string::npos) {
      // z-ordered sequential index (TMS_Geom::GetPlaneNumberForCurrentNode(), built from real
      // node z-positions), not the raw ROOT copy number -- raw node numbers reset per parent
      // volume and are NOT globally unique across planes (confirmed via the geometry survey:
      // one real geometry has 82 physically distinct planes but only 52 distinct raw plane-node
      // numbers). PlaneNumber is used as a literal z-ordered coordinate elsewhere in
      // reconstruction (e.g. TMS_Reco.cpp's A* merge candidate search assumes monotonically
      // increasing z order to early-exit its scan), so it must be collision-free and monotonic
      // with z, which only the lookup-table scheme guarantees. CheckBar() validates this
      // against TMS_Geom::GetNumberOfScintillatorPlanes() (also lookup-table-based), not the
      // raw-node-based geometry survey set, to keep both sides consistent.
      PlaneNumber = TMS_Geom::GetInstance().GetPlaneNumberForCurrentNode();
      // There are two rotations of bars, and their names are literally "modulelayervol1" and "modulelayervol2"
      if (NodeName.find(ModuleLayerNameU) != std::string::npos) BarOrient = kUBar;       // +3 degrees tilt from pure y orientation
      else if (NodeName.find(ModuleLayerNameV) != std::string::npos) BarOrient = kVBar;  // -3 degrees rotated/tilted from kYBar orientation
      else if (NodeName.find(ModuleLayerNameX) != std::string::npos) BarOrient = kXBar;  // 90 degrees rotated/tiled from vertical
      else if (NodeName.find(ModuleLayerNameY) != std::string::npos) BarOrient = kYBar;  // not yet implemented! Pure y orientation
      else {
        std::cout<<"Error: Do not understand ModuleLayer volume name '"<<NodeName<<"'. This bar will have type kError."<<std::endl;
        BarOrient = kError;
      }
    }

    // This is the furthest down hit we have: scintillator bar
    else if (NodeName.find(ScintLayerName) != std::string::npos || NodeName.find(ScintLayerOrthoName) != std::string::npos || NodeName.find(ScintLayerParallelName) != std::string::npos) {
      BarNumber = TMS_Geom::GetInstance().GetBarNumberForCurrentNode();

      // Get the width
      TGeoBBox *box = dynamic_cast<TGeoBBox*>(geom->GetCurrentVolume()->GetShape());
      // ROOT saves half of the width of a shape
      xw = 2*box->GetDX();
      yw = 2*box->GetDY();
      zw = 2*box->GetDZ();
      TMS_Geom::GetInstance().Scale(xw, yw, zw);

      BarWidth  = (xw>yw) ? yw : xw; // Set smaller value as width
      BarLength = (xw>yw) ? xw : yw; // Set larger value as length :)

      Double_t position[3];
      const Double_t local[3] = {0, 0, 0};
      geom->GetCurrentMatrix()->LocalToMaster(local, position);
      TMS_Geom::GetInstance().Scale(position[0], position[1], position[2]);
      // Update the hit value to be the bar, not the exact hit position
      x = position[0];
      y = position[1];
      z = position[2];

    }

    else if (NodeName.find(ModuleName) != std::string::npos) {
      GlobalBarNumber = geom->GetCurrentNode()->GetNumber();
    }
   
    geom->CdUp();
    NodeName = std::string(geom->GetCurrentNode()->GetName());
  }

  // If we've reached the world volume we don't have a scintillator hit -> return some mad bad value
  if (BarNumber == -1 || PlaneNumber == -1 || GlobalBarNumber == -1) {
    // Since the bar has already been created as a "error" in the above empty constructor we can just return
    return false;
  }

  // Set the width and orientation of the bar
  /*
  if (PlaneNumber % 2 == 0) BarOrient = kXBar;
  else BarOrient = kYBar;
  */

  // If this is a y-bar, remove the y coordinate
  if (BarOrient == kXBar) {
    AxisReadoutCenter = x;
    x = -99999000;
    // Flip the widths
    double tempyw = yw;
    yw = xw;
    xw = tempyw;
  } else if (BarOrient == kUBar || BarOrient == kVBar || BarOrient == kYBar) {
    AxisReadoutCenter = y;
    y = -99999000;
    // Don't need to flip the widths because they're already correct (yw = large, xw = 4cm)
  } else {
    AxisReadoutCenter = -99999000;
    x = -99999000;
    y = -99999000;
    z = -99999000;
  }

  // Reset the geom navigator node level in case it's used again
  TMS_Geom::GetInstance().FindNode(xval,yval,zval);

  CheckBar();

  return true;
}

// Find which bar a given x,y,z position corresponds to
// Maybe this function should be moved to the singleton instead
int TMS_Bar::FindBar(double x_pos, double y_pos, double z_pos) {

  // Use the ROOT geometry to figure it out if available
  TGeoManager *geom = TMS_Geom::GetInstance().GetGeometry();

  TVector3 vec_tmp = TVector3(x_pos,y_pos,z_pos);
  if (TMS_Geom::GetInstance().IsInsideTMS(vec_tmp)) {
    // Since the bar has already been created as a "error" in the above empty constructor we can just return
    std::cout << "position given to TMS_Bar::FindBar is outside detector!" << std::endl;
    return -1;
  }

  TMS_Manager &manager = TMS_Manager::GetInstance();
  const std::string &ModuleLayerName = manager.Get_GEOMETRY_VOLUME_ModuleLayer();
  const std::string &DetEnclosureName = manager.Get_GEOMETRY_VOLUME_DetEnclosure();

  // Find which node this position is equivalent too
  TMS_Geom::GetInstance().FindNode(x_pos,y_pos,z_pos);
  TGeoNavigator *nav = geom->GetCurrentNavigator();
  std::string NodeName = std::string(nav->GetCurrentNode()->GetName());
  // cd up in the geometry to find the right name
  while (NodeName.find(ModuleLayerName) == std::string::npos &&
      NodeName.find(DetEnclosureName) == std::string::npos) {
    nav->CdUp();
    NodeName = std::string(nav->GetCurrentNode()->GetName());
  }

  // If we've reached the world volume we don't have a scintillator hit -> return some mad bad value
  if (NodeName.find(DetEnclosureName) != std::string::npos) {
    // Since the bar has already been created as a "error" in the above empty constructor we can just return
    std::cout << "Bar position not found in TMS_Bar::FindBar!" << std::endl;
    return -1;
  }

  int Number = TMS_Geom::GetInstance().GetPlaneNumberForCurrentNode();

  return Number;
}

std::string TMS_Bar::BarType_ToString(BarType bar) const {
  if (bar == kXBar) return std::string("X-bar");
  else if (bar == kYBar) return std::string("Y-bar");
  else if (bar == kUBar) return std::string("U-bar");
  else if (bar == kVBar) return std::string("V-bar");
  else if (bar == kXBar) return std::string("X-bar");
  return std::string("ERROR");
}

void TMS_Bar::Print() const {
  std::cout << "Printing TMS bar: " << std::endl;
  std::cout << "(x, y, z) = " << "(" << x << ", " << y << ", " << z << ")" << std::endl;
  std::cout << "(xw, yw, zw) = " << "(" << xw << ", " << yw << ", " << zw << ")" << std::endl;
  std::cout << "BarOrient: " << BarType_ToString(BarOrient) << " (" << BarOrient << ")" << std::endl;
  std::cout << "PlaneNumber: " << PlaneNumber << std::endl;
  std::cout << "BarNumber: " << BarNumber << std::endl;
  std::cout << "GlobalBarNumber: " << GlobalBarNumber << std::endl;
}

double TMS_Bar::FindYbar(double yval) {
  const double ymin = TMS_Const::TMS_End[2];

  // The total range of y
  // Splits into how many 40mm slices
  // Change to 36 mm slices?
  int bin = (yval-ymin)/40;
  // Return the center of the bin
  double val = ymin+bin*40+20;

  return val;
}

bool TMS_Bar::CheckBar() {

  // Sanity check the global bar number
  if (BarNumber < 0) {
    std::cerr << "Bar number was not found in the geometry bar lookup." << std::endl;
    std::cerr << "Has the geometry been updated without updating the geometry constants in TMS_Constants.h?" << std::endl;
    std::cout << "Bar number: " << BarNumber << std::endl;
    throw std::runtime_error("TMS_Bar::CheckBar(): bar number not found in the geometry bar lookup");
  }

  // This check validates GlobalBarNumber against what the geometry survey actually found
  // in the loaded GDML (TMS_Geom::SurveyGeometry()) rather than hardcoded expected counts,
  // so it can't drift out of sync with a geometry change. If no survey is available (e.g.
  // SetGeometry() wasn't called), TMS_Geom falls back to the TMS_Constants.h range instead
  // of failing every bar.
  TMS_Geom &tmsGeom = TMS_Geom::GetInstance();
  if (tmsGeom.HasGeometrySurvey()) {
    if (!tmsGeom.IsSurveyedModuleNumber(GlobalBarNumber)) {
      std::cerr << "Global bar number " << GlobalBarNumber << " was not among the "
        << tmsGeom.GetNDistinctModuleNumbers() << " distinct module numbers found by the geometry survey." << std::endl;
      throw std::runtime_error("TMS_Bar::CheckBar(): global bar number not among the surveyed module numbers");
    }
  } else {
    if (GlobalBarNumber >= TMS_Const::nModules || GlobalBarNumber < 0) {
      std::cerr << "Global bar number does not agree with expectation of between 0 to " << TMS_Const::nModules << std::endl;
      std::cerr << "Has the geometry been updated without updating the geometry constants in TMS_Constants.h?" << std::endl;
      throw std::runtime_error("TMS_Bar::CheckBar(): global bar number out of the TMS_Constants.h expected range");
    }
  }

  // PlaneNumber comes from TMS_Geom's z-ordered lookup table (GetPlaneNumberForCurrentNode()),
  // not the geometry survey's raw-node-number set -- those two schemes are not interchangeable
  // (see the comment in FindModules() above). The lookup table assigns a contiguous 0..count-1
  // index by construction, whether or not a bbox survey has run, so a simple range check against
  // GetNumberOfScintillatorPlanes() (itself lookup-table-based, with a config fallback) applies
  // uniformly here -- no HasGeometrySurvey() split needed, unlike GlobalBarNumber above.
  int expected_nplanes = TMS_Geom::GetInstance().GetNumberOfScintillatorPlanes();
  if (PlaneNumber >= expected_nplanes || PlaneNumber < 0) {
    std::cerr << "Plane number does not agree with expectation of between 0 to " << expected_nplanes << std::endl;
    std::cerr << "Has the geometry been updated without updating the geometry constants in TMS_Constants.h?" << std::endl;
    throw std::runtime_error("TMS_Bar::CheckBar(): plane number out of the expected range");
  }

  return true;
}
