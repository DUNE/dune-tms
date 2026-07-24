#include "TMS_Geom.h"

// Safety valve: abort the survey rather than hang if the geometry tree turns out
// to have far more nodes than expected (e.g. a finely-segmented LAr TPC hierarchy
// under the same top volume). Ten million is comfortably above the node count of
// the TMS-only region (order 10^4-10^5), so this should only ever trip on a
// genuinely pathological tree -- if it does, SurveyGeometry() prints a warning
// and TMS_Geom falls back to the compiled TMS_Constants.h values, same as if no
// survey had run at all.
static const long kMaxNodesToVisit = 10000000;

void TMS_Geom::SurveyNodeRecursive(TMS_GeometryLayout &layout, const std::string &ModuleLayerName,
    const std::string &ModuleName, const std::string &ScintLayerName,
    const std::string &ScintLayerOrthoName, const std::string &ScintLayerParallelName,
    long &nodesVisited) {
  if (nodesVisited >= kMaxNodesToVisit) {
    layout.truncated = true;
    return;
  }
  nodesVisited += 1;

  TGeoNode *node = geom->GetCurrentNode();
  std::string name(node->GetName());

  if (name.find(ModuleLayerName) != std::string::npos) {
    layout.planeNumbers.insert(node->GetNumber());
    layout.nPlaneNodesVisited += 1;
  } else if (name.find(ModuleName) != std::string::npos) {
    layout.moduleNumbers.insert(node->GetNumber());
    layout.nModuleNodesVisited += 1;
  } else if (name.find(ScintLayerName) != std::string::npos ||
             name.find(ScintLayerOrthoName) != std::string::npos ||
             name.find(ScintLayerParallelName) != std::string::npos) {
    layout.nBars += 1;

    // Same local-to-master pattern already used in TMS_Bar::FindModules()
    TGeoBBox *box = dynamic_cast<TGeoBBox*>(node->GetVolume()->GetShape());
    if (box) {
      const double local[3] = {0, 0, 0};
      double master[3];
      geom->GetCurrentMatrix()->LocalToMaster(local, master);
      double x = master[0], y = master[1], z = master[2];
      Scale(x, y, z);
      layout.ExpandToInclude(x, y, z);
    }
    // A scintillator bar node has no daughters worth descending into
    return;
  }

  int nDaughters = node->GetNdaughters();
  for (int i = 0; i < nDaughters && nodesVisited < kMaxNodesToVisit; ++i) {
    geom->CdDown(i);
    SurveyNodeRecursive(layout, ModuleLayerName, ModuleName, ScintLayerName,
        ScintLayerOrthoName, ScintLayerParallelName, nodesVisited);
    geom->CdUp();
  }
}

void TMS_Geom::SurveyGeometry() {
  TMS_Manager &manager = TMS_Manager::GetInstance();
  const std::string ModuleLayerName = manager.Get_GEOMETRY_VOLUME_ModuleLayer();
  const std::string ModuleName = manager.Get_GEOMETRY_VOLUME_Module();
  const std::string ScintLayerName = manager.Get_GEOMETRY_VOLUME_ScintLayer();
  const std::string ScintLayerOrthoName = manager.Get_GEOMETRY_VOLUME_ScintLayerOrtho();
  const std::string ScintLayerParallelName = manager.Get_GEOMETRY_VOLUME_ScintLayerParallel();

  fLayout = TMS_GeometryLayout();
  long nodesVisited = 0;

  geom->CdTop();
  SurveyNodeRecursive(fLayout, ModuleLayerName, ModuleName, ScintLayerName,
      ScintLayerOrthoName, ScintLayerParallelName, nodesVisited);
  // Leave the navigator back where SetGeometry() found it
  geom->CdTop();

  std::cout << "[TMS_Geom] Geometry survey visited " << nodesVisited << " nodes and found "
    << fLayout.nPlaneNodesVisited << " planes, "
    << fLayout.nModuleNodesVisited << " module-name matches, "
    << fLayout.nBars << " scintillator-bar nodes." << std::endl;
  std::cout << "[TMS_Geom]   Of those, " << fLayout.planeNumbers.size() << " distinct raw plane-node numbers "
    << "and " << fLayout.moduleNumbers.size() << " distinct raw module-node numbers appear. Neither of these "
    << "affects TMS_Bar::CheckBar()'s pass/fail logic (that only checks set membership), but if the "
    << "module-name match count isn't close to the distinct module-number count, it usually means the "
    << "'Module' volume nests inside 'ModuleLayer' and gets revisited once per plane -- treat the distinct "
    << "count as the physically meaningful one in that case, not the match count above." << std::endl;
  if (fLayout.nModuleNodesVisited > 0 && !fLayout.moduleNumbers.empty()
      && fLayout.nModuleNodesVisited % (int)fLayout.moduleNumbers.size() == 0
      && fLayout.nModuleNodesVisited / (int)fLayout.moduleNumbers.size() == fLayout.nPlaneNodesVisited) {
    std::cout << "[TMS_Geom]   NOTE: module-name match count (" << fLayout.nModuleNodesVisited
      << ") is exactly plane count (" << fLayout.nPlaneNodesVisited << ") times distinct module count ("
      << fLayout.moduleNumbers.size() << ") -- this is the nesting case above. This geometry most likely "
      << "has " << fLayout.moduleNumbers.size() << " modules per plane, not " << fLayout.nModuleNodesVisited
      << "." << std::endl;
  }
  if (fLayout.valid) {
    std::cout << "[TMS_Geom]   Bar-region bounding box: x [" << fLayout.xMin << ", " << fLayout.xMax
      << "], y [" << fLayout.yMin << ", " << fLayout.yMax
      << "], z [" << fLayout.zMin << ", " << fLayout.zMax << "] mm" << std::endl;
  }
  std::cout << "[TMS_Geom]   For comparison, TMS_Constants.h expects "
    << TMS_Const::nPlanes << " planes and " << TMS_Const::nModules << " modules, "
    << "and a bar region of x [" << TMS_Const::TMS_Start_Bars_Only[0] << ", " << TMS_Const::TMS_End_Bars_Only[0]
    << "], y [" << TMS_Const::TMS_Start_Bars_Only[1] << ", " << TMS_Const::TMS_End_Bars_Only[1]
    << "], z [" << TMS_Const::TMS_Start_Bars_Only[2] << ", " << TMS_Const::TMS_End_Bars_Only[2] << "] mm."
    << std::endl;

  if (fLayout.truncated) {
    std::cerr << "[TMS_Geom] WARNING: Geometry survey hit its " << kMaxNodesToVisit
      << "-node safety cap before finishing. Results above are incomplete -- "
      << "falling back to TMS_Constants.h for anything the survey would otherwise provide."
      << std::endl;
    fLayout.valid = false;
  } else if (fLayout.nBars == 0) {
    std::cerr << "[TMS_Geom] WARNING: Geometry survey found no scintillator-bar nodes. Check "
      << "config/TMS_Default_Config.toml [Geometry.VolumeNames] against this geometry's actual "
      << "volume names." << std::endl;
  }
}
