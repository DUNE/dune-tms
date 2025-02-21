#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TGeoManager.h"
#include "TStopwatch.h"

// EDepSim includes
#include "EDepSim/TG4Event.h"
#include "EDepSim/TG4PrimaryVertex.h"

// TMS includes
// Geometry singleton
#include "TMS_Geom.h"
// Event class
#include "TMS_Event.h"
// Event viewer singleton
#include "TMS_EventViewer.h"
// Reconstructor
#include "TMS_Reco.h"
// TTree writer
#include "TMS_TreeWriter.h"
// General manager
#include "TMS_Manager.h"

int main(int argc, char **argv) {
  if (argc != 8) {
    std::cerr << "Need eight arguments: edep_sim_file startx starty starz endx endy endz" << std::endl;
    return -1;
  }

  std::string filename = std::string(argv[1]);
  TMS_Manager::GetInstance().SetFileName(filename);

  double startx = std::atof(argv[2]);
  double starty = std::atof(argv[3]);
  double startz = std::atof(argv[4]);

  double endx = std::atof(argv[5]);
  double endy = std::atof(argv[6]);
  double endz = std::atof(argv[7]);

  // The input file
  TFile *input = new TFile(filename.c_str(), "open");
  // Get the detector geometry
  TGeoManager *geom = (TGeoManager*)input->Get("EDepSimGeometry");

  // Load up the geometry
  TMS_Geom::GetInstance().SetGeometry(geom);
  TVector3 start(startx, starty, startz);
  TVector3 stop(endx, endy, endz);

  std::cout << "Shooting ray between " << std::endl;
  start.Print();
  std::cout << "and " << std::endl;
  stop.Print();

  // Get the position of materials
  std::vector<std::pair<TGeoMaterial*, TVector3> > Materials = TMS_Geom::GetInstance().GetMaterialsPos(start, stop);
  int counter = 0;
  std::cout << "Found " << Materials.size() << " materials: " << std::endl;
  for (auto it = Materials.begin(); it != Materials.end(); ++it, counter++) {
    std::cout << "***" << std::endl;
    std::cout << "Material " << counter+1 << "/" << Materials.size() << std::endl;
    (*it).first->Print();
    std::cout << "Position: " << std::endl;
    (*it).second.Print();
    std::cout << std::endl;
  }

  // Get the width of materials
  std::vector<std::pair<TGeoMaterial*, double> > MaterialsW = TMS_Geom::GetInstance().GetMaterials(start, stop);
  counter = 0;
  std::cout << "Found " << MaterialsW.size() << " materials: " << std::endl;
  for (auto it = MaterialsW.begin(); it != MaterialsW.end(); ++it, counter++) {
    std::cout << "***" << std::endl;
    std::cout << "Material " << counter+1 << "/" << MaterialsW.size() << std::endl;
    (*it).first->Print();
    std::cout << "Width: " << (*it).second << std::endl;
    std::cout << std::endl;
  }

  std::vector<std::pair<std::string, TVector3> > Nodes = TMS_Geom::GetInstance().GetNodes(start, stop);
  std::cout << "Found " << Nodes.size() << " nodes: " << std::endl;
  counter = 0;
  for (auto it = Nodes.begin(); it != Nodes.end(); ++it, counter++) {
    std::cout << "***" << std::endl;
    std::cout << "Node " << counter+1 << "/" << Nodes.size() << std::endl;
    std::cout << "   " << (*it).first << std::endl;
    std::cout << "   position: " << std::endl;
    (*it).second.Print();
    std::cout << std::endl;
  }

  std::vector<std::pair<int*, TVector3> > Planes = TMS_Geom::GetInstance().GetUniquePlaneBarIdent(start, stop);
  std::cout << "Found " << Planes.size() << " TMS modules: " << std::endl;
  counter = 0;
  for (std::vector<std::pair<int*, TVector3> >::iterator it = Planes.begin(); it != Planes.end(); ++it, counter++) {
    std::cout << "***" << std::endl;
    std::cout << "Module " << counter+1 << "/" << Planes.size() << std::endl;
    std::cout << "Plane enumeration: ";
    for (int i = 0; i < 3; ++i) {
      std::cout << (*it).first[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "At position: " << std::endl;
    (*it).second.Print();
    std::cout << std::endl;
  }

  return 0;
}
