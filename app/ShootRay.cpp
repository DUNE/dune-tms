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
  std::vector<std::pair<TGeoMaterial*, double> > Materials = TMS_Geom::GetInstance().GetMaterials(start, stop);

  std::cout << "Shooting ray between " << std::endl;
  start.Print();
  std::cout << "and " << std::endl;
  stop.Print();
  int counter = 0;
  std::cout << "Found " << Materials.size() << " materials: " << std::endl;
  for (auto it = Materials.begin(); it != Materials.end(); ++it, counter++) {
    std::cout << "Material " << counter << "/" << Materials.size() << std::endl;
    (*it).first->Print();
    std::cout << "Thickness: " << (*it).second << std::endl;
  }

  std::vector<std::pair<std::string, const double*> > Nodes = TMS_Geom::GetInstance().GetNodes(start, stop);
  std::cout << "Found " << Nodes.size() << " nodes: " << std::endl;
  counter = 0;
  for (auto it = Nodes.begin(); it != Nodes.end(); ++it, counter++) {
    std::cout << "Node " << counter << "/" << Nodes.size() << std::endl;
    std::cout << "   " << (*it).first << std::endl;
    std::cout << "   position: ";
    for (int i = 0; i < 3; ++i) std::cout << (*it).second[i] << " ";
    std::cout << std::endl;
  }

  std::vector<std::pair<int*, const double*> > Planes = TMS_Geom::GetInstance().GetUniquePlaneBarIdent(start, stop);
  std::cout << "Found " << Planes.size() << " nodes: " << std::endl;
  counter = 0;
  for (auto it = Planes.begin(); it != Planes.end(); ++it, counter++) {
    std::cout << "Node " << counter << "/" << Planes.size() << std::endl;
    std::cout << "position: " << std::endl;
    for (int i = 0; i < 3; ++i) {
      std::cout << (*it).first[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < 3; ++i) {
      std::cout << (*it).second[i] << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}
