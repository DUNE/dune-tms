// python code was very slow...
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
  double TrackLength = TMS_Geom::GetInstance().GetTrackLength(start, stop);
  std::cout << "Track-length between (x,y,z)=(" << startx << "," << starty << "," << startz << ") and (x,y,z)=(" << endx << "," << endy << "," << endz << ") = " << TrackLength << " g/cm3" << std::endl;

  return 0;
}

