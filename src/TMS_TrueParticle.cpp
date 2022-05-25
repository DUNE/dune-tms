#include "TMS_TrueParticle.h"

void TMS_TrueParticle::Print() {
  std::cout << "Printing TMS_TrueParticle class: " << std::endl;
  std::cout << "  Parent: " << Parent << std::endl;
  std::cout << "  TrackId: " << TrackId << std::endl;
  std::cout << "  PDG: " << PDG << std::endl;
  std::cout << "  Vertex ID: " << VertexID << std::endl;
  std::cout << "  Birth position: ";
  BirthPosition.Print();
  std::cout << "  Birth momentum: " << BirthMomentum.Mag() << std::endl;
  BirthMomentum.Print();

  std::cout << "  Death position: ";
  DeathPosition.Print();
  std::cout << "  Death momentum: " << DeathMomentum.Mag() << std::endl;
  DeathMomentum.Print();

  std::cout << "  Size of trajectory points vector: " << PositionPoints.size() << std::endl;
  std::cout << "  Size of momentum at trajectory points vector: " << MomentumPoints.size() << std::endl;
  std::cout << "  Position of trajectory: " << std::endl;
  int it = 0;
  for (auto i = PositionPoints.begin(); i != PositionPoints.end(); ++i,++it) {
    std::cout << "  Point " << it+1 << "/" << PositionPoints.size() << std::endl;
    std::cout << "  Position: " << std::endl;
    (*i).Print();
    std::cout << "  Momentum: " << std::endl;
    MomentumPoints[it].Print();
    std::cout << "  GEANT4 Process, Subprocess: " << Process[it] << ", " << Subprocess[it] << std::endl;
  }
}
