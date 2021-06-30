#include "TMS_TrueParticle.h"

void TMS_TrueParticle::Print() {
  std::cout << "Printing TMS_TrueParticle class: " << std::endl;
  std::cout << "  Parent: " << Parent << std::endl;
  std::cout << "  TrackId: " << TrackId << std::endl;
  std::cout << "  PDG: " << PDG << std::endl;
  std::cout << "  Vertex ID: " << VertexID << std::endl;
  std::cout << "  Birth position: ";
  Position.Print();
  std::cout << "  Birth momentum: " << Birth_Momentum.Mag() << std::endl;
  Birth_Momentum.Print();

  std::cout << "  TMS momentum: " << TMS_Momentum.Mag() << std::endl;
  TMS_Momentum.Print();
  std::cout << "  TMS entry point: " << std::endl;
  TMS_EntryPoint.Print();

  std::cout << "  Size of points vector: " << PositionPoints.size() << std::endl;
  std::cout << "  Position of hits: " << std::endl;
  int it = 0;
  for (auto i = PositionPoints.begin(); i != PositionPoints.end(); ++i) {
    std::cout << "  Point " << it << "/" << PositionPoints.size() << std::endl;
    (*i).Print();
  }

  std::cout << "  Momentum at trajectory points: " << MomentumPoints.size() << std::endl;
  it = 0;
  for (auto i = MomentumPoints.begin(); i != MomentumPoints.end(); ++i) {
    std::cout << "  Point " << it << "/" << MomentumPoints.size() << std::endl;
    (*i).Print();
  }

}
