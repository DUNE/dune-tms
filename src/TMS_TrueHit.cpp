#include "TMS_TrueHit.h"

/*
TMS_TrueHit::TMS_TrueHit() :
  x(-999.99),
  y(-999.99),
  z(-999.99),
  t(-999.99),
  EnergyDeposit(-999.99)
{};

TMS_TrueHit::TMS_TrueHit(double x, double y, double z, double t, double E) {
  SetX(x);
  SetY(y);
  SetZ(z);
  SetT(t);
  SetE(E);
}
*/

TMS_TrueHit::TMS_TrueHit(TG4HitSegment &edep_seg, int vertex_id) : VertexId(vertex_id) {

  // Set the energy
  SetE(edep_seg.GetEnergyDeposit());

  // Set the primary contributor

  // Set the true x,y,z,t of the hit (not just the bar)
  TLorentzVector avg = (edep_seg.GetStart()+edep_seg.GetStop());
  avg *= 0.5;
  SetX(avg.X());
  SetY(avg.Y());
  SetZ(avg.Z());
  SetT(avg.T());

  PrimaryId = edep_seg.GetPrimaryId();
}

void TMS_TrueHit::Print() const {
  std::cout << "TMS_TrueHit: " << std::endl;
  std::cout << "(x,y,z,t,E): (" << GetX() << ", " << GetY() << ", " << GetZ() << ", " << GetT() << ", " << GetE() << ")" << std::endl;
  std::cout << "PrimaryId: " << PrimaryId  << std::endl;
}
