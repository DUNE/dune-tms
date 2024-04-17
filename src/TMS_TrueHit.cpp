#include "TMS_TrueHit.h"
#include "TMS_Readout_Manager.h"

TMS_TrueHit::TMS_TrueHit() :
  x(-99999999),
  y(-99999999),
  z(-99999999),
  t(-99999999),
  EnergyDeposit(-99999999)
{};

/*
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
  TLorentzVector diff = (edep_seg.GetStop() - edep_seg.GetStart());
  SetdX(diff.P());
  SetPE(GetE() * TMS_Readout_Manager::GetInstance().Get_Sim_Optical_LightYield());
  SetPEAfterFibers(GetPE());
  SetPEAfterFibersLongPath(0);
  SetPEAfterFibersShortPath(GetPE());

  PrimaryId = edep_seg.GetPrimaryId();
}

void TMS_TrueHit::Print() const {
  std::cout << "TMS_TrueHit: " << std::endl;
  std::cout << "(x,y,z,t,E): (" << GetX() << ", " << GetY() << ", " << GetZ() << ", " << GetT() << ", " << GetE() << ")" << std::endl;
  std::cout << "PrimaryId: " << PrimaryId  << std::endl;
}

void TMS_TrueHit::MergeWith(TMS_TrueHit& hit) {
  double new_true_e = GetE() + hit.GetE();
  double new_true_pe = GetPE() + hit.GetPE();
  double new_true_short_path_pe = GetPEAfterFibersShortPath() + hit.GetPEAfterFibersShortPath();
  double new_true_long_path_pe = GetPEAfterFibersLongPath() + hit.GetPEAfterFibersLongPath();
  double new_true_t = std::min(GetT(), hit.GetT());
  SetE(new_true_e);
  SetT(new_true_t);
  SetPE(new_true_pe);
  SetPEAfterFibersShortPath(new_true_short_path_pe);
  SetPEAfterFibersLongPath(new_true_long_path_pe);
  // TODO how do we handle dedx? and x/y positions?
}




