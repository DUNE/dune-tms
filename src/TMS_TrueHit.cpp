#include "TMS_TrueHit.h"
#include "TMS_Readout_Manager.h"

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

TMS_TrueHit::TMS_TrueHit(TG4HitSegment &edep_seg, int vertex_id) {

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

  PrimaryIds.push_back(edep_seg.GetPrimaryId());
  VertexIds.push_back(vertex_id);
  EnergyShare.push_back(GetE());
}

void TMS_TrueHit::Print() const {
  std::cout << "TMS_TrueHit: " << std::endl;
  std::cout << "(x,y,z,t,E): (" << GetX() << ", " << GetY() << ", " << GetZ() << ", " << GetT() << ", " << GetE() << ")" << std::endl;
  std::cout << "PrimaryId: " << GetPrimaryId()  << std::endl;
}

void TMS_TrueHit::MergeWith(TMS_TrueHit& hit) {
  // All these variables are simply sums
  double new_true_e = GetE() + hit.GetE();
  double new_true_pe = GetPE() + hit.GetPE();
  double new_true_short_path_pe = GetPEAfterFibersShortPath() + hit.GetPEAfterFibersShortPath();
  double new_true_long_path_pe = GetPEAfterFibersLongPath() + hit.GetPEAfterFibersLongPath();
  // T is defined as the first hit, so take the minimum
  double new_true_t = std::min(GetT(), hit.GetT());
  
  // Take energy weighted average of positions
  // Doesn't make sense in some cases:
  // like two ~same E hits on opposite lengths of scint bar would have an average in the middle
  double new_x = (GetX() * GetE() + hit.GetX() * hit.GetE()) / new_true_e;
  double new_y = (GetY() * GetE() + hit.GetY() * hit.GetE()) / new_true_e;
  double new_z = (GetZ() * GetE() + hit.GetZ() * hit.GetE()) / new_true_e;
  
  // Is dx just a sum of all the parts? For a single particle that would make sense. Not sure about multiple particles
  double new_dx = GetdX() + hit.GetdX();
  
  // Now set the new variables
  // Make sure to calculate above this line or else values get adjusted
  SetX(new_x);
  SetY(new_y);
  SetZ(new_z);
  SetdX(new_dx);
  SetE(new_true_e);
  SetT(new_true_t);
  SetPE(new_true_pe);
  SetPEAfterFibersShortPath(new_true_short_path_pe);
  SetPEAfterFibersLongPath(new_true_long_path_pe);
  
  // Add to the pid vectors
  for (size_t i = 0; i < hit.PrimaryIds.size(); i++) {
    PrimaryIds.push_back(hit.PrimaryIds[i]);
    VertexIds.push_back(hit.VertexIds[i]);
    EnergyShare.push_back(hit.EnergyShare[i]);
  }
}




