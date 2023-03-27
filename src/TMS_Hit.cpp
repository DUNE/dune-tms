#include "TMS_Hit.h"

// The constructor for a hit in the TMS, from a edep-sim hit

// Set the bar and truehit
TMS_Hit::TMS_Hit(TG4HitSegment &edep_seg, int vertex_id) : 
  TrueHit(edep_seg, vertex_id), 
  Bar(edep_seg),
  EnergyDeposit(edep_seg.GetEnergyDeposit()),
  // Define time as the average between start and stop of hit
  Time((edep_seg.GetStop().T()+edep_seg.GetStart().T())/2), 
  Slice(0), 
  PedSuppressed(false), 
  PE(edep_seg.GetEnergyDeposit() * TMS_Const::TMS_EtoPE) {

  // The true particle
  //TrueParticle = TMS_TrueParticle(edep_seg);

}

bool TMS_Hit::NextToGap() {
  // There are no gaps in xy view, only xz view
  if (Bar.GetBarType() != TMS_Bar::kYBar) return false;

  double pos = GetNotZ() + GetNotZw();
  double neg = GetNotZ() - GetNotZw();
  // Check the top
  if ((pos > TMS_Const::TMS_Dead_Top[0] && pos < TMS_Const::TMS_Dead_Top[1]) ||
      (neg < TMS_Const::TMS_Dead_Top[1] && neg > TMS_Const::TMS_Dead_Top[0])) return true;
  // Check the center
  else if ((pos > TMS_Const::TMS_Dead_Center[0] && pos < TMS_Const::TMS_Dead_Center[1]) ||
      (neg < TMS_Const::TMS_Dead_Center[1] && neg > TMS_Const::TMS_Dead_Center[0])) return true;
  // Check the bottom
  else if ((pos > TMS_Const::TMS_Dead_Bottom[0] && pos < TMS_Const::TMS_Dead_Bottom[1]) ||
      (neg < TMS_Const::TMS_Dead_Bottom[1] && neg > TMS_Const::TMS_Dead_Bottom[0])) return true;

  else return false;

  return false;
}

void TMS_Hit::Print() const {
  std::cout << "Printing TMS hit" << std::endl;
  std::cout << "EnergyDeposit: " << EnergyDeposit << std::endl;
  std::cout << "Time: " << Time << std::endl;
  std::cout << "Bar: " << std::endl;
  Bar.Print();

  std::cout << "TrueHit: " << std::endl;
  TrueHit.Print();
}
