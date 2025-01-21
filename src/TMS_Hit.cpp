#include "TMS_Hit.h"
#include "TMS_Readout_Manager.h"

// The constructor for a hit in the TMS, from a edep-sim hit

// Set the bar and truehit
TMS_Hit::TMS_Hit(TG4HitSegment &edep_seg, int vertex_id) : 
  TrueHit(edep_seg, vertex_id), 
  Bar(edep_seg),
  EnergyDeposit(edep_seg.GetEnergyDeposit()),
  // Define time as the average between start and stop of hit
  Time(std::fmod((edep_seg.GetStop().T()+edep_seg.GetStart().T())/2, TMS_Manager::GetInstance().Get_Nersc_Spill_Period())), 
  Slice(0), 
  #ifdef RECORD_HIT_DEADTIME
  DeadtimeStart(-999.0),
  DeadtimeStop(-999.0),
  #endif
  PedSuppressed(false), 
  PE(edep_seg.GetEnergyDeposit() * TMS_Readout_Manager::GetInstance().Get_Sim_Optical_LightYield()) {

  // The true particle
  //TrueParticle = TMS_TrueParticle(edep_seg);

}

bool TMS_Hit::NextToGap() {
  // There are no gaps in xy view, only xz view
  if (Bar.GetBarType() != TMS_Bar::kUBar && Bar.GetBarType() != TMS_Bar::kVBar) return false;

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

void TMS_Hit::MergeWith(TMS_Hit& hit) {
  SetE(GetE() + hit.GetE());
  SetPE(GetPE() + hit.GetPE());
  SetT(std::min(GetT(), hit.GetT()));
  
  // And merge truth info
  GetAdjustableTrueHit().MergeWith(hit.GetAdjustableTrueHit());
}

double TMS_Hit::GetTrueDistanceFromReadout() {
  // Note that you want to do always do more positive - less positive, or else you get a sign error
  if (GetBar().GetBarType() == TMS_Bar::kXBar) {
    // Readout from sides
    if (GetTrueHit().GetX() < 0) return GetTrueHit().GetX() - TMS_Geom::GetInstance().XBarNegReadoutLocation();
    else return TMS_Geom::GetInstance().XBarPosReadoutLocation() - GetTrueHit().GetX();
  }
  else {
    // Readout from top. Assuming U ~ V ~ Y for now
    return TMS_Geom::GetInstance().YBarReadoutLocation() - GetTrueHit().GetY();
  }
}

double TMS_Hit::GetTrueLongDistanceFromReadout() {
  double additional_length;
  if (GetBar().GetBarType() == TMS_Bar::kXBar) {
    // Readout from sides
    additional_length = 2 * TMS_Geom::GetInstance().XBarLength();
  }
  else {
    // Readout from top. Assuming U ~ V ~ Y for now
    additional_length = 2 * TMS_Geom::GetInstance().YBarLength();
  }
  return additional_length - GetTrueDistanceFromReadout();
}

double TMS_Hit::GetTrueDistanceFromMiddle() {
  // Subtract out half a bar length
  double additional_length;
  if (GetBar().GetBarType() == TMS_Bar::kXBar) {
    additional_length = 0.5 * TMS_Geom::GetInstance().XBarLength();
  }
  else {
    additional_length = 0.5 * TMS_Geom::GetInstance().YBarLength();
  }
  return GetTrueDistanceFromReadout() - additional_length;
}

double TMS_Hit::GetTrueLongDistanceFromMiddle() {
  // Subtract out half a bar length
  double additional_length;
  if (GetBar().GetBarType() == TMS_Bar::kXBar) {
    additional_length = 0.5 * TMS_Geom::GetInstance().XBarLength();
  }
  else {
    additional_length = 0.5 * TMS_Geom::GetInstance().YBarLength();
  }
  return GetTrueLongDistanceFromReadout() - additional_length;
}








