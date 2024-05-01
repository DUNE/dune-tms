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

TLorentzVector TMS_TrueParticle::GetMomentumAtZ(double z, double max_z_dist) {
  TVector3 out(-99999999, -99999999, -99999999);
  double z_dist_found = 999999;
  for (size_t i = 0; i < GetPositionPoints().size(); i++) {
    double distance = abs(GetPositionPoints()[i].Z() - z);
    if (distance <= max_z_dist) {
      // Found a candidate
      if (distance < z_dist_found) {
        out = GetMomentumPoints()[i];
        z_dist_found = distance;
      }
    }
  } 
  double energy = GetEnergyFromMomentum(out);
  return TLorentzVector(out.Px(), out.Py(), out.Pz(), energy);
}

std::vector<TVector3> TMS_TrueParticle::GetPositionPoints(double z_start, double z_end, bool onlyInsideTMS) {
  std::vector<TVector3> out;
  for (size_t i = 0; i < GetPositionPoints().size(); i++) {
    double z = GetPositionPoints()[i].Z();
    if (z >= z_start && z <= z_end) {
      double x = GetPositionPoints()[i].X();
      double y = GetPositionPoints()[i].Y();
      TVector3 point(x, y, z);
      bool add = false;
      if (onlyInsideTMS) {
        add = TMS_Geom::StaticIsInsideTMS(point);
      }
      else add = TMS_Geom::StaticIsInsideReasonableSize(point);
      if (add) out.push_back(point);
    }
  }
  return out;
}

TLorentzVector TMS_TrueParticle::GetPositionAtZ(double z, double max_z_dist) {
  TLorentzVector out(-99999999, -99999999, -99999999, -99999999);
  double z_dist_found = 999999;
  for (size_t i = 0; i < GetPositionPoints().size(); i++) {
    double distance = abs(GetPositionPoints()[i].Z() - z);
    if (distance <= max_z_dist) {
      // Found a candidate
      if (distance < z_dist_found) {
        out = GetPositionPoints()[i];
        z_dist_found = distance;
      }
    }
  } 
  return out;
}

TLorentzVector TMS_TrueParticle::GetPositionEntering(IsInsideFunctionType isInside) {
  TLorentzVector out(-99999999, -99999999, -99999999, -99999999);
  for (size_t i = 0; i < GetPositionPoints().size(); i++) {
    // First time this is true means we are inside the volume
    if (isInside(GetPositionPoints()[i].Vect())) {
      out = GetPositionPoints()[i];
      break;
    }
  } 
  return out;
}

TLorentzVector TMS_TrueParticle::GetPositionLeaving(IsInsideFunctionType isInside) {
  TLorentzVector out(-99999999, -99999999, -99999999, -99999999);
  bool areInside = false;
  for (size_t i = 0; i < GetPositionPoints().size(); i++) {
    // First time this is true means we are inside the volume
    // but then the first time it's false means we left the volume
    if (isInside(GetPositionPoints()[i].Vect())) {
      areInside = true;
    }
    else {
      if (areInside) {
        // We were inside but are no longer inside
        areInside = false;
        // At this point we should break since we want the first exiting point
        break;
      }
    }
    // Update the position as long as we're inside the volume
    if (areInside) out = GetPositionPoints()[i];
  } 
  return out;
}

TLorentzVector TMS_TrueParticle::GetMomentumEntering(IsInsideFunctionType isInside) {
  TVector3 out(-99999999, -99999999, -99999999);
  for (size_t i = 0; i < GetPositionPoints().size(); i++) {
    // First time this is true means we are inside the volume
    if (isInside(GetPositionPoints()[i].Vect())) {
      out = GetMomentumPoints()[i];
      break;
    }
  } 
  double energy = GetEnergyFromMomentum(out);
  return TLorentzVector(out.Px(), out.Py(), out.Pz(), energy);
}

TLorentzVector TMS_TrueParticle::GetMomentumLeaving(IsInsideFunctionType isInside) {
  TVector3 out(-99999999, -99999999, -99999999);
  bool areInside = false;
  for (size_t i = 0; i < GetPositionPoints().size(); i++) {
    // First time this is true means we are inside the volume
    // but then the first time it's false means we left the volume
    if (isInside(GetPositionPoints()[i].Vect())) {
      areInside = true;
    }
    else {
      // Were we inside the volume yet?
      if (areInside) {
        // We were inside but are no longer inside
        areInside = false;
        // At this point we should break since we want the first exiting point
        break;
      }
    }
    // Update the momentum as long as we're inside the volume
    if (areInside) out = GetMomentumPoints()[i];
  } 
  double energy = GetEnergyFromMomentum(out);
  return TLorentzVector(out.Px(), out.Py(), out.Pz(), energy);
}
