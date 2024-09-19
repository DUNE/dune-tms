#include "TMS_TrueParticle.h"

#include <algorithm>    // std::swap

TMS_TrueParticle::~TMS_TrueParticle() {
  // Don't delete the elements since root will do that
  PositionPoints.clear();
  MomentumPoints.clear();
  
  DeathMomentum = TVector3();
  BirthMomentum = TVector3();
  BirthPosition = TLorentzVector();
  DeathPosition = TLorentzVector();
}

void TMS_TrueParticle::Print(bool small) {
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
  if (!small) {
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
}

TLorentzVector TMS_TrueParticle::GetMomentumAtZ(double z, double max_z_dist) {
  // Finds the true momentum of the particle at z.
  // If z is outside the true range of the particle, will still return edges if within range + max_z_dist.
  // This is useful because reco only measures z in the middle of the scintillator, but the true particle
  // can start after that. 
  // If within the range of the particle, it does it lerp (linear interpolation) between 
  // the two closest points so that you get the exact x,y, and t position where the z was hit
  
  bool found_out = false;
  TVector3 out(-99999999, -99999999, -99999999);
  
  // Need at least one position point to check
  if (GetPositionPoints().size() > 0) {
    // Then check for cases where position is outside the range
    double z_start = GetPositionPoints()[0].Z();
    double z_end = GetPositionPoints()[GetPositionPoints().size()-1].Z();
    // It's possible that the track is going backwards, so check for that
    bool is_backwards = false;
    if (z_start > z_end) is_backwards = true;
    if (is_backwards) {
      // In backwards case, z_start and z_end are reversed
      std::swap(z_start, z_end);   
    }
    if (z < z_start || z > z_end || GetPositionPoints().size() == 1) {
      // Case where Z is outside range, but it's possible that it's still within max_z_dist
      // Or there's only one point so that z_start == z_end
      if (z < z_start && z >= z_start - max_z_dist) {
        if (!is_backwards) out = GetMomentumPoints()[0];
        else out = GetMomentumPoints()[GetMomentumPoints().size()-1];
        found_out = true;
        //std::cout<<"For z="<<z<<",\tmax_z_dist="<<max_z_dist<<",\tfound out: "<<out.X()<<",\t"<<out.Y()<<",\t"<<out.Z()<<",\t"<<out.T()<<",\tz<z_start case,\t"<<(z_start - max_z_dist)<<std::endl;
      }
      if (z > z_end && z <= z_end + max_z_dist) {
        if (is_backwards) out = GetMomentumPoints()[0];
        else out = GetMomentumPoints()[GetMomentumPoints().size()-1];
        found_out = true;
        //std::cout<<"For z="<<z<<",\tmax_z_dist="<<max_z_dist<<",\tfound out: "<<out.X()<<",\t"<<out.Y()<<",\t"<<out.Z()<<",\t"<<out.T()<<",\tz>z_end case,\t"<<(z_end + max_z_dist)<<std::endl;
      }
    }
    else {
      // Main cases where z is inside range of true particle
      for (size_t i = 0; i < GetPositionPoints().size() - 1; i++) {
        auto upstream_point = GetPositionPoints()[i].Z();
        auto downstream_point = GetPositionPoints()[i+1].Z();
        if (is_backwards) std::swap(upstream_point, downstream_point);
        if (upstream_point <= z && z <= downstream_point) {
          // Found case where z is between upstream and downstream points
          // Find lerping constant
          // We want a = 1 if z == upstream_point, and a = 0 if z == downstream_point
          double a = 1;
          // Make sure to not get nan
          if (upstream_point != downstream_point) a = (z - downstream_point) / (upstream_point - downstream_point);
          // By definition, this will give the reverse answer
          if (is_backwards) a = 1 - a;
          out = GetMomentumPoints()[i] * a + GetMomentumPoints()[i+1] * (1 - a);
          //std::cout<<"For z="<<z<<",\tup_z="<<upstream_point<<",\tdown_z="<<downstream_point<<",\tnpoints="<<GetPositionPoints().size()<<",\tfound out: "<<out.X()<<",\t"<<out.Y()<<",\t"<<out.Z()<<",\ta="<<a<<std::endl;
          found_out = true;
        }
      }
    }
  }
  else { std::cout<<"Found GetPositionPoints().size()==0 case"<<std::endl; }
  
  double energy = -99999999; 
  if (found_out) energy = GetEnergyFromMomentum(out);
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
  // Finds the true position of the particle at z.
  // If z is outside the true range of the particle, will still return edges if within range + max_z_dist.
  // This is useful because reco only measures z in the middle of the scintillator, but the true particle
  // can start after that. 
  // If within the range of the particle, it does it lerp (linear interpolation) between 
  // the two closest points so that you get the exact x,y, and t position where the z was hit
  
  TLorentzVector out(-99999999, -99999999, -99999999, -99999999);
  
  // Need at least one position point to check
  if (GetPositionPoints().size() > 0) {
    // Then check for cases where position is outside the range
    double z_start = GetPositionPoints()[0].Z();
    double z_end = GetPositionPoints()[GetPositionPoints().size()-1].Z();
    // It's possible that the track is going backwards, so check for that
    bool is_backwards = false;
    if (z_start > z_end) is_backwards = true;
    if (is_backwards) {
      // In backwards case, z_start and z_end are reversed
      std::swap(z_start, z_end);   
    }
    if (z < z_start || z > z_end || GetPositionPoints().size() == 1) {
      // Case where Z is outside range, but it's possible that it's still within max_z_dist
      // Or there's only one point so that z_start == z_end
      if (z < z_start && z >= z_start - max_z_dist) {
        if (!is_backwards) out = GetPositionPoints()[0];
        else out = GetPositionPoints()[GetPositionPoints().size()-1];
        //std::cout<<"For z="<<z<<",\tmax_z_dist="<<max_z_dist<<",\tfound out: "<<out.X()<<",\t"<<out.Y()<<",\t"<<out.Z()<<",\t"<<out.T()<<",\tz<z_start case,\t"<<(z_start - max_z_dist)<<std::endl;
      }
      if (z > z_end && z <= z_end + max_z_dist) {
        if (is_backwards) out = GetPositionPoints()[0];
        else out = GetPositionPoints()[GetPositionPoints().size()-1];
        //std::cout<<"For z="<<z<<",\tmax_z_dist="<<max_z_dist<<",\tfound out: "<<out.X()<<",\t"<<out.Y()<<",\t"<<out.Z()<<",\t"<<out.T()<<",\tz>z_end case,\t"<<(z_end + max_z_dist)<<std::endl;
      }
    }
    else {
      // Main cases where z is inside range of true particle
      for (size_t i = 0; i < GetPositionPoints().size() - 1; i++) {
        auto upstream_point = GetPositionPoints()[i].Z();
        auto downstream_point = GetPositionPoints()[i+1].Z();
        if (is_backwards) std::swap(upstream_point, downstream_point);
        if (upstream_point <= z && z <= downstream_point) {
          // Found case where z is between upstream and downstream points
          // Find lerping constant
          // We want a = 1 if z == upstream_point, and a = 0 if z == downstream_point
          double a = 1;
          // Make sure to not get nan
          if (upstream_point != downstream_point) a = (z - downstream_point) / (upstream_point - downstream_point);
          // By definition, this will give the reverse answer
          if (is_backwards) a = 1 - a;
          out = GetPositionPoints()[i] * a + GetPositionPoints()[i+1] * (1 - a);
          // std::cout<<"For z="<<z<<",\tfound out: "<<out.X()<<",\t"<<out.Y()<<",\t"<<out.Z()<<",\t"<<out.T()<<",\ta="<<a<<std::endl;
        }
      }
    }
  }
  else { std::cout<<"Found GetPositionPoints().size()==0 case"<<std::endl; }
  
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

bool TMS_TrueParticle::EntersVolume(IsInsideFunctionType isInside) {
  bool out = false;
  for (size_t i = 0; i < GetPositionPoints().size(); i++) {
    if (isInside(GetPositionPoints()[i].Vect())) {
      // Found a single example, we can end our search
      out = true;
      break;
    }
  }
  return out;
}
