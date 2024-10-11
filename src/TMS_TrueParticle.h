#ifndef _TMS_TRUEPARTICLE_H_SEEN_
#define _TMS_TRUEPARTICLE_H_SEEN_

// edep-sim classes for TG4PrimaryParticle
#include "EDepSim/TG4PrimaryVertex.h"
#include "EDepSim/TG4HitSegment.h"

#include "TLorentzVector.h"

#include <iostream>
#include "TMS_Constants.h"

#include "TMS_Geom.h"

#include <functional>

using IsInsideFunctionType = std::function<bool(TVector3)>;
//using IsInsideFunctionType = bool (TMS_Geom::*InInsideFunction)(Vector3);

class TMS_TrueParticle {
  public:

    TMS_TrueParticle() :
      VertexID(-999), 
      Parent(-999), 
      TrackId(-999), 
      PDG(-999),
      TrueVisibleEnergy(-999) {
    }

    // Copy over the edep-sim info
    TMS_TrueParticle(TG4PrimaryParticle &edep_part) :
      VertexID(-999),
      Parent(-999),
      TrackId(edep_part.GetTrackId()),
      PDG(edep_part.GetPDGCode()),
      TrueVisibleEnergy(-999),
      BirthMomentum(edep_part.GetMomentum().Vect()) {
    }

    TMS_TrueParticle(TG4PrimaryParticle &edep_part, TG4PrimaryVertex &vtx) : 
      VertexID(vtx.GetInteractionNumber()),
      Parent(-999),
      TrackId(edep_part.GetTrackId()),
      PDG(edep_part.GetPDGCode()),
      TrueVisibleEnergy(-999),
      BirthMomentum(edep_part.GetMomentum().Vect()),
      BirthPosition(vtx.GetPosition()) {
      if (VertexID < 0) {
        std::cout<<"Fatal in TMS_TrueParticle: Get a vertex id < 0: "<<VertexID<<std::endl;
        throw std::runtime_error("Fatal in TMS_TrueParticle: Get a vertex id < 0");
      }
    }

    // Construct directly from edep-sim
    TMS_TrueParticle(int ParentVal, int TrackVal, int PDGVal, int VertexIDVal) :
      VertexID(VertexIDVal),
      Parent(ParentVal), 
      TrackId(TrackVal), 
      PDG(PDGVal), 
      TrueVisibleEnergy(-999) {
      if (VertexID < 0) {
        std::cout<<"Fatal: Get a vertex id < 0: "<<VertexID<<std::endl;
        throw std::runtime_error("Fatal: Get a vertex id < 0");
      }
    }

    // Print
    void Print(bool small = false);

    // Add a point in for this true particle
    void AddPoint(TLorentzVector &Position, TVector3 &Momentum) {
      PositionPoints.emplace_back(Position);
      MomentumPoints.emplace_back(Momentum);
    }

    void AddPoint(TLorentzVector &Position) {
      PositionPoints.emplace_back(Position);
    }

    void AddPoint(TLorentzVector &Position, TVector3 &Momentum, int &G4Process, int &G4Subprocess)  {
      PositionPoints.emplace_back(Position);
      MomentumPoints.emplace_back(Momentum);
      Process.emplace_back(G4Process);
      Subprocess.emplace_back(G4Subprocess);
    }

    void SetParent(int num) { Parent = num; };
    void SetTrackId(int num) { TrackId = num; };
    void SetPDG(int num) { PDG = num; };

    int GetPDG() const { return PDG; };
    int GetParent() const { return Parent; };
    bool IsPrimary() const { return Parent < 0; };
    int GetTrackId() const { return TrackId; };
    int GetVertexID() const { return VertexID; };
    double GetTrueVisibleEnergy(bool slice = false) const { if (slice) { return TrueVisibleEnergySlice; } else { return TrueVisibleEnergy; } };
    void SetTrueVisibleEnergy(double energy, bool slice) { if (slice) { TrueVisibleEnergySlice = energy; } else { TrueVisibleEnergy = energy; } };
    int GetNTrueHits(bool slice) const { if (slice) { return NTrueHitsSlice; } else { return NTrueHits; } };
    void SetNTrueHits(int n, bool slice) { if (slice) { NTrueHitsSlice = n; } else { NTrueHits = n; } };

    std::vector<TLorentzVector> &GetPositionPoints() { return PositionPoints; };
    std::vector<TVector3> &GetMomentumPoints() { return MomentumPoints; };
    std::vector<int> &GetProcessPoints() { return Process; };
    std::vector<int> &GetSubprocessPoints() { return Subprocess; };

    void SetVertexID(int vtx) {
      VertexID = vtx;
    }

    void SetBirthMomentum(TVector3 &birthmom) { BirthMomentum = birthmom; };
    void SetBirthPosition(TLorentzVector &birthpos) { BirthPosition = birthpos; };

    void SetDeathMomentum(TVector3 &deathmom) { DeathMomentum = deathmom; };
    void SetDeathPosition(TLorentzVector &deathpos) { DeathPosition = deathpos; };

    TVector3       &GetBirthMomentum() { return BirthMomentum; };
    TLorentzVector GetBirthMomentumAsLorentz() 
      { auto mom = GetBirthMomentum(); auto en = GetBirthEnergy(); return TLorentzVector(mom.X(), mom.Y(), mom.Z(), en); };
    TLorentzVector &GetBirthPosition() { return BirthPosition; };

    TVector3       &GetDeathMomentum() { return DeathMomentum; };
    TLorentzVector &GetDeathPosition() { return DeathPosition; };

    TVector3       &GetInitialMomentum() { return MomentumPoints.back(); };
    TLorentzVector &GetInitialPoint() { return PositionPoints.back(); };
    
    TLorentzVector GetPositionAtZ(double z, double max_z_dist = 220); // About 2 planes in either direction is the max z distance we'll tolerate, 110mm / thick plane
    TLorentzVector GetPositionZIsLArEnd() { return GetPositionAtZ(TMS_Geom::GetInstance().GetZEndOfLAr()); };
    TLorentzVector GetPositionZIsTMSStart() { return GetPositionAtZ(TMS_Geom::GetInstance().GetZStartOfTMS()); };
    TLorentzVector GetPositionZIsTMSEnd() { return GetPositionAtZ(TMS_Geom::GetInstance().GetZEndOfTMS()); };
    TLorentzVector GetPositionEntering(IsInsideFunctionType isInside);
    TLorentzVector GetPositionLeaving(IsInsideFunctionType isInside);
    TLorentzVector GetPositionEnteringTMS() { return GetPositionEntering(TMS_Geom::StaticIsInsideTMS); };
    TLorentzVector GetPositionLeavingTMS() { return GetPositionLeaving(TMS_Geom::StaticIsInsideTMS); };
    TLorentzVector GetPositionEnteringTMSThin() { return GetPositionEntering(TMS_Geom::StaticIsInsideTMSThin); };
    TLorentzVector GetPositionLeavingTMSThin() { return GetPositionLeaving(TMS_Geom::StaticIsInsideTMSThin); };
    TLorentzVector GetPositionEnteringTMSFirstTwoModules() { return GetPositionEntering(TMS_Geom::StaticIsInsideTMSFirstTwoModules); };
    TLorentzVector GetPositionLeavingTMSFirstTwoModules() { return GetPositionLeaving(TMS_Geom::StaticIsInsideTMSFirstTwoModules); };
    TLorentzVector GetPositionEnteringLAr() { return GetPositionEntering(TMS_Geom::StaticIsInsideLAr); };
    TLorentzVector GetPositionLeavingLAr() { return GetPositionLeaving(TMS_Geom::StaticIsInsideLAr); };
    
    std::vector<TVector3> GetPositionPoints(double z_start, double z_end, bool onlyInsideTMS = false);
    
    TLorentzVector GetMomentumAtZ(double z, double max_z_dist = 220); // About 2 planes in either direction is the max z distance we'll tolerate, 110mm / thick plane
    TLorentzVector GetMomentumZIsLArEnd() { return GetMomentumAtZ(TMS_Geom::GetInstance().GetZEndOfLAr()); };
    TLorentzVector GetMomentumZIsTMSStart() { return GetMomentumAtZ(TMS_Geom::GetInstance().GetZStartOfTMS()); };
    TLorentzVector GetMomentumZIsTMSEnd() { return GetMomentumAtZ(TMS_Geom::GetInstance().GetZEndOfTMS()); };
    
    TLorentzVector GetMomentumEntering(IsInsideFunctionType isInside);
    TLorentzVector GetMomentumLeaving(IsInsideFunctionType isInside);
    TLorentzVector GetMomentumEnteringTMS() { return GetMomentumEntering(TMS_Geom::StaticIsInsideTMS); };
    TLorentzVector GetMomentumLeavingTMS() { return GetMomentumLeaving(TMS_Geom::StaticIsInsideTMS); };
    TLorentzVector GetMomentumEnteringTMSThin() { return GetMomentumEntering(TMS_Geom::StaticIsInsideTMSThin); };
    TLorentzVector GetMomentumLeavingTMSThin() { return GetMomentumLeaving(TMS_Geom::StaticIsInsideTMSThin); };
    TLorentzVector GetMomentumEnteringTMSFirstTwoModules() { return GetMomentumEntering(TMS_Geom::StaticIsInsideTMSFirstTwoModules); };
    TLorentzVector GetMomentumLeavingTMSFirstTwoModules() { return GetMomentumLeaving(TMS_Geom::StaticIsInsideTMSFirstTwoModules); };
    TLorentzVector GetMomentumEnteringLAr() { return GetMomentumEntering(TMS_Geom::StaticIsInsideLAr); };
    TLorentzVector GetMomentumLeavingLAr() { return GetMomentumLeaving(TMS_Geom::StaticIsInsideLAr); };
    
    bool EntersVolume(IsInsideFunctionType isInside);
    
    double GetEnergyFromMomentum(TVector3 momentum) {
      double mass = TMS_KinConst::GetMass(PDG);
      return sqrt(momentum.Mag2()+mass*mass);
    }

    double GetBirthEnergy() { 
      return GetEnergyFromMomentum(BirthMomentum);
    }

    double GetDeathEnergy() { 
      return GetEnergyFromMomentum(DeathMomentum);
    }

  private:
    int VertexID;
    int Parent;
    int TrackId;
    int PDG;
    double TrueVisibleEnergy;
    int NTrueHits;
    double TrueVisibleEnergySlice;
    int NTrueHitsSlice;

    std::vector<TLorentzVector> PositionPoints;
    std::vector<TVector3> MomentumPoints;
    std::vector<int> Process;
    std::vector<int> Subprocess;

    TVector3 BirthMomentum;
    TLorentzVector BirthPosition;

    TVector3 DeathMomentum;
    TLorentzVector DeathPosition;
};

#endif
