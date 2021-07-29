#ifndef _TMS_TRUEPARTICLE_H_SEEN_
#define _TMS_TRUEPARTICLE_H_SEEN_

// edep-sim classes for TG4PrimaryParticle
#include "EDepSim/TG4PrimaryVertex.h"
#include "EDepSim/TG4HitSegment.h"

#include "TLorentzVector.h"

#include <iostream>
#include "TMS_Constants.h"

class TMS_TrueParticle {
  public:

    TMS_TrueParticle() :
      VertexID(-999), 
      Parent(-999), 
      TrackId(-999), 
      PDG(-999) {
    }

    // Copy over the edep-sim info
    TMS_TrueParticle(TG4PrimaryParticle &edep_part) :
      VertexID(-999),
      Parent(-999),
      TrackId(edep_part.GetTrackId()),
      PDG(edep_part.GetPDGCode()),
      BirthMomentum(edep_part.GetMomentum().Vect()) {
    }

    TMS_TrueParticle(TG4PrimaryParticle &edep_part, TG4PrimaryVertex &vtx) : 
      VertexID(-999),
      Parent(-999),
      TrackId(edep_part.GetTrackId()),
      PDG(edep_part.GetPDGCode()),
      BirthMomentum(edep_part.GetMomentum().Vect()),
      BirthPosition(vtx.GetPosition()) {
    }

    // Construct directly from edep-sim
    TMS_TrueParticle(int ParentVal, int TrackVal, int PDGVal) :
      VertexID(-999),
      Parent(ParentVal), 
      TrackId(TrackVal), 
      PDG(PDGVal) {
    }

    // Print
    void Print();

    // Add a point in for this true particle
    void AddPoint(TLorentzVector &Position, TVector3 &Momentum) {
      PositionPoints.emplace_back(Position);
      MomentumPoints.emplace_back(Momentum);
    }

    void AddPoint(TLorentzVector &Position) {
      PositionPoints.emplace_back(Position);
    }

    void SetParent(int num) { Parent = num; };
    void SetTrackId(int num) { TrackId = num; };
    void SetPDG(int num) { PDG = num; };

    int GetPDG() { return PDG; };
    int GetParent() { return Parent; };
    int GetTrackId() { return TrackId; };

    std::vector<TLorentzVector> &GetPositionPoints() { return PositionPoints; };
    std::vector<TVector3> &GetMomentumPoints() { return MomentumPoints; };

    void SetVertexID(int vtx) {
      VertexID = vtx;
    }

    void SetBirthMomentum(TVector3 &birthmom) { BirthMomentum = birthmom; };
    void SetBirthPosition(TLorentzVector &birthpos) { BirthPosition = birthpos; };

    void SetDeathMomentum(TVector3 &deathmom) { DeathMomentum = deathmom; };
    void SetDeathPosition(TLorentzVector &deathpos) { DeathPosition = deathpos; };

    TVector3       &GetBirthMomentum() { return BirthMomentum; };
    TLorentzVector &GetBirthPosition() { return BirthPosition; };

    TVector3       &GetDeathMomentum() { return DeathMomentum; };
    TLorentzVector &GetDeathPosition() { return DeathPosition; };

    TVector3       &GetInitialMomentum() { return MomentumPoints.back(); };
    TLorentzVector &GetInitialPoint() { return PositionPoints.back(); };

    double GetBirthEnergy() { 
      double mass = 0;
      if (abs(PDG) == 13) mass = TMS_KinConst::mass_mu;
      else if (abs(PDG) == 11) mass = TMS_KinConst::mass_e;
      else if (abs(PDG) == 15) mass = TMS_KinConst::mass_tau;
      else if (abs(PDG) == 211) mass = TMS_KinConst::mass_pic;
      else if (abs(PDG) == 111) mass = TMS_KinConst::mass_pi0;
      else if (abs(PDG) == 2212) mass = TMS_KinConst::mass_proton;
      else if (abs(PDG) == 2112) mass = TMS_KinConst::mass_neutron;
      return sqrt(BirthMomentum.Mag2()+mass*mass);
    }

  private:
    int VertexID;
    int Parent;
    int TrackId;
    int PDG;

    std::vector<TLorentzVector> PositionPoints;
    std::vector<TVector3> MomentumPoints;

    TVector3 BirthMomentum;
    TLorentzVector BirthPosition;

    TVector3 DeathMomentum;
    TLorentzVector DeathPosition;
};

#endif
