#ifndef _TMS_EVENT_H_SEEN_
#define _TMS_EVENT_H_SEEN_

#include <string>
#include <iostream>

// Include the constants
#include "TMS_Constants.h"
#include "TMS_Hit.h"
#include "TMS_TrueParticle.h"
#include "TMS_Geom.h"
#include <random>
#include "TMS_Track.h"

// The edep-sim event class
#include "EDepSim/TG4Event.h"

// The general event class
class TMS_Event {
  public:
    TMS_Event(TG4Event &event, bool FillEvent = true);
    TMS_Event(TMS_Event &event, int slice);
    TMS_Event();
    //~TMS_Event();

    void AddEvent(TMS_Event &);

    // The getters once the class is completed
    const std::vector<TMS_Hit> GetHits(int slice = -1, bool include_ped_sup = false);
    std::vector<TMS_Hit> GetHitsRaw() { return TMS_Hits; };
    void SetHitsRaw(std::vector<TMS_Hit> hits) { TMS_Hits = hits; };
    // Reconstructed tracks
    std::vector<TMS_Track> GetTracks() {return TMS_Tracks;}; // Needs filled
    // The true particles
    const std::vector<TMS_TrueParticle> &GetTrueParticles() const { return TMS_TrueParticles; };

    double GetMuonTrueKE();
    double GetMuonTrueTrackLength();

    void Print();

    // Does event contain any TMS hits? (Can be skipped for TMS reco)
    bool IsEmpty() { 
      return (TMS_Hits.size() == 0 ? true : false);
    }
    
    int GetNHits() { return TMS_Hits.size(); };

    int GetEventNumber() { return EventNumber; };
    //void SetEventNumber(int num) { EventNumber = num; };
    std::string GetReaction() { return Reaction; };
    
    // Include some truth metadata, like process, energy, lepton momentum
    void FillTruthFromGRooTracker(int pdg[100], double p4[100][4], double vtx[100][4]);

    int GetNeutrinoPDG() { return TrueNeutrino.second; };
    TLorentzVector GetNeutrinoP4() { return TrueNeutrino.first; };
    TLorentzVector GetNeutrinoVtx() { return TrueNeutrinoPosition; };
    int GetLeptonPDG() { return TrueLeptonPDG; };
    TLorentzVector GetLeptonX4() { return TrueLeptonPosition; };
    TLorentzVector GetLeptonP4() { return TrueLeptonMomentum; };
    
    void FillTrueLeptonInfo(int pdg, TLorentzVector position, TLorentzVector momentum);
    
    int GetNSlices() { return NSlices; }; 
    void SetNSlices(int n) { NSlices = n; };
    
    int GetSliceNumber() { return SliceNumber; };
    void SetSliceNumber(int slice) { SliceNumber = slice; };
    
    int GetSpillNumber() { return SpillNumber; };
    void SetSpillNumber(int spill) { SpillNumber = spill; };
    
    void SortHits(bool(*comp)(TMS_Hit& a, TMS_Hit& b)) { std::sort(TMS_Hits.begin(), TMS_Hits.end(), comp); };
    
    std::pair<double, double> GetEventTimeRange();
    
    std::map<int, double>& GetTrueVisibleEnergyPerVertex() { return TrueVisibleEnergyPerVertex; };
    
    void SetTotalVisibleEnergyFromVertex(double energy) { TotalVisibleEnergyFromVertex = energy; };
    int GetVertexIdOfMostVisibleEnergy();
    double GetVisibleEnergyFromVertexInSlice() { return VisibleEnergyFromVertexInSlice; };
    double GetTotalVisibleEnergyFromVertex() { return TotalVisibleEnergyFromVertex; };

    int GetNVertices() { return nVertices; };
    int GetNTrueForgottenParticles() { return nTrueForgottenParticles; };
    double GetVisibleEnergyFromOtherVerticesInSlice() { return VisibleEnergyFromOtherVerticesInSlice; };

    
    std::vector<std::pair<float, float>> GetDeadChannelPositions() { return ChannelPositions; };
    std::vector<std::pair<float, float>> GetDeadChannelTimes() { return DeadChannelTimes; };
    std::vector<std::pair<float, float>> GetReadChannelPositions() { return ChannelPositions; };
    std::vector<std::pair<float, float>> GetReadChannelTimes() { return ReadChannelTimes; };
    
    int GetTrueParticleIndex(int trackid);

  private:
    bool LightWeight; // Don't save all true trajectories; only save significant ones

    // Hits
    std::vector<TMS_Hit> TMS_Hits;
    
    void ApplyReconstructionEffects();
    void MergeCoincidentHits();
    void SimulateOpticalModel();
    void SimulateDeadtime();
    void SimulatePedestalSubtraction();
    void SimulateTimingModel();
    void SimulateDarkCount();
    void SimulateReadoutNoise();

    int GetUniqIDForDeadtime(const TMS_Hit& hit) const;

    // True particles that create trajectories in TMS or LAr; after G4 is run
    std::vector<TMS_TrueParticle> TMS_TrueParticles;
    int nTrueForgottenParticles;

    // Primary particles from neutrino event; before G4 is run
    std::vector<TMS_TrueParticle> TMS_TruePrimaryParticles;

    // Reconstructed tracks
    std::vector<TMS_Track> TMS_Tracks;

    // The number of true trajectories right out of edep-sim
    // No energy cuts, or number of deposits etc checked
    int nTrueTrajectories;
    int nVertices;

    std::string Reaction;
 
    // Counts how many times constructor has been called
    static int EventCounter;

    // Saves the event number for a constructed event
    int EventNumber;
    int SliceNumber;
    
    int NSlices;
    int SpillNumber;

    // Saves down the true Neutrino information from the gRooTracker passthrough (not available in edep-sim or G4)
    std::pair<TLorentzVector,int> TrueNeutrino;
    TLorentzVector TrueNeutrinoPosition;
    int TrueLeptonPDG;
    TLorentzVector TrueLeptonPosition;
    TLorentzVector TrueLeptonMomentum;
    std::map<int, double> TrueVisibleEnergyPerVertex;
    std::map<int, double> TrueVisibleEnergyPerParticle;
    
    int VertexIdOfMostEnergyInEvent;
    double VisibleEnergyFromVertexInSlice;
    double TotalVisibleEnergyFromVertex;
    double VisibleEnergyFromOtherVerticesInSlice;

    std::vector<std::pair<float, float>> ChannelPositions;
    std::vector<std::pair<float, float>> DeadChannelTimes;
    std::vector<std::pair<float, float>> ReadChannelTimes;

    std::default_random_engine generator;
    
};

#endif
