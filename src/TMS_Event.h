#ifndef _TMS_EVENT_H_SEEN_
#define _TMS_EVENT_H_SEEN_

#include <string>
#include <iostream>

// Include the constants
#include "TMS_Constants.h"
#include "TMS_Hit.h"
#include "TMS_TrueParticle.h"
#include "TMS_Geom.h"
//#include "TMS_Tracks.h"

// The edep-sim event class
#include "EDepSim/TG4Event.h"

// The general event class
class TMS_Event {
  public:
    TMS_Event(TG4Event &event, bool FillEvent = true);
    TMS_Event();
    //~TMS_Event();

    void AddEvent(TMS_Event &);

    // The getters once the class is completed
    const std::vector<TMS_Hit> GetHits(int slice = -1, bool include_ped_sup = false);
    // Reconstructed tracks
    //std::vector<TMS_Track> GetTracks() {return TMS_Tracks;};
    // The true particles
    const std::vector<TMS_TrueParticle> &GetTrueParticles() const { return TMS_TrueParticles; };

    double GetMuonTrueKE();
    double GetMuonTrueTrackLength();

    void Print();

    // Does event contain any TMS hits? (Can be skipped for TMS reco)
    bool IsEmpty() { 
      return (TMS_Hits.size() == 0 ? true : false);
    }

    int GetEventNumber() { return EventNumber; };
    std::string GetReaction() { return Reaction; };
    
    // Include some truth metadata, like process, energy, lepton momentum
    void FillTruthFromGRooTracker(int pdg[100], double p4[100][4]);

    int GetNeutrinoPDG() { return TrueNeutrino.second; };
    TLorentzVector GetNeutrinoP4() { return TrueNeutrino.first; };
    
    int GetNSlices() { return NSlices; }; 
    void SetNSlices(int n) { NSlices = n; };

  private:
    bool LightWeight; // Don't save all true trajectories; only save significant ones

    // Hits
    std::vector<TMS_Hit> TMS_Hits;
    
    void ApplyReconstructionEffects();
    void MergeCoincidentHits();
    void SimulateOpticalModel();
    void SimulatePedestalSubtraction();
    void SimulateTimingModel();

    // True particles that create trajectories in TMS or LAr; after G4 is run
    std::vector<TMS_TrueParticle> TMS_TrueParticles;

    // Primary particles from neutrino event; before G4 is run
    std::vector<TMS_TrueParticle> TMS_TruePrimaryParticles;

    // Reconstructed tracks
    //std::vector<TMS_Track> TMS_Tracks;

    // The number of true trajectories right out of edep-sim
    // No energy cuts, or number of deposits etc checked
    int nTrueTrajectories;

    std::string Reaction;
 
    // Counts how many times constructor has been called
    static int EventCounter;

    // Saves the event number for a constructed event
    int EventNumber;
    
    int NSlices;

    // Saves down the true Neutrino information from the gRooTracker passthrough (not available in edep-sim or G4)
    std::pair<TLorentzVector,int> TrueNeutrino;
};

#endif
