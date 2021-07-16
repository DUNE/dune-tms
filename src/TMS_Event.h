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
    TMS_Event(TG4Event &event);
    //~TMS_Event();

    // The getters once the class is completed
    const std::vector<TMS_Hit> &GetHits() const {return TMS_Hits;};
    // Reconstructed tracks
    //std::vector<TMS_Track> GetTracks() {return TMS_Tracks;};
    // The true particles
    const std::vector<TMS_TrueParticle> &GetTrueParticles() const { return TMS_TrueParticles; };

    void Print();

    // Does event contain any TMS hits? (Can be skipped for TMS reco)
    bool IsEmpty() { 
      return (TMS_Hits.size() == 0 ? true : false);
    }

    int GetEventNumber() { return EventNumber; };
    std::string GetReaction() { return Reaction; };
    
    // Include some truth metadata, like process, energy, lepton momentum?
    void FillTruthFromGRooTracker(int pdg[100], double p4[100][4]);

    int GetNeutrinoPDG() { return TrueNeutrino.second; };
    TLorentzVector GetNeutrinoP4() { return TrueNeutrino.first; };

  private:
    // Hits
    std::vector<TMS_Hit> TMS_Hits;
    // True particles
    std::vector<TMS_TrueParticle> TMS_TrueParticles;
    // Reconstructed tracks
    //std::vector<TMS_Track> TMS_Tracks;

    // Spill number (can have many events in a spill)?
    //int SpillNumber;
    std::string Reaction;
 
    // Counts how many times constructor has been called
    static int EventCounter;
    // Saves the event number for a constructed event
    int EventNumber;

    std::pair<TLorentzVector,int> TrueNeutrino;
};

#endif
