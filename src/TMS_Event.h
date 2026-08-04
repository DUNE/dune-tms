#ifndef _TMS_EVENT_H_SEEN_
#define _TMS_EVENT_H_SEEN_

#include <string>
#include <iostream>

// Include the constants
#include "TMS_Constants.h"
#include "TMS_Hit.h"
#include "TMS_TrueParticle.h"
#include "TMS_Geom.h"
#include <map>
#include <random>
#include "TMS_Track.h"

// The edep-sim event class
#include "EDepSim/TG4Event.h"

struct Vtx_Info {
  TLorentzVector vtx;
  TLorentzVector p4;
  
  int pdg { -999 };
  int run_id { -999 };
  int vtx_id { -999 };
  std::string reaction;
  double hadronic_energy_lar_shell {};
  double hadronic_energy_lar {};
  double hadronic_energy_tms {};
  double hadronic_energy_total {};
  double true_visible_energy_tms {};
  double true_visible_energy_lar {};
  double true_visible_energy_total {};
  
  bool fiducial_cut {};
  // Default to true since hadronic_energy_lar_shell defaults to zero
  bool shell_energy_cut { true };
  bool nd_physics_cut {};
  
  void AddEnergyFromHit(const TMS_TrueHit& hit, int index);
  void UpdateShellEnergyCut() {
    shell_energy_cut = hadronic_energy_lar_shell <= TMS_Manager::GetInstance().Get_ND_PHYSICS_MUON_LAR_SHELL_CUT_ENERGY();
    UpdateNDPhysicsCut(); // Update dirty cut
  };
  void UpdateNDPhysicsCut() {
    nd_physics_cut = fiducial_cut && shell_energy_cut;
  };
  void SetVtx(TLorentzVector pos) {
    vtx = pos;
    fiducial_cut = TMS_Geom::GetInstance().IsInsideLarFiducial(vtx.Vect());
    UpdateNDPhysicsCut(); // Update dirty cut
  };
};

using VertexInfoMap = std::map<long long, Vtx_Info>;

// The general event class
class TMS_Event {
  public:
    TMS_Event(TG4Event event, bool FillEvent = true);
    TMS_Event(TMS_Event &event, int slice);
    TMS_Event();
    //~TMS_Event();
    
    void ProcessTG4Event(TG4Event &event, bool FillEvent = true);

    void AddEvent(TMS_Event &);
    void OverlayEvents(std::vector<TMS_Event>& overlay_events);
    void FinalizeEvent();

    // The getters once the class is completed
    const std::vector<TMS_Hit> GetHits(int slice = -1, bool include_ped_sup = false);
    std::vector<TMS_Hit> GetHitsRaw() { return TMS_Hits; };
    void SetHitsRaw(std::vector<TMS_Hit> hits) { TMS_Hits = hits; };

    // Phase III truth/reco decoupling: per-hit truth, keyed by TMS_Hit::GetHitId() rather than
    // vector position, since hits get sorted/merged/sliced and a positional index wouldn't
    // survive that. Absent for a hit with no truth (e.g. real data) -- callers must null-check.
    int NextHitId() { return HitIdCounter++; };
    const TMS_TrueHit* GetTrueHit(int hitId) const {
      auto it = TrueHitByHitId.find(hitId);
      return it == TrueHitByHitId.end() ? nullptr : &it->second;
    };
    TMS_TrueHit* GetAdjustableTrueHit(int hitId) {
      auto it = TrueHitByHitId.find(hitId);
      return it == TrueHitByHitId.end() ? nullptr : &it->second;
    };
    // Not std::map::operator[] -- TMS_TrueHit has no default constructor, so operator[]'s
    // "default-construct then assign" on a missing key won't compile. find()+emplace/assign
    // only ever uses the copy constructor/assignment, both of which TMS_TrueHit does have.
    void SetTrueHit(int hitId, const TMS_TrueHit& hit) {
      auto it = TrueHitByHitId.find(hitId);
      if (it != TrueHitByHitId.end()) it->second = hit;
      else TrueHitByHitId.emplace(hitId, hit);
    };
    void EraseTrueHit(int hitId) { TrueHitByHitId.erase(hitId); };
    // For TMS_SignalProcessing::MergeCoincidentHits(): merge mergedAwayHitId's truth into
    // survivingHitId's (TMS_TrueHit::MergeWith), then drop the merged-away entry. No-op if
    // either hit has no truth (e.g. real data).
    void MergeTrueHit(int survivingHitId, int mergedAwayHitId) {
      auto itSurvive = TrueHitByHitId.find(survivingHitId);
      auto itMerged = TrueHitByHitId.find(mergedAwayHitId);
      if (itSurvive != TrueHitByHitId.end() && itMerged != TrueHitByHitId.end()) {
        itSurvive->second.MergeWith(itMerged->second);
      }
      if (itMerged != TrueHitByHitId.end()) TrueHitByHitId.erase(itMerged);
    };
    // Reference access for TMS_DetectorSimulation/TMS_SignalProcessing to mutate hits in place
    // without the copy cost of GetHitsRaw()/SetHitsRaw().
    std::vector<TMS_Hit>& GetHitsRawRef() { return TMS_Hits; };
    // Reconstructed tracks
    std::vector<TMS_Track> GetTracks() {return TMS_Tracks;}; // Needs filled
    // The true particles
    const std::vector<TMS_TrueParticle> &GetTrueParticles() const { return TMS_TrueParticles; };
    
    bool IsInTimeSlice(double t) const;

    double GetMuonTrueKE();
    double GetMuonTrueTrackLength();

    void Print();

    // Does event contain any TMS hits? (Can be skipped for TMS reco)
    bool IsEmpty() { 
      return (TMS_Hits.size() == 0 ? true : false);
    }
    
    int GetNHits() { return TMS_Hits.size(); };

    int GetRunNumber() { return RunNumber; };
    int GetEventNumber() { return EventNumber; };
    //void SetEventNumber(int num) { EventNumber = num; };
    std::string GetReaction() { return Reaction; };
    void SetReaction(std::string reaction) { Reaction = reaction; };
    
    // Include some truth metadata, like process, energy, lepton momentum
    void FillTruthFromGRooTracker(int pdg[100], double p4[100][4], double vtx[100][4]);

    int GetNeutrinoPDG() { return TrueNeutrino.second; };
    TLorentzVector GetNeutrinoP4() { return TrueNeutrino.first; };
    TLorentzVector GetNeutrinoVtx() { return TrueNeutrinoPosition; };
    int GetLeptonPDG() { return TrueLeptonPDG; };
    TLorentzVector GetLeptonX4() { return TrueLeptonPosition; };
    TLorentzVector GetLeptonP4() { return TrueLeptonMomentum; };
    
    void FillTrueLeptonInfo(int pdg, TLorentzVector position, TLorentzVector momentum, int vertexid);
    
    int GetNSlices() const { return NSlices; }; 
    void SetNSlices(int n) { NSlices = n; };
    
    int GetSliceNumber() const { return SliceNumber; };
    void SetSliceNumber(int slice) { SliceNumber = slice; };
    
    int GetSpillNumber() const { return SpillNumber; };
    void SetSpillNumber(int spill) { SpillNumber = spill; };
    
    void SortHits(bool(*comp)(TMS_Hit& a, TMS_Hit& b)) { std::sort(TMS_Hits.begin(), TMS_Hits.end(), comp); };
    
    std::pair<double, double> GetEventTimeRange();
    
    std::map<long long, double>& GetTrueVisibleEnergyPerVertex() { return TrueVisibleEnergyPerVertex; };
    
    void SetTotalVisibleEnergyFromVertex(double energy) { TotalVisibleEnergyFromVertex = energy; };
    int GetVertexIdOfMostVisibleEnergy();
    long long GetVertexGlobalIdOfMostVisibleEnergy() { return VertexGlobalIdOfMostEnergyInEvent; };
    double GetVisibleEnergyFromVertexInSlice() { return VisibleEnergyFromVertexInSlice; };
    double GetTotalVisibleEnergyFromVertex() { return TotalVisibleEnergyFromVertex; };

    int GetNVertices() { return nVertices; };
    int GetNTrueForgottenParticles() { return nTrueForgottenParticles; };
    double GetVisibleEnergyFromOtherVerticesInSlice() { return VisibleEnergyFromOtherVerticesInSlice; };

    
    std::vector<std::pair<float, float>> GetDeadChannelPositions() { return ChannelPositions; };
    std::vector<std::pair<float, float>> GetDeadChannelTimes() { return DeadChannelTimes; };
    std::vector<std::pair<float, float>> GetReadChannelPositions() { return ChannelPositions; };
    std::vector<std::pair<float, float>> GetReadChannelTimes() { return ReadChannelTimes; };
    // For TMS_DetectorSimulation::SimulateDeadtime() to record a channel's readout/deadtime window
    void AddDeadtimeChannelRecord(std::pair<float, float> position, std::pair<float, float> deadtime_range, std::pair<float, float> readout_range) {
      ChannelPositions.push_back(position);
      DeadChannelTimes.push_back(deadtime_range);
      ReadChannelTimes.push_back(readout_range);
    };
    
    int GetTrueParticleIndex(long long vertexglobalid, int trackid);
    
    void ApplyReconstructionEffects();
    
    void SetLeptonInfoUsingGlobalVertexID(long long vertexglobalid);
    
    void AddTimeSliceInformation(std::vector<std::pair<double, double>> time_slice_bounds) 
         { TimeSliceBounds.insert(TimeSliceBounds.end(), time_slice_bounds.begin(), time_slice_bounds.end()); };
     std::pair<double, double> GetTimeSliceBounds(int slice = -1);
         
     double CalculateEnergyInLArOuterShell(double thickness, long long vertexglobalid = -1);
     double CalculateEnergyInLAr(long long vertexglobalid = -1);
     double CalculateTotalNonTMSEnergy(long long vertexglobalid = -1);
     
     void ConnectTrueHitWithTrueParticle(bool slide);
     
     VertexInfoMap GetVertexInfo() { return info_about_vtx; };
     Vtx_Info* GetVertexInfo(int run_id, int vertex_id);
     Vtx_Info* GetVertexInfoByGlobalID(long long vertex_global_id);

  private:
    bool LightWeight; // Don't save all true trajectories; only save significant ones

    // Hits
    std::vector<TMS_Hit> TMS_Hits;
    std::vector<TMS_TrueHit> NonTMS_Hits;

    // See NextHitId()/GetTrueHit(int)/SetTrueHit(int, ...) above.
    int HitIdCounter = 0;
    std::map<int, TMS_TrueHit> TrueHitByHitId;

    int GetPrimaryLeptonOfGlobalVertexID(long long vertexglobalid);
    
    void SaveKeyVertexInfo(const TMS_TrueHit& hit);

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
    std::map<long long, std::string> Reactions;
 
    // Counts how many times constructor has been called
    static int EventCounter;

    // Saves the event number for a constructed event
    int RunNumber;
    int EventNumber;
    int SliceNumber;
    
    int NSlices;
    int SpillNumber;

    // Saves down the true Neutrino information from the gRooTracker passthrough (not available in edep-sim or G4)
    std::pair<TLorentzVector,int> TrueNeutrino;
    TLorentzVector TrueNeutrinoPosition;
    int TrueLeptonPDG;
    int TrueLeptonVertexID;
    TLorentzVector TrueLeptonPosition;
    TLorentzVector TrueLeptonMomentum;
    std::map<long long, double> TrueVisibleEnergyPerVertex;
    std::map<std::pair<long long, int>, double> TrueVisibleEnergyPerParticle;
    
    int VertexIdOfMostEnergyInEvent;
    long long VertexGlobalIdOfMostEnergyInEvent;
    double VisibleEnergyFromVertexInSlice;
    double TotalVisibleEnergyFromVertex;
    double VisibleEnergyFromOtherVerticesInSlice;

    std::vector<std::pair<float, float>> ChannelPositions;
    std::vector<std::pair<float, float>> DeadChannelTimes;
    std::vector<std::pair<float, float>> ReadChannelTimes;
    std::vector<std::pair<double, double>> TimeSliceBounds;
    
    VertexInfoMap info_about_vtx;

    std::default_random_engine generator;
    
};

#endif
