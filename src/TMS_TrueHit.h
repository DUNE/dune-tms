#ifndef _TMS_TRUEHIT_H_
#define _TMS_TRUEHIT_H_

#include <vector>
#include <iostream>

// Include the constants
#include "TMS_Constants.h"

#include "EDepSim/TG4HitSegment.h"

// Essentially a copy of the edep-sim THit
class TMS_TrueHit {
  public:
    TMS_TrueHit(TG4HitSegment &edep_seg, long long vertex_global_id);
    
    TMS_TrueHit() = delete;
    
    // Explicitly create copy construct and assignment operators
    // Or else true particle information is lost 
    // Copy constructor
    TMS_TrueHit(const TMS_TrueHit& other) : PrimaryIds(other.PrimaryIds),
      VertexGlobalIds(other.VertexGlobalIds), EnergyShare(other.EnergyShare), EnergyShareIsLeptonic(other.EnergyShareIsLeptonic),
      LightPrimaryIds(other.LightPrimaryIds), LightVertexGlobalIds(other.LightVertexGlobalIds), LightPhotonCounts(other.LightPhotonCounts)
    {
        if (this != &other) {
          x = other.x;
          y = other.y;
          z = other.z;
          t = other.t;
          dx = other.dx;
          EnergyDeposit = other.EnergyDeposit;
          pe = other.pe;
          peAfterFibers = other.peAfterFibers;
          peAfterFibersLongPath = other.peAfterFibersLongPath;
          peAfterFibersShortPath = other.peAfterFibersShortPath;
          NPhotons = other.NPhotons;
          PrimaryIdByLight = other.PrimaryIdByLight;
          VertexGlobalIdByLight = other.VertexGlobalIdByLight;
          LightShare = other.LightShare;
          FirstPhotonPrimaryId = other.FirstPhotonPrimaryId;
          FirstPhotonVertexGlobalId = other.FirstPhotonVertexGlobalId;
        }
    }
    
    // Copy assignment operator
    TMS_TrueHit& operator=(const TMS_TrueHit& other) {
        if (this != &other) {
          PrimaryIds = other.PrimaryIds;
          VertexGlobalIds = other.VertexGlobalIds;
          EnergyShare = other.EnergyShare;
          EnergyShareIsLeptonic = other.EnergyShareIsLeptonic;
          LightPrimaryIds = other.LightPrimaryIds;
          LightVertexGlobalIds = other.LightVertexGlobalIds;
          LightPhotonCounts = other.LightPhotonCounts;
          x = other.x;
          y = other.y;
          z = other.z;
          t = other.t;
          dx = other.dx;
          EnergyDeposit = other.EnergyDeposit;
          pe = other.pe;
          peAfterFibers = other.peAfterFibers;
          peAfterFibersLongPath = other.peAfterFibersLongPath;
          peAfterFibersShortPath = other.peAfterFibersShortPath;
          NPhotons = other.NPhotons;
          PrimaryIdByLight = other.PrimaryIdByLight;
          VertexGlobalIdByLight = other.VertexGlobalIdByLight;
          LightShare = other.LightShare;
          FirstPhotonPrimaryId = other.FirstPhotonPrimaryId;
          FirstPhotonVertexGlobalId = other.FirstPhotonVertexGlobalId;
        }
        return *this;
    }
    
    // Move assignment operator
    TMS_TrueHit& operator=(TMS_TrueHit&& other) noexcept {
        if (this != &other) {
          PrimaryIds = std::move(other.PrimaryIds);
          VertexGlobalIds = std::move(other.VertexGlobalIds);
          EnergyShare = std::move(other.EnergyShare);
          EnergyShareIsLeptonic = std::move(other.EnergyShareIsLeptonic);
          LightPrimaryIds = std::move(other.LightPrimaryIds);
          LightVertexGlobalIds = std::move(other.LightVertexGlobalIds);
          LightPhotonCounts = std::move(other.LightPhotonCounts);
          x = other.x;
          y = other.y;
          z = other.z;
          t = other.t;
          dx = other.dx;
          EnergyDeposit = other.EnergyDeposit;
          pe = other.pe;
          peAfterFibers = other.peAfterFibers;
          peAfterFibersLongPath = other.peAfterFibersLongPath;
          peAfterFibersShortPath = other.peAfterFibersShortPath;
          NPhotons = other.NPhotons;
          PrimaryIdByLight = other.PrimaryIdByLight;
          VertexGlobalIdByLight = other.VertexGlobalIdByLight;
          LightShare = other.LightShare;
          FirstPhotonPrimaryId = other.FirstPhotonPrimaryId;
          FirstPhotonVertexGlobalId = other.FirstPhotonVertexGlobalId;
        }
        return *this;
    }

    double GetX() const {return x;};
    double GetY() const {return y;};
    double GetZ() const {return z;};
    double GetT() const {return t;};
    double GetdX() const {return dx;};
    double GetE() const {return EnergyDeposit; };
    double GetdEdx() const {return GetE()/GetdX(); };
    double GetPE() const {return pe; };
    double GetPEAfterFibers() const {return peAfterFibers; };
    double GetPEAfterFibersLongPath() const {return peAfterFibersLongPath; };
    double GetPEAfterFibersShortPath() const {return peAfterFibersShortPath; };
    
    int GetPrimaryId() const { return PrimaryIds.at(0); };
    int GetPrimaryIds(int index) const { return PrimaryIds.at(index); };
    long long GetVertexGlobalIds(int index) const { return VertexGlobalIds.at(index); };
    double GetEnergyShare(int index) const { return EnergyShare.at(index); };
    // GetPrimaryId()/GetVertexGlobalIds(0) return whichever contributor happened to be
    // pushed first (construction order, or first-in-merge-order after MergeWith()) --
    // that can differ between otherwise-identical runs whenever merge order differs, even
    // though the hit's own PE/energy is unaffected. These two instead return the highest
    // *energy-share* contributor, which is order-independent and matches the convention
    // already used for RecoHitPrimary* branches via TMS_Utils::GetPrimaryIdsByEnergy().
    // Ties on EnergyShare (exact ties are physically possible, e.g. two contributors of equal
    // deposited energy) break on (VertexGlobalIds, PrimaryIds) ascending -- both are intrinsic
    // identifiers of the contributor, not merge-order artifacts, so this stays order-independent.
    // Same tie preference GetSumAndHighest() gets for free from its std::map<pair<...>> key order.
    size_t IndexOfHighestEnergyContributor() const {
      size_t best = 0;
      for (size_t i = 1; i < EnergyShare.size(); i++) {
        if (EnergyShare[i] > EnergyShare[best] ||
            (EnergyShare[i] == EnergyShare[best] &&
             (VertexGlobalIds[i] < VertexGlobalIds[best] ||
              (VertexGlobalIds[i] == VertexGlobalIds[best] && PrimaryIds[i] < PrimaryIds[best])))) {
          best = i;
        }
      }
      return best;
    };
    // Light provenance of the readout this true hit belongs to (response-element pipeline only,
    // filled by TMS_Event::FillLightProvenance(); -999 otherwise): number of detected photons,
    // the particle (trajectory, vertex) that produced the most of them and its fraction, and the
    // particle that produced the first one (which sets the hit time).
    void SetLightProvenance(int n_photons, int primary_id, long long vertex_id, double share,
                            int first_primary_id, long long first_vertex_id) {
      NPhotons = n_photons; PrimaryIdByLight = primary_id; VertexGlobalIdByLight = vertex_id;
      LightShare = share; FirstPhotonPrimaryId = first_primary_id; FirstPhotonVertexGlobalId = first_vertex_id;
    };
    // The full light breakdown of the same readout: one entry per (trajectory, vertex) that
    // produced at least one detected photon, ascending in (vertex, trajectory), with its photon
    // count (fraction = count / NPhotons). Empty unless FillLightProvenance() ran.
    void SetLightContributions(std::vector<int> primary_ids, std::vector<long long> vertex_ids, std::vector<int> photon_counts) {
      LightPrimaryIds = std::move(primary_ids); LightVertexGlobalIds = std::move(vertex_ids); LightPhotonCounts = std::move(photon_counts);
    };
    size_t GetNLightContributions() const { return LightPhotonCounts.size(); };
    int GetLightContributionPrimaryId(size_t i) const { return LightPrimaryIds.at(i); };
    long long GetLightContributionVertexGlobalId(size_t i) const { return LightVertexGlobalIds.at(i); };
    int GetLightContributionPhotons(size_t i) const { return LightPhotonCounts.at(i); };
    int GetNPhotons() const { return NPhotons; };
    int GetPrimaryIdByLight() const { return PrimaryIdByLight; };
    long long GetVertexGlobalIdByLight() const { return VertexGlobalIdByLight; };
    double GetLightShare() const { return LightShare; };
    int GetFirstPhotonPrimaryId() const { return FirstPhotonPrimaryId; };
    long long GetFirstPhotonVertexGlobalId() const { return FirstPhotonVertexGlobalId; };
    int GetPrimaryIdByEnergy() const { return PrimaryIds.at(IndexOfHighestEnergyContributor()); };
    long long GetVertexGlobalIdByEnergy() const { return VertexGlobalIds.at(IndexOfHighestEnergyContributor()); };
    double GetEnergySharePortion(int index) const { return EnergyShare.at(index) / GetE(); };
    //void SetVertexId(int id) { VertexId = id; };
    size_t GetNTrueParticles() const { return EnergyShare.size(); };
    double GetLeptonicEnergy() const;
    double GetHadronicEnergy() const { return GetE() - GetLeptonicEnergy(); };

    void SetX(double pos) {x = pos;};
    void SetY(double pos) {y = pos;};
    void SetZ(double pos) {z = pos;};
    void SetT(double pos) {t = pos;};
    void SetdX(double dX) {dx = dX;};
    void SetE(double E) {EnergyDeposit = E;};
    void SetPE(double PE) {pe = PE;};
    void SetPEAfterFibers(double PE) {peAfterFibers = PE;};
    void SetPEAfterFibersLongPath(double PE) {peAfterFibersLongPath = PE;};
    void SetPEAfterFibersShortPath(double PE) {peAfterFibersShortPath = PE;};
    void SetEnergyLeptonic(int index, bool value = true) { EnergyShareIsLeptonic[index] = value; };

    void Print() const;
    
    void MergeWith(TMS_TrueHit& hit);

  private:
    double x;
    double y;
    double z;
    double t;
    double dx;
    double EnergyDeposit;
    double pe;
    double peAfterFibers;
    double peAfterFibersLongPath;
    double peAfterFibersShortPath;

    // See SetLightProvenance()
    int NPhotons = -999;
    int PrimaryIdByLight = -999;
    long long VertexGlobalIdByLight = -999;
    double LightShare = -999;
    int FirstPhotonPrimaryId = -999;
    long long FirstPhotonVertexGlobalId = -999;
    // See SetLightContributions()
    std::vector<int> LightPrimaryIds;
    std::vector<long long> LightVertexGlobalIds;
    std::vector<int> LightPhotonCounts;
    
    // Store individual particles for later particle identication
    std::vector<int> PrimaryIds;
    std::vector<long long> VertexGlobalIds;
    std::vector<double> EnergyShare;
    std::vector<bool> EnergyShareIsLeptonic;
};

#endif
