#ifndef _TMS_TRUEHIT_H_
#define _TMS_TRUEHIT_H_

#include <vector>
#include <iostream>

// Include the constants
#include "TMS_Constants.h"

#include "TG4HitSegment.h"

// Essentially a copy of the edep-sim THit
class TMS_TrueHit {
  public:
    TMS_TrueHit(TG4HitSegment &edep_seg, int vertex_id);
    
    TMS_TrueHit() = delete;
    
    // Explicitly create copy construct and assignment operators
    // Or else true particle information is lost 
    // Copy constructor
    TMS_TrueHit(const TMS_TrueHit& other) : PrimaryIds(other.PrimaryIds),
      VertexIds(other.VertexIds), EnergyShare(other.EnergyShare)
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
        }
    }
    
    // Copy assignment operator
    TMS_TrueHit& operator=(const TMS_TrueHit& other) {
        if (this != &other) {
          PrimaryIds = other.PrimaryIds;
          VertexIds = other.VertexIds;
          EnergyShare = other.EnergyShare;
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
        }
        return *this;
    }
    
    // Move assignment operator
    TMS_TrueHit& operator=(TMS_TrueHit&& other) noexcept {
        if (this != &other) {
          PrimaryIds = std::move(other.PrimaryIds);
          VertexIds = std::move(other.VertexIds);
          EnergyShare = std::move(other.EnergyShare);
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
        }
        return *this;
    }

    double GetX() const {return x;};
    double GetY() const {return y;};
    double GetZ() const {return z;};
    double GetT() const {return t;};
    double GetdX() const {return dx;};
    double GetE() const {return EnergyDeposit; };
    double GetPE() const {return pe; };
    double GetPEAfterFibers() const {return peAfterFibers; };
    double GetPEAfterFibersLongPath() const {return peAfterFibersLongPath; };
    double GetPEAfterFibersShortPath() const {return peAfterFibersShortPath; };
    
    int GetPrimaryId() const { return PrimaryIds.at(0); };
    int GetVertexId() const { 
      if (VertexIds.size() == 1) return VertexIds.at(0); 
      else if (VertexIds.size() == 0) return -999;
      std::cout<<"Fatal: Using GetVertexId with > 1 possible VertexId. Use GetVertexIds instead."<<std::endl;
      exit(1);
    };
    int GetPrimaryIds(int index) const { return PrimaryIds.at(index); };
    int GetVertexIds(int index) const { return VertexIds.at(index); };
    double GetEnergyShare(int index) const { return EnergyShare.at(index); };
    //void SetVertexId(int id) { VertexId = id; };
    size_t GetNTrueParticles() const { return EnergyShare.size(); };

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
    
    // Store individual particles for later particle identication
    std::vector<int> PrimaryIds;
    std::vector<int> VertexIds;
    std::vector<double> EnergyShare;
};

#endif
