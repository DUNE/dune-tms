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
    TMS_TrueHit(TG4HitSegment &edep_seg, int vertex_id);
    TMS_TrueHit();

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
    
    int GetPrimaryId() const { return PrimaryId; };
    int GetVertexId() const { return VertexId; };
    void SetVertexId(int id) { VertexId = id; };

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
    int PrimaryId;
    int VertexId;
};

#endif
