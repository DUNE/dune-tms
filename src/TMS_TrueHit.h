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

    double GetX() const {return x;};
    double GetY() const {return y;};
    double GetZ() const {return z;};
    double GetT() const {return t;};
    double GetE() const {return EnergyDeposit; };
    
    int GetPrimaryId() const { return PrimaryId; };
    int GetVertexId() const { return VertexId; };
    void SetVertexId(int id) { VertexId = id; };

    void SetX(double pos) {x = pos;};
    void SetY(double pos) {y = pos;};
    void SetZ(double pos) {z = pos;};
    void SetT(double pos) {t = pos;};
    void SetE(double E) {EnergyDeposit = E;};

    void Print() const;

  private:
    double x;
    double y;
    double z;
    double t;
    double EnergyDeposit;
    int PrimaryId;
    int VertexId;
};

#endif
