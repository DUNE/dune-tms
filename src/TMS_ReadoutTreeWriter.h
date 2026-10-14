#ifndef __TMS_ReadoutTreeWriter_H__
#define __TMS_ReadoutTreeWriter_H__
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include "TMS_Manager.h"
#include "TMS_Reco.h"
#include "TMS_Event.h"

// Just a simple tree writer for the output tree
class TMS_ReadoutTreeWriter {

  public:
    static TMS_ReadoutTreeWriter& GetWriter() {
      static TMS_ReadoutTreeWriter Instance;
      return Instance;
    }

    void Fill(TMS_Event &event);

    void Write() {
      Output->cd();
      TMS_Readout->Write();
      std::cout << "TMS_ReadoutTreeWriter wrote output to " << Output->GetName() << std::endl;
      Output->Close();
    }

  private:
    TMS_ReadoutTreeWriter();
    TMS_ReadoutTreeWriter(TMS_ReadoutTreeWriter const &) = delete;
    void operator=(TMS_ReadoutTreeWriter const &) = delete;
    ~TMS_ReadoutTreeWriter() {};


    TFile *Output; // The output TFile
    TTree *TMS_Readout; 

    void Clear();
    void MakeBranches(); // Make the output branches
    
    // Want ability to think only about reco when dealing with real data
    bool hasTruth = true;
    // Per-hit RecoHitTrue* branches in the Readout file, set from [Truth] ReadoutHitTruth (default off). Needs hasTruth.
    bool writeHitTruth = false;
    
    // Ideally I'd use std::vector<Float_t> but they don't fill with the right info
    #define __MAX_READOUT_TREE_ARRAY_LENGTH__ 100000
    #define VAR(x) float x[__MAX_READOUT_TREE_ARRAY_LENGTH__]
    #define INTVAR(x) int x[__MAX_READOUT_TREE_ARRAY_LENGTH__]
    
    // True hit branches
    int NTrueHits;
    VAR(TrueHitX);
    VAR(TrueHitY);
    VAR(TrueHitZ);
    VAR(TrueHitT);
    VAR(TrueHitE);
    VAR(TrueHitPE);
    VAR(TrueHitPEAfterFibers);
    VAR(TrueHitPEAfterFibersLongPath);
    VAR(TrueHitPEAfterFibersShortPath);
    
    // Reco hit branches
    int NRecoHits;
    VAR(RecoHitX);
    VAR(RecoHitY);
    VAR(RecoHitZ);
    VAR(RecoHitNotZ);
    VAR(RecoHitT);
    VAR(RecoHitE);
    VAR(RecoHitPE);
    bool RecoHitIsPedSupped[__MAX_READOUT_TREE_ARRAY_LENGTH__];
    INTVAR(RecoHitBar);
    INTVAR(RecoHitPlane);
    // Ground truth per reco hit (same TMS_Utils::GetPrimaryIdsByEnergy lookup as the SpacePointHitTrue* branches in
    // TMS_TreeWriter), written for every hit of the spill, not only those in space points. Vertex id -1 / track id -999 =
    // no truth (noise hit). The first contributor is the one with the most energy, the second the next largest.
    long long RecoHitTrueVertexGlobalId[__MAX_READOUT_TREE_ARRAY_LENGTH__];
    INTVAR(RecoHitTrueTrackId);
    VAR(RecoHitTrueEnergyFrac);
    long long RecoHitTrue2VertexGlobalId[__MAX_READOUT_TREE_ARRAY_LENGTH__];
    INTVAR(RecoHitTrue2TrackId);
    VAR(RecoHitTrue2EnergyFrac);
    INTVAR(RecoHitNTrueParticles);
    #ifdef RECORD_HIT_DEADTIME
    VAR(RecoHitDeadtimeStart);
    VAR(RecoHitDeadtimeStop);
    #endif
};


#endif
