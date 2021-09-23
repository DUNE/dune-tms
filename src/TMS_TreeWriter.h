#ifndef __TMS_TREEWRITER_H__
#define __TMS_TREEWRITER_H__
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include "TMS_Manager.h"
#include "TMS_Reco.h"
#include "TMS_Event.h"

#define __TMS_MAX_LINES__ 20

#define __TMS_AUTOSAVE__ 1000

// Just a simple tree writer for the output tree
class TMS_TreeWriter {

  public:
    static TMS_TreeWriter& GetWriter() {
      static TMS_TreeWriter Instance;
      return Instance;
    }

    void Fill(TMS_Event &event);

    void Write() {
      Output->cd();
      Branch_Lines->Write();
      Truth_Info->Write();
      std::cout << "TMS_TreeWriter wrote output to " << Output->GetName() << std::endl;
      Output->Close();
    }

  private:
    TMS_TreeWriter();
    TMS_TreeWriter(TMS_TreeWriter const &) = delete;
    void operator=(TMS_TreeWriter const &) = delete;
    ~TMS_TreeWriter() {};


    TFile *Output; // The output TFile
    TTree *Branch_Lines; // The TTree
    TTree *Truth_Info; // Truth info

    void Clear();
    void MakeBranches(); // Make the output branches

    // The variables
    int EventNo;
    double Slope[__TMS_MAX_LINES__];
    double Intercept[__TMS_MAX_LINES__];
    double DirectionZ[__TMS_MAX_LINES__];
    double DirectionX[__TMS_MAX_LINES__];
    int nLines;
    double FirstHit[__TMS_MAX_LINES__][2];
    double LastHit[__TMS_MAX_LINES__][2];
    int FirstPlane[__TMS_MAX_LINES__];
    int LastPlane[__TMS_MAX_LINES__];
    bool TMSStart;
    double Occupancy[__TMS_MAX_LINES__];
    double TrackLength[__TMS_MAX_LINES__];

    double MuonP4[4];
    double Muon_Vertex[4];
    double Muon_TrueKE;
    int nParticles;
    std::string Reaction;
    int NeutrinoPDG;
    double NeutrinoP4[4];
};


#endif
