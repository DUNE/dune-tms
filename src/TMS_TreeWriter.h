#ifndef __TMS_TREEWRITER_H__
#define __TMS_TREEWRITER_H__
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include "TMS_Manager.h"
#include "TMS_Reco.h"

#define __TMS_MAX_LINES__ 20

// Just a simple tree writer for the output tree
class TMS_TreeWriter {

  public:
    static TMS_TreeWriter& GetWriter() {
      static TMS_TreeWriter Instance;
      return Instance;
    }


    void Fill(int i); // Fill the variables

  private:
    TMS_TreeWriter();
    TMS_TreeWriter(TMS_TreeWriter const &) = delete;
    void operator=(TMS_TreeWriter const &) = delete;
    ~TMS_TreeWriter() {
      Output->cd();
      Branch_Lines->Write();
      std::cout << "TMS_TreeWriter wrote output to " << Output->GetName() << std::endl;
      Branch_Lines->Delete();
      Output->Close();
    }

    TFile *Output; // The output TFile
    TTree *Branch_Lines; // The TTree

    void MakeBranches(); // Make the output branches

    // The variables
    int EventNo;
    double Slope[__TMS_MAX_LINES__];
    double Intercept[__TMS_MAX_LINES__];
    double DirectionZ[__TMS_MAX_LINES__];
    double DirectionX[__TMS_MAX_LINES__];
    int nLines;
    double FirstHit[2];
    int FirstPlane;
    bool TMSStart;
};


#endif
