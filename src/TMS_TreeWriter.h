#ifndef __TMS_TREEWRITER_H__
#define __TMS_TREEWRITER_H__
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include "TMS_Manager.h"
#include "TMS_Reco.h"

// Just a simple tree writer for the output tree
class TMS_TreeWriter {

  public:
    static TMS_TreeWriter& GetWriter() {
      static TMS_TreeWriter Instance;
      return Instance;
    }

    // Destructor writes the TTree to file and closes the file
    ~TMS_TreeWriter() {
      Output->cd();
      Branch_Lines->Write();
      std::cout << "TMS_TreeWriter wrote output to " << Output->GetName() << std::endl;
      Output->Close();
    }

    void Fill(int i); // Fill the variables

  private:
    TMS_TreeWriter();
    TMS_TreeWriter(TMS_TreeWriter const &) = delete;
    void operator=(TMS_TreeWriter const &) = delete;

    TFile *Output; // The output TFile
    TTree *Branch_Lines; // The TTree

    void MakeBranches(); // Make the output branches

    // The variables
    int EventNo;
    double Slope[20];
    double Intercept[20];
    double DirectionZ[20];
    double DirectionX[20];
    int nLines;
    double FirstHit[2];
    int FirstPlane;
    bool TMSStart;

};


#endif
