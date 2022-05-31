#ifndef __TMS_TREEWRITER_H__
#define __TMS_TREEWRITER_H__
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include "TMS_Manager.h"
#include "TMS_Reco.h"
#include "TMS_Event.h"

// Only hard-coded for constant ROOT strucutre
// Could (and probably should) be replaced by vectors
#define __TMS_MAX_LINES__ 10 // Maximum number of lines in an event
#define __TMS_MAX_HITS__ 2000 // Maximum number of hits in an event
#define __TMS_MAX_LINE_HITS__ 200 // Maximum number of hits in a track
#define __TMS_MAX_CLUSTERS__ 50 // Maximum number of clusters in an event
#define __TMS_AUTOSAVE__ 1000 // Auto save to root file

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
    float Slope[__TMS_MAX_LINES__];
    float Intercept[__TMS_MAX_LINES__];
    float DirectionZ[__TMS_MAX_LINES__];
    float DirectionX[__TMS_MAX_LINES__];
    int nLines;
    float FirstHit[__TMS_MAX_LINES__][2]; // [0] is Z, [1] is NotZ
    float LastHit[__TMS_MAX_LINES__][2]; // [0] is Z, [1] is NotZ
    int FirstPlane[__TMS_MAX_LINES__];
    int LastPlane[__TMS_MAX_LINES__];
    bool TMSStart;
    float Occupancy[__TMS_MAX_LINES__];
    float TrackLength[__TMS_MAX_LINES__];
    float TotalTrackEnergy[__TMS_MAX_LINES__];
    float TrackHitEnergy[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__]; // Energy per track hit
    float TrackHitPos[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__][2]; // [0] is Z, [1] is NotZ
    int nHitsInTrack[__TMS_MAX_LINES__];

    // Cluster information
    int nClusters; // Number of clusters
    float ClusterEnergy[__TMS_MAX_CLUSTERS__]; // Energy in cluster
    int nHitsInCluster[__TMS_MAX_CLUSTERS__]; // Number of hits in cluster
    float ClusterPosMean[__TMS_MAX_CLUSTERS__][2]; // Mean cluster position, [0] is Z, [1] is NotZ
    float ClusterPosStdDev[__TMS_MAX_CLUSTERS__][2]; // Cluster standard deviation, [0] is Z, [1] is NotZ
    float ClusterHitPos[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__][2]; // Cluster hit position
    float ClusterHitEnergy[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__]; // Cluster hit energy

    // Reco information
    int nHits; // How many hits in event
    float RecoHitPos[__TMS_MAX_HITS__][4]; // Position of hit; [0] is x, [1] is y, [2] is z, [3] is time
    float RecoHitEnergy[__TMS_MAX_HITS__]; // Energy in hit

    // Truth information
    float MuonP4[4];
    float Muon_Vertex[4];
    float Muon_Death[4];
    float Muon_TrueKE;
    int nParticles;
    std::string Reaction;
    int NeutrinoPDG;
    float NeutrinoP4[4];
    float Muon_TrueTrackLength;
    bool IsCC;
};


#endif
