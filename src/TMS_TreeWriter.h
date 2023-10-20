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
#define __TMS_MAX_LINES__ 100 // Maximum number of lines in an event
#define __TMS_MAX_HITS__ 20000 // Maximum number of hits in an event
#define __TMS_MAX_LINE_HITS__ 200 // Maximum number of hits in a track
#define __TMS_MAX_CLUSTERS__ 500 // Maximum number of clusters in an event
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
    int nLines;
    
    int SliceNo;
    int SpillNo;
    
    int VertexIdOfMostEnergyInEvent;
    float VisibleEnergyFromVertexInSlice;
    float TotalVisibleEnergyFromVertex;
    float VisibleEnergyFromOtherVerticesInSlice;
    float VertexVisibleEnergyFractionInSlice;
    float PrimaryVertexVisibleEnergyFraction;

    float Slope[__TMS_MAX_LINES__];
    float Intercept[__TMS_MAX_LINES__];

    float Slope_Downstream[__TMS_MAX_LINES__];
    float Intercept_Downstream[__TMS_MAX_LINES__];

    float Slope_Upstream[__TMS_MAX_LINES__];
    float Intercept_Upstream[__TMS_MAX_LINES__];

    float DirectionZOne[__TMS_MAX_LINES__];
    float DirectionZOther[__TMS_MAX_LINES__];
    float DirectionXOne[__TMS_MAX_LINES__];
    float DirectionXOther[__TMS_MAX_LINES__];

    float DirectionZOne_Upstream[__TMS_MAX_LINES__];
    float DirectionZOther_Upstream[__TMS_MAX_LINES__];
    float DirectionXOne_Upstream[__TMS_MAX_LINES__];
    float DirectionXOther_Upstream[__TMS_MAX_LINES__];

    float DirectionZOne_Downstream[__TMS_MAX_LINES__];
    float DirectionZOther_Downstream[__TMS_MAX_LINES__];
    float DirectionXOne_Downstream[__TMS_MAX_LINES__];
    float DirectionXOther_Downstream[__TMS_MAX_LINES__];

    float FirstHitOne[__TMS_MAX_LINES__][2]; // [0] is Z, [1] is NotZ
    float FirstHitOther[__TMS_MAX_LINES__][2];
    float LastHitOne[__TMS_MAX_LINES__][2]; // [0] is Z, [1] is NotZ
    float LastHitOther[__TMS_MAX_LINES__][2];
    float FirstHitTimeOne[__TMS_MAX_LINES__]; 
    float FirstHitTimeOther[__TMS_MAX_LINES__];
    float LastHitTimeOne[__TMS_MAX_LINES__];
    float LastHitTimeOther[__TMS_MAX_LINES__];
    float EarliestHitTimeOne[__TMS_MAX_LINES__]; 
    float EarliestHitTimeOther[__TMS_MAX_LINES__];
    float LatestHitTimeOne[__TMS_MAX_LINES__];
    float LatestHitTimeOther[__TMS_MAX_LINES__];
    int FirstPlaneOne[__TMS_MAX_LINES__];
    int FirstPlaneOther[__TMS_MAX_LINES__];
    int LastPlaneOne[__TMS_MAX_LINES__];
    int LastPlaneOther[__TMS_MAX_LINES__];
    bool TMSStart;
    float TMSStartTime;
    float OccupancyOne[__TMS_MAX_LINES__];
    float OccupancyOther[__TMS_MAX_LINES__];
    float TrackLengthOne[__TMS_MAX_LINES__];
    float TrackLengthOther[__TMS_MAX_LINES__];
    float TotalTrackEnergyOne[__TMS_MAX_LINES__];
    float TotalTrackEnergyOther[__TMS_MAX_LINES__];
    bool TrackStoppingOne[__TMS_MAX_LINES__];
    bool TrackStoppingOther[__TMS_MAX_LINES__];

    float TrackHitEnergyOne[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__]; // Energy per track hit
    float TrackHitEnergyOther[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__];
    float TrackHitTimeOne[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__];
    float TrackHitTimeOther[__TMS_MAX_LINES][__TMS_MAX_LINE_HITS__];
    float TrackHitPosOne[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__][2]; // [0] is Z, [1] is NotZ
    float TrackHitPosOther[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__][2];
    int nHitsInTrack[__TMS_MAX_LINES__];

    // Cluster information
    int nClusters; // Number of clusters
    float ClusterEnergy[__TMS_MAX_CLUSTERS__]; // Energy in cluster
    float ClusterTime[__TMS_MAX_CLUSTERS__]; // Energy in cluster
    int nHitsInCluster[__TMS_MAX_CLUSTERS__]; // Number of hits in cluster
    float ClusterPosMean[__TMS_MAX_CLUSTERS__][2]; // Mean cluster position, [0] is Z, [1] is NotZ
    float ClusterPosStdDev[__TMS_MAX_CLUSTERS__][2]; // Cluster standard deviation, [0] is Z, [1] is NotZ
    float ClusterHitPos[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__][2]; // Cluster hit position
    float ClusterHitEnergy[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__]; // Cluster hit energy
    float ClusterHitTime[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__]; // Cluster hit energy
    int ClusterHitSlice[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__]; // Cluster hit slice

    // Reco information
    int nHits; // How many hits in event
    float RecoHitPos[__TMS_MAX_HITS__][4]; // Position of hit; [0] is x, [1] is y, [2] is z, [3] is time
    float RecoHitEnergy[__TMS_MAX_HITS__]; // Energy in hit
    int RecoHitSlice[__TMS_MAX_HITS__];

    // Truth information
    float MuonP4[4];
    float Muon_Vertex[4];
    float Muon_Death[4];
    float Muon_TrueKE;
    int nParticles;
    std::string Reaction;
    int NeutrinoPDG;
    float NeutrinoP4[4];
    float NeutrinoX4[4];
    int LeptonPDG;
    float LeptonP4[4];
    float LeptonX4[4];
    float Muon_TrueTrackLength;
    bool IsCC;
};


#endif
