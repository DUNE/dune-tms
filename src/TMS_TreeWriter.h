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
#define __TMS_MAX_TRACKS__ 100 // Maximum number of lines in an event
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
      Reco_Tree->Write();
      Truth_Info->Write();
      std::cout << "TMS_TreeWriter wrote output to " << Output->GetName() << std::endl;
      Output->Close();
    }

    // 3D Track Object Info
    int nTracks;
    int nHitsIn3DTrack[__TMS_MAX_TRACKS__];
    float RecoTrackHitPos[__TMS_MAX_TRACKS__][__TMS_MAX_LINE_HITS__][3]; // Due to a lack of variables, but as this is taken from line hits, it would make sense (maybe times 2?)
    float RecoTrackStartPos[__TMS_MAX_TRACKS__][3];
    float RecoTrackDirection[__TMS_MAX_TRACKS__][3];
    float RecoTrackEndPos[__TMS_MAX_TRACKS__][3];
    float RecoTrackEnergy[__TMS_MAX_TRACKS__];
    float RecoTrackEnergyDeposit[__TMS_MAX_TRACKS__];
    float RecoTrackLength[__TMS_MAX_TRACKS__];

  private:
    TMS_TreeWriter();
    TMS_TreeWriter(TMS_TreeWriter const &) = delete;
    void operator=(TMS_TreeWriter const &) = delete;
    ~TMS_TreeWriter() {};


    TFile* Output; // The output TFile
    TTree* Branch_Lines; // The TTree
    TTree* Reco_Tree; // The TTree 
    TTree* Truth_Info; // Truth info

    void Clear();
    void MakeBranches(); // Make the output branches

    // The variables
    int EventNo;
    int nLinesOne;
    int nLinesOther;

    int nLines3D;
    
    int SliceNo;
    int SpillNo;
    
    int VertexIdOfMostEnergyInEvent;
    float VisibleEnergyFromVertexInSlice;
    float TotalVisibleEnergyFromVertex;
    float VisibleEnergyFromOtherVerticesInSlice;
    float VertexVisibleEnergyFractionInSlice;
    float PrimaryVertexVisibleEnergyFraction;

    float SlopeOne[__TMS_MAX_LINES__];
    float SlopeOther[__TMS_MAX_LINES__];
    float InterceptOne[__TMS_MAX_LINES__];
    float InterceptOther[__TMS_MAX_LINES__];

    float Slope_DownstreamOne[__TMS_MAX_LINES__];
    float Slope_DownstreamOther[__TMS_MAX_LINES__];
    float Intercept_DownstreamOne[__TMS_MAX_LINES__];
    float Intercept_DownstreamOther[__TMS_MAX_LINES__];

    float Slope_UpstreamOne[__TMS_MAX_LINES__];
    float Slope_UpstreamOther[__TMS_MAX_LINES__];
    float Intercept_UpstreamOne[__TMS_MAX_LINES__];
    float Intercept_UpstreamOther[__TMS_MAX_LINES__];

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

    float FirstHit3D[__TMS_MAX_LINES__][3]; // [0] is Z, [1] is 'X', [2] is Y
    float LastHit3D[__TMS_MAX_LINES__][3];
    float FirstHitTime3D[__TMS_MAX_LINES__];
    float LastHitTime3D[__TMS_MAX_LINES__];
    float EarliestHitTime3D[__TMS_MAX_LINES__];
    float LatestHitTime3D[__TMS_MAX_LINES__];
    float FirstPlane3D[__TMS_MAX_LINES__];
    float LastPlane3D[__TMS_MAX_LINES__];
    float Occupancy3D[__TMS_MAX_LINES__];
    float TrackLength3D[__TMS_MAX_LINES__];
    float TotalTrackEnergy3D[__TMS_MAX_LINES__];
    float TrackStopping3D[__TMS_MAX_LINES__];
    float TrackHitEnergy3D[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__];
    float TrackHitTime3D[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__];
    float TrackHitPos3D[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__][3];
    float nHitsInTrack3D[__TMS_MAX_LINES__];

    float TrackHitEnergyOne[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__]; // Energy per track hit
    float TrackHitEnergyOther[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__];
    float TrackHitTimeOne[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__];
    float TrackHitTimeOther[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__];
    float TrackHitPosOne[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__][2]; // [0] is Z, [1] is NotZ
    float TrackHitPosOther[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__][2];
    int nHitsInTrackOne[__TMS_MAX_LINES__];
    int nHitsInTrackOther[__TMS_MAX_LINES__];

    // Cluster information
    int nClustersOne; // Number of clusters
    int nClustersOther;
    float ClusterEnergyOne[__TMS_MAX_CLUSTERS__]; // Energy in cluster
    float ClusterEnergyOther[__TMS_MAX_CLUSTERS__];
    float ClusterTimeOne[__TMS_MAX_CLUSTERS__]; // Time in cluster
    float ClusterTimeOther[__TMS_MAX_CLUSTERS__];
    int nHitsInClusterOne[__TMS_MAX_CLUSTERS__]; // Number of hits in cluster
    int nHitsInClusterOther[__TMS_MAX_CLUSTERS__];
    float ClusterPosMeanOne[__TMS_MAX_CLUSTERS__][2]; // Mean cluster position, [0] is Z, [1] is NotZ
    float ClusterPosMeanOther[__TMS_MAX_CLUSTERS__][2];
    float ClusterPosStdDevOne[__TMS_MAX_CLUSTERS__][2]; // Cluster standard deviation, [0] is Z, [1] is NotZ
    float ClusterPosStdDevOther[__TMS_MAX_CLUSTERS__][2];
    float ClusterHitPosOne[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__][2]; // Cluster hit position
    float ClusterHitPosOther[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__][2];
    float ClusterHitEnergyOne[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__]; // Cluster hit energy
    float ClusterHitEnergyOther[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__];
    float ClusterHitTimeOne[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__]; // Cluster hit energy
    float ClusterHitTimeOther[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__];
    int ClusterHitSliceOne[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__]; // Cluster hit slice
    int ClusterHitSliceOther[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__];

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
