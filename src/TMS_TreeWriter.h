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
    int nLinesU;
    int nLinesV;

    int nLines3D;
    
    int SliceNo;
    int SpillNo;
    
    int VertexIdOfMostEnergyInEvent;
    float VisibleEnergyFromUVertexInSlice;
    float TotalVisibleEnergyFromVertex;
    float VisibleEnergyFromVVerticesInSlice;
    float VertexVisibleEnergyFractionInSlice;
    float PrimaryVertexVisibleEnergyFraction;

    float SlopeU[__TMS_MAX_LINES__];
    float SlopeV[__TMS_MAX_LINES__];
    float InterceptU[__TMS_MAX_LINES__];
    float InterceptV[__TMS_MAX_LINES__];

    float Slope_DownstreamU[__TMS_MAX_LINES__];
    float Slope_DownstreamV[__TMS_MAX_LINES__];
    float Intercept_DownstreamU[__TMS_MAX_LINES__];
    float Intercept_DownstreamV[__TMS_MAX_LINES__];

    float Slope_UpstreamU[__TMS_MAX_LINES__];
    float Slope_UpstreamV[__TMS_MAX_LINES__];
    float Intercept_UpstreamU[__TMS_MAX_LINES__];
    float Intercept_UpstreamV[__TMS_MAX_LINES__];

    float DirectionZU[__TMS_MAX_LINES__];
    float DirectionZV[__TMS_MAX_LINES__];
    float DirectionXU[__TMS_MAX_LINES__];
    float DirectionXV[__TMS_MAX_LINES__];

    float DirectionZU_Upstream[__TMS_MAX_LINES__];
    float DirectionZV_Upstream[__TMS_MAX_LINES__];
    float DirectionXU_Upstream[__TMS_MAX_LINES__];
    float DirectionXV_Upstream[__TMS_MAX_LINES__];

    float DirectionZU_Downstream[__TMS_MAX_LINES__];
    float DirectionZV_Downstream[__TMS_MAX_LINES__];
    float DirectionXU_Downstream[__TMS_MAX_LINES__];
    float DirectionXV_Downstream[__TMS_MAX_LINES__];

    float FirstHitU[__TMS_MAX_LINES__][2]; // [0] is Z, [1] is NotZ
    float FirstHitV[__TMS_MAX_LINES__][2];
    float LastHitU[__TMS_MAX_LINES__][2]; // [0] is Z, [1] is NotZ
    float LastHitV[__TMS_MAX_LINES__][2];
    float FirstHitTimeU[__TMS_MAX_LINES__]; 
    float FirstHitTimeV[__TMS_MAX_LINES__];
    float LastHitTimeU[__TMS_MAX_LINES__];
    float LastHitTimeV[__TMS_MAX_LINES__];
    float EarliestHitTimeU[__TMS_MAX_LINES__]; 
    float EarliestHitTimeV[__TMS_MAX_LINES__];
    float LatestHitTimeU[__TMS_MAX_LINES__];
    float LatestHitTimeV[__TMS_MAX_LINES__];
    int FirstPlaneU[__TMS_MAX_LINES__];
    int FirstPlaneV[__TMS_MAX_LINES__];
    int LastPlaneU[__TMS_MAX_LINES__];
    int LastPlaneV[__TMS_MAX_LINES__];
    bool TMSStart;
    float TMSStartTime;
    float OccupancyU[__TMS_MAX_LINES__];
    float OccupancyV[__TMS_MAX_LINES__];
    float TrackLengthU[__TMS_MAX_LINES__];
    float TrackLengthV[__TMS_MAX_LINES__];
    float TotalTrackEnergyU[__TMS_MAX_LINES__];
    float TotalTrackEnergyV[__TMS_MAX_LINES__];
    bool TrackStoppingU[__TMS_MAX_LINES__];
    bool TrackStoppingV[__TMS_MAX_LINES__];

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

    float TrackHitEnergyU[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__]; // Energy per track hit
    float TrackHitEnergyV[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__];
    float TrackHitTimeU[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__];
    float TrackHitTimeV[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__];
    float TrackHitPosU[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__][2]; // [0] is Z, [1] is NotZ
    float TrackHitPosV[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__][2];
    int nHitsInTrackU[__TMS_MAX_LINES__];
    int nHitsInTrackV[__TMS_MAX_LINES__];

    // Cluster information
    int nClustersU; // Number of clusters
    int nClustersV;
    float ClusterEnergyU[__TMS_MAX_CLUSTERS__]; // Energy in cluster
    float ClusterEnergyV[__TMS_MAX_CLUSTERS__];
    float ClusterTimeU[__TMS_MAX_CLUSTERS__]; // Time in cluster
    float ClusterTimeV[__TMS_MAX_CLUSTERS__];
    int nHitsInClusterU[__TMS_MAX_CLUSTERS__]; // Number of hits in cluster
    int nHitsInClusterV[__TMS_MAX_CLUSTERS__];
    float ClusterPosMeanU[__TMS_MAX_CLUSTERS__][2]; // Mean cluster position, [0] is Z, [1] is NotZ
    float ClusterPosMeanV[__TMS_MAX_CLUSTERS__][2];
    float ClusterPosStdDevU[__TMS_MAX_CLUSTERS__][2]; // Cluster standard deviation, [0] is Z, [1] is NotZ
    float ClusterPosStdDevV[__TMS_MAX_CLUSTERS__][2];
    float ClusterHitPosU[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__][2]; // Cluster hit position
    float ClusterHitPosV[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__][2];
    float ClusterHitEnergyU[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__]; // Cluster hit energy
    float ClusterHitEnergyV[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__];
    float ClusterHitTimeU[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__]; // Cluster hit energy
    float ClusterHitTimeV[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__];
    int ClusterHitSliceU[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__]; // Cluster hit slice
    int ClusterHitSliceV[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__];

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
