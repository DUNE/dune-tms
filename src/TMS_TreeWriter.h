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
#define __TMS_MAX_TRUE_PARTICLES__ 20000 // Maximum number of true particles to save info about
#define __TMS_MAX_TRUE_NONTMS_HITS__ 100000 // Maximum number of true particles to save info about

// Just a simple tree writer for the output tree
class TMS_TreeWriter {

  public:
    static TMS_TreeWriter& GetWriter() {
      static TMS_TreeWriter Instance;
      return Instance;
    }

    void Fill(TMS_Event &event);
    void FillSpill(TMS_Event &event, int truth_info_entry_number, int truth_info_n_slices);
    void FillTruthInfo(TMS_Event &event);

    void Write() {
      Output->cd();
      Branch_Lines->Write();
      Reco_Tree->Write();
      Truth_Info->Write();
      Truth_Spill->Write();
      std::cout << "TMS_TreeWriter wrote output to " << Output->GetName() << std::endl;
      Output->Close();
    }

    // 3D Track Object Info
    int nTracks;
    int nHitsIn3DTrack[__TMS_MAX_TRACKS__];
    int nKalmanNodes[__TMS_MAX_TRACKS__];
    float RecoTrackKalmanPos[__TMS_MAX_TRACKS__][__TMS_MAX_LINE_HITS__][3];
    int RecoTrackKalmanFirstPlaneBarView[__TMS_MAX_TRACKS__][3];
    int RecoTrackKalmanLastPlaneBarView[__TMS_MAX_TRACKS__][3];
    int RecoTrackKalmanPlaneBarView[__TMS_MAX_TRACKS__][__TMS_MAX_LINE_HITS__][3];
    float RecoTrackKalmanTruePos[__TMS_MAX_TRACKS__][__TMS_MAX_LINE_HITS__][3];
    int RecoTrackKalmanFirstPlaneBarViewTrue[__TMS_MAX_TRACKS__][3];
    int RecoTrackKalmanLastPlaneBarViewTrue[__TMS_MAX_TRACKS__][3];
    int RecoTrackKalmanPlaneBarViewTrue[__TMS_MAX_TRACKS__][__TMS_MAX_LINE_HITS__][3];
    float RecoTrackHitPos[__TMS_MAX_TRACKS__][__TMS_MAX_LINE_HITS__][3]; // Due to a lack of variables, but as this is taken from line hits, it would make sense (maybe times 2?)
    float RecoTrackHitEnergies[__TMS_MAX_TRACKS__][__TMS_MAX_LINE_HITS__]; // Due to a lack of variables, but as this is taken from line hits, it would make sense (maybe times 2?)
    float RecoTrackStartPos[__TMS_MAX_TRACKS__][3];
    float RecoTrackStartDirection[__TMS_MAX_TRACKS__][3];
    float RecoTrackEndPos[__TMS_MAX_TRACKS__][3];
    float RecoTrackEndDirection[__TMS_MAX_TRACKS__][3];
    float RecoTrackEnergyRange[__TMS_MAX_TRACKS__];
    float RecoTrackEnergyDeposit[__TMS_MAX_TRACKS__];
    float RecoTrackMomentum[__TMS_MAX_TRACKS__];
    float RecoTrackTrueMomentum[__TMS_MAX_TRACKS__];
    float RecoTrackLength[__TMS_MAX_TRACKS__];
    float RecoTrackChi2[__TMS_MAX_TRACKS__];
    int RecoTrackCharge[__TMS_MAX_TRACKS__];
    
  private:
    TMS_TreeWriter();
    TMS_TreeWriter(TMS_TreeWriter const &) = delete;
    void operator=(TMS_TreeWriter const &) = delete;
    ~TMS_TreeWriter() {};


    TFile* Output; // The output TFile
    TTree* Branch_Lines; // The TTree
    TTree* Reco_Tree; // The TTree 
    TTree* Truth_Info; // Truth info
    TTree* Truth_Spill; // Truth spill

    void Clear();
    void MakeBranches(); // Make the output branches
    void MakeTruthBranches(TTree* truth); // Make the output branches
    
    float TimeSliceStartTime;
    float TimeSliceEndTime;

    // The variables
    int EventNo;
    int RunNo;
    int nLinesU;
    int nLinesV;
    int nLinesX;

    int nLines3D;
    
    int SliceNo;
    int SpillNo;
    
    int VertexIdOfMostEnergyInEvent;
    float VisibleEnergyFromVertexInSlice;
    float TotalVisibleEnergyFromVertex;
    float VisibleEnergyFromOtherVerticesInSlice;
    float VertexVisibleEnergyFractionInSlice;
    float PrimaryVertexVisibleEnergyFraction;

    float LArOuterShellEnergy;
    float LArOuterShellEnergyFromVertex;
    float LArTotalEnergy;
    float LArTotalEnergyFromVertex;
    float TotalNonTMSEnergy;
    float TotalNonTMSEnergyFromVertex;

    float SlopeU[__TMS_MAX_LINES__];
    float SlopeV[__TMS_MAX_LINES__];
    float SlopeX[__TMS_MAX_LINES__];
    float InterceptU[__TMS_MAX_LINES__];
    float InterceptV[__TMS_MAX_LINES__];
    float InterceptX[__TMS_MAX_LINES__];

    float Slope_DownstreamU[__TMS_MAX_LINES__];
    float Slope_DownstreamV[__TMS_MAX_LINES__];
    float Slope_DownstreamX[__TMS_MAX_LINES__];
    float Intercept_DownstreamU[__TMS_MAX_LINES__];
    float Intercept_DownstreamV[__TMS_MAX_LINES__];
    float Intercept_DownstreamX[__TMS_MAX_LINES__];

    float Slope_UpstreamU[__TMS_MAX_LINES__];
    float Slope_UpstreamV[__TMS_MAX_LINES__];
    float Slope_UpstreamX[__TMS_MAX_LINES__];
    float Intercept_UpstreamU[__TMS_MAX_LINES__];
    float Intercept_UpstreamV[__TMS_MAX_LINES__];
    float Intercept_UpstreamX[__TMS_MAX_LINES__];

    float DirectionZU[__TMS_MAX_LINES__];
    float DirectionZV[__TMS_MAX_LINES__];
    float DirectionZX[__TMS_MAX_LINES__];
    float DirectionXU[__TMS_MAX_LINES__];
    float DirectionXV[__TMS_MAX_LINES__];
    float DirectionYX[__TMS_MAX_LINES__];

    float DirectionZU_Upstream[__TMS_MAX_LINES__];
    float DirectionZV_Upstream[__TMS_MAX_LINES__];
    float DirectionZX_Upstream[__TMS_MAX_LINES__];
    float DirectionXU_Upstream[__TMS_MAX_LINES__];
    float DirectionXV_Upstream[__TMS_MAX_LINES__];
    float DirectionYX_Upstream[__TMS_MAX_LINES__];

    float DirectionZU_Downstream[__TMS_MAX_LINES__];
    float DirectionZV_Downstream[__TMS_MAX_LINES__];
    float DirectionZX_Downstream[__TMS_MAX_LINES__];
    float DirectionXU_Downstream[__TMS_MAX_LINES__];
    float DirectionXV_Downstream[__TMS_MAX_LINES__];
    float DirectionYX_Downstream[__TMS_MAX_LINES__];

    float FirstHitU[__TMS_MAX_LINES__][2]; // [0] is Z, [1] is NotZ
    float FirstHitV[__TMS_MAX_LINES__][2];
    float FirstHitX[__TMS_MAX_LINES__][2];
    float LastHitU[__TMS_MAX_LINES__][2]; // [0] is Z, [1] is NotZ
    float LastHitV[__TMS_MAX_LINES__][2];
    float LastHitX[__TMS_MAX_LINES__][2];
    float FirstHitTimeU[__TMS_MAX_LINES__]; 
    float FirstHitTimeV[__TMS_MAX_LINES__];
    float FirstHitTimeX[__TMS_MAX_LINES__];
    float LastHitTimeU[__TMS_MAX_LINES__];
    float LastHitTimeV[__TMS_MAX_LINES__];
    float LastHitTimeX[__TMS_MAX_LINES__];
    float EarliestHitTimeU[__TMS_MAX_LINES__]; 
    float EarliestHitTimeV[__TMS_MAX_LINES__];
    float EarliestHitTimeX[__TMS_MAX_LINES__];
    float LatestHitTimeU[__TMS_MAX_LINES__];
    float LatestHitTimeV[__TMS_MAX_LINES__];
    float LatestHitTimeX[__TMS_MAX_LINES__];
    int FirstPlaneU[__TMS_MAX_LINES__];
    int FirstPlaneV[__TMS_MAX_LINES__];
    int FirstPlaneX[__TMS_MAX_LINES__];
    int LastPlaneU[__TMS_MAX_LINES__];
    int LastPlaneV[__TMS_MAX_LINES__];
    int LastPlaneX[__TMS_MAX_LINES__];
    bool TMSStart;
    float TMSStartTime;
    float OccupancyU[__TMS_MAX_LINES__];
    float OccupancyV[__TMS_MAX_LINES__];
    float OccupancyX[__TMS_MAX_LINES__];
    float TrackLengthU[__TMS_MAX_LINES__];
    float TrackLengthV[__TMS_MAX_LINES__];
    float TrackLengthX[__TMS_MAX_LINES__];
    float TotalTrackEnergyU[__TMS_MAX_LINES__];
    float TotalTrackEnergyV[__TMS_MAX_LINES__];
    float TotalTrackEnergyX[__TMS_MAX_LINES__];
    bool TrackStoppingU[__TMS_MAX_LINES__];
    bool TrackStoppingV[__TMS_MAX_LINES__];
    bool TrackStoppingX[__TMS_MAX_LINES__];

    /*float FirstHit3D[__TMS_MAX_LINES__][3]; // [0] is Z, [1] is 'X', [2] is Y
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
    float nHitsInTrack3D[__TMS_MAX_LINES__];*/

    float TrackHitEnergyU[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__]; // Energy per track hit
    float TrackHitEnergyV[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__];
    float TrackHitEnergyX[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__];
    float TrackHitTimeU[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__];
    float TrackHitTimeV[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__];
    float TrackHitTimeX[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__];
    float TrackHitPosU[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__][2]; // [0] is Z, [1] is NotZ
    float TrackHitPosV[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__][2];
    float TrackHitPosX[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__][2];
    int nHitsInTrackU[__TMS_MAX_LINES__];
    int nHitsInTrackV[__TMS_MAX_LINES__];
    int nHitsInTrackX[__TMS_MAX_LINES__];

    // Cluster information
    int nClustersU; // Number of clusters
    int nClustersV;
    int nClustersX;
    float ClusterEnergyU[__TMS_MAX_CLUSTERS__]; // Energy in cluster
    float ClusterEnergyV[__TMS_MAX_CLUSTERS__];
    float ClusterEnergyX[__TMS_MAX_CLUSTERS__];

    float ClusterTimeU[__TMS_MAX_CLUSTERS__]; // Time in cluster
    float ClusterTimeV[__TMS_MAX_CLUSTERS__];
    float ClusterTimeX[__TMS_MAX_CLUSTERS__];
    int nHitsInClusterU[__TMS_MAX_CLUSTERS__]; // Number of hits in cluster
    int nHitsInClusterV[__TMS_MAX_CLUSTERS__];
    int nHitsInClusterX[__TMS_MAX_CLUSTERS__];
    float ClusterPosMeanU[__TMS_MAX_CLUSTERS__][2]; // Mean cluster position, [0] is Z, [1] is NotZ
    float ClusterPosMeanV[__TMS_MAX_CLUSTERS__][2];
    float ClusterPosMeanX[__TMS_MAX_CLUSTERS__][2];
    float ClusterPosStdDevU[__TMS_MAX_CLUSTERS__][2]; // Cluster standard deviation, [0] is Z, [1] is NotZ
    float ClusterPosStdDevV[__TMS_MAX_CLUSTERS__][2];
    float ClusterPosStdDevX[__TMS_MAX_CLUSTERS__][2];
    float ClusterHitPosU[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__][2]; // Cluster hit position
    float ClusterHitPosV[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__][2];
    float ClusterHitPosX[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__][2];
    float ClusterHitEnergyU[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__]; // Cluster hit energy
    float ClusterHitEnergyV[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__];
    float ClusterHitEnergyX[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__];
    float ClusterHitTimeU[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__]; // Cluster hit energy
    float ClusterHitTimeV[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__];
    float ClusterHitTimeX[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__];
    int ClusterHitSliceU[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__]; // Cluster hit slice
    int ClusterHitSliceV[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__];
    int ClusterHitSliceX[__TMS_MAX_CLUSTERS__][__TMS_MAX_LINE_HITS__];

    // Reco information
    int nHits; // How many hits in event
    float RecoHitPos[__TMS_MAX_HITS__][4]; // Position of hit; [0] is x, [1] is y, [2] is z, [3] is time
    float RecoHitEnergy[__TMS_MAX_HITS__]; // Energy in hit
    float RecoHitPE[__TMS_MAX_HITS__];
    int RecoHitBar[__TMS_MAX_HITS__];
    int RecoHitPlane[__TMS_MAX_HITS__];
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
    
    int nTrueParticles;
    int nTruePrimaryParticles;
    int nTrueForgottenParticles;
    
    int VertexID[__TMS_MAX_TRUE_PARTICLES__];
    int Parent[__TMS_MAX_TRUE_PARTICLES__];
    int TrackId[__TMS_MAX_TRUE_PARTICLES__];
    int PDG[__TMS_MAX_TRUE_PARTICLES__];
    bool IsPrimary[__TMS_MAX_TRUE_PARTICLES__];
    float TrueVisibleEnergy[__TMS_MAX_TRUE_PARTICLES__];
    int TrueNHits[__TMS_MAX_TRUE_PARTICLES__];
    float TrueVisibleEnergyInSlice[__TMS_MAX_TRUE_PARTICLES__];
    int TrueNHitsInSlice[__TMS_MAX_TRUE_PARTICLES__];
    float TruePathLength[__TMS_MAX_TRUE_PARTICLES__];
    float TruePathLengthIgnoreY[__TMS_MAX_TRUE_PARTICLES__];
    float TruePathLengthInTMS[__TMS_MAX_TRUE_PARTICLES__];
    float TruePathLengthInTMSIgnoreY[__TMS_MAX_TRUE_PARTICLES__];

    int TrueNonTMSNHits;
    float TrueNonTMSHitPos[__TMS_MAX_TRUE_NONTMS_HITS__][4];
    float TrueNonTMSHitEnergy[__TMS_MAX_TRUE_NONTMS_HITS__];
    float TrueNonTMSHitHadronicEnergy[__TMS_MAX_TRUE_NONTMS_HITS__];
    float TrueNonTMSHitDx[__TMS_MAX_TRUE_NONTMS_HITS__];
    float TrueNonTMSHitdEdx[__TMS_MAX_TRUE_NONTMS_HITS__];
    int TrueNonTMSHitVertexID[__TMS_MAX_TRUE_NONTMS_HITS__];

    // Flags for easy use
    bool InteractionTMSFiducial;
    bool InteractionTMSFirstTwoModules;
    bool InteractionTMSThin;
    bool InteractionLArFiducial;
    bool TMSFiducialStart[__TMS_MAX_TRUE_PARTICLES__];
    bool TMSFiducialTouch[__TMS_MAX_TRUE_PARTICLES__];
    bool TMSFiducialEnd[__TMS_MAX_TRUE_PARTICLES__];
    bool LArFiducialStart[__TMS_MAX_TRUE_PARTICLES__];
    bool LArFiducialTouch[__TMS_MAX_TRUE_PARTICLES__];
    bool LArFiducialEnd[__TMS_MAX_TRUE_PARTICLES__];

    float BirthMomentum[__TMS_MAX_TRUE_PARTICLES__][4];
    float BirthPosition[__TMS_MAX_TRUE_PARTICLES__][4];

    float DeathMomentum[__TMS_MAX_TRUE_PARTICLES__][4];
    float DeathPosition[__TMS_MAX_TRUE_PARTICLES__][4];

    float MomentumZIsLArEnd[__TMS_MAX_TRUE_PARTICLES__][4];
    float PositionZIsLArEnd[__TMS_MAX_TRUE_PARTICLES__][4];

    float MomentumLArStart[__TMS_MAX_TRUE_PARTICLES__][4];
    float PositionLArStart[__TMS_MAX_TRUE_PARTICLES__][4];

    float MomentumLArEnd[__TMS_MAX_TRUE_PARTICLES__][4];
    float PositionLArEnd[__TMS_MAX_TRUE_PARTICLES__][4];

    float MomentumZIsTMSStart[__TMS_MAX_TRUE_PARTICLES__][4];
    float PositionZIsTMSStart[__TMS_MAX_TRUE_PARTICLES__][4];

    float MomentumZIsTMSEnd[__TMS_MAX_TRUE_PARTICLES__][4];
    float PositionZIsTMSEnd[__TMS_MAX_TRUE_PARTICLES__][4];

    float MomentumTMSStart[__TMS_MAX_TRUE_PARTICLES__][4];
    float PositionTMSStart[__TMS_MAX_TRUE_PARTICLES__][4];

    float MomentumTMSFirstTwoModulesEnd[__TMS_MAX_TRUE_PARTICLES__][4];
    float PositionTMSFirstTwoModulesEnd[__TMS_MAX_TRUE_PARTICLES__][4];

    float MomentumTMSThinEnd[__TMS_MAX_TRUE_PARTICLES__][4];
    float PositionTMSThinEnd[__TMS_MAX_TRUE_PARTICLES__][4];

    float MomentumTMSEnd[__TMS_MAX_TRUE_PARTICLES__][4];
    float PositionTMSEnd[__TMS_MAX_TRUE_PARTICLES__][4];
    
    int RecoTrackN;
    float RecoTrackTrueVisibleEnergy[__TMS_MAX_LINES__];
    // deprecated, with pileup we can't guarentee a 1-1 relationship
    int RecoTrackPrimaryParticleIndex[__TMS_MAX_LINES__];
    float RecoTrackPrimaryParticleTrueVisibleEnergy[__TMS_MAX_LINES__];
    int RecoTrackPrimaryParticleTrueNHits[__TMS_MAX_LINES__];
    // deprecated, with pileup we can't guarentee a 1-1 relationship
    int RecoTrackSecondaryParticleIndex[__TMS_MAX_LINES__];
    float RecoTrackSecondaryParticleTrueVisibleEnergy[__TMS_MAX_LINES__]; 
    int RecoTrackSecondaryParticleTrueNHits[__TMS_MAX_LINES__]; 
    
    // Save truth info about primary particle
    int RecoTrackPrimaryParticlePDG[__TMS_MAX_LINES__];
    int RecoTrackPrimaryParticleIsPrimary[__TMS_MAX_LINES__];
    float RecoTrackPrimaryParticleTrueMomentum[__TMS_MAX_LINES__][4];
    float RecoTrackPrimaryParticleTruePositionStart[__TMS_MAX_LINES__][4];
    float RecoTrackPrimaryParticleTruePositionEnd[__TMS_MAX_LINES__][4];
    float RecoTrackPrimaryParticleTrueTrackLength[__TMS_MAX_LINES__];
    float RecoTrackPrimaryParticleTrueTrackLengthIgnoreY[__TMS_MAX_LINES__];
    float RecoTrackPrimaryParticleTrueMomentumEnteringTMS[__TMS_MAX_LINES__][4];
    float RecoTrackPrimaryParticleTruePositionEnteringTMS[__TMS_MAX_LINES__][4];
    float RecoTrackPrimaryParticleTrueMomentumLeavingTMS[__TMS_MAX_LINES__][4];
    float RecoTrackPrimaryParticleTruePositionLeavingTMS[__TMS_MAX_LINES__][4];
    float RecoTrackPrimaryParticleTrueMomentumLeavingLAr[__TMS_MAX_LINES__][4];
    float RecoTrackPrimaryParticleTruePositionLeavingLAr[__TMS_MAX_LINES__][4];
    
    int RecoTrackNHits[__TMS_MAX_LINES__];
    float RecoTrackTrueHitPosition[__TMS_MAX_LINES__][__TMS_MAX_LINE_HITS__][4];
    
    bool RecoTrackPrimaryParticleTMSFiducialStart[__TMS_MAX_LINES__];
    bool RecoTrackPrimaryParticleTMSFiducialTouch[__TMS_MAX_LINES__];
    bool RecoTrackPrimaryParticleTMSFiducialEnd[__TMS_MAX_LINES__];
    bool RecoTrackPrimaryParticleLArFiducialStart[__TMS_MAX_LINES__];
    bool RecoTrackPrimaryParticleLArFiducialTouch[__TMS_MAX_LINES__];
    bool RecoTrackPrimaryParticleLArFiducialEnd[__TMS_MAX_LINES__];
    
    // Save truth info about secondary particle
    int RecoTrackSecondaryParticlePDG[__TMS_MAX_LINES__];
    int RecoTrackSecondaryParticleIsPrimary[__TMS_MAX_LINES__];
    float RecoTrackSecondaryParticleTrueMomentum[__TMS_MAX_LINES__][4];
    float RecoTrackSecondaryParticleTruePositionStart[__TMS_MAX_LINES__][4];
    float RecoTrackSecondaryParticleTruePositionEnd[__TMS_MAX_LINES__][4];

    float RecoTrackPrimaryParticleTrueMomentumTrackStart[__TMS_MAX_LINES__][4];
    float RecoTrackPrimaryParticleTruePositionTrackStart[__TMS_MAX_LINES__][4];
    float RecoTrackPrimaryParticleTrueMomentumTrackEnd[__TMS_MAX_LINES__][4];
    float RecoTrackPrimaryParticleTruePositionTrackEnd[__TMS_MAX_LINES__][4];

    float RecoTrackPrimaryParticleTrueTrackLengthAsMeasured[__TMS_MAX_LINES__];
    float RecoTrackPrimaryParticleTrueTrackLengthAsMeasuredIgnoreY[__TMS_MAX_LINES__];
    float RecoTrackPrimaryParticleTrueTrackLengthRecoStart[__TMS_MAX_LINES__];
    float RecoTrackPrimaryParticleTrueTrackLengthRecoStartIgnoreY[__TMS_MAX_LINES__];
    float RecoTrackPrimaryParticleTrueTrackLengthInTMS[__TMS_MAX_LINES__];
    float RecoTrackPrimaryParticleTrueTrackLengthInTMSIgnoreY[__TMS_MAX_LINES__];
    
    int TruthInfoIndex;
    int TruthInfoNSlices;
    
    int nPrimaryVertices;
    bool HasPileup;
    
    // Ideally I'd use std::vector<Float_t> but they don't fill with the right info
    #define __MAX_TRUE_TREE_ARRAY_LENGTH__ 100000
    #define MYVAR(x) float x[__MAX_TRUE_TREE_ARRAY_LENGTH__]
    #define INTMYVAR(x) int x[__MAX_TRUE_TREE_ARRAY_LENGTH__]
    
    // True hit branches
    int NTrueHits;
    MYVAR(TrueHitX);
    MYVAR(TrueHitY);
    MYVAR(TrueHitZ);
    MYVAR(TrueHitT);
    MYVAR(TrueHitE);
    MYVAR(TrueHitPE);
    MYVAR(TrueHitPEAfterFibers);
    MYVAR(TrueHitPEAfterFibersLongPath);
    MYVAR(TrueHitPEAfterFibersShortPath);
    INTMYVAR(TrueHitBar);
    INTMYVAR(TrueHitPlane);
    INTMYVAR(TrueNTrueParticles);
    MYVAR(TrueLeptonicEnergy);
    MYVAR(TrueHadronicEnergy);
    
    
    // Reco info saved in truth tree to make comparisons easier
    MYVAR(TrueRecoHitX);
    MYVAR(TrueRecoHitY);
    MYVAR(TrueRecoHitZ);
    MYVAR(TrueRecoHitTrackX);
    MYVAR(TrueRecoHitTrackY);
    MYVAR(TrueRecoHitTrackXUncertainty);
    MYVAR(TrueRecoHitTrackYUncertainty);
    MYVAR(TrueRecoHitNotZ);
    MYVAR(TrueRecoHitT);
    MYVAR(TrueRecoHitE);
    MYVAR(TrueRecoHitPE);
    MYVAR(TrueRecoHitEVis);
    bool TrueRecoHitIsPedSupped[__MAX_TRUE_TREE_ARRAY_LENGTH__];
};


#endif
