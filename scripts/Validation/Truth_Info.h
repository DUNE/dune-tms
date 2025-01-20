//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May  2 20:30:27 2024 by ROOT version 6.22/08
// from TTree Truth_Info/Truth_Info
// found on file: ../../neutrino.0_1698176146.edep_TMS_RecoCandidates_Hough_Cluster1.root
//////////////////////////////////////////////////////////

#ifndef Truth_Info_h
#define Truth_Info_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "string"

class Truth_Info {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           EventNo;
   Int_t           SpillNo;
   Int_t           RunNo;
   Bool_t          IsCC;
   std::string    *Interaction;
   Int_t           TruthInfoIndex;
   Int_t           TruthInfoNSlices;
   Int_t           nPrimaryVertices;
   Bool_t          HasPileup;
   Int_t           NeutrinoPDG;
   Float_t         NeutrinoP4[4];
   Float_t         NeutrinoX4[4];
   Int_t           nTrueParticles;
   Int_t           nTruePrimaryParticles;
   Int_t           nTrueForgottenParticles;
   Int_t           VertexID[20000];   //[nTrueParticles]
   Int_t           Parent[20000];   //[nTrueParticles]
   Int_t           TrackId[20000];   //[nTrueParticles]
   Int_t           PDG[20000];   //[nTrueParticles]
   Bool_t          IsPrimary[20000];   //[nTrueParticles]
   Float_t         TrueVisibleEnergy[20000];   //[nTrueParticles]
   Int_t           TrueNHits[20000];   //[nTrueParticles]
   Float_t         TrueVisibleEnergyInSlice[20000];   //[nTrueParticles]
   Int_t           TrueNHitsInSlice[20000];   //[nTrueParticles]
   Float_t         TruePathLength[20000];   //[nTrueParticles]
   Float_t         TruePathLengthIgnoreY[20000];   //[nTrueParticles]
   Float_t         TruePathLengthInTMS[20000];   //[nTrueParticles]
   Float_t         TruePathLengthInTMSIgnoreY[20000];   //[nTrueParticles]
   Bool_t          InteractionTMSFiducial;
   Bool_t          InteractionTMSFirstTwoModules;
   Bool_t          InteractionTMSThin;
   Bool_t          InteractionLArFiducial;
   Bool_t          TMSFiducialStart[20000];   //[nTrueParticles]
   Bool_t          TMSFiducialTouch[20000];   //[nTrueParticles]
   Bool_t          TMSFiducialEnd[20000];   //[nTrueParticles]
   Bool_t          LArFiducialStart[20000];   //[nTrueParticles]
   Bool_t          LArFiducialTouch[20000];   //[nTrueParticles]
   Bool_t          LArFiducialEnd[20000];   //[nTrueParticles]
   Float_t         BirthMomentum[20000][4];   //[nTrueParticles]
   Float_t         BirthPosition[20000][4];   //[nTrueParticles]
   Float_t         DeathMomentum[20000][4];   //[nTrueParticles]
   Float_t         DeathPosition[20000][4];   //[nTrueParticles]
   Float_t         MomentumLArStart[20000][4];   //[nTrueParticles]
   Float_t         PositionLArStart[20000][4];   //[nTrueParticles]
   Float_t         MomentumLArEnd[20000][4];   //[nTrueParticles]
   Float_t         PositionLArEnd[20000][4];   //[nTrueParticles]
   Float_t         MomentumTMSStart[20000][4];   //[nTrueParticles]
   Float_t         PositionTMSStart[20000][4];   //[nTrueParticles]
   Float_t         MomentumTMSFirstTwoModulesEnd[20000][4];   //[nTrueParticles]
   Float_t         PositionTMSFirstTwoModulesEnd[20000][4];   //[nTrueParticles]
   Float_t         MomentumTMSThinEnd[20000][4];   //[nTrueParticles]
   Float_t         PositionTMSThinEnd[20000][4];   //[nTrueParticles]
   Float_t         MomentumTMSEnd[20000][4];   //[nTrueParticles]
   Float_t         PositionTMSEnd[20000][4];   //[nTrueParticles]
   Float_t         MomentumZIsLArEnd[20000][4];   //[nTrueParticles]
   Float_t         PositionZIsLArEnd[20000][4];   //[nTrueParticles]
   Float_t         MomentumZIsTMSStart[20000][4];   //[nTrueParticles]
   Float_t         PositionZIsTMSStart[20000][4];   //[nTrueParticles]
   Float_t         MomentumZIsTMSEnd[20000][4];   //[nTrueParticles]
   Float_t         PositionZIsTMSEnd[20000][4];   //[nTrueParticles]
   Int_t           nParticles;
   Int_t           LeptonPDG;
   Float_t         LeptonP4[4];
   Float_t         LeptonX4[4];
   Float_t         MuonP4[4];
   Float_t         Muon_Vertex[4];
   Float_t         Muon_Death[4];
   Float_t         Muon_TrueKE;
   Float_t         Muon_TrueTrackLength;
   Float_t         LArOuterShellEnergy;
   Float_t         LArOuterShellEnergyFromVertex;
   Float_t         LArTotalEnergy;
   Float_t         LArTotalEnergyFromVertex;
   Float_t         TotalNonTMSEnergy;
   Float_t         TotalNonTMSEnergyFromVertex;
   Int_t           VertexIdOfMostEnergyInEvent;
   Float_t         VisibleEnergyFromUVertexInSlice;
   Float_t         TotalVisibleEnergyFromVertex;
   Float_t         VisibleEnergyFromVVerticesInSlice;
   Float_t         VertexVisibleEnergyFractionInSlice;
   Float_t         PrimaryVertexVisibleEnergyFraction;
   Int_t           RecoTrackN;
   Float_t         RecoTrackTrueVisibleEnergy[10];   //[RecoTrackN]
   Int_t           RecoTrackPrimaryParticleIndex[10];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueVisibleEnergy[10];   //[RecoTrackN]
   Int_t           RecoTrackPrimaryParticleTrueNHits[10]; //[RecoTrackN]
   Int_t           RecoTrackSecondaryParticleIndex[10];   //[RecoTrackN]
   Float_t         RecoTrackSecondaryParticleTrueVisibleEnergy[10];   //[RecoTrackN]
   Int_t           RecoTrackSecondaryParticleTrueNHits[10];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueMomentumTrackStart[10][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTruePositionTrackStart[10][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueMomentumTrackEnd[10][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTruePositionTrackEnd[10][4];   //[RecoTrackN]
   Int_t           RecoTrackNHits[10];   //[RecoTrackN]
   Float_t         RecoTrackTrueHitPosition[10][300][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueTrackLengthAsMeasured[10];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueTrackLengthAsMeasuredIgnoreY[10];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueTrackLengthRecoStart[10];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueTrackLengthRecoStartIgnoreY[10];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueTrackLengthInTMSIgnoreY[10];   //[RecoTrackN]
   Int_t           RecoTrackPrimaryParticlePDG[10];   //[RecoTrackN]
   Bool_t          RecoTrackPrimaryParticleIsPrimary[10];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueMomentum[10][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTruePositionStart[10][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTruePositionEnd[10][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueTrackLength[10];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueTrackLengthIgnoreY[10];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueTrackLengthInTMS[10];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTruePositionEnteringTMS[10][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueMomentumEnteringTMS[10][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueMomentumLeavingTMS[10][4];   //[RecoTrackN]
   Bool_t          RecoTrackPrimaryParticleTMSFiducialStart[10];   //[RecoTrackN]
   Bool_t          RecoTrackPrimaryParticleTMSFiducialTouch[10];   //[RecoTrackN]
   Bool_t          RecoTrackPrimaryParticleTMSFiducialEnd[10];   //[RecoTrackN]
   Bool_t          RecoTrackPrimaryParticleLArFiducialStart[10];   //[RecoTrackN]
   Bool_t          RecoTrackPrimaryParticleLArFiducialTouch[10];   //[RecoTrackN]
   Bool_t          RecoTrackPrimaryParticleLArFiducialEnd[10];   //[RecoTrackN]
   Int_t           RecoTrackSecondaryParticlePDG[10];   //[RecoTrackN]
   Bool_t          RecoTrackSecondaryParticleIsPrimary[10];   //[RecoTrackN]
   Float_t         RecoTrackSecondaryParticleTrueMomentum[10][4];   //[RecoTrackN]
   Float_t         RecoTrackSecondaryParticleTruePositionStart[10][4];   //[RecoTrackN]
   Float_t         RecoTrackSecondaryParticleTruePositionEnd[10][4];   //[RecoTrackN]
   Int_t TrueNonTMSNHits;
   Float_t TrueNonTMSHitPos[100000][4]; //[TrueNonTMSNHits]
   Float_t TrueNonTMSHitEnergy[100000]; //[TrueNonTMSNHits]
   Float_t TrueNonTMSHitHadronicEnergy[100000]; //[TrueNonTMSNHits]
   Float_t TrueNonTMSHitDx[100000]; //[TrueNonTMSNHits]
   Float_t TrueNonTMSHitdEdx[100000]; //[TrueNonTMSNHits]
   Int_t TrueNonTMSHitVertexID[100000]; //[TrueNonTMSNHits]
   
   Int_t           NTrueHits;
   Float_t         TrueHitX[100000];   //[NTrueHits]
   Float_t         TrueHitY[100000];   //[NTrueHits]
   Float_t         TrueHitZ[100000];   //[NTrueHits]
   Float_t         TrueHitT[100000];   //[NTrueHits]
   Float_t         TrueHitE[100000];   //[NTrueHits]
   Float_t         TrueHitPE[100000];   //[NTrueHits]
   Float_t         TrueHitPEAfterFibers[100000];   //[NTrueHits]
   Float_t         TrueHitPEAfterFibersLongPath[100000];   //[NTrueHits]
   Float_t         TrueHitPEAfterFibersShortPath[100000];   //[NTrueHits]
   Int_t           TrueNTrueParticles[100000];   //[NTrueHits]
   Float_t         TrueLeptonicEnergy[100000];   //[NTrueHits]
   Float_t         TrueHadronicEnergy[100000];   //[NTrueHits]
   Float_t         TrueRecoHitX[100000];   //[NTrueHits]
   Float_t         TrueRecoHitY[100000];   //[NTrueHits]
   Float_t         TrueRecoHitZ[100000];   //[NTrueHits]
   Float_t         TrueRecoHitTrackX[100000];   //[NTrueHits]
   Float_t         TrueRecoHitTrackY[100000];   //[NTrueHits]
   Float_t         TrueRecoHitTrackXUncertainty[100000];   //[NTrueHits]
   Float_t         TrueRecoHitTrackYUncertainty[100000];   //[NTrueHits]
   Float_t         TrueRecoHitNotZ[100000];   //[NTrueHits]
   Float_t         TrueRecoHitT[100000];   //[NTrueHits]
   Float_t         TrueRecoHitE[100000];   //[NTrueHits]
   Float_t         TrueRecoHitPE[100000];   //[NTrueHits]
   Float_t         TrueRecoHitEVis[100000];   //[NTrueHits]
   Bool_t          TrueRecoHitIsPedSupped[100000];   //[NTrueHits]
   Int_t           TrueHitBar[100000];   //[NTrueHits]
   Int_t           TrueHitView[100000];   //[NTrueHits]
   Int_t           TrueHitPlane[100000];   //[NTrueHits]

   // List of branches
   TBranch        *b_EventNo;   //!
   TBranch        *b_SpillNo;   //!
   TBranch        *b_RunNo;   //!
   TBranch        *b_IsCC;   //!
   TBranch        *b_Interaction;   //!
   TBranch        *b_TruthInfoIndex;   //!
   TBranch        *b_TruthInfoNSlices;   //!
   TBranch        *b_nPrimaryVertices;   //!
   TBranch        *b_HasPileup;   //!
   TBranch        *b_NeutrinoPDG;   //!
   TBranch        *b_NeutrinoP4;   //!
   TBranch        *b_NeutrinoX4;   //!
   TBranch        *b_nTrueParticles;   //!
   TBranch        *b_nTruePrimaryParticles;   //!
   TBranch        *b_nTrueForgottenParticles;   //!
   TBranch        *b_VertexID;   //!
   TBranch        *b_Parent;   //!
   TBranch        *b_TrackId;   //!
   TBranch        *b_PDG;   //!
   TBranch        *b_IsPrimary;   //!
   TBranch        *b_TrueVisibleEnergy;   //!
   TBranch        *b_TrueNHits;   //!
   TBranch        *b_TrueVisibleEnergyInSlice;   //!
   TBranch        *b_TrueNHitsInSlice;   //!
   TBranch        *b_TruePathLength;   //!
   TBranch        *b_TruePathLengthIgnoreY;   //!
   TBranch        *b_TruePathLengthInTMS;   //!
   TBranch        *b_TruePathLengthInTMSIgnoreY;   //!
   TBranch        *b_InteractionTMSFiducial;   //!
   TBranch        *b_InteractionTMSFirstTwoModules;   //!
   TBranch        *b_InteractionTMSThin;   //!
   TBranch        *b_InteractionLArFiducial;   //!
   TBranch        *b_TMSFiducialStart;   //!
   TBranch        *b_TMSFiducialTouch;   //!
   TBranch        *b_TMSFiducialEnd;   //!
   TBranch        *b_LArFiducialStart;   //!
   TBranch        *b_LArFiducialTouch;   //!
   TBranch        *b_LArFiducialEnd;   //!
   TBranch        *b_BirthMomentum;   //!
   TBranch        *b_BirthPosition;   //!
   TBranch        *b_DeathMomentum;   //!
   TBranch        *b_DeathPosition;   //!
   TBranch        *b_MomentumLArStart;   //!
   TBranch        *b_PositionLArStart;   //!
   TBranch        *b_MomentumLArEnd;   //!
   TBranch        *b_PositionLArEnd;   //!
   TBranch        *b_MomentumTMSStart;   //!
   TBranch        *b_PositionTMSStart;   //!
   TBranch        *b_MomentumTMSFirstTwoModulesEnd;   //!
   TBranch        *b_PositionTMSFirstTwoModulesEnd;   //!
   TBranch        *b_MomentumTMSThinEnd;   //!
   TBranch        *b_PositionTMSThinEnd;   //!
   TBranch        *b_MomentumTMSEnd;   //!
   TBranch        *b_PositionTMSEnd;   //!
   TBranch        *b_MomentumZIsLArEnd;   //!
   TBranch        *b_PositionZIsLArEnd;   //!
   TBranch        *b_MomentumZIsTMSStart;   //!
   TBranch        *b_PositionZIsTMSStart;   //!
   TBranch        *b_MomentumZIsTMSEnd;   //!
   TBranch        *b_PositionZIsTMSEnd;   //!
   TBranch        *b_nParticles;   //!
   TBranch        *b_LeptonPDG;   //!
   TBranch        *b_LeptonP4;   //!
   TBranch        *b_LeptonX4;   //!
   TBranch        *b_MuonP4;   //!
   TBranch        *b_Muon_Vertex;   //!
   TBranch        *b_Muon_Death;   //!
   TBranch        *b_Muon_TrueKE;   //!
   TBranch        *b_Muon_TrueTrackLength;   //!
   TBranch        *b_LArOuterShellEnergy;   //!
   TBranch        *b_LArOuterShellEnergyFromVertex;   //!
   TBranch        *b_LArTotalEnergy;   //!
   TBranch        *b_LArTotalEnergyFromVertex;   //!
   TBranch        *b_TotalNonTMSEnergy;   //!
   TBranch        *b_TotalNonTMSEnergyFromVertex;   //!
   TBranch        *b_VertexIdOfMostEnergyInEvent;   //!
   TBranch        *b_VisibleEnergyFromVertexInSlice;   //!
   TBranch        *b_TotalVisibleEnergyFromVertex;   //!
   TBranch        *b_VisibleEnergyFromOtherVerticesInSlice;   //!
   TBranch        *b_VertexVisibleEnergyFractionInSlice;   //!
   TBranch        *b_PrimaryVertexVisibleEnergyFraction;   //!
   TBranch        *b_VisibleEnergyFromUVertexInSlice;   //!
   TBranch        *b_VisibleEnergyFromVVerticesInSlice;   //!
   TBranch        *b_RecoTrackN;   //!
   TBranch        *b_RecoTrackTrueVisibleEnergy;   //!
   TBranch        *b_RecoTrackPrimaryParticleIndex;   //!
   TBranch        *b_RecoTrackPrimaryParticleTrueVisibleEnergy;   //!
   TBranch        *b_RecoTrackPrimaryParticleTrueNHits;   //!
   TBranch        *b_RecoTrackSecondaryParticleIndex;   //!
   TBranch        *b_RecoTrackSecondaryParticleTrueVisibleEnergy;   //!
   TBranch        *b_RecoTrackSecondaryParticleTrueNHits;   //!
   TBranch        *b_RecoTrackPrimaryParticleTrueMomentumTrackStart;   //!
   TBranch        *b_RecoTrackPrimaryParticleTruePositionTrackStart;   //!
   TBranch        *b_RecoTrackPrimaryParticleTrueMomentumTrackEnd;   //!
   TBranch        *b_RecoTrackPrimaryParticleTruePositionTrackEnd;   //!
   TBranch        *b_RecoTrackNHits;   //!
   TBranch        *b_RecoTrackTrueHitPosition;   //!
   TBranch        *b_RecoTrackPrimaryParticleTrueTrackLengthAsMeasured;   //!
   TBranch        *b_RecoTrackPrimaryParticleTrueTrackLengthAsMeasuredIgnoreY;   //!
   TBranch        *b_RecoTrackPrimaryParticleTrueTrackLengthRecoStart;   //!
   TBranch        *b_RecoTrackPrimaryParticleTrueTrackLengthRecoStartIgnoreY;   //!
   TBranch        *b_RecoTrackPrimaryParticleTrueTrackLengthInTMSIgnoreY;   //!
   TBranch        *b_RecoTrackPrimaryParticlePDG;   //!
   TBranch        *b_RecoTrackPrimaryParticleIsPrimary;   //!
   TBranch        *b_RecoTrackPrimaryParticleTrueMomentum;   //!
   TBranch        *b_RecoTrackPrimaryParticleTruePositionStart;   //!
   TBranch        *b_RecoTrackPrimaryParticleTruePositionEnd;   //!
   TBranch        *b_RecoTrackPrimaryParticleTrueTrackLength;   //!
   TBranch        *b_RecoTrackPrimaryParticleTrueTrackLengthIgnoreY;   //!
   TBranch        *b_RecoTrackPrimaryParticleTrueTrackLengthInTMS;   //!
   TBranch        *b_RecoTrackPrimaryParticleTruePositionEnteringTMS;   //!
   TBranch        *b_RecoTrackPrimaryParticleTrueMomentumEnteringTMS;   //!
   TBranch        *b_RecoTrackPrimaryParticleTrueMomentumLeavingTMS;   //!
   TBranch        *b_RecoTrackPrimaryParticleTMSFiducialStart;   //!
   TBranch        *b_RecoTrackPrimaryParticleTMSFiducialTouch;   //!
   TBranch        *b_RecoTrackPrimaryParticleTMSFiducialEnd;   //!
   TBranch        *b_RecoTrackPrimaryParticleLArFiducialStart;   //!
   TBranch        *b_RecoTrackPrimaryParticleLArFiducialTouch;   //!
   TBranch        *b_RecoTrackPrimaryParticleLArFiducialEnd;   //!
   TBranch        *b_RecoTrackSecondaryParticlePDG;   //!
   TBranch        *b_RecoTrackSecondaryParticleIsPrimary;   //!
   TBranch        *b_RecoTrackSecondaryParticleTrueMomentum;   //!
   TBranch        *b_RecoTrackSecondaryParticleTruePositionStart;   //!
   TBranch        *b_RecoTrackSecondaryParticleTruePositionEnd;   //!
   TBranch        *b_TrueNonTMSNHits;   //!
   TBranch        *b_TrueNonTMSHitPos;   //!
   TBranch        *b_TrueNonTMSHitEnergy;   //!
   TBranch        *b_TrueNonTMSHitHadronicEnergy;   //!
   TBranch        *b_TrueNonTMSHitDx;   //!
   TBranch        *b_TrueNonTMSHitdEdx;   //!
   TBranch        *b_TrueNonTMSHitVertexID;   //!
   TBranch        *b_NTrueHits;   //!
   TBranch        *b_TrueHitX;   //!
   TBranch        *b_TrueHitY;   //!
   TBranch        *b_TrueHitZ;   //!
   TBranch        *b_TrueHitT;   //!
   TBranch        *b_TrueHitE;   //!
   TBranch        *b_TrueHitPE;   //!
   TBranch        *b_TrueHitPEAfterFibers;   //!
   TBranch        *b_TrueHitPEAfterFibersLongPath;   //!
   TBranch        *b_TrueHitPEAfterFibersShortPath;   //!
   TBranch        *b_TrueNTrueParticles;   //!
   TBranch        *b_TrueLeptonicEnergy;   //!
   TBranch        *b_TrueHadronicEnergy;   //!
   TBranch        *b_TrueRecoHitX;   //!
   TBranch        *b_TrueRecoHitY;   //!
   TBranch        *b_TrueRecoHitZ;   //!
   TBranch        *b_TrueRecoHitTrackX;   //!
   TBranch        *b_TrueRecoHitTrackY;   //!
   TBranch        *b_TrueRecoHitTrackXUncertainty;   //!
   TBranch        *b_TrueRecoHitTrackYUncertainty;   //!
   TBranch        *b_TrueRecoHitNotZ;   //!
   TBranch        *b_TrueRecoHitT;   //!
   TBranch        *b_TrueRecoHitE;   //!
   TBranch        *b_TrueRecoHitPE;   //!
   TBranch        *b_TrueRecoHitEVis;   //!
   TBranch        *b_TrueRecoHitIsPedSupped;   //!
   TBranch        *b_TrueHitBar;   //!
   TBranch        *b_TrueHitView;   //!
   TBranch        *b_TrueHitPlane;   //!

   Truth_Info(TTree *tree=0);
   virtual ~Truth_Info();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Long64_t GetEntriesFast() { return fChain->GetEntriesFast(); };
   virtual bool HasBranch(const char* branch) { return fChain->GetBranch(branch) != NULL; };
};

#endif

#define Truth_Info_cxx
#ifdef Truth_Info_cxx
Truth_Info::Truth_Info(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../neutrino.0_1698176146.edep_TMS_RecoCandidates_Hough_Cluster1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../neutrino.0_1698176146.edep_TMS_RecoCandidates_Hough_Cluster1.root");
      }
      f->GetObject("Truth_Info",tree);

   }
   Init(tree);
}

Truth_Info::~Truth_Info()
{
   if (!fChain) return;
   //delete fChain->GetCurrentFile();
}

Int_t Truth_Info::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Truth_Info::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Truth_Info::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Interaction = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EventNo", &EventNo, &b_EventNo);
   fChain->SetBranchAddress("SpillNo", &SpillNo, &b_SpillNo);
   fChain->SetBranchAddress("RunNo", &RunNo, &b_RunNo);
   fChain->SetBranchAddress("IsCC", &IsCC, &b_IsCC);
   fChain->SetBranchAddress("Interaction", &Interaction, &b_Interaction);
   fChain->SetBranchAddress("TruthInfoIndex", &TruthInfoIndex, &b_TruthInfoIndex);
   fChain->SetBranchAddress("TruthInfoNSlices", &TruthInfoNSlices, &b_TruthInfoNSlices);
   fChain->SetBranchAddress("nPrimaryVertices", &nPrimaryVertices, &b_nPrimaryVertices);
   fChain->SetBranchAddress("HasPileup", &HasPileup, &b_HasPileup);
   fChain->SetBranchAddress("NeutrinoPDG", &NeutrinoPDG, &b_NeutrinoPDG);
   fChain->SetBranchAddress("NeutrinoP4", NeutrinoP4, &b_NeutrinoP4);
   fChain->SetBranchAddress("NeutrinoX4", NeutrinoX4, &b_NeutrinoX4);
   fChain->SetBranchAddress("nTrueParticles", &nTrueParticles, &b_nTrueParticles);
   fChain->SetBranchAddress("nTruePrimaryParticles", &nTruePrimaryParticles, &b_nTruePrimaryParticles);
   fChain->SetBranchAddress("nTrueForgottenParticles", &nTrueForgottenParticles, &b_nTrueForgottenParticles);
   fChain->SetBranchAddress("VertexID", VertexID, &b_VertexID);
   fChain->SetBranchAddress("Parent", Parent, &b_Parent);
   fChain->SetBranchAddress("TrackId", TrackId, &b_TrackId);
   fChain->SetBranchAddress("PDG", PDG, &b_PDG);
   fChain->SetBranchAddress("IsPrimary", IsPrimary, &b_IsPrimary);
   fChain->SetBranchAddress("TrueVisibleEnergy", TrueVisibleEnergy, &b_TrueVisibleEnergy);
   fChain->SetBranchAddress("TrueNHits", TrueNHits, &b_TrueNHits);
   fChain->SetBranchAddress("TrueVisibleEnergyInSlice", TrueVisibleEnergyInSlice, &b_TrueVisibleEnergyInSlice);
   fChain->SetBranchAddress("TrueNHitsInSlice", TrueNHitsInSlice, &b_TrueNHitsInSlice);
   fChain->SetBranchAddress("TruePathLength", TruePathLength, &b_TruePathLength);
   fChain->SetBranchAddress("TruePathLengthIgnoreY", TruePathLengthIgnoreY, &b_TruePathLengthIgnoreY);
   fChain->SetBranchAddress("TruePathLengthInTMS", TruePathLengthInTMS, &b_TruePathLengthInTMS);
   fChain->SetBranchAddress("TruePathLengthInTMSIgnoreY", TruePathLengthInTMSIgnoreY, &b_TruePathLengthInTMSIgnoreY);
   fChain->SetBranchAddress("InteractionTMSFiducial", &InteractionTMSFiducial, &b_InteractionTMSFiducial);
   fChain->SetBranchAddress("InteractionTMSFirstTwoModules", &InteractionTMSFirstTwoModules, &b_InteractionTMSFirstTwoModules);
   fChain->SetBranchAddress("InteractionTMSThin", &InteractionTMSThin, &b_InteractionTMSThin);
   fChain->SetBranchAddress("InteractionLArFiducial", &InteractionLArFiducial, &b_InteractionLArFiducial);
   fChain->SetBranchAddress("TMSFiducialStart", TMSFiducialStart, &b_TMSFiducialStart);
   fChain->SetBranchAddress("TMSFiducialTouch", TMSFiducialTouch, &b_TMSFiducialTouch);
   fChain->SetBranchAddress("TMSFiducialEnd", TMSFiducialEnd, &b_TMSFiducialEnd);
   fChain->SetBranchAddress("LArFiducialStart", LArFiducialStart, &b_LArFiducialStart);
   fChain->SetBranchAddress("LArFiducialTouch", LArFiducialTouch, &b_LArFiducialTouch);
   fChain->SetBranchAddress("LArFiducialEnd", LArFiducialEnd, &b_LArFiducialEnd);
   fChain->SetBranchAddress("BirthMomentum", BirthMomentum, &b_BirthMomentum);
   fChain->SetBranchAddress("BirthPosition", BirthPosition, &b_BirthPosition);
   fChain->SetBranchAddress("DeathMomentum", DeathMomentum, &b_DeathMomentum);
   fChain->SetBranchAddress("DeathPosition", DeathPosition, &b_DeathPosition);
   fChain->SetBranchAddress("MomentumLArStart", MomentumLArStart, &b_MomentumLArStart);
   fChain->SetBranchAddress("PositionLArStart", PositionLArStart, &b_PositionLArStart);
   fChain->SetBranchAddress("MomentumLArEnd", MomentumLArEnd, &b_MomentumLArEnd);
   fChain->SetBranchAddress("PositionLArEnd", PositionLArEnd, &b_PositionLArEnd);
   fChain->SetBranchAddress("MomentumTMSStart", MomentumTMSStart, &b_MomentumTMSStart);
   fChain->SetBranchAddress("PositionTMSStart", PositionTMSStart, &b_PositionTMSStart);
   fChain->SetBranchAddress("MomentumTMSFirstTwoModulesEnd", MomentumTMSFirstTwoModulesEnd, &b_MomentumTMSFirstTwoModulesEnd);
   fChain->SetBranchAddress("PositionTMSFirstTwoModulesEnd", PositionTMSFirstTwoModulesEnd, &b_PositionTMSFirstTwoModulesEnd);
   fChain->SetBranchAddress("MomentumTMSThinEnd", MomentumTMSThinEnd, &b_MomentumTMSThinEnd);
   fChain->SetBranchAddress("PositionTMSThinEnd", PositionTMSThinEnd, &b_PositionTMSThinEnd);
   fChain->SetBranchAddress("MomentumTMSEnd", MomentumTMSEnd, &b_MomentumTMSEnd);
   fChain->SetBranchAddress("PositionTMSEnd", PositionTMSEnd, &b_PositionTMSEnd);
   fChain->SetBranchAddress("MomentumZIsLArEnd", MomentumZIsLArEnd, &b_MomentumZIsLArEnd);
   fChain->SetBranchAddress("PositionZIsLArEnd", PositionZIsLArEnd, &b_PositionZIsLArEnd);
   fChain->SetBranchAddress("MomentumZIsTMSStart", MomentumZIsTMSStart, &b_MomentumZIsTMSStart);
   fChain->SetBranchAddress("PositionZIsTMSStart", PositionZIsTMSStart, &b_PositionZIsTMSStart);
   fChain->SetBranchAddress("MomentumZIsTMSEnd", MomentumZIsTMSEnd, &b_MomentumZIsTMSEnd);
   fChain->SetBranchAddress("PositionZIsTMSEnd", PositionZIsTMSEnd, &b_PositionZIsTMSEnd);
   fChain->SetBranchAddress("nParticles", &nParticles, &b_nParticles);
   fChain->SetBranchAddress("LeptonPDG", &LeptonPDG, &b_LeptonPDG);
   fChain->SetBranchAddress("LeptonP4", LeptonP4, &b_LeptonP4);
   fChain->SetBranchAddress("LeptonX4", LeptonX4, &b_LeptonX4);
   fChain->SetBranchAddress("MuonP4", MuonP4, &b_MuonP4);
   fChain->SetBranchAddress("Muon_Vertex", Muon_Vertex, &b_Muon_Vertex);
   fChain->SetBranchAddress("Muon_Death", Muon_Death, &b_Muon_Death);
   fChain->SetBranchAddress("Muon_TrueKE", &Muon_TrueKE, &b_Muon_TrueKE);
   fChain->SetBranchAddress("Muon_TrueTrackLength", &Muon_TrueTrackLength, &b_Muon_TrueTrackLength);
   fChain->SetBranchAddress("LArOuterShellEnergy", &LArOuterShellEnergy, &b_LArOuterShellEnergy);
   fChain->SetBranchAddress("LArOuterShellEnergyFromVertex", &LArOuterShellEnergyFromVertex, &b_LArOuterShellEnergyFromVertex);
   fChain->SetBranchAddress("LArTotalEnergy", &LArTotalEnergy, &b_LArTotalEnergy);
   fChain->SetBranchAddress("LArTotalEnergyFromVertex", &LArTotalEnergyFromVertex, &b_LArTotalEnergyFromVertex);
   fChain->SetBranchAddress("TotalNonTMSEnergy", &TotalNonTMSEnergy, &b_TotalNonTMSEnergy);
   fChain->SetBranchAddress("TotalNonTMSEnergyFromVertex", &TotalNonTMSEnergyFromVertex, &b_TotalNonTMSEnergyFromVertex);
   fChain->SetBranchAddress("VertexIdOfMostEnergyInEvent", &VertexIdOfMostEnergyInEvent, &b_VertexIdOfMostEnergyInEvent);
   fChain->SetBranchAddress("VisibleEnergyFromUVertexInSlice", &VisibleEnergyFromUVertexInSlice, &b_VisibleEnergyFromUVertexInSlice);
   fChain->SetBranchAddress("TotalVisibleEnergyFromVertex", &TotalVisibleEnergyFromVertex, &b_TotalVisibleEnergyFromVertex);
   fChain->SetBranchAddress("VisibleEnergyFromVVerticesInSlice", &VisibleEnergyFromVVerticesInSlice, &b_VisibleEnergyFromVVerticesInSlice);
   fChain->SetBranchAddress("VertexVisibleEnergyFractionInSlice", &VertexVisibleEnergyFractionInSlice, &b_VertexVisibleEnergyFractionInSlice);
   fChain->SetBranchAddress("PrimaryVertexVisibleEnergyFraction", &PrimaryVertexVisibleEnergyFraction, &b_PrimaryVertexVisibleEnergyFraction);
   fChain->SetBranchAddress("RecoTrackN", &RecoTrackN, &b_RecoTrackN);
   fChain->SetBranchAddress("RecoTrackTrueVisibleEnergy", RecoTrackTrueVisibleEnergy, &b_RecoTrackTrueVisibleEnergy);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleIndex", RecoTrackPrimaryParticleIndex, &b_RecoTrackPrimaryParticleIndex);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueVisibleEnergy", RecoTrackPrimaryParticleTrueVisibleEnergy, &b_RecoTrackPrimaryParticleTrueVisibleEnergy);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueNHits", RecoTrackPrimaryParticleTrueNHits, &b_RecoTrackPrimaryParticleTrueNHits);
   fChain->SetBranchAddress("RecoTrackSecondaryParticleIndex", RecoTrackSecondaryParticleIndex, &b_RecoTrackSecondaryParticleIndex);
   fChain->SetBranchAddress("RecoTrackSecondaryParticleTrueVisibleEnergy", RecoTrackSecondaryParticleTrueVisibleEnergy, &b_RecoTrackSecondaryParticleTrueVisibleEnergy);
   fChain->SetBranchAddress("RecoTrackSecondaryParticleTrueNHits", RecoTrackSecondaryParticleTrueNHits, &b_RecoTrackSecondaryParticleTrueNHits);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueMomentumTrackStart", RecoTrackPrimaryParticleTrueMomentumTrackStart, &b_RecoTrackPrimaryParticleTrueMomentumTrackStart);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTruePositionTrackStart", RecoTrackPrimaryParticleTruePositionTrackStart, &b_RecoTrackPrimaryParticleTruePositionTrackStart);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueMomentumTrackEnd", RecoTrackPrimaryParticleTrueMomentumTrackEnd, &b_RecoTrackPrimaryParticleTrueMomentumTrackEnd);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTruePositionTrackEnd", RecoTrackPrimaryParticleTruePositionTrackEnd, &b_RecoTrackPrimaryParticleTruePositionTrackEnd);
   fChain->SetBranchAddress("RecoTrackNHits", RecoTrackNHits, &b_RecoTrackNHits);
   fChain->SetBranchAddress("RecoTrackTrueHitPosition", RecoTrackTrueHitPosition, &b_RecoTrackTrueHitPosition);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueTrackLengthAsMeasured", RecoTrackPrimaryParticleTrueTrackLengthAsMeasured, &b_RecoTrackPrimaryParticleTrueTrackLengthAsMeasured);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueTrackLengthAsMeasuredIgnoreY", RecoTrackPrimaryParticleTrueTrackLengthAsMeasuredIgnoreY, &b_RecoTrackPrimaryParticleTrueTrackLengthAsMeasuredIgnoreY);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueTrackLengthRecoStart", RecoTrackPrimaryParticleTrueTrackLengthRecoStart, &b_RecoTrackPrimaryParticleTrueTrackLengthRecoStart);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueTrackLengthRecoStartIgnoreY", RecoTrackPrimaryParticleTrueTrackLengthRecoStartIgnoreY, &b_RecoTrackPrimaryParticleTrueTrackLengthRecoStartIgnoreY);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueTrackLengthInTMSIgnoreY", RecoTrackPrimaryParticleTrueTrackLengthInTMSIgnoreY, &b_RecoTrackPrimaryParticleTrueTrackLengthInTMSIgnoreY);
   fChain->SetBranchAddress("RecoTrackPrimaryParticlePDG", RecoTrackPrimaryParticlePDG, &b_RecoTrackPrimaryParticlePDG);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleIsPrimary", RecoTrackPrimaryParticleIsPrimary, &b_RecoTrackPrimaryParticleIsPrimary);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueMomentum", RecoTrackPrimaryParticleTrueMomentum, &b_RecoTrackPrimaryParticleTrueMomentum);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTruePositionStart", RecoTrackPrimaryParticleTruePositionStart, &b_RecoTrackPrimaryParticleTruePositionStart);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTruePositionEnd", RecoTrackPrimaryParticleTruePositionEnd, &b_RecoTrackPrimaryParticleTruePositionEnd);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueTrackLength", RecoTrackPrimaryParticleTrueTrackLength, &b_RecoTrackPrimaryParticleTrueTrackLength);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueTrackLengthIgnoreY", RecoTrackPrimaryParticleTrueTrackLengthIgnoreY, &b_RecoTrackPrimaryParticleTrueTrackLengthIgnoreY);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueTrackLengthInTMS", RecoTrackPrimaryParticleTrueTrackLengthInTMS, &b_RecoTrackPrimaryParticleTrueTrackLengthInTMS);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTruePositionEnteringTMS", RecoTrackPrimaryParticleTruePositionEnteringTMS, &b_RecoTrackPrimaryParticleTruePositionEnteringTMS);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueMomentumEnteringTMS", RecoTrackPrimaryParticleTrueMomentumEnteringTMS, &b_RecoTrackPrimaryParticleTrueMomentumEnteringTMS);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueMomentumLeavingTMS", RecoTrackPrimaryParticleTrueMomentumLeavingTMS, &b_RecoTrackPrimaryParticleTrueMomentumLeavingTMS);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTMSFiducialStart", RecoTrackPrimaryParticleTMSFiducialStart, &b_RecoTrackPrimaryParticleTMSFiducialStart);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTMSFiducialTouch", RecoTrackPrimaryParticleTMSFiducialTouch, &b_RecoTrackPrimaryParticleTMSFiducialTouch);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTMSFiducialEnd", RecoTrackPrimaryParticleTMSFiducialEnd, &b_RecoTrackPrimaryParticleTMSFiducialEnd);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleLArFiducialStart", RecoTrackPrimaryParticleLArFiducialStart, &b_RecoTrackPrimaryParticleLArFiducialStart);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleLArFiducialTouch", RecoTrackPrimaryParticleLArFiducialTouch, &b_RecoTrackPrimaryParticleLArFiducialTouch);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleLArFiducialEnd", RecoTrackPrimaryParticleLArFiducialEnd, &b_RecoTrackPrimaryParticleLArFiducialEnd);
   fChain->SetBranchAddress("RecoTrackSecondaryParticlePDG", RecoTrackSecondaryParticlePDG, &b_RecoTrackSecondaryParticlePDG);
   fChain->SetBranchAddress("RecoTrackSecondaryParticleIsPrimary", RecoTrackSecondaryParticleIsPrimary, &b_RecoTrackSecondaryParticleIsPrimary);
   fChain->SetBranchAddress("RecoTrackSecondaryParticleTrueMomentum", RecoTrackSecondaryParticleTrueMomentum, &b_RecoTrackSecondaryParticleTrueMomentum);
   fChain->SetBranchAddress("RecoTrackSecondaryParticleTruePositionStart", RecoTrackSecondaryParticleTruePositionStart, &b_RecoTrackSecondaryParticleTruePositionStart);
   fChain->SetBranchAddress("RecoTrackSecondaryParticleTruePositionEnd", RecoTrackSecondaryParticleTruePositionEnd, &b_RecoTrackSecondaryParticleTruePositionEnd);
   fChain->SetBranchAddress("TrueNonTMSNHits", &TrueNonTMSNHits, &b_TrueNonTMSNHits);
   fChain->SetBranchAddress("TrueNonTMSHitPos", TrueNonTMSHitPos, &b_TrueNonTMSHitPos);
   fChain->SetBranchAddress("TrueNonTMSHitEnergy", TrueNonTMSHitEnergy, &b_TrueNonTMSHitEnergy);
   fChain->SetBranchAddress("TrueNonTMSHitHadronicEnergy", TrueNonTMSHitHadronicEnergy, &b_TrueNonTMSHitHadronicEnergy);
   fChain->SetBranchAddress("TrueNonTMSHitDx", TrueNonTMSHitDx, &b_TrueNonTMSHitDx);
   fChain->SetBranchAddress("TrueNonTMSHitdEdx", TrueNonTMSHitdEdx, &b_TrueNonTMSHitdEdx);
   fChain->SetBranchAddress("TrueNonTMSHitVertexID", TrueNonTMSHitVertexID, &b_TrueNonTMSHitVertexID);
   fChain->SetBranchAddress("NTrueHits", &NTrueHits, &b_NTrueHits);
   fChain->SetBranchAddress("TrueHitX", TrueHitX, &b_TrueHitX);
   fChain->SetBranchAddress("TrueHitY", TrueHitY, &b_TrueHitY);
   fChain->SetBranchAddress("TrueHitZ", TrueHitZ, &b_TrueHitZ);
   fChain->SetBranchAddress("TrueHitT", TrueHitT, &b_TrueHitT);
   fChain->SetBranchAddress("TrueHitE", TrueHitE, &b_TrueHitE);
   fChain->SetBranchAddress("TrueHitPE", TrueHitPE, &b_TrueHitPE);
   fChain->SetBranchAddress("TrueHitPEAfterFibers", TrueHitPEAfterFibers, &b_TrueHitPEAfterFibers);
   fChain->SetBranchAddress("TrueHitPEAfterFibersLongPath", TrueHitPEAfterFibersLongPath, &b_TrueHitPEAfterFibersLongPath);
   fChain->SetBranchAddress("TrueHitPEAfterFibersShortPath", TrueHitPEAfterFibersShortPath, &b_TrueHitPEAfterFibersShortPath);
   fChain->SetBranchAddress("TrueNTrueParticles", TrueNTrueParticles, &b_TrueNTrueParticles);
   fChain->SetBranchAddress("TrueLeptonicEnergy", TrueLeptonicEnergy, &b_TrueLeptonicEnergy);
   fChain->SetBranchAddress("TrueHadronicEnergy", TrueHadronicEnergy, &b_TrueHadronicEnergy);
   fChain->SetBranchAddress("TrueRecoHitX", TrueRecoHitX, &b_TrueRecoHitX);
   fChain->SetBranchAddress("TrueRecoHitY", TrueRecoHitY, &b_TrueRecoHitY);
   fChain->SetBranchAddress("TrueRecoHitZ", TrueRecoHitZ, &b_TrueRecoHitZ);
   fChain->SetBranchAddress("TrueRecoHitTrackX", TrueRecoHitTrackX, &b_TrueRecoHitTrackX);
   fChain->SetBranchAddress("TrueRecoHitTrackY", TrueRecoHitTrackY, &b_TrueRecoHitTrackY);
   fChain->SetBranchAddress("TrueRecoHitTrackXUncertainty", TrueRecoHitTrackXUncertainty, &b_TrueRecoHitTrackXUncertainty);
   fChain->SetBranchAddress("TrueRecoHitTrackYUncertainty", TrueRecoHitTrackYUncertainty, &b_TrueRecoHitTrackYUncertainty);
   fChain->SetBranchAddress("TrueRecoHitNotZ", TrueRecoHitNotZ, &b_TrueRecoHitNotZ);
   fChain->SetBranchAddress("TrueRecoHitT", TrueRecoHitT, &b_TrueRecoHitT);
   fChain->SetBranchAddress("TrueRecoHitE", TrueRecoHitE, &b_TrueRecoHitE);
   fChain->SetBranchAddress("TrueRecoHitPE", TrueRecoHitPE, &b_TrueRecoHitPE);
   fChain->SetBranchAddress("TrueRecoHitEVis", TrueRecoHitEVis, &b_TrueRecoHitEVis);
   fChain->SetBranchAddress("TrueRecoHitIsPedSupped", TrueRecoHitIsPedSupped, &b_TrueRecoHitIsPedSupped);
   fChain->SetBranchAddress("TrueHitBar", TrueHitBar, &b_TrueHitBar);
   fChain->SetBranchAddress("TrueHitView", TrueHitView, &b_TrueHitView);
   fChain->SetBranchAddress("TrueHitPlane", TrueHitPlane, &b_TrueHitPlane);
   Notify();
}

Bool_t Truth_Info::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Truth_Info::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Truth_Info::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   entry += 1;
   return 1;
}


void Truth_Info::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Truth_Info.C
//      root> Truth_Info t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
#endif // #ifdef Truth_Info_cxx
