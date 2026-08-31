//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 30 11:57:19 2026 by ROOT version 6.28/12
// from TTree Truth_Info/Truth_Info
// found on file: /exp/dune/app/users/kleykamp/2026_dune-tms/dune-tms/MicroProdN4p1_NDComplex_FHC.spill.full.0002459.EDEPSIM_SPILLS_TMS_RecoCandidates_Hough_Cluster1.root
//////////////////////////////////////////////////////////

#ifndef Truth_Info_h
#define Truth_Info_h
#define Truth_Info_cxx

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "string"

using std::string;

class Truth_Info {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   enum {
      kMaxRecoTracks = 100,
      kMaxRecoTrackHits = 2000,
      kMaxTrueHits = 100000
   };

   // Declaration of leaf types
   Int_t           EventNo;
   Int_t           SpillNo;
   Int_t           RunNo;
   Bool_t          IsCC;
   string          *Interaction;
   Int_t           TruthInfoIndex;
   Int_t           TruthInfoNSlices;
   Int_t           nPrimaryVertices;
   Bool_t          HasPileup;
   Int_t           NeutrinoPDG;
   Float_t         NeutrinoP4[4];
   Float_t         NeutrinoX4[4];
   Int_t           nParticles;
   Int_t           LeptonPDG;
   Float_t         LeptonP4[4];
   Float_t         LeptonX4[4];
   Float_t         MuonP4[4];
   Float_t         Muon_Vertex[4];
   Float_t         Muon_Death[4];
   Float_t         Muon_TrueKE;
   Float_t         Muon_TrueTrackLength;
   Int_t           VertexIdOfMostEnergyInEvent;
   Float_t         VisibleEnergyFromVertexInSlice;
   Float_t         TotalVisibleEnergyFromVertex;
   Float_t         VisibleEnergyFromOtherVerticesInSlice;
   Float_t         VertexVisibleEnergyFractionInSlice;
   Float_t         PrimaryVertexVisibleEnergyFraction;
   Float_t         LArOuterShellEnergy;
   Float_t         LArTotalEnergy;
   Float_t         TotalNonTMSEnergy;
   Float_t         LArOuterShellEnergyFromVertex;
   Float_t         LArTotalEnergyFromVertex;
   Float_t         TotalNonTMSEnergyFromVertex;
   Int_t           RecoTrackN;
   Float_t         RecoTrackTrueVisibleEnergy[kMaxRecoTracks];   //[RecoTrackN]
   Int_t           RecoTrackPrimaryParticleIndex[kMaxRecoTracks];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueVisibleEnergy[kMaxRecoTracks];   //[RecoTrackN]
   Int_t           RecoTrackPrimaryParticleTrueNHits[kMaxRecoTracks];   //[RecoTrackN]
   Int_t           RecoTrackSecondaryParticleIndex[kMaxRecoTracks];   //[RecoTrackN]
   Float_t         RecoTrackSecondaryParticleTrueVisibleEnergy[kMaxRecoTracks];   //[RecoTrackN]
   Int_t           RecoTrackSecondaryParticleTrueNHits[kMaxRecoTracks];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueMomentumTrackStart[kMaxRecoTracks][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTruePositionTrackStart[kMaxRecoTracks][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueMomentumTrackEnd[kMaxRecoTracks][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTruePositionTrackEnd[kMaxRecoTracks][4];   //[RecoTrackN]
   Int_t           RecoTrackNHits[kMaxRecoTracks];   //[RecoTrackN]
   Float_t         RecoTrackTrueHitPosition[kMaxRecoTracks][kMaxRecoTrackHits][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueTrackLengthAsMeasured[kMaxRecoTracks];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueTrackLengthAsMeasuredIgnoreY[kMaxRecoTracks];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueTrackLengthRecoStart[kMaxRecoTracks];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueTrackLengthRecoStartIgnoreY[kMaxRecoTracks];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueTrackLengthInTMSIgnoreY[kMaxRecoTracks];   //[RecoTrackN]
   Int_t           RecoTrackPrimaryParticlePDG[kMaxRecoTracks];   //[RecoTrackN]
   Bool_t          RecoTrackPrimaryParticleIsPrimary[kMaxRecoTracks];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueMomentum[kMaxRecoTracks][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTruePositionStart[kMaxRecoTracks][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTruePositionEnd[kMaxRecoTracks][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueTrackLength[kMaxRecoTracks];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueTrackLengthIgnoreY[kMaxRecoTracks];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueTrackLengthInTMS[kMaxRecoTracks];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueMomentumEnteringTMS[kMaxRecoTracks][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueMomentumLeavingTMS[kMaxRecoTracks][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTruePositionEnteringTMS[kMaxRecoTracks][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTruePositionLeavingTMS[kMaxRecoTracks][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTruePositionLeavingLAr[kMaxRecoTracks][4];   //[RecoTrackN]
   Float_t         RecoTrackPrimaryParticleTrueMomentumLeavingLAr[kMaxRecoTracks][4];   //[RecoTrackN]
   Bool_t          RecoTrackPrimaryParticleTMSFiducialStart[kMaxRecoTracks];   //[RecoTrackN]
   Bool_t          RecoTrackPrimaryParticleTMSFiducialTouch[kMaxRecoTracks];   //[RecoTrackN]
   Bool_t          RecoTrackPrimaryParticleTMSFiducialEnd[kMaxRecoTracks];   //[RecoTrackN]
   Bool_t          RecoTrackPrimaryParticleLArFiducialStart[kMaxRecoTracks];   //[RecoTrackN]
   Bool_t          RecoTrackPrimaryParticleLArFiducialTouch[kMaxRecoTracks];   //[RecoTrackN]
   Bool_t          RecoTrackPrimaryParticleLArFiducialEnd[kMaxRecoTracks];   //[RecoTrackN]
   Int_t           RecoTrackPrimaryParticleVtxId[kMaxRecoTracks];   //[RecoTrackN]
   Bool_t          RecoTrackPrimaryParticleVtxFiducialCut[kMaxRecoTracks];   //[RecoTrackN]
   Bool_t          RecoTrackPrimaryParticleVtxShellEnergyCut[kMaxRecoTracks];   //[RecoTrackN]
   Bool_t          RecoTrackPrimaryParticleVtxNDPhysicsCut[kMaxRecoTracks];   //[RecoTrackN]
   Int_t           RecoTrackSecondaryParticlePDG[kMaxRecoTracks];   //[RecoTrackN]
   Bool_t          RecoTrackSecondaryParticleIsPrimary[kMaxRecoTracks];   //[RecoTrackN]
   Float_t         RecoTrackSecondaryParticleTrueMomentum[kMaxRecoTracks][4];   //[RecoTrackN]
   Float_t         RecoTrackSecondaryParticleTruePositionStart[kMaxRecoTracks][4];   //[RecoTrackN]
   Float_t         RecoTrackSecondaryParticleTruePositionEnd[kMaxRecoTracks][4];   //[RecoTrackN]
   Int_t           NTrueHits;
   Float_t         TrueHitX[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueHitY[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueHitZ[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueHitT[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueHitE[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueHitPE[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueHitPEAfterFibers[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueHitPEAfterFibersLongPath[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueHitPEAfterFibersShortPath[kMaxTrueHits];   //[NTrueHits]
   Int_t           TrueNTrueParticles[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueLeptonicEnergy[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueHadronicEnergy[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueRecoHitX[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueRecoHitY[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueRecoHitZ[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueRecoHitTrackX[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueRecoHitTrackY[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueRecoHitTrackXUncertainty[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueRecoHitTrackYUncertainty[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueRecoHitNotZ[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueRecoHitT[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueRecoHitE[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueRecoHitPE[kMaxTrueHits];   //[NTrueHits]
   Float_t         TrueRecoHitEVis[kMaxTrueHits];   //[NTrueHits]
   Bool_t          TrueRecoHitIsPedSupped[kMaxTrueHits];   //[NTrueHits]
   Int_t           TrueHitBar[kMaxTrueHits];   //[NTrueHits]
   Int_t           TrueHitPlane[kMaxTrueHits];   //[NTrueHits]
   Int_t           TrueHitView[kMaxTrueHits];   //[NTrueHits]

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
   TBranch        *b_nParticles;   //!
   TBranch        *b_LeptonPDG;   //!
   TBranch        *b_LeptonP4;   //!
   TBranch        *b_LeptonX4;   //!
   TBranch        *b_MuonP4;   //!
   TBranch        *b_Muon_Vertex;   //!
   TBranch        *b_Muon_Death;   //!
   TBranch        *b_Muon_TrueKE;   //!
   TBranch        *b_Muon_TrueTrackLength;   //!
   TBranch        *b_VertexIdOfMostEnergyInEvent;   //!
   TBranch        *b_VisibleEnergyFromVertexInSlice;   //!
   TBranch        *b_TotalVisibleEnergyFromVertex;   //!
   TBranch        *b_VisibleEnergyFromOtherVerticesInSlice;   //!
   TBranch        *b_VertexVisibleEnergyFractionInSlice;   //!
   TBranch        *b_PrimaryVertexVisibleEnergyFraction;   //!
   TBranch        *b_LArOuterShellEnergy;   //!
   TBranch        *b_LArTotalEnergy;   //!
   TBranch        *b_TotalNonTMSEnergy;   //!
   TBranch        *b_LArOuterShellEnergyFromVertex;   //!
   TBranch        *b_LArTotalEnergyFromVertex;   //!
   TBranch        *b_TotalNonTMSEnergyFromVertex;   //!
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
   TBranch        *b_RecoTrackPrimaryParticleTrueMomentumEnteringTMS;   //!
   TBranch        *b_RecoTrackPrimaryParticleTrueMomentumLeavingTMS;   //!
   TBranch        *b_RecoTrackPrimaryParticleTruePositionEnteringTMS;   //!
   TBranch        *b_RecoTrackPrimaryParticleTruePositionLeavingTMS;   //!
   TBranch        *b_RecoTrackPrimaryParticleTruePositionLeavingLAr;   //!
   TBranch        *b_RecoTrackPrimaryParticleTrueMomentumLeavingLAr;   //!
   TBranch        *b_RecoTrackPrimaryParticleTMSFiducialStart;   //!
   TBranch        *b_RecoTrackPrimaryParticleTMSFiducialTouch;   //!
   TBranch        *b_RecoTrackPrimaryParticleTMSFiducialEnd;   //!
   TBranch        *b_RecoTrackPrimaryParticleLArFiducialStart;   //!
   TBranch        *b_RecoTrackPrimaryParticleLArFiducialTouch;   //!
   TBranch        *b_RecoTrackPrimaryParticleLArFiducialEnd;   //!
   TBranch        *b_RecoTrackPrimaryParticleVtxId;   //!
   TBranch        *b_RecoTrackPrimaryParticleVtxFiducialCut;   //!
   TBranch        *b_RecoTrackPrimaryParticleVtxShellEnergyCut;   //!
   TBranch        *b_RecoTrackPrimaryParticleVtxNDPhysicsCut;   //!
   TBranch        *b_RecoTrackSecondaryParticlePDG;   //!
   TBranch        *b_RecoTrackSecondaryParticleIsPrimary;   //!
   TBranch        *b_RecoTrackSecondaryParticleTrueMomentum;   //!
   TBranch        *b_RecoTrackSecondaryParticleTruePositionStart;   //!
   TBranch        *b_RecoTrackSecondaryParticleTruePositionEnd;   //!
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
   TBranch        *b_TrueHitPlane;   //!
   TBranch        *b_TrueHitView;   //!

   Truth_Info(TTree *tree=0);
   virtual ~Truth_Info();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual Long64_t GetEntriesFast() { return fChain->GetEntriesFast(); };
   virtual bool HasBranch(const char *branch) {
      return fChain && fChain->GetBranch(branch);
   };
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Truth_Info_cxx
Truth_Info::Truth_Info(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   Init(tree);
}

Truth_Info::~Truth_Info()
{
   if (!fChain) return;
   // The TChain owns its files in validation jobs.
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
   fChain->SetBranchAddress("nParticles", &nParticles, &b_nParticles);
   fChain->SetBranchAddress("LeptonPDG", &LeptonPDG, &b_LeptonPDG);
   fChain->SetBranchAddress("LeptonP4", LeptonP4, &b_LeptonP4);
   fChain->SetBranchAddress("LeptonX4", LeptonX4, &b_LeptonX4);
   fChain->SetBranchAddress("MuonP4", MuonP4, &b_MuonP4);
   fChain->SetBranchAddress("Muon_Vertex", Muon_Vertex, &b_Muon_Vertex);
   fChain->SetBranchAddress("Muon_Death", Muon_Death, &b_Muon_Death);
   fChain->SetBranchAddress("Muon_TrueKE", &Muon_TrueKE, &b_Muon_TrueKE);
   fChain->SetBranchAddress("Muon_TrueTrackLength", &Muon_TrueTrackLength, &b_Muon_TrueTrackLength);
   fChain->SetBranchAddress("VertexIdOfMostEnergyInEvent", &VertexIdOfMostEnergyInEvent, &b_VertexIdOfMostEnergyInEvent);
   fChain->SetBranchAddress("VisibleEnergyFromVertexInSlice", &VisibleEnergyFromVertexInSlice, &b_VisibleEnergyFromVertexInSlice);
   fChain->SetBranchAddress("TotalVisibleEnergyFromVertex", &TotalVisibleEnergyFromVertex, &b_TotalVisibleEnergyFromVertex);
   fChain->SetBranchAddress("VisibleEnergyFromOtherVerticesInSlice", &VisibleEnergyFromOtherVerticesInSlice, &b_VisibleEnergyFromOtherVerticesInSlice);
   fChain->SetBranchAddress("VertexVisibleEnergyFractionInSlice", &VertexVisibleEnergyFractionInSlice, &b_VertexVisibleEnergyFractionInSlice);
   fChain->SetBranchAddress("PrimaryVertexVisibleEnergyFraction", &PrimaryVertexVisibleEnergyFraction, &b_PrimaryVertexVisibleEnergyFraction);
   fChain->SetBranchAddress("LArOuterShellEnergy", &LArOuterShellEnergy, &b_LArOuterShellEnergy);
   fChain->SetBranchAddress("LArTotalEnergy", &LArTotalEnergy, &b_LArTotalEnergy);
   fChain->SetBranchAddress("TotalNonTMSEnergy", &TotalNonTMSEnergy, &b_TotalNonTMSEnergy);
   fChain->SetBranchAddress("LArOuterShellEnergyFromVertex", &LArOuterShellEnergyFromVertex, &b_LArOuterShellEnergyFromVertex);
   fChain->SetBranchAddress("LArTotalEnergyFromVertex", &LArTotalEnergyFromVertex, &b_LArTotalEnergyFromVertex);
   fChain->SetBranchAddress("TotalNonTMSEnergyFromVertex", &TotalNonTMSEnergyFromVertex, &b_TotalNonTMSEnergyFromVertex);
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
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueMomentumEnteringTMS", RecoTrackPrimaryParticleTrueMomentumEnteringTMS, &b_RecoTrackPrimaryParticleTrueMomentumEnteringTMS);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueMomentumLeavingTMS", RecoTrackPrimaryParticleTrueMomentumLeavingTMS, &b_RecoTrackPrimaryParticleTrueMomentumLeavingTMS);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTruePositionEnteringTMS", RecoTrackPrimaryParticleTruePositionEnteringTMS, &b_RecoTrackPrimaryParticleTruePositionEnteringTMS);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTruePositionLeavingTMS", RecoTrackPrimaryParticleTruePositionLeavingTMS, &b_RecoTrackPrimaryParticleTruePositionLeavingTMS);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTruePositionLeavingLAr", RecoTrackPrimaryParticleTruePositionLeavingLAr, &b_RecoTrackPrimaryParticleTruePositionLeavingLAr);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTrueMomentumLeavingLAr", RecoTrackPrimaryParticleTrueMomentumLeavingLAr, &b_RecoTrackPrimaryParticleTrueMomentumLeavingLAr);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTMSFiducialStart", RecoTrackPrimaryParticleTMSFiducialStart, &b_RecoTrackPrimaryParticleTMSFiducialStart);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTMSFiducialTouch", RecoTrackPrimaryParticleTMSFiducialTouch, &b_RecoTrackPrimaryParticleTMSFiducialTouch);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleTMSFiducialEnd", RecoTrackPrimaryParticleTMSFiducialEnd, &b_RecoTrackPrimaryParticleTMSFiducialEnd);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleLArFiducialStart", RecoTrackPrimaryParticleLArFiducialStart, &b_RecoTrackPrimaryParticleLArFiducialStart);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleLArFiducialTouch", RecoTrackPrimaryParticleLArFiducialTouch, &b_RecoTrackPrimaryParticleLArFiducialTouch);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleLArFiducialEnd", RecoTrackPrimaryParticleLArFiducialEnd, &b_RecoTrackPrimaryParticleLArFiducialEnd);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleVtxId", RecoTrackPrimaryParticleVtxId, &b_RecoTrackPrimaryParticleVtxId);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleVtxFiducialCut", RecoTrackPrimaryParticleVtxFiducialCut, &b_RecoTrackPrimaryParticleVtxFiducialCut);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleVtxShellEnergyCut", RecoTrackPrimaryParticleVtxShellEnergyCut, &b_RecoTrackPrimaryParticleVtxShellEnergyCut);
   fChain->SetBranchAddress("RecoTrackPrimaryParticleVtxNDPhysicsCut", RecoTrackPrimaryParticleVtxNDPhysicsCut, &b_RecoTrackPrimaryParticleVtxNDPhysicsCut);
   fChain->SetBranchAddress("RecoTrackSecondaryParticlePDG", RecoTrackSecondaryParticlePDG, &b_RecoTrackSecondaryParticlePDG);
   fChain->SetBranchAddress("RecoTrackSecondaryParticleIsPrimary", RecoTrackSecondaryParticleIsPrimary, &b_RecoTrackSecondaryParticleIsPrimary);
   fChain->SetBranchAddress("RecoTrackSecondaryParticleTrueMomentum", RecoTrackSecondaryParticleTrueMomentum, &b_RecoTrackSecondaryParticleTrueMomentum);
   fChain->SetBranchAddress("RecoTrackSecondaryParticleTruePositionStart", RecoTrackSecondaryParticleTruePositionStart, &b_RecoTrackSecondaryParticleTruePositionStart);
   fChain->SetBranchAddress("RecoTrackSecondaryParticleTruePositionEnd", RecoTrackSecondaryParticleTruePositionEnd, &b_RecoTrackSecondaryParticleTruePositionEnd);
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
   fChain->SetBranchAddress("TrueHitPlane", TrueHitPlane, &b_TrueHitPlane);
   fChain->SetBranchAddress("TrueHitView", TrueHitView, &b_TrueHitView);
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
void Truth_Info::Loop()
{
}
Int_t Truth_Info::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Truth_Info_cxx
