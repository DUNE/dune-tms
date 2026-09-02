//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 30 11:57:24 2026 by ROOT version 6.28/12
// from TTree Truth_Spill/Truth_Spill
// found on file: /exp/dune/app/users/kleykamp/2026_dune-tms/dune-tms/MicroProdN4p1_NDComplex_FHC.spill.full.0002459.EDEPSIM_SPILLS_TMS_RecoCandidates_Hough_Cluster1.root
//////////////////////////////////////////////////////////

#ifndef Truth_Spill_h
#define Truth_Spill_h
#define Truth_Spill_cxx

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "string"
#include "vector"

using std::string;
using std::vector;

class Truth_Spill {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   enum {
      kMaxTrueParticles = 50000,
      kMaxTrueVertices = 5000
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
   Int_t           nTrueParticles;
   Int_t           nTruePrimaryParticles;
   Int_t           nTrueForgottenParticles;
   Int_t           VertexID[kMaxTrueParticles];   //[nTrueParticles]
   Int_t           Parent[kMaxTrueParticles];   //[nTrueParticles]
   Int_t           TrackId[kMaxTrueParticles];   //[nTrueParticles]
   Int_t           PDG[kMaxTrueParticles];   //[nTrueParticles]
   Bool_t          IsPrimary[kMaxTrueParticles];   //[nTrueParticles]
   Float_t         TrueVisibleEnergy[kMaxTrueParticles];   //[nTrueParticles]
   Int_t           TrueNHits[kMaxTrueParticles];   //[nTrueParticles]
   Float_t         TruePathLength[kMaxTrueParticles];   //[nTrueParticles]
   Float_t         TruePathLengthIgnoreY[kMaxTrueParticles];   //[nTrueParticles]
   Float_t         TruePathLengthInTMS[kMaxTrueParticles];   //[nTrueParticles]
   Float_t         TruePathLengthInTMSIgnoreY[kMaxTrueParticles];   //[nTrueParticles]
   Bool_t          InteractionTMSFiducial;
   Bool_t          InteractionTMSFirstTwoModules;
   Bool_t          InteractionTMSThin;
   Bool_t          InteractionLArFiducial;
   Bool_t          TMSFiducialStart[kMaxTrueParticles];   //[nTrueParticles]
   Bool_t          TMSFiducialTouch[kMaxTrueParticles];   //[nTrueParticles]
   Bool_t          TMSFiducialEnd[kMaxTrueParticles];   //[nTrueParticles]
   Bool_t          LArFiducialStart[kMaxTrueParticles];   //[nTrueParticles]
   Bool_t          LArFiducialTouch[kMaxTrueParticles];   //[nTrueParticles]
   Bool_t          LArFiducialEnd[kMaxTrueParticles];   //[nTrueParticles]
   Float_t         BirthMomentum[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         BirthPosition[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         DeathMomentum[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         DeathPosition[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         MomentumLArStart[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         PositionLArStart[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         MomentumLArEnd[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         PositionLArEnd[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         MomentumTMSStart[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         PositionTMSStart[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         MomentumTMSFirstTwoModulesEnd[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         PositionTMSFirstTwoModulesEnd[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         MomentumTMSThinEnd[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         PositionTMSThinEnd[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         MomentumTMSEnd[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         PositionTMSEnd[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         MomentumZIsLArEnd[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         PositionZIsLArEnd[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         MomentumZIsTMSStart[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         PositionZIsTMSStart[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         MomentumZIsTMSEnd[kMaxTrueParticles][4];   //[nTrueParticles]
   Float_t         PositionZIsTMSEnd[kMaxTrueParticles][4];   //[nTrueParticles]
   Int_t           TrueVtxN;
   Float_t         TrueVtxX[kMaxTrueVertices];   //[TrueVtxN]
   Float_t         TrueVtxY[kMaxTrueVertices];   //[TrueVtxN]
   Float_t         TrueVtxZ[kMaxTrueVertices];   //[TrueVtxN]
   Float_t         TrueVtxT[kMaxTrueVertices];   //[TrueVtxN]
   Float_t         TrueVtxPx[kMaxTrueVertices];   //[TrueVtxN]
   Float_t         TrueVtxPy[kMaxTrueVertices];   //[TrueVtxN]
   Float_t         TrueVtxPz[kMaxTrueVertices];   //[TrueVtxN]
   Float_t         TrueVtxE[kMaxTrueVertices];   //[TrueVtxN]
   Int_t           TrueVtxPDG[kMaxTrueVertices];   //[TrueVtxN]
   Int_t           TrueVtxID[kMaxTrueVertices];   //[TrueVtxN]
   vector<string>  *TrueVtxReaction;
   Float_t         TrueVtxHadronicELarShell[kMaxTrueVertices];   //[TrueVtxN]
   Float_t         TrueVtxHadronicELAr[kMaxTrueVertices];   //[TrueVtxN]
   Float_t         TrueVtxHadronicETMS[kMaxTrueVertices];   //[TrueVtxN]
   Float_t         TrueVtxHadronicE[kMaxTrueVertices];   //[TrueVtxN]
   Float_t         TrueVtxVisibleETMS[kMaxTrueVertices];   //[TrueVtxN]
   Float_t         TrueVtxVisibleELAr[kMaxTrueVertices];   //[TrueVtxN]
   Float_t         TrueVtxVisibleE[kMaxTrueVertices];   //[TrueVtxN]
   Bool_t          TrueVtxFiducialCut[kMaxTrueVertices];   //[TrueVtxN]
   Bool_t          TrueVtxShellEnergyCut[kMaxTrueVertices];   //[TrueVtxN]
   Bool_t          TrueVtxNDPhysicsCut[kMaxTrueVertices];   //[TrueVtxN]

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
   TBranch        *b_TrueVtxN;   //!
   TBranch        *b_TrueVtxX;   //!
   TBranch        *b_TrueVtxY;   //!
   TBranch        *b_TrueVtxZ;   //!
   TBranch        *b_TrueVtxT;   //!
   TBranch        *b_TrueVtxPx;   //!
   TBranch        *b_TrueVtxPy;   //!
   TBranch        *b_TrueVtxPz;   //!
   TBranch        *b_TrueVtxE;   //!
   TBranch        *b_TrueVtxPDG;   //!
   TBranch        *b_TrueVtxID;   //!
   TBranch        *b_TrueVtxReaction;   //!
   TBranch        *b_TrueVtxHadronicELarShell;   //!
   TBranch        *b_TrueVtxHadronicELAr;   //!
   TBranch        *b_TrueVtxHadronicETMS;   //!
   TBranch        *b_TrueVtxHadronicE;   //!
   TBranch        *b_TrueVtxVisibleETMS;   //!
   TBranch        *b_TrueVtxVisibleELAr;   //!
   TBranch        *b_TrueVtxVisibleE;   //!
   TBranch        *b_TrueVtxFiducialCut;   //!
   TBranch        *b_TrueVtxShellEnergyCut;   //!
   TBranch        *b_TrueVtxNDPhysicsCut;   //!

   Truth_Spill(TTree *tree=0);
   virtual ~Truth_Spill();
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

#ifdef Truth_Spill_cxx
Truth_Spill::Truth_Spill(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   Init(tree);
}

Truth_Spill::~Truth_Spill()
{
   if (!fChain) return;
   // The TChain owns its files in validation jobs.
}

Int_t Truth_Spill::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Truth_Spill::LoadTree(Long64_t entry)
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

void Truth_Spill::Init(TTree *tree)
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
   TrueVtxReaction = 0;
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
   fChain->SetBranchAddress("TrueVtxN", &TrueVtxN, &b_TrueVtxN);
   fChain->SetBranchAddress("TrueVtxX", TrueVtxX, &b_TrueVtxX);
   fChain->SetBranchAddress("TrueVtxY", TrueVtxY, &b_TrueVtxY);
   fChain->SetBranchAddress("TrueVtxZ", TrueVtxZ, &b_TrueVtxZ);
   fChain->SetBranchAddress("TrueVtxT", TrueVtxT, &b_TrueVtxT);
   fChain->SetBranchAddress("TrueVtxPx", TrueVtxPx, &b_TrueVtxPx);
   fChain->SetBranchAddress("TrueVtxPy", TrueVtxPy, &b_TrueVtxPy);
   fChain->SetBranchAddress("TrueVtxPz", TrueVtxPz, &b_TrueVtxPz);
   fChain->SetBranchAddress("TrueVtxE", TrueVtxE, &b_TrueVtxE);
   fChain->SetBranchAddress("TrueVtxPDG", TrueVtxPDG, &b_TrueVtxPDG);
   fChain->SetBranchAddress("TrueVtxID", TrueVtxID, &b_TrueVtxID);
   fChain->SetBranchAddress("TrueVtxReaction", &TrueVtxReaction, &b_TrueVtxReaction);
   fChain->SetBranchAddress("TrueVtxHadronicELarShell", TrueVtxHadronicELarShell, &b_TrueVtxHadronicELarShell);
   fChain->SetBranchAddress("TrueVtxHadronicELAr", TrueVtxHadronicELAr, &b_TrueVtxHadronicELAr);
   fChain->SetBranchAddress("TrueVtxHadronicETMS", TrueVtxHadronicETMS, &b_TrueVtxHadronicETMS);
   fChain->SetBranchAddress("TrueVtxHadronicE", TrueVtxHadronicE, &b_TrueVtxHadronicE);
   fChain->SetBranchAddress("TrueVtxVisibleETMS", TrueVtxVisibleETMS, &b_TrueVtxVisibleETMS);
   fChain->SetBranchAddress("TrueVtxVisibleELAr", TrueVtxVisibleELAr, &b_TrueVtxVisibleELAr);
   fChain->SetBranchAddress("TrueVtxVisibleE", TrueVtxVisibleE, &b_TrueVtxVisibleE);
   fChain->SetBranchAddress("TrueVtxFiducialCut", TrueVtxFiducialCut, &b_TrueVtxFiducialCut);
   fChain->SetBranchAddress("TrueVtxShellEnergyCut", TrueVtxShellEnergyCut, &b_TrueVtxShellEnergyCut);
   fChain->SetBranchAddress("TrueVtxNDPhysicsCut", TrueVtxNDPhysicsCut, &b_TrueVtxNDPhysicsCut);
   Notify();
}

Bool_t Truth_Spill::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Truth_Spill::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
void Truth_Spill::Loop()
{
}
Int_t Truth_Spill::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Truth_Spill_cxx
