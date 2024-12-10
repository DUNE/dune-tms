//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May  9 20:29:33 2024 by ROOT version 6.22/08
// from TTree Reco_Tree/Reco_Tree
// found on file: /pnfs/dune/persistent/users/kleykamp/tmsreco_combined_files/2024-04-19_bfield_0p0T.tmsreco.root
//////////////////////////////////////////////////////////

#ifndef Reco_Tree_h
#define Reco_Tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Reco_Tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
    #define __TMS_MAX_TRACKS__ 100 

   // Declaration of leaf types
   Int_t           EventNo;
   Int_t           SliceNo;
   Int_t           SpillNo;
   Int_t           nTracks;
   Int_t           nHits[__TMS_MAX_TRACKS__];   //[nTracks]
   Float_t         TrackHitPos[__TMS_MAX_TRACKS__][200][3];   //[nTracks]
   Int_t           nKalmanNodes[__TMS_MAX_TRACKS__];   //[nTracks]
   Float_t         KalmanPos[__TMS_MAX_TRACKS__][200][3];   //[nTracks]
   Float_t         KalmanTruePos[__TMS_MAX_TRACKS__][200][3];   //[nTracks]
   Float_t         StartDirection[__TMS_MAX_TRACKS__][3];   //[nTracks]
   Float_t         Direction[__TMS_MAX_TRACKS__][3];   //[nTracks]
   Float_t         EndDirection[__TMS_MAX_TRACKS__][3];   //[nTracks]
   Float_t         StartPos[__TMS_MAX_TRACKS__][3];   //[nTracks]
   Float_t         EndPos[__TMS_MAX_TRACKS__][3];   //[nTracks]
   Float_t         EnergyRange[__TMS_MAX_TRACKS__];   //[nTracks]
   Float_t         EnergyDeposit[__TMS_MAX_TRACKS__];   //[nTracks]
   Float_t         Momentum[__TMS_MAX_TRACKS__];   //[nTracks]
   Float_t         Length[__TMS_MAX_TRACKS__];   //[nTracks]
   Int_t           Charge[__TMS_MAX_TRACKS__];   //[nTracks]
   Float_t         TrackHitEnergies[__TMS_MAX_TRACKS__][200];

   // List of branches
   TBranch        *b_EventNo;   //!
   TBranch        *b_SliceNo;   //!
   TBranch        *b_SpillNo;   //!
   TBranch        *b_nTracks;   //!
   TBranch        *b_nHits;   //!
   TBranch        *b_TrackHitPos;   //!
   TBranch        *b_nKalmanNodes;   //!
   TBranch        *b_KalmanPos;   //!
   TBranch        *b_KalmanTruePos;   //!
   TBranch        *b_StartDirection;   //!
   TBranch        *b_Direction;   //!
   TBranch        *b_EndDirection;   //!
   TBranch        *b_StartPos;   //!
   TBranch        *b_EndPos;   //!
   TBranch        *b_EnergyRange;   //!
   TBranch        *b_EnergyDeposit;   //!
   TBranch        *b_Momentum;   //!
   TBranch        *b_Length;   //!
   TBranch        *b_Charge;   //!
   TBranch        *b_RecoTrackHitEnergies;   //!

   Reco_Tree(TTree *tree=0);
   virtual ~Reco_Tree();
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

#define Reco_Tree_cxx
#ifdef Reco_Tree_cxx
Reco_Tree::Reco_Tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/pnfs/dune/persistent/users/kleykamp/tmsreco_combined_files/2024-04-19_bfield_0p0T.tmsreco.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/pnfs/dune/persistent/users/kleykamp/tmsreco_combined_files/2024-04-19_bfield_0p0T.tmsreco.root");
      }
      f->GetObject("Reco_Tree",tree);

   }
   Init(tree);
}

Reco_Tree::~Reco_Tree()
{
   if (!fChain) return;
   //delete fChain->GetCurrentFile();
}

Int_t Reco_Tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Reco_Tree::LoadTree(Long64_t entry)
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

void Reco_Tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EventNo", &EventNo, &b_EventNo);
   fChain->SetBranchAddress("SliceNo", &SliceNo, &b_SliceNo);
   fChain->SetBranchAddress("SpillNo", &SpillNo, &b_SpillNo);
   fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
   fChain->SetBranchAddress("nHits", nHits, &b_nHits);
   fChain->SetBranchAddress("TrackHitPos", TrackHitPos, &b_TrackHitPos);
   fChain->SetBranchAddress("nKalmanNodes", nKalmanNodes, &b_nKalmanNodes);
   fChain->SetBranchAddress("KalmanPos", KalmanPos, &b_KalmanPos);
   fChain->SetBranchAddress("KalmanTruePos", KalmanTruePos, &b_KalmanTruePos);
   // Used to be Direction, now is StartDirection, check for both options depending on the file
   //if (HasBranch("Direction")) fChain->SetBranchAddress("Direction", Direction, &b_Direction);
   //else fChain->SetBranchAddress("StartDirection", Direction, &b_Direction);
   fChain->SetBranchAddress("StartDirection", StartDirection, &b_StartDirection);
   fChain->SetBranchAddress("EndDirection", EndDirection, &b_EndDirection);
   fChain->SetBranchAddress("StartPos", StartPos, &b_StartPos);
   fChain->SetBranchAddress("EndPos", EndPos, &b_EndPos);
   fChain->SetBranchAddress("EnergyRange", EnergyRange, &b_EnergyRange);
   fChain->SetBranchAddress("EnergyDeposit", EnergyDeposit, &b_EnergyDeposit);
   fChain->SetBranchAddress("Momentum", Momentum, &b_Momentum);
   fChain->SetBranchAddress("Length", Length, &b_Length);
   fChain->SetBranchAddress("Charge", Charge, &b_Charge);
   fChain->SetBranchAddress("TrackHitEnergies", TrackHitEnergies, &b_RecoTrackHitEnergies);
   Notify();
}

Bool_t Reco_Tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Reco_Tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Reco_Tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   entry += 1;
   return 1;
}

void Reco_Tree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Reco_Tree.C
//      root> Reco_Tree t
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
#endif // #ifdef Reco_Tree_cxx
