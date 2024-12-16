//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 28 22:24:07 2024 by ROOT version 6.22/08
// from TTree Line_Candidates/Line_Candidates
// found on file: neutrino.0_1698176146.edep_TMS_RecoCandidates_Hough_Cluster1.root
//////////////////////////////////////////////////////////

#ifndef Line_Candidates_h
#define Line_Candidates_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Line_Candidates {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           EventNo;
   Int_t           SliceNo;
   Int_t           SpillNo;
   Int_t           nLinesU;
   Int_t           nLinesV;
   Int_t           nLinesX;
   Int_t           nLines3D;
   Float_t         SlopeU[20];   //[nLinesU]
   Float_t         InterceptU[20];   //[nLinesU]
   Float_t         Slope_DownstreamU[20];   //[nLinesU]
   Float_t         Intercept_DownstreamU[20];   //[nLinesU]
   Float_t         Slope_UpstreamU[20];   //[nLinesU]
   Float_t         Intercept_UpstreamU[20];   //[nLinesU]
   Float_t         SlopeV[20];   //[nLinesV]
   Float_t         InterceptV[20];   //[nLinesV]
   Float_t         Slope_DownstreamV[20];   //[nLinesV]
   Float_t         Intercept_DownstreamV[20];   //[nLinesV]
   Float_t         Slope_UpstreamV[20];   //[nLinesV]
   Float_t         Intercept_UpstreamV[20];   //[nLinesV]
   Float_t         SlopeX[20];   //[nLinesX]
   Float_t         InterceptX[20];   //[nLinesX]
   Float_t         Slope_DownstreamX[20];   //[nLinesX]
   Float_t         Intercept_DownstreamX[20];   //[nLinesX]
   Float_t         Slope_UpstreamX[20];   //[nLinesX]
   Float_t         Intercept_UpstreamX[20];   //[nLinesX]
   /*Float_t         DirectionZU[20];   //[nLinesU]
   Float_t         DirectionZV[20];   //[nLinesV]
   Float_t         DirectionZX[20];   //[nLinesX]
   Float_t         DirectionXU[20];   //[nLinesU]
   Float_t         DirectionXV[20];   //[nLinesV]
   Float_t         DirectionYX[20];   //[nLinesX]
   Float_t         DirectionZU_Downstream[20];   //[nLinesU]
   Float_t         DirectionXU_Downstream[20];   //[nLinesU]
   Float_t         DirectionZU_Upstream[20];   //[nLinesU]
   Float_t         DirectionXU_Upstream[20];   //[nLinesU]
   Float_t         DirectionZV_Downstream[20];   //[nLinesV]
   Float_t         DirectionXV_Downstream[20];   //[nLinesV]
   Float_t         DirectionZV_Upstream[20];   //[nLinesV]
   Float_t         DirectionXV_Upstream[20];   //[nLinesV]
   Float_t         DirectionZX_Downstream[20];   //[nLinesX]
   Float_t         DirectionYX_Downstream[20];   //[nLinesX]
   Float_t         DirectionZX_Upstream[20];   //[nLinesX]
   Float_t         DirectionYX_Upstream[20];   //[nLinesX]
   Float_t         FirstHoughHitU[25][2];   //[nLinesU]
   Float_t         FirstHoughHitV[25][2];   //[nLinesV]
   Float_t         FirstHoughHitX[25][2];   //[nLinesX]
   Float_t         LastHoughHitU[25][2];   //[nLinesU]
   Float_t         LastHoughHitV[25][2];   //[nLinesV]
   Float_t         LastHoughHitX[25][2];   //[nLinesX]
   Float_t         FirstHoughHitTimeU[20];   //[nLinesU]
   Float_t         FirstHoughHitTimeV[20];   //[nLinesV]
   Float_t         FirstHoughHitTimeX[20];   //[nLinesX]
   Float_t         LastHoughHitTimeU[20];   //[nLinesU]
   Float_t         LastHoughHitTimeV[20];   //[nLinesV]
   Float_t         LastHoughHitTImeX[20];   //[nLinesX]
   Float_t         HoughEarliestHitTimeU[20];   //[nLinesU]
   Float_t         HoughEarliestHitTimeV[20];   //[nLinesV]
   Float_t         HoughEarliestHitTimeX[20];   //[nLinesX]
   Float_t         HoughLatestHitTimeU[20];   //[nLinesU]
   Float_t         HoughLatestHitTimeV[20];   //[nLinesV]
   Float_t         HoughLatestHitTimeX[20];   //[nLinesX]
   Int_t           FirstHoughPlaneU[20];   //[nLinesU]
   Int_t           FirstHoughPlaneV[20];   //[nLinesV]
   Int_t           FirstHoughPlaneX[20];   //[nLinesX]
   Int_t           LastHoughPlaneU[20];   //[nLinesU]
   Int_t           LastHoughPlaneV[20];   //[nLinesV]
   Int_t           LastHoughPlaneX[20];   //[nLinesX] */
   Bool_t          TMSStart;
   Float_t         TMSStartTime;
   Float_t         OccupancyU[20];   //[nLinesU]
   Float_t         OccupancyV[20];   //[nLinesV]
   Float_t         OccupancyX[20];   //[nLinesX] 
   Float_t         TrackLengthU[20];   //[nLinesU]
   Float_t         TrackLengthV[20];   //[nLinesV]
   Float_t         TrackLengthX[20];   //[nLinesX]
   Float_t         TotalTrackEnergyU[20];   //[nLinesU]
   Float_t         TotalTrackEnergyV[20];   //[nLinesV]
   Float_t         TotalTrackEnergyX[20];   //[nLinesX]
   Bool_t          TrackStoppingU[20];   //[nLinesU]
   Bool_t          TrackStoppingV[20];   //[nLinesV]
   Bool_t          TrackStoppingX[20];   //[nLinesX]
   Int_t           nHitsInTrackU[20];   //[nLinesU]
   Int_t           nHitsInTrackV[20];   //[nLinesV]
   Int_t           nHitsInTrackX[20];   //[nLinesX]
   Float_t         TrackHitEnergyU[10][500];
   Float_t         TrackHitEnergyV[10][500];
   Float_t         TrackHitEnergyX[10][500]; 
   Float_t         TrackHitPosU[10][500][2]; 
   Float_t         TrackHitPosV[10][500][2]; 
   Float_t         TrackHitPosX[10][500][2]; 
   Float_t         TrackHitTimeU[10][500];
   Float_t         TrackHitTimeV[10][500];
   Float_t         TrackHitTimeX[10][500]; 
   /*Int_t           nClustersU;
   Int_t           nClustersV;
   Int_t           nClusterX;
   Float_t         ClusterEnergyU[20];   //[nClustersU]
   Float_t         ClusterEnergyV[20];   //[nClustersV]
   Float_t         ClusterEnergyX[20];   //[nClustersX]
   Float_t         ClusterTimeU[20];   //[nClustersU]
   Float_t         ClusterTimeV[20];   //[nClustersV]
   Float_t         ClusterTimeX[20];   //[nClustersX]
   Float_t         ClusterPosMeanU[25][2];
   Float_t         ClusterPosMeanV[25][2];
   Float_t         ClusterPosMeanX[25][2];
   Float_t         ClusterPosStdDevU[25][2];
   Float_t         ClusterPosStdDevV[25][2];
   Float_t         ClusterPosStdDevX[25][2];
   Int_t           nHitsInClusterU[20];   //[nClustersU]
   Int_t           nHitsInClusterV[20];   //[nClustersV]
   Int_t           nHitsInClusterX[20];   //[nClustersX]
   Float_t         ClusterHitPosU[25][200][2];
   Float_t         ClusterHitPosV[25][200][2];
   Float_t         ClusterHitPosX[25][200][2];
   Float_t         ClusterHitEnergyU[25][200];
   Float_t         ClusterHitEnergyV[25][200];
   Float_t         ClusterHitEnergyX[25][200];
   Float_t         ClusterHitTimeU[25][200];
   Float_t         ClusterHitTimeV[25][200];
   Float_t         ClusterHitTimeX[25][200];
   Int_t           ClusterHitSliceU[25][200];
   Int_t           ClusterHitSliceV[25][200];
   Int_t           ClusterHitSliceX[25][200]; */
   Int_t           nHits;
   Float_t         RecoHitPos[2000][4];   //[nHits]
   Float_t         RecoHitEnergy[2000];   //[nHits]
   Int_t           RecoHitSlice[2000];   //[nHits] 
   Float_t         RecoHitPE[2000];   //[nHits]
   Int_t           RecoHitBar[2000];   //[nHits]
   Int_t           RecoHitPlane[2000];   //[nHits]

   // List of branches
   TBranch        *b_EventNo;   //!
   TBranch        *b_SliceNo;   //!
   TBranch        *b_SpillNo;   //!
   TBranch        *b_nLinesU;   //!
   TBranch        *b_nLinesV;   //!
   TBranch        *b_nLinesX;   //!
   TBranch        *b_nLines3D;   //!
   TBranch        *b_SlopeU;   //!
   TBranch        *b_InterceptU;   //!
   TBranch        *b_Slope_DownstreamU;   //!
   TBranch        *b_Intercept_DownstreamU;   //!
   TBranch        *b_Slope_UpstreamU;   //!
   TBranch        *b_Intercept_UpstreamU;   //!
   TBranch        *b_SlopeV;   //!
   TBranch        *b_InterceptV;   //!
   TBranch        *b_Slope_DownstreamV;   //!
   TBranch        *b_Intercept_DownstreamV;   //!
   TBranch        *b_Slope_UpstreamV;   //!
   TBranch        *b_Intercept_UpstreamV;   //!
   TBranch        *b_SlopeX;   //!
   TBranch        *b_InterceptX;   //!
   TBranch        *b_Slope_DownstreamX;   //!
   TBranch        *b_Intercept_DownstreamX;   //!
   TBranch        *b_Slope_UpstreamX;   //!
   TBranch        *b_Intercept_UpstreamX;   //!
   /*TBranch        *b_DirectionZU;   //!
   TBranch        *b_DirectionZV;   //!
   TBranch        *b_DirectionZX;   //!
   TBranch        *b_DirectionXU;   //!
   TBranch        *b_DirectionXV;   //!
   TBranch        *b_DirectionYX;   //!
   TBranch        *b_DirectionZU_Downstream;   //!
   TBranch        *b_DirectionXU_Downstream;   //!
   TBranch        *b_DirectionZU_Upstream;   //!
   TBranch        *b_DirectionXU_Upstream;   //!
   TBranch        *b_DirectionZV_Downstream;   //!
   TBranch        *b_DirectionXV_Downstream;   //!
   TBranch        *b_DirectionZV_Upstream;   //!
   TBranch        *b_DirectionXV_Upstream;   //!
   TBranch        *b_DirectionZX_Downstream;   //!
   TBranch        *b_DirectionYX_Downstream;   //!
   TBranch        *b_DirectionZX_Upstream;   //!
   TBranch        *b_DirectionYX_Upstream;   //!
   TBranch        *b_FirstHoughHitU;   //!
   TBranch        *b_FirstHoughHitV;   //!
   TBranch        *b_FirstHoughHitX;   //!
   TBranch        *b_LastHoughHitU;   //!
   TBranch        *b_LastHoughHitV;   //!
   TBranch        *b_LastHoughHitX;   //!
   TBranch        *b_FirstHoughHitTimeU;   //!
   TBranch        *b_FirstHoughHitTimeV;   //!
   TBranch        *b_FirstHoughHitTimeX;   //!
   TBranch        *b_LastHoughHitTimeU;   //!
   TBranch        *b_LastHoughHitTimeV;   //!
   TBranch        *b_LastHoughHitTImeX;   //!
   TBranch        *b_HoughEarliestHitTimeU;   //!
   TBranch        *b_HoughEarliestHitTimeV;   //!
   TBranch        *b_HoughEarliestHitTimeX;   //!
   TBranch        *b_HoughLatestHitTimeU;   //!
   TBranch        *b_HoughLatestHitTimeV;   //!
   TBranch        *b_HoughLatestHitTimeX;   //!
   TBranch        *b_FirstHoughPlaneU;   //!
   TBranch        *b_FirstHoughPlaneV;   //!
   TBranch        *b_FirstHoughPlaneX;   //!
   TBranch        *b_LastHoughPlaneU;   //!
   TBranch        *b_LastHoughPlaneV;   //!
   TBranch        *b_LastHoughPlaneX;   //! */
   TBranch        *b_TMSStart;   //!
   TBranch        *b_TMSStartTime;   //!
   TBranch        *b_OccupancyU;   //!
   TBranch        *b_OccupancyV;   //!
   TBranch        *b_OccupancyX;   //! 
   TBranch        *b_TrackLengthU;   //!
   TBranch        *b_TrackLengthV;   //!
   TBranch        *b_TrackLengthX;   //!
   TBranch        *b_TotalTrackEnergyU;   //!
   TBranch        *b_TotalTrackEnergyV;   //!
   TBranch        *b_TotalTrackEnergyX;   //!
   TBranch        *b_TrackStoppingU;   //!
   TBranch        *b_TrackStoppingV;   //!
   TBranch        *b_TrackStoppingX;   //!
   TBranch        *b_nHitsInTrackU;   //!
   TBranch        *b_nHitsInTrackV;   //!
   TBranch        *b_nHitsInTrackX;   //!
   TBranch        *b_TrackHitEnergyU;   //!
   TBranch        *b_TrackHitEnergyV;   //!
   TBranch        *b_TrackHitEnergyX;   //!
   TBranch        *b_TrackHitPosU;   //!
   TBranch        *b_TrackHitPosV;   //!
   TBranch        *b_TrackHitPosX;   //!
   TBranch        *b_TrackHitTimeU;   //!
   TBranch        *b_TrackHitTimeV;   //!
   TBranch        *b_TrackHitTimeX;   //!
   /*TBranch        *b_nClustersU;   //!
   TBranch        *b_nClustersV;   //!
   TBranch        *b_nClustersX;   //!
   TBranch        *b_ClusterEnergyU;   //!
   TBranch        *b_ClusterEnergyV;   //!
   TBranch        *b_ClusterEnergyX;   //!
   TBranch        *b_ClusterTimeU;   //!
   TBranch        *b_ClusterTimeV;   //!
   TBranch        *b_ClusterTimeX;   //!
   TBranch        *b_ClusterPosMeanU;   //!
   TBranch        *b_ClusterPosMeanV;   //!
   TBranch        *b_ClusterPosMeanX;   //!
   TBranch        *b_ClusterPosStdDevU;   //!
   TBranch        *b_ClusterPosStdDevV;   //!
   TBranch        *b_ClusterPosStdDevX;   //!
   TBranch        *b_nHitsInClusterU;   //!
   TBranch        *b_nHitsInClusterV;   //!
   TBranch        *b_nHitsInClusterX;   //!
   TBranch        *b_ClusterHitPosU;   //!
   TBranch        *b_ClusterHitPosV;   //!
   TBranch        *b_ClusterHitPosX;   //!
   TBranch        *b_ClusterHitEnergyU;   //!
   TBranch        *b_ClusterHitEnergyV;   //!
   TBranch        *b_ClusterHitEnergyX;   //!
   TBranch        *b_ClusterHitTimeU;   //!
   TBranch        *b_ClusterHitTimeV;   //!
   TBranch        *b_ClusterHitTimeX;   //!
   TBranch        *b_ClusterHitSliceU;   //!
   TBranch        *b_ClusterHitSliceV;   //!
   TBranch        *b_ClusterHitSliceX;   //! */
   TBranch        *b_nHits;   //!
   TBranch        *b_RecoHitPos;   //!
   TBranch        *b_RecoHitEnergy;   //!
   TBranch        *b_RecoHitPE;   //!
   TBranch        *b_RecoHitBar;   //!
   TBranch        *b_RecoHitPlane;   //!
   TBranch        *b_RecoHitSlice;   //! 

   Line_Candidates(TTree *tree=0);
   virtual ~Line_Candidates();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Long64_t GetEntriesFast() { return fChain->GetEntriesFast(); };
};

#endif

/*

Line_Candidates::~Line_Candidates()
{
   if (!fChain) return;
}
Line_Candidates::Line_Candidates(TTree *tree) : fChain(0) 
{

}

Int_t Line_Candidates::Cut(Long64_t entry) { return 0; }
Int_t    Line_Candidates::GetEntry(Long64_t entry) { return 0; }
Long64_t Line_Candidates::LoadTree(Long64_t entry) { return 0; }
void     Line_Candidates::Init(TTree *tree) { }
void     Line_Candidates::Loop() {}
Bool_t   Line_Candidates::Notify() { return true; }
void     Line_Candidates::Show(Long64_t entry)  {}*/

#define Line_Candidates_cxx
#ifdef Line_Candidates_cxx
Line_Candidates::Line_Candidates(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("neutrino.0_1698176146.edep_TMS_RecoCandidates_Hough_Cluster1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("neutrino.0_1698176146.edep_TMS_RecoCandidates_Hough_Cluster1.root");
      }
      f->GetObject("Line_Candidates",tree);

   }
   Init(tree);
}

Line_Candidates::~Line_Candidates()
{
   if (!fChain) return;
}

Int_t Line_Candidates::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Line_Candidates::LoadTree(Long64_t entry)
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

void Line_Candidates::Init(TTree *tree)
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
   fChain->SetBranchAddress("nLinesU", &nLinesU, &b_nLinesU);
   fChain->SetBranchAddress("nLinesV", &nLinesV, &b_nLinesV);
   fChain->SetBranchAddress("nLinesX", &nLinesX, &b_nLinesX);
   fChain->SetBranchAddress("nLines3D", &nLines3D, &b_nLines3D);
   fChain->SetBranchAddress("SlopeU", SlopeU, &b_SlopeU);
   fChain->SetBranchAddress("InterceptU", InterceptU, &b_InterceptU);
   fChain->SetBranchAddress("Slope_DownstreamU", Slope_DownstreamU, &b_Slope_DownstreamU);
   fChain->SetBranchAddress("Intercept_DownstreamU", Intercept_DownstreamU, &b_Intercept_DownstreamU);
   fChain->SetBranchAddress("Slope_UpstreamU", Slope_UpstreamU, &b_Slope_UpstreamU);
   fChain->SetBranchAddress("Intercept_UpstreamU", Intercept_UpstreamU, &b_Intercept_UpstreamU);
   fChain->SetBranchAddress("SlopeV", SlopeV, &b_SlopeV);
   fChain->SetBranchAddress("InterceptV", InterceptV, &b_InterceptV);
   fChain->SetBranchAddress("Slope_DownstreamV", Slope_DownstreamV, &b_Slope_DownstreamV);
   fChain->SetBranchAddress("Intercept_DownstreamV", Intercept_DownstreamV, &b_Intercept_DownstreamV);
   fChain->SetBranchAddress("Slope_UpstreamV", Slope_UpstreamV, &b_Slope_UpstreamV);
   fChain->SetBranchAddress("Intercept_UpstreamV", Intercept_UpstreamV, &b_Intercept_UpstreamV);
   fChain->SetBranchAddress("SlopeX", &SlopeX, &b_SlopeX);
   fChain->SetBranchAddress("InterceptX", &InterceptX, &b_InterceptX);
   fChain->SetBranchAddress("Slope_DownstreamX", &Slope_DownstreamX, &b_Slope_DownstreamX);
   fChain->SetBranchAddress("Intercept_DownstreamX", &Intercept_DownstreamX, &b_Intercept_DownstreamX);
   fChain->SetBranchAddress("Slope_UpstreamX", &Slope_UpstreamX, &b_Slope_UpstreamX);
   fChain->SetBranchAddress("Intercept_UpstreamX", &Intercept_UpstreamX, &b_Intercept_UpstreamX);
   /*fChain->SetBranchAddress("DirectionZU", DirectionZU, &b_DirectionZU);
   fChain->SetBranchAddress("DirectionZV", DirectionZV, &b_DirectionZV);
   fChain->SetBranchAddress("DirectionZX", &DirectionZX, &b_DirectionZX);
   fChain->SetBranchAddress("DirectionXU", DirectionXU, &b_DirectionXU);
   fChain->SetBranchAddress("DirectionXV", DirectionXV, &b_DirectionXV);
   fChain->SetBranchAddress("DirectionYX", &DirectionYX, &b_DirectionYX);
   fChain->SetBranchAddress("DirectionZU_Downstream", DirectionZU_Downstream, &b_DirectionZU_Downstream);
   fChain->SetBranchAddress("DirectionXU_Downstream", DirectionXU_Downstream, &b_DirectionXU_Downstream);
   fChain->SetBranchAddress("DirectionZU_Upstream", DirectionZU_Upstream, &b_DirectionZU_Upstream);
   fChain->SetBranchAddress("DirectionXU_Upstream", DirectionXU_Upstream, &b_DirectionXU_Upstream);
   fChain->SetBranchAddress("DirectionZV_Downstream", DirectionZV_Downstream, &b_DirectionZV_Downstream);
   fChain->SetBranchAddress("DirectionXV_Downstream", DirectionXV_Downstream, &b_DirectionXV_Downstream);
   fChain->SetBranchAddress("DirectionZV_Upstream", DirectionZV_Upstream, &b_DirectionZV_Upstream);
   fChain->SetBranchAddress("DirectionXV_Upstream", DirectionXV_Upstream, &b_DirectionXV_Upstream);
   fChain->SetBranchAddress("DirectionZX_Downstream", &DirectionZX_Downstream, &b_DirectionZX_Downstream);
   fChain->SetBranchAddress("DirectionYX_Downstream", &DirectionYX_Downstream, &b_DirectionYX_Downstream);
   fChain->SetBranchAddress("DirectionZX_Upstream", &DirectionZX_Upstream, &b_DirectionZX_Upstream);
   fChain->SetBranchAddress("DirectionYX_Upstream", &DirectionYX_Upstream, &b_DirectionYX_Upstream);
   fChain->SetBranchAddress("FirstHoughHitU", FirstHoughHitU, &b_FirstHoughHitU);
   fChain->SetBranchAddress("FirstHoughHitV", FirstHoughHitV, &b_FirstHoughHitV);
   fChain->SetBranchAddress("FirstHoughHitX", &FirstHoughHitX, &b_FirstHoughHitX);
   fChain->SetBranchAddress("LastHoughHitU", LastHoughHitU, &b_LastHoughHitU);
   fChain->SetBranchAddress("LastHoughHitV", LastHoughHitV, &b_LastHoughHitV);
   fChain->SetBranchAddress("LastHoughHitX", &LastHoughHitX, &b_LastHoughHitX);
   fChain->SetBranchAddress("FirstHoughHitTimeU", FirstHoughHitTimeU, &b_FirstHoughHitTimeU);
   fChain->SetBranchAddress("FirstHoughHitTimeV", FirstHoughHitTimeV, &b_FirstHoughHitTimeV);
   fChain->SetBranchAddress("FirstHoughHitTimeX", &FirstHoughHitTimeX, &b_FirstHoughHitTimeX);
   fChain->SetBranchAddress("LastHoughHitTimeU", LastHoughHitTimeU, &b_LastHoughHitTimeU);
   fChain->SetBranchAddress("LastHoughHitTimeV", LastHoughHitTimeV, &b_LastHoughHitTimeV);
   fChain->SetBranchAddress("LastHoughHitTImeX", &LastHoughHitTImeX, &b_LastHoughHitTImeX);
   fChain->SetBranchAddress("HoughEarliestHitTimeU", HoughEarliestHitTimeU, &b_HoughEarliestHitTimeU);
   fChain->SetBranchAddress("HoughEarliestHitTimeV", HoughEarliestHitTimeV, &b_HoughEarliestHitTimeV);
   fChain->SetBranchAddress("HoughEarliestHitTimeX", &HoughEarliestHitTimeX, &b_HoughEarliestHitTimeX);
   fChain->SetBranchAddress("HoughLatestHitTimeU", HoughLatestHitTimeU, &b_HoughLatestHitTimeU);
   fChain->SetBranchAddress("HoughLatestHitTimeV", HoughLatestHitTimeV, &b_HoughLatestHitTimeV);
   fChain->SetBranchAddress("HoughLatestHitTimeX", &HoughLatestHitTimeX, &b_HoughLatestHitTimeX);
   fChain->SetBranchAddress("FirstHoughPlaneU", FirstHoughPlaneU, &b_FirstHoughPlaneU);
   fChain->SetBranchAddress("FirstHoughPlaneV", FirstHoughPlaneV, &b_FirstHoughPlaneV);
   fChain->SetBranchAddress("FirstHoughPlaneX", &FirstHoughPlaneX, &b_FirstHoughPlaneX);
   fChain->SetBranchAddress("LastHoughPlaneU", LastHoughPlaneU, &b_LastHoughPlaneU);
   fChain->SetBranchAddress("LastHoughPlaneV", LastHoughPlaneV, &b_LastHoughPlaneV);
   fChain->SetBranchAddress("LastHoughPlaneX", &LastHoughPlaneX, &b_LastHoughPlaneX); 
   fChain->SetBranchAddress("TMSStart", &TMSStart, &b_TMSStart);
   fChain->SetBranchAddress("TMSStartTime", &TMSStartTime, &b_TMSStartTime);
   fChain->SetBranchAddress("OccupancyU", OccupancyU, &b_OccupancyU);
   fChain->SetBranchAddress("OccupancyV", OccupancyV, &b_OccupancyV);
   fChain->SetBranchAddress("OccupancyX", &OccupancyX, &b_OccupancyX); 
   fChain->SetBranchAddress("TrackLengthU", TrackLengthU, &b_TrackLengthU);
   fChain->SetBranchAddress("TrackLengthV", TrackLengthV, &b_TrackLengthV);
   fChain->SetBranchAddress("TrackLengthX", &TrackLengthX, &b_TrackLengthX);
   fChain->SetBranchAddress("TotalTrackEnergyU", TotalTrackEnergyU, &b_TotalTrackEnergyU);
   fChain->SetBranchAddress("TotalTrackEnergyV", TotalTrackEnergyV, &b_TotalTrackEnergyV);
   fChain->SetBranchAddress("TotalTrackEnergyX", &TotalTrackEnergyX, &b_TotalTrackEnergyX);
   fChain->SetBranchAddress("TrackStoppingU", TrackStoppingU, &b_TrackStoppingU);
   fChain->SetBranchAddress("TrackStoppingV", TrackStoppingV, &b_TrackStoppingV);
   fChain->SetBranchAddress("TrackStoppingX", &TrackStoppingX, &b_TrackStoppingX);
   fChain->SetBranchAddress("nHitsInTrackU", nHitsInTrackU, &b_nHitsInTrackU);
   fChain->SetBranchAddress("nHitsInTrackV", nHitsInTrackV, &b_nHitsInTrackV);
   fChain->SetBranchAddress("nHitsInTrackX", &nHitsInTrackX, &b_nHitsInTrackX);
   fChain->SetBranchAddress("TrackHitEnergyU", TrackHitEnergyU, &b_TrackHitEnergyU);
   fChain->SetBranchAddress("TrackHitEnergyV", TrackHitEnergyV, &b_TrackHitEnergyV);
   fChain->SetBranchAddress("TrackHitEnergyX", TrackHitEnergyX, &b_TrackHitEnergyX);
   fChain->SetBranchAddress("TrackHitPosU", TrackHitPosU, &b_TrackHitPosU);
   fChain->SetBranchAddress("TrackHitPosV", TrackHitPosV, &b_TrackHitPosV);
   fChain->SetBranchAddress("TrackHitPosX", TrackHitPosX, &b_TrackHitPosX);
   fChain->SetBranchAddress("TrackHitTimeU", TrackHitTimeU, &b_TrackHitTimeU);
   fChain->SetBranchAddress("TrackHitTimeV", TrackHitTimeV, &b_TrackHitTimeV);
   fChain->SetBranchAddress("TrackHitTimeX", TrackHitTimeX, &b_TrackHitTimeX);
   fChain->SetBranchAddress("nClustersU", &nClustersU, &b_nClustersU);
   fChain->SetBranchAddress("nClustersV", &nClustersV, &b_nClustersV);
   fChain->SetBranchAddress("nClusterX", &nClusterX, &b_nClustersX);
   fChain->SetBranchAddress("ClusterEnergyU", ClusterEnergyU, &b_ClusterEnergyU);
   fChain->SetBranchAddress("ClusterEnergyV", ClusterEnergyV, &b_ClusterEnergyV);
   fChain->SetBranchAddress("ClusterEnergyX", &ClusterEnergyX, &b_ClusterEnergyX);
   fChain->SetBranchAddress("ClusterTimeU", ClusterTimeU, &b_ClusterTimeU);
   fChain->SetBranchAddress("ClusterTimeV", ClusterTimeV, &b_ClusterTimeV);
   fChain->SetBranchAddress("ClusterTimeX", &ClusterTimeX, &b_ClusterTimeX);
   fChain->SetBranchAddress("ClusterPosMeanU", ClusterPosMeanU, &b_ClusterPosMeanU);
   fChain->SetBranchAddress("ClusterPosMeanV", ClusterPosMeanV, &b_ClusterPosMeanV);
   fChain->SetBranchAddress("ClusterPosMeanX", ClusterPosMeanX, &b_ClusterPosMeanX);
   fChain->SetBranchAddress("ClusterPosStdDevU", ClusterPosStdDevU, &b_ClusterPosStdDevU);
   fChain->SetBranchAddress("ClusterPosStdDevV", ClusterPosStdDevV, &b_ClusterPosStdDevV);
   fChain->SetBranchAddress("ClusterPosStdDevX", ClusterPosStdDevX, &b_ClusterPosStdDevX);
   fChain->SetBranchAddress("nHitsInClusterU", nHitsInClusterU, &b_nHitsInClusterU);
   fChain->SetBranchAddress("nHitsInClusterV", nHitsInClusterV, &b_nHitsInClusterV);
   fChain->SetBranchAddress("nHitsInClusterX", &nHitsInClusterX, &b_nHitsInClusterX);
   fChain->SetBranchAddress("ClusterHitPosU", ClusterHitPosU, &b_ClusterHitPosU);
   fChain->SetBranchAddress("ClusterHitPosV", ClusterHitPosV, &b_ClusterHitPosV);
   fChain->SetBranchAddress("ClusterHitPosX", ClusterHitPosX, &b_ClusterHitPosX);
   fChain->SetBranchAddress("ClusterHitEnergyU", ClusterHitEnergyU, &b_ClusterHitEnergyU);
   fChain->SetBranchAddress("ClusterHitEnergyV", ClusterHitEnergyV, &b_ClusterHitEnergyV);
   fChain->SetBranchAddress("ClusterHitEnergyX", ClusterHitEnergyX, &b_ClusterHitEnergyX);
   fChain->SetBranchAddress("ClusterHitTimeU", ClusterHitTimeU, &b_ClusterHitTimeU);
   fChain->SetBranchAddress("ClusterHitTimeV", ClusterHitTimeV, &b_ClusterHitTimeV);
   fChain->SetBranchAddress("ClusterHitTimeX", ClusterHitTimeX, &b_ClusterHitTimeX);
   fChain->SetBranchAddress("ClusterHitSliceU", ClusterHitSliceU, &b_ClusterHitSliceU);
   fChain->SetBranchAddress("ClusterHitSliceV", ClusterHitSliceV, &b_ClusterHitSliceV);
   fChain->SetBranchAddress("ClusterHitSliceX", ClusterHitSliceX, &b_ClusterHitSliceX); */
   fChain->SetBranchAddress("nHits", &nHits, &b_nHits);
   fChain->SetBranchAddress("RecoHitPos", RecoHitPos, &b_RecoHitPos);
   fChain->SetBranchAddress("RecoHitEnergy", RecoHitEnergy, &b_RecoHitEnergy);
   fChain->SetBranchAddress("RecoHitPE", RecoHitPE, &b_RecoHitPE);
   fChain->SetBranchAddress("RecoHitBar", RecoHitBar, &b_RecoHitBar);
   fChain->SetBranchAddress("RecoHitPlane", RecoHitPlane, &b_RecoHitPlane);
   fChain->SetBranchAddress("RecoHitSlice", RecoHitSlice, &b_RecoHitSlice); 
   Notify();
}

Bool_t Line_Candidates::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Line_Candidates::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Line_Candidates::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   entry += 1;
   return 1;
}


void Line_Candidates::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Line_Candidates.C
//      root> Line_Candidates t
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
#endif // #ifdef Line_Candidates_cxx
