#include "TMS_ReadoutTreeWriter.h"
#include "TMS_Hit.h"

// Saves every n events
#define __TMS_READOUT_AUTOSAVE__ 1000

TMS_ReadoutTreeWriter::TMS_ReadoutTreeWriter() {
  
  // Make the output file
  std::string filename = TMS_Manager::GetInstance().GetFileName();
  // Removes directory path
  while (filename.find("/") != std::string::npos) {
    filename = filename.substr(filename.find("/")+1, filename.size());
  }
  
  TString Outputname = filename.c_str();
  Outputname.ReplaceAll(".root", "_TMS_Readout.root");
  // Make an output file
  Output = new TFile(Outputname, "recreate");
  
  TMS_Readout = new TTree("TMS", "TMS");
  TMS_Readout->SetDirectory(Output);
  TMS_Readout->SetAutoSave(__TMS_READOUT_AUTOSAVE__);
  
  // Now make the branches
  MakeBranches();
}

void TMS_ReadoutTreeWriter::MakeBranches() {
  // Makes the branches
  
  // Truth branches
  if (hasTruth) {
    TMS_Readout->Branch("NTrueHits", &NTrueHits);
    TMS_Readout->Branch("TrueHitX", &TrueHitX, "TrueHitX[NTrueHits]/F");
    TMS_Readout->Branch("TrueHitY", &TrueHitY, "TrueHitY[NTrueHits]/F");
    TMS_Readout->Branch("TrueHitZ", &TrueHitZ, "TrueHitZ[NTrueHits]/F");
    TMS_Readout->Branch("TrueHitT", &TrueHitT, "TrueHitT[NTrueHits]/F");
    TMS_Readout->Branch("TrueHitE", &TrueHitE, "TrueHitE[NTrueHits]/F");
    TMS_Readout->Branch("TrueHitPE", &TrueHitPE, "TrueHitPE[NTrueHits]/F");
    TMS_Readout->Branch("TrueHitPEAfterFibers", &TrueHitPEAfterFibers, "TrueHitPEAfterFibers[NTrueHits]/F");
    TMS_Readout->Branch("TrueHitPEAfterFibersLongPath", &TrueHitPEAfterFibersLongPath, "TrueHitPEAfterFibersLongPath[NTrueHits]/F");
    TMS_Readout->Branch("TrueHitPEAfterFibersShortPath", &TrueHitPEAfterFibersShortPath, "TrueHitPEAfterFibersShortPath[NTrueHits]/F");
  }
  
  // Reco branches
  TMS_Readout->Branch("NRecoHits", &NRecoHits);
  TMS_Readout->Branch("RecoHitX", &RecoHitX, "RecoHitX[NRecoHits]/F");
  TMS_Readout->Branch("RecoHitY", &RecoHitY, "RecoHitY[NRecoHits]/F");
  TMS_Readout->Branch("RecoHitZ", &RecoHitZ, "RecoHitZ[NRecoHits]/F");
  TMS_Readout->Branch("RecoHitNotZ", &RecoHitNotZ, "RecoHitNotZ[NRecoHits]/F");
  TMS_Readout->Branch("RecoHitT", &RecoHitT, "RecoHitT[NRecoHits]/F");
  TMS_Readout->Branch("RecoHitE", &RecoHitE, "RecoHitE[NRecoHits]/F");
  TMS_Readout->Branch("RecoHitPE", &RecoHitPE, "RecoHitPE[NRecoHits]/F");
  TMS_Readout->Branch("RecoHitIsPedSupped", &RecoHitIsPedSupped, "RecoHitIsPedSupped[NRecoHits]/O");
  TMS_Readout->Branch("RecoHitBar", &RecoHitBar, "RecoHitBar[NRecoHits]/I");
  TMS_Readout->Branch("RecoHitPlane", &RecoHitPlane, "RecoHitPlane[NRecoHits]/I");
  #ifdef RECORD_HIT_DEADTIME
  TMS_Readout->Branch("RecoHitDeadtimeStart", &RecoHitDeadtimeStart, "RecoHitDeadtimeStart[NRecoHits]/F");
  TMS_Readout->Branch("RecoHitDeadtimeStop", &RecoHitDeadtimeStop, "RecoHitDeadtimeStop[NRecoHits]/F");
  #endif

}

void TMS_ReadoutTreeWriter::Fill(TMS_Event &event) {
  // Clear branches
  NRecoHits = 0;
  NTrueHits = 0;

  // Now fill branches
  // Get all hits including ped supped ones
  bool include_ped_sup = true;
  const int slice = -1; // Want all slices
  int index = 0;
  for (auto hit : event.GetHits(slice, include_ped_sup)) {
    if (index >= __MAX_READOUT_TREE_ARRAY_LENGTH__) {
      std::cout<<"TMS_ReadoutTreeWriter WARNING: Too many hits in event. Increase __MAX_READOUT_TREE_ARRAY_LENGTH__. Saving partial event"<<std::endl;
      break;
    }
    if (index < __MAX_READOUT_TREE_ARRAY_LENGTH__) {
      // In theory a reco hit should have many true hits based on how the merging worked
      // Plus true hits should have noise hits which don't have any parent info
      auto true_hit = hit.GetTrueHit();
      
      if (hasTruth) {
        // True info
        NTrueHits += 1;
        TrueHitX[index] = true_hit.GetX();
        TrueHitY[index] = true_hit.GetY();
        TrueHitZ[index] = true_hit.GetZ();
        TrueHitT[index] = true_hit.GetT();
        TrueHitE[index] = true_hit.GetE();
        TrueHitPE[index] = true_hit.GetPE();
        TrueHitPEAfterFibers[index] = true_hit.GetPEAfterFibers();
        TrueHitPEAfterFibersLongPath[index] = true_hit.GetPEAfterFibersLongPath();
        TrueHitPEAfterFibersShortPath[index] = true_hit.GetPEAfterFibersShortPath();
      }
      
      // Reco info
      NRecoHits += 1;
      RecoHitX[index] = hit.GetX();
      RecoHitY[index] = hit.GetY();
      RecoHitZ[index] = hit.GetZ();
      RecoHitNotZ[index] = hit.GetNotZ();
      RecoHitT[index] = hit.GetT();
      RecoHitE[index] = hit.GetE();
      RecoHitPE[index] = hit.GetPE();
      RecoHitIsPedSupped[index] = hit.GetPedSup();
      RecoHitBar[index] = hit.GetBarNumber();
      RecoHitPlane[index] = hit.GetPlaneNumber();
      #ifdef RECORD_HIT_DEADTIME
      RecoHitDeadtimeStart[index] = hit.GetDeadtimeStart();
      RecoHitDeadtimeStop[index] = hit.GetDeadtimeStop();
      #endif
    }
    
    index += 1;
  }
  
  // Finally fill
  TMS_Readout->Fill();
}
