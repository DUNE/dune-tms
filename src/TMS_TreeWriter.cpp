#include "TMS_TreeWriter.h"

TMS_TreeWriter::TMS_TreeWriter() {
  // Make the output file
  std::string filename = TMS_Manager::GetInstance().GetFileName();

  while (filename.find("/") != std::string::npos) {
    filename = filename.substr(filename.find("/")+1, filename.size());
  }
  TString Outputname = filename.c_str();
  Outputname.ReplaceAll(".root", "_LineCandidates.root");
  // Make an output file
  Output = new TFile(Outputname, "recreate");
  if (!Output->IsOpen()) {
    std::cerr << "Could not write to file " << Outputname << std::endl;
    std::cerr << "Are you sure you have write access to the directory?" << std::endl;
    throw;
  }
  Output->cd();
  Branch_Lines = new TTree("Line_Candidates", "Line_Candidates");
  Branch_Lines->SetAutoSave(1E3); // Every 1000 events (negative for MB)

  MakeBranches();
}

void TMS_TreeWriter::MakeBranches() {
  Output->cd();
  Branch_Lines->Branch("nLines", &nLines, "nLines/I");
  Branch_Lines->Branch("EventNo", &EventNo, "EventNo/I");
  Branch_Lines->Branch("Slope", Slope, "Slope[nLines]/D");
  Branch_Lines->Branch("Intercept", Intercept, "Intercept[nLines]/D");
  Branch_Lines->Branch("DirectionZ", DirectionZ, "DirectionZ[nLines]/D");
  Branch_Lines->Branch("DirectionX", DirectionX, "DirectionX[nLines]/D");
  Branch_Lines->Branch("FirstHoughHit", FirstHit, "FirstHoughHit[2]/D");
  Branch_Lines->Branch("FirstHoughPlane", &FirstPlane, "FirstHoughPlane/I");
  Branch_Lines->Branch("TMSStart", &TMSStart, "TMSStart/O");
}

void TMS_TreeWriter::Fill(int i) {
  EventNo = i;
  std::vector<std::pair<bool, TF1*>> HoughLines = TMS_TrackFinder::GetFinder().GetHoughLines();
  nLines = HoughLines.size();
  // Skip the event if there aren't any Hough Lines
  if (nLines == 0) return;
  int it = 0;
  for (auto &Lines: HoughLines) {
    Intercept[it] = Lines.second->GetParameter(0);
    Slope[it] = Lines.second->GetParameter(1);
    // Calculate the z and x vectors by evaling the TF1 in thin and thick target
    double xlow = TMS_Const::TMS_Thin_Start;
    double xhi = TMS_Const::TMS_Thick_Start;
    double ylow = Lines.second->Eval(xlow);
    double yhi = Lines.second->Eval(xhi);

    double xlen = xhi-xlow;
    double ylen = yhi-ylow;
    double len = sqrt(xlen*xlen+ylen*ylen);
    xlen = xlen/len;
    ylen = ylen/len;
    DirectionZ[it] = xlen;
    DirectionX[it] = ylen;

    it++;
  }

  // Find where the first hough hit is
  TMSStart = false;

  std::vector<std::vector<TMS_Hit> > HoughCands = TMS_TrackFinder::GetFinder().GetHoughCandidates();
  TMS_Hit *FirstTrack = NULL;
  for (auto &Candidates: HoughCands) {
    for (auto &hit: Candidates) {
      if (FirstTrack == NULL) {
        FirstTrack = &hit;
      } else if (hit.GetZ() < FirstTrack->GetZ()) {
        FirstTrack = &hit;
      }
    }
  }
  // Was the first hit within the first 3 layers?
  if (FirstTrack->GetPlaneNumber() < 4) TMSStart = false;
  else TMSStart = true;
  // Then save the hit info
  FirstHit[0] = FirstTrack->GetZ();
  FirstHit[1] = FirstTrack->GetNotZ();
  FirstPlane = FirstTrack->GetPlaneNumber();

  Branch_Lines->Fill();
}


