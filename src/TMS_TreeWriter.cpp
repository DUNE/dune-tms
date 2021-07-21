#include "TMS_TreeWriter.h"

TMS_TreeWriter::TMS_TreeWriter() {
  // Make the output file
  std::string filename = TMS_Manager::GetInstance().GetFileName();

  while (filename.find("/") != std::string::npos) {
    filename = filename.substr(filename.find("/")+1, filename.size());
  }
  TString Outputname = filename.c_str();
  Outputname.ReplaceAll(".root", "_LineCandidates");
  Outputname += "_"+TMS_Manager::GetInstance().Get_Reco_TrackMethod();
  Outputname += Form("_Cluster%i", TMS_Manager::GetInstance().Get_Reco_Clustering());
  Outputname += ".root";
  // Make an output file
  Output = new TFile(Outputname, "recreate");
  if (!Output->IsOpen()) {
    std::cerr << "Could not write to file " << Outputname << std::endl;
    std::cerr << "Are you sure you have write access to the directory?" << std::endl;
    throw;
  }
  Output->cd();

  Branch_Lines = new TTree("Line_Candidates", "Line_Candidates");
  Branch_Lines->SetDirectory(Output);
  Branch_Lines->SetAutoSave(__TMS_AUTOSAVE__); // Every 1000 events (negative for MB)

  Truth_Info = new TTree("Truth_Info", "Truth_Info");
  Truth_Info->SetDirectory(Output);
  Truth_Info->SetAutoSave(__TMS_AUTOSAVE__);

  MakeBranches();
}

void TMS_TreeWriter::MakeBranches() {
  Branch_Lines->Branch("nLines", &nLines, "nLines/I");
  Branch_Lines->Branch("EventNo", &EventNo, "EventNo/I");
  Branch_Lines->Branch("Slope", Slope, "Slope[nLines]/D");
  Branch_Lines->Branch("Intercept", Intercept, "Intercept[nLines]/D");
  Branch_Lines->Branch("DirectionZ", DirectionZ, "DirectionZ[nLines]/D");
  Branch_Lines->Branch("DirectionX", DirectionX, "DirectionX[nLines]/D");
  Branch_Lines->Branch("FirstHoughHit", FirstHit, "FirstHoughHit[nLines][2]/D");
  Branch_Lines->Branch("LastHoughHit", LastHit, "LastHoughHit[nLines][2]/D");
  Branch_Lines->Branch("FirstHoughPlane", FirstPlane, "FirstHoughPlane[nLines]/I");
  Branch_Lines->Branch("LastHoughPlane", LastPlane, "LastHoughPlane[nLines]/I");
  Branch_Lines->Branch("TMSStart", &TMSStart, "TMSStart/O");
  Branch_Lines->Branch("Occupancy", Occupancy, "Occupancy[nLines]/D");

  Truth_Info->Branch("MuonP4", MuonP4, "MuonP4[4]/D");
  Truth_Info->Branch("Muon_Vertex", Muon_Vertex, "Muon_Vertex[4]/D");
  Truth_Info->Branch("nParticles", &nParticles, "nParticles/I");
  Truth_Info->Branch("Interaction", &Reaction);
  Truth_Info->Branch("EventNo", &EventNo, "EventNo/I");
  Truth_Info->Branch("NeutrinoPDG", &NeutrinoPDG, "NeutrinoPDG/I");
  Truth_Info->Branch("NeutrinoP4", NeutrinoP4, "NeutrinoP4[4]/D");
}

void TMS_TreeWriter::Fill(TMS_Event &event) {

  // Clear old info
  Clear();

  // Fill the truth info
  EventNo = event.GetEventNumber();
  Reaction = event.GetReaction();

  NeutrinoPDG = event.GetNeutrinoPDG();
  NeutrinoP4[0] = event.GetNeutrinoP4().X();
  NeutrinoP4[1] = event.GetNeutrinoP4().Y();
  NeutrinoP4[2] = event.GetNeutrinoP4().Z();
  NeutrinoP4[3] = event.GetNeutrinoP4().T();

  // Get the truth info
  std::vector<TMS_TrueParticle> TrueParticles = event.GetTrueParticles();
  nParticles = TrueParticles.size();
  for (auto it = TrueParticles.begin(); it != TrueParticles.end(); ++it) {

    // Only save muon info for now
    if (abs((*it).GetPDG()) != 13) continue;

    MuonP4[0] = (*it).GetMomentum().Px();
    MuonP4[1] = (*it).GetMomentum().Py();
    MuonP4[2] = (*it).GetMomentum().Pz();
    MuonP4[3] = (*it).GetMomentum().E();
    Muon_Vertex[0] = (*it).GetPosition().X();
    Muon_Vertex[1] = (*it).GetPosition().Y();
    Muon_Vertex[2] = (*it).GetPosition().Z();
    Muon_Vertex[3] = (*it).GetPosition().T();
  }

  Truth_Info->Fill();

  // Fill the reco info
  std::vector<std::pair<bool, TF1*>> HoughLines = TMS_TrackFinder::GetFinder().GetHoughLines();
  nLines = HoughLines.size();
  // Also get the size of the hits to get a measure of relative goodness
  std::vector<std::vector<TMS_Hit> > HoughCandidates = TMS_TrackFinder::GetFinder().GetHoughCandidates();
  // Get the cleaned hits for a reference of the "total"
  int TotalHits = TMS_TrackFinder::GetFinder().GetCleanedHits().size();
 
  // Skip the event if there aren't any Hough Lines
  if (nLines == 0) return;
  if (nLines > __TMS_MAX_LINES__) {
    std::cerr << "Exceeded max number of HoughLines to write to file" << std::endl;
    std::cerr << "Not writing event" << std::endl;
    return;
  }

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
    Occupancy[it] = (double)HoughCandidates[it].size()/TotalHits;

    it++;
  }

  // Find where the first hough hit is
  TMSStart = false;

  std::vector<std::vector<TMS_Hit> > HoughCands = TMS_TrackFinder::GetFinder().GetHoughCandidates();
  TMS_Hit *FirstTrack = NULL;

  it = 0;
  for (auto &Candidates: HoughCands) {

    for (auto &hit: Candidates) {
      if (FirstTrack == NULL) {
        FirstTrack = &hit;
      } else if (hit.GetZ() < FirstTrack->GetZ()) {
        FirstTrack = &hit;
      }
    }

    // Then save the hit info
    FirstPlane[it] = Candidates.front().GetPlaneNumber();
    FirstHit[it][0] = Candidates.front().GetZ();
    FirstHit[it][1] = Candidates.front().GetNotZ();

    LastPlane[it] = Candidates.back().GetPlaneNumber();
    LastHit[it][0] = Candidates.back().GetZ();
    LastHit[it][1] = Candidates.back().GetNotZ();

    it++;
  }

  // Was the first hit within the first 3 layers?
  if (FirstTrack->GetPlaneNumber() < 4) TMSStart = false;
  else TMSStart = true;

  Branch_Lines->Fill();
}

void TMS_TreeWriter::Clear() {

  // The variables
  EventNo = nLines = nParticles = NeutrinoPDG = -999;
  TMSStart = false;
  Reaction = "";


  for (int i = 0; i < 4; ++i) {
    MuonP4[i]=-999;
    Muon_Vertex[i]=-999;
    NeutrinoP4[i]=-999;
  }

  for (int i = 0; i < __TMS_MAX_LINES__; ++i) {
    Slope[i]=-999;
    Intercept[i]=-999;
    DirectionZ[i]=-999;
    DirectionX[i]=-999;
    Occupancy[i]=-999;
    FirstPlane[i]=-999;
    LastPlane[i]=-999;
    for (int j = 0; j < 2; ++j) {
      FirstHit[i][j] = -999;
      LastHit[i][j] = -999;
    }
  }

}

