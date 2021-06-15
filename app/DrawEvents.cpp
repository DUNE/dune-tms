// python code was very slow...
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TGeoManager.h"
#include "TStopwatch.h"
#include "TLorentzVector.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGeoManager.h"
#include "TStyle.h"
#include "TPaveText.h"

// EDepSim includes
#include "EDepSim/TG4Event.h"

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cout << "Need one argument, the input filename" << std::endl;
    return -1;
  }

  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(255);
  gStyle->SetPalette(89);

  std::string filename = std::string(argv[1]);

  std::cout << "Got " << filename << std::endl;

  // The input file
  TFile *input = new TFile(filename.c_str(), "open");
  // The EDepSim events
  TTree *events = (TTree*)input->Get("EDepSimEvents");
  // The generator pass-through information
  TTree *gRoo = (TTree*)input->Get("DetSimPassThru/gRooTracker");
  // Get the detector geometry
  //TGeoManager *geom = (TGeoManager*)input->Get("EDepSimGeometry");

  // Get the event
  TG4Event *event = NULL;
  events->SetBranchAddress("Event", &event);

  int N_entries = events->GetEntries();
  //int N_entries = 100;
  const double zmin = 1000;
  const double zmax = 20000;
  const double xmin = -4000;
  const double xmax = 4000;

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->SetLeftMargin(canv->GetLeftMargin()*1.4);
  canv->SetRightMargin(canv->GetRightMargin()*1.4);
  // Make the output name
  std::string outputname = filename;
  while (outputname.find("/") != std::string::npos) {
    outputname = outputname.substr(outputname.find("/")+1, outputname.size());
  }
  TString Outputname = outputname.c_str();
  Outputname.ReplaceAll(".root", "_Deposits.pdf");
  canv->Print(Outputname+"[");

  TPaveText *text = new TPaveText(0,0,1,1,"NDC");
  text->AddText(Form("Exec: %s", argv[0]));
  text->AddText(Form("Input: %s", filename.c_str()));
  text->AddText(Form("Output: %s", Outputname.Data()));
  text->SetBorderSize(0);
  text->SetFillColor(0);
  text->Draw();
  canv->Print(Outputname);

  gErrorIgnoreLevel = kWarning;
  int i = 0;
  TH2D *plot = new TH2D("xz", "xz-view;z (mm); x (mm); Energy deposit (MeV)", 1000, zmin, zmax, 1000, xmin, xmax);
  for (; i < N_entries; ++i) {
    events->GetEntry(i);
    gRoo->GetEntry(i);

    plot->Reset();
    plot->SetTitle(Form("Event Number %i", i));
    if (i % (N_entries/10) == 0) {
      std::cout << "Processed " << i << "/" << N_entries << " (" << double(i)*100./N_entries << "%)" << std::endl;
    }

    // Just extract the hits from the edep-sim event

    // Loop over the primary vertices
    /*
    for (TG4PrimaryVertexContainer::iterator it = event->Primaries.begin(); it != event->Primaries.end(); ++it) {
      TG4PrimaryVertex vtx = *it;
    }
    */

    // Loop over the
    for (TG4HitSegmentDetectors::iterator jt = event->SegmentDetectors.begin(); jt != event->SegmentDetectors.end(); ++jt) {
      TG4HitSegmentContainer tms_hits = (*jt).second;
      for (TG4HitSegmentContainer::iterator kt = tms_hits.begin(); kt != tms_hits.end(); ++kt) {
        TG4HitSegment edep_hit = *kt;
        TLorentzVector Position = (edep_hit.GetStop()+edep_hit.GetStart());
        Position *= 0.5;
        plot->Fill(Position.Z(), Position.X(), edep_hit.GetEnergyDeposit());
      }
    }

    if (plot->Integral() == 0) continue;

    /*
    // Loop over the trajectories
    for (TG4TrajectoryContainer::iterator jt = event->Trajectories.begin(); jt != event->Trajectories.end(); ++jt) {
      TG4Trajectory traj = *jt;
      //int PDGcode = traj.GetPDGCode();
      //int TrackId = traj.GetTrackId();

      // Loop over each point (hit) in the trajectory
      for (std::vector<TG4TrajectoryPoint>::iterator kt = traj.Points.begin(); kt != traj.Points.end(); kt++) {
        TG4TrajectoryPoint pt = *kt;
        TLorentzVector Position = pt.GetPosition();
        TVector3 Momentum = pt.GetMomentum();
      }
    }
    */
    plot->SetMinimum(0);
    plot->SetMaximum(5);
    plot->Draw("colz");
    canv->Print(Outputname);

  }

  canv->Print(Outputname+"]");
}
