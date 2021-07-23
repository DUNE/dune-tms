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
#include "TBox.h"

// EDepSim includes
#include "EDepSim/TG4Event.h"

#include "TMS_Event.h"
#include "TMS_Constants.h"

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
  // Get the true neutrino vector from the gRooTracker object
  int StdHepPdg[__EDEP_SIM_MAX_PART__];
  double StdHepP4[__EDEP_SIM_MAX_PART__][4];
  gRoo->SetBranchStatus("*", false);
  gRoo->SetBranchStatus("StdHepPdg", true);
  gRoo->SetBranchStatus("StdHepP4", true);
  gRoo->SetBranchAddress("StdHepPdg", StdHepPdg);
  gRoo->SetBranchAddress("StdHepP4", StdHepP4);


  // Get the event
  TG4Event *event = NULL;
  events->SetBranchAddress("Event", &event);

  int N_entries = events->GetEntries();
  N_entries = 1000;
  const double zmin = 3000;
  const double zmax = 19000;
  const double xmin = -4000;
  const double xmax = 4000;

  // Make some boxes for the fiducial volumes
  // Full view from inspecting all hits
  TBox *xz_box_Full = new TBox(TMS_Const::TMS_Thin_Start,
      -3485,
      TMS_Const::TMS_Thick_End,
      3485);
  xz_box_Full->SetLineColor(kMagenta+2);
  xz_box_Full->SetFillStyle(0);

  TBox *xz_box_LAr_Full = new TBox(TMS_Const::LAr_Start_Exact[2],
      TMS_Const::LAr_Start_Exact[0],
      TMS_Const::LAr_End_Exact[2],
      TMS_Const::LAr_End_Exact[0]);
  xz_box_LAr_Full->SetLineColor(kMagenta+2);
  xz_box_LAr_Full->SetFillStyle(0);

  // FV just taking 50 cm in from the full
  TBox *xz_box_FV = new TBox(xz_box_Full->GetX1(),
      xz_box_Full->GetY1()+500,
      xz_box_Full->GetX2()-500,
      xz_box_Full->GetY2()-500);
  xz_box_FV->SetLineColor(kRed);
  xz_box_FV->SetLineStyle(kDashed);
  xz_box_FV->SetFillStyle(0);

  TBox *xz_box_LAr_FV = new TBox(xz_box_LAr_Full->GetX1()+500,
      xz_box_LAr_Full->GetY1()+500,
      xz_box_LAr_Full->GetX2()-500,
      xz_box_LAr_Full->GetY2()-500);
  xz_box_LAr_FV->SetLineColor(kRed);
  xz_box_LAr_FV->SetLineStyle(kDashed);
  xz_box_LAr_FV->SetFillStyle(0);


  // Include the dead region boxes
  TBox *xz_dead_top = new TBox(TMS_Const::TMS_Thin_Start,
      TMS_Const::TMS_Dead_Top[0],
      TMS_Const::TMS_Thick_End,
      TMS_Const::TMS_Dead_Top[1]);
  TBox *xz_dead_center = new TBox(TMS_Const::TMS_Thin_Start,
      TMS_Const::TMS_Dead_Center[0],
      TMS_Const::TMS_Thick_End,
      TMS_Const::TMS_Dead_Center[0]);
  TBox *xz_dead_bottom = new TBox(TMS_Const::TMS_Thin_Start,
      TMS_Const::TMS_Dead_Bottom[0],
      TMS_Const::TMS_Thick_End,
      TMS_Const::TMS_Dead_Bottom[1]);
  xz_dead_top->SetFillStyle(3003);
  xz_dead_center->SetFillStyle(3003);
  xz_dead_bottom->SetFillStyle(3003);
  xz_dead_top->SetFillColor(kGray);
  xz_dead_center->SetFillColor(kGray);
  xz_dead_bottom->SetFillColor(kGray);
  
  // And a line at the thin/thick divide
  TLine *xz_Thin_Thick = new TLine(TMS_Const::TMS_Thick_Start,
      TMS_Const::TMS_Start_Exact[0],
      TMS_Const::TMS_Thick_Start,
      TMS_Const::TMS_End_Exact[0]);
  xz_Thin_Thick->SetLineColor(kGray);
  xz_Thin_Thick->SetLineStyle(kDashed);


  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->SetLeftMargin(canv->GetLeftMargin()*1.4);
  canv->SetRightMargin(canv->GetRightMargin()*1.5);
  canv->SetTopMargin(canv->GetTopMargin()*1.2);
  // Make the output name
  std::string outputname = filename;
  while (outputname.find("/") != std::string::npos) {
    outputname = outputname.substr(outputname.find("/")+1, outputname.size());
  }
  TString Outputname = outputname.c_str();
  Outputname.ReplaceAll(".root", "_Deposits.pdf");
  canv->Print(Outputname+"[");

  TString Outputname_root = outputname.c_str();
  Outputname_root.ReplaceAll(".root", "_Deposits.root");
  TFile *outfile = new TFile(Outputname_root, "recreate");
  TTree *outtree = new TTree("energy_dep", "energy_dep");
  float pos[4];
  float energy;
  outtree->Branch("Position", pos, "Position[4]/F");
  outtree->Branch("Energy", &energy, "Energy/F");

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
  plot->GetZaxis()->SetTitleOffset(plot->GetZaxis()->GetTitleOffset()*1.3);

  for (; i < N_entries; ++i) {
    events->GetEntry(i);
    gRoo->GetEntry(i);

    plot->Reset();
    if (i % (N_entries/10) == 0) {
      std::cout << "Processed " << i << "/" << N_entries << " (" << double(i)*100./N_entries << "%)" << std::endl;
    }

    TMS_Event tms_event = TMS_Event(*event, false);
    tms_event.FillTruthFromGRooTracker(StdHepPdg, StdHepP4);

    int pdg = tms_event.GetNeutrinoPDG();
    double enu = tms_event.GetNeutrinoP4().E();
    TString reaction = tms_event.GetReaction().c_str();
    reaction.ReplaceAll(";",",");
    plot->SetTitle(Form("#splitline{Event %i, #nu PDG: %i, E_{#nu}=%.2f GeV}{%s}", i, pdg, enu, reaction.Data()));

    // Just extract the hits from the edep-sim event
    for (TG4HitSegmentDetectors::iterator jt = event->SegmentDetectors.begin(); jt != event->SegmentDetectors.end(); ++jt) {
      TG4HitSegmentContainer tms_hits = (*jt).second;
      for (TG4HitSegmentContainer::iterator kt = tms_hits.begin(); kt != tms_hits.end(); ++kt) {
        TG4HitSegment edep_hit = *kt;
        TLorentzVector Position = (edep_hit.GetStop()+edep_hit.GetStart());
        Position *= 0.5;
        plot->Fill(Position.Z(), Position.X(), edep_hit.GetEnergyDeposit());
        pos[0] = Position.X();
        pos[1] = Position.Y();
        pos[2] = Position.Z();
        pos[3] = Position.T();
        energy = edep_hit.GetEnergyDeposit();
        outtree->Fill();
      }
    }

    if (plot->Integral() == 0) continue;

    plot->SetMinimum(0);
    plot->SetMaximum(5);

    plot->Draw("colz");
    xz_box_FV->Draw("same");
    xz_box_LAr_FV->Draw("same");
    xz_box_Full->Draw("same");
    xz_box_LAr_Full->Draw("same");
    xz_dead_top->Draw("same");
    xz_dead_center->Draw("same");
    xz_dead_bottom->Draw("same");
    xz_Thin_Thick->Draw("same");

    canv->Print(Outputname);

  }

  outfile->cd();
  outtree->Write();
  outfile->Close();

  input->Close();

  canv->Print(Outputname+"]");
}
