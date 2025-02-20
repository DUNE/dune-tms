#include "paul_tol_colors.hpp"
void efficiency() {

  TFile *f = new TFile("neutrino_merge_0p1Ecut_updtrue.root");
  TTree *t = (TTree*)f->Get("Truth_Info");
  TTree *r = (TTree*)f->Get("Line_Candidates");
  t->SetBranchStatus("*", false);
  r->SetBranchStatus("*", false);

  std::string *Interaction;
  float Muon_Vertex[4];
  float Muon_Death[4];
  float Muon_TrueKE;
  t->SetBranchStatus("Interaction", true);
  t->SetBranchAddress("Interaction", &Interaction);
  t->SetBranchStatus("Muon_Vertex", true);
  t->SetBranchAddress("Muon_Vertex", Muon_Vertex);
  t->SetBranchStatus("Muon_Death", true);
  t->SetBranchAddress("Muon_Death", Muon_Death);
  t->SetBranchStatus("Muon_TrueKE", true);
  t->SetBranchAddress("Muon_TrueKE", &Muon_TrueKE);

  int nLines;
  int nClusters;
  r->SetBranchStatus("nLines", true);
  r->SetBranchAddress("nLines", &nLines);
  r->SetBranchStatus("nClusters", true);
  r->SetBranchAddress("nClusters", &nClusters);

  // Just get the interaction string
  int nCC = 0;
  int nNC = 0;

  int nCCwTrueMuon = 0;

  int nCCwTrack = 0;
  int nCCwClust = 0;
  int nCCwTrackOrClust = 0;

  int nNCwTrack = 0;
  int nNCwClust = 0;
  int nNCwTrackOrClust = 0;

  int nEntries = t->GetEntries();

  // let's make a TH1D too in muon true KE
  TH1D *AllCC = new TH1D("CC", "All CC", 25, 0, 5);
  TH1D *AllTrack = new TH1D("AllTrack", "All CC with track", 25, 0, 5);
  TH1D *AllClust = new TH1D("AllClust", "All CC with cluster", 25, 0, 5);
  TH1D *AllTrackOrClust = new TH1D("AllTrackOrClust", "All CC with track or cluster", 25, 0, 5);

  for (int i = 0; i < nEntries; ++i) {
    t->GetEntry(i);
    r->GetEntry(i);

    bool isCC = true;
    if ((*Interaction).find("[CC]") != std::string::npos) isCC = true;
    else if ((*Interaction).find("[NC]") != std::string::npos) isCC = false;
    Muon_TrueKE *= 1.E-3;
    if (Muon_TrueKE < 0) Muon_TrueKE = 0;

    if (isCC) {
      nCC++;
      if (Muon_Vertex[2] > 0) {
        nCCwTrueMuon++;
        AllCC->Fill(Muon_TrueKE);
      }

      if (nLines > 0) {
        nCCwTrack++;
        AllTrack->Fill(Muon_TrueKE);
      }

      if (nClusters > 0) {
        nCCwClust++;
        AllClust->Fill(Muon_TrueKE);
      }

      if (nLines > 0 || nClusters > 0) {
        nCCwTrackOrClust++;
        AllTrackOrClust->Fill(Muon_TrueKE);
      }

    } else {
      nNC++;
      if (nLines > 0) nNCwTrack++;
      if (nClusters > 0) nNCwClust++;
      if (nLines > 0 || nClusters > 0) nNCwTrackOrClust++;
    }

  }

  std::cout << "Truth: " << std::endl;
  std::cout << nCC << "/" << nEntries << " CC events (" << double(nCC)*100./nEntries << "%)" << std::endl;
  std::cout << nNC << "/" << nEntries << " NC events (" << double(nNC)*100./nEntries << "%)" << std::endl;

  std::cout << nCCwTrueMuon << "/" << nEntries << " CC events with true muon (" << double(nCCwTrueMuon)*100./nEntries << "%)" << std::endl;

  std::cout << std::endl;
  std::cout << "Reco: " << std::endl;
  std::cout << nCCwTrack << "/" << nCC << " CC events with track (" << double(nCCwTrack)*100./nCC << "%)" << std::endl;
  std::cout << nCCwClust << "/" << nCC << " CC events with cluster (" << double(nCCwClust)*100./nCC << "%)" << std::endl;
  std::cout << nCCwTrackOrClust << "/" << nCC << " CC events with track or cluster (" << double(nCCwTrackOrClust)*100./nCC << "%)" << std::endl;

  std::cout << std::endl;
  std::cout << nNCwTrack << "/" << nNC << " NC events with track (" << double(nNCwTrack)*100./nNC << "%)" << std::endl;
  std::cout << nNCwClust << "/" << nNC << " NC events with cluster (" << double(nNCwClust)*100./nNC << "%)" << std::endl;
  std::cout << nNCwTrackOrClust << "/" << nNC << " NC events with track or cluster (" << double(nNCwTrackOrClust)*100./nNC << "%)" << std::endl;

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->SetLeftMargin(canv->GetLeftMargin()*1.4);
  canv->Print("tracking.pdf[");

  AllCC->GetYaxis()->SetTitleOffset(AllCC->GetYaxis()->GetTitleOffset()*1.4);
  AllCC->GetYaxis()->SetRangeUser(0, AllCC->GetMaximum()*1.2);
  AllCC->Draw();
  AllCC->GetXaxis()->SetTitle("Muon true KE (GeV)");
  AllCC->GetYaxis()->SetTitle("Number of events");
  AllTrack->Draw("same");
  AllClust->Draw("same");
  AllTrackOrClust->Draw("same");
  TLegend *leg = canv->BuildLegend();
  leg->SetX1(0.55);
  leg->SetX2(0.9);
  leg->SetY1(0.4);
  leg->SetY2(0.9);

  AllTrack->SetLineColor(tolcols::kTBriBlue);
  AllClust->SetLineColor(tolcols::kTBriGreen);
  AllTrackOrClust->SetLineColor(tolcols::kTBriRed);

  AllCC->SetMarkerSize(0);
  AllTrack->SetMarkerSize(0);
  AllClust->SetMarkerSize(0);
  AllTrackOrClust->SetMarkerSize(0);

  AllCC->SetFillStyle(0);
  AllTrack->SetFillStyle(0);
  AllClust->SetFillStyle(0);
  AllTrackOrClust->SetFillStyle(0);

  AllCC->SetLineWidth(2);
  AllTrack->SetLineWidth(2);
  AllClust->SetLineWidth(2);
  AllTrackOrClust->SetLineWidth(2);

  AllTrack->SetLineStyle(7);
  AllClust->SetLineStyle(5);
  AllTrackOrClust->SetLineStyle(2);

  canv->Print("tracking.pdf");

  AllTrack->Divide(AllCC);
  AllClust->Divide(AllCC);
  AllTrackOrClust->Divide(AllCC);

  AllTrack->GetXaxis()->SetTitle(AllCC->GetXaxis()->GetTitle());
  AllTrack->GetYaxis()->SetTitle("Ratio relative all CC");
  AllTrack->GetYaxis()->SetRangeUser(0, 1.1);

  AllTrack->Draw();
  AllClust->Draw("same");
  AllTrackOrClust->Draw("same");
  TLine *line = new TLine(0, 1, 5, 1);
  line->SetLineColor(kBlack);
  line->SetLineWidth(3);
  line->SetLineStyle(kDashed);

  leg = canv->BuildLegend();
  leg->SetX1(0.25);
  leg->SetX2(0.9);
  leg->SetY1(0.15);
  leg->SetY2(0.45);

  leg->Draw("same");
  line->Draw("same");
  canv->Print("tracking.pdf");

  canv->Print("tracking.pdf]");
}

int main() {
  efficiency();
  return 0;
}
