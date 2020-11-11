#include "TMS_EventViewer.h"

TMS_EventViewer::TMS_EventViewer() :
DrawTrackFinding(false)
{

  nDraws = 0;

  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(255);

  const double zmin = (700+TMS_Const::TMS_Det_Offset[2])*10;
  const double zmax = (1500+TMS_Const::TMS_Det_Offset[2])*10;
  const double xmin = (-400+TMS_Const::TMS_Det_Offset[0])*10;
  const double xmax = (400+TMS_Const::TMS_Det_Offset[0])*10;
  const double ymin = (-250+TMS_Const::TMS_Det_Offset[1])*10;
  const double ymax = (100+TMS_Const::TMS_Det_Offset[1])*10;
  // Scint bars are 4 by 1 cm
  const int nbinsz = ((zmax-zmin)/10)/6;
  const int nbinsx = ((xmax-xmin)/10)/6;
  const int nbinsy = ((ymax-ymin)/10)/6;

  // The 2D views
  xz_view = new TH2D("TMS_Viewer_xz", "TMS viewer xz;z (mm); x (mm); Energy Deposit (MeV)", nbinsz, zmin, zmax, nbinsx, xmin, xmax);
  yz_view = new TH2D("TMS_Viewer_yz", "TMS viewer yz;z (mm); y (mm); Energy Deposit (MeV)", nbinsz, zmin, zmax, nbinsy, ymin, ymax);

  xz_view->SetMinimum(-0.01);
  yz_view->SetMinimum(-0.01);
  xz_view->SetMaximum(4);
  yz_view->SetMaximum(4);

  yz_view->GetZaxis()->SetTitleOffset(yz_view->GetZaxis()->GetTitleOffset()*1.2);
  xz_view->GetZaxis()->SetTitleOffset(xz_view->GetZaxis()->GetTitleOffset()*1.2);

  yz_view->GetYaxis()->SetTitleOffset(1.8);
  xz_view->GetYaxis()->SetTitleOffset(1.8);

  yz_view->GetYaxis()->SetMaxDigits(3);
  yz_view->GetXaxis()->SetMaxDigits(3);
  xz_view->GetYaxis()->SetMaxDigits(3);
  xz_view->GetXaxis()->SetMaxDigits(3);

  TGaxis::SetExponentOffset(0, -0.031, "x");

  // The canvas
  Canvas = new TCanvas("TMS_EventViewer", "TMS_EventViewer", 1024, 1024);
  Canvas->Divide(2);
  Canvas->cd(1)->SetLeftMargin(Canvas->GetLeftMargin()*1.2);
  Canvas->cd(2)->SetLeftMargin(Canvas->GetLeftMargin()*1.2);
  Canvas->cd(1)->SetRightMargin(Canvas->GetRightMargin()*1.5);
  Canvas->cd(2)->SetRightMargin(Canvas->GetRightMargin()*1.5);

  // Full view from inspecting all hits
  xz_box_Full = new TBox((730+TMS_Const::TMS_Det_Offset[2])*10,
      (-348.5+TMS_Const::TMS_Det_Offset[0])*10, 
      (1415+TMS_Const::TMS_Det_Offset[2])*10, 
      (348.5+TMS_Const::TMS_Det_Offset[0])*10);
  xz_box_Full->SetLineColor(kGreen);
  xz_box_Full->SetFillStyle(0);

  yz_box_Full = new TBox((730+TMS_Const::TMS_Det_Offset[2])*10, 
      (-234+TMS_Const::TMS_Det_Offset[1])*10, 
      (1415+TMS_Const::TMS_Det_Offset[2])*10, 
      (87+TMS_Const::TMS_Det_Offset[1])*10);
  yz_box_Full->SetLineColor(kGreen);
  yz_box_Full->SetFillStyle(0);

  // FV just taking 50 cm in from the full
  xz_box_FV = new TBox(xz_box_Full->GetX1(),
      xz_box_Full->GetY1()+500,
      xz_box_Full->GetX2()-500,
      xz_box_Full->GetY2()-500);
  xz_box_FV->SetLineColor(kRed);
  xz_box_FV->SetLineStyle(kDashed);
  xz_box_FV->SetFillStyle(0);

  yz_box_FV = new TBox(yz_box_Full->GetX1(),
      yz_box_Full->GetY1()+500,
      yz_box_Full->GetX2()-500,
      yz_box_Full->GetY2()-500);
  yz_box_FV->SetLineColor(kRed);
  yz_box_FV->SetLineStyle(kDashed);
  yz_box_FV->SetFillStyle(0);

  // Include the dead region boxes
  xz_dead_top = new TBox((730+TMS_Const::TMS_Det_Offset[2])*10, 
      (171.7+TMS_Const::TMS_Det_Offset[0])*10, 
      (1415+TMS_Const::TMS_Det_Offset[2])*10, 
      (180.4+TMS_Const::TMS_Det_Offset[0])*10);
  xz_dead_center = new TBox((730+TMS_Const::TMS_Det_Offset[2])*10, 
      (-3.3+TMS_Const::TMS_Det_Offset[0])*10, 
      (1415+TMS_Const::TMS_Det_Offset[2])*10, 
      (3.3+TMS_Const::TMS_Det_Offset[0])*10);
  xz_dead_bottom = new TBox((730+TMS_Const::TMS_Det_Offset[2])*10, 
      (-180.4+TMS_Const::TMS_Det_Offset[0])*10, 
      (1415+TMS_Const::TMS_Det_Offset[2])*10, 
      (-171.7+TMS_Const::TMS_Det_Offset[0])*10);
  xz_dead_top->SetFillStyle(3003);
  xz_dead_center->SetFillStyle(3003);
  xz_dead_bottom->SetFillStyle(3003);
  xz_dead_top->SetFillColor(kGray);
  xz_dead_center->SetFillColor(kGray);
  xz_dead_bottom->SetFillColor(kGray);
  
  // And a line at the thin/thick divide
  xz_Thin_Thick = new TLine((TMS_Const::TMS_Trans_Start+TMS_Const::TMS_Det_Offset[2])*10, 
      (-348.5+TMS_Const::TMS_Det_Offset[0])*10, 
      (TMS_Const::TMS_Trans_Start+TMS_Const::TMS_Det_Offset[2])*10, 
      (348.5+TMS_Const::TMS_Det_Offset[0])*10);
  xz_Thin_Thick->SetLineColor(kGray);
  xz_Thin_Thick->SetLineStyle(kDashed);

  yz_Thin_Thick = new TLine((TMS_Const::TMS_Trans_Start+TMS_Const::TMS_Det_Offset[2])*10, 
      (-234+TMS_Const::TMS_Det_Offset[1])*10, 
      (TMS_Const::TMS_Trans_Start+TMS_Const::TMS_Det_Offset[2])*10, 
      (87+TMS_Const::TMS_Det_Offset[1])*10);
  yz_Thin_Thick->SetLineColor(kGray);
  yz_Thin_Thick->SetLineStyle(kDashed);

  // Open up the pdf
  CanvasName = "TMS_EventViewer_Collection_Hough";
  Canvas->Print(CanvasName+".pdf[");
}

// Draw the finished event
void TMS_EventViewer::Draw(TMS_Event &event) {

  xz_view->Reset();
  yz_view->Reset();

  int EventNumber = event.GetEventNumber();
  xz_view->SetTitle(Form("TMS viewer xz, Event %i", EventNumber));
  yz_view->SetTitle(Form("TMS viewer yz, Event %i", EventNumber));

  std::vector<TMS_Hit> TMS_Hits = event.GetHits();

  // Check that there are hits
  if (TMS_Hits.size() < 50) {
    //std::cerr << "Trying to draw an event that has no hits in the TMS, returning..." << std::endl;
    return;
  }

  // Loop over the hits and add them
  for (std::vector<TMS_Hit>::iterator it = TMS_Hits.begin(); it != TMS_Hits.end(); ++it) {
    TMS_Bar bar = (*it).GetBar();  
    double x = bar.GetX();
    double y = bar.GetY();
    double z = bar.GetZ();
    int BarType = bar.GetBarType();
    double e = (*it).GetE();

    // Bar along y (no x info)
    if (BarType == TMS_Bar::BarType::kYBar) {
      xz_view->Fill(z, x, e);
    }
    // Bar along x (no y info)
    else if (BarType == TMS_Bar::BarType::kXBar) {
      yz_view->Fill(z, y, e);
    }
  }

  // Loop over the reconstructed tracks
  std::vector<TMS_Hit> Candidates = TMS_TrackFinder::GetFinder().GetCandidates();
  //std::vector<std::vector<TMS_Hit> > TotalCandidates = TMS_TrackFinder::GetFinder().GetTotalCandidates();

  //int iter = 0;
  // Loop over each total candidates
  //for (std::vector<std::vector<TMS_Hit> >::iterator it = TotalCandidates.begin(); it != TotalCandidates.end(); ++it) {
    //std::vector<TMS_Hit> Candidates = (*it);

    //double e = 3;
    //if (iter == 0) e = 100;
    //else e = 2;
    double e = 4;
    for (std::vector<TMS_Hit>::iterator jt = Candidates.begin(); jt != Candidates.end(); ++jt) {

      TMS_Bar bar = (*jt).GetBar();  
      double x = bar.GetX();
      double y = bar.GetY();
      double z = bar.GetZ();
      int BarType = bar.GetBarType();

      // Bar along y (no x info)
      if (BarType == TMS_Bar::BarType::kYBar) {
        xz_view->Fill(z, x, e);
      }
      // Bar along x (no y info)
      else if (BarType == TMS_Bar::BarType::kXBar) {
        yz_view->Fill(z, y, e);
      }

    }
    //iter++;
  //}


  gStyle->SetPalette(kBird);
  Canvas->cd(1); 
  xz_view->Draw("colz");
  xz_box_FV->Draw("same");
  xz_box_Full->Draw("same");
  xz_dead_top->Draw("same");
  xz_dead_center->Draw("same");
  xz_dead_bottom->Draw("same");
  xz_Thin_Thick->Draw("same");
  std::vector<TF1*> HoughLines_zx = TMS_TrackFinder::GetFinder().GetHoughLines_zx();
  for (std::vector<TF1*>::iterator it = HoughLines_zx.begin(); it != HoughLines_zx.end(); ++it) (*it)->Draw("same");

  Canvas->cd(2); 
  yz_view->Draw("colz");
  yz_box_FV->Draw("same");
  yz_box_Full->Draw("same");
  yz_Thin_Thick->Draw("same");
  std::vector<TF1*> HoughLines_zy = TMS_TrackFinder::GetFinder().GetHoughLines_zy();
  for (std::vector<TF1*>::iterator it = HoughLines_zy.begin(); it != HoughLines_zy.end(); ++it) (*it)->Draw("same");

  Canvas->Print(CanvasName+".pdf");

  if (DrawTrackFinding) {
    gStyle->SetPalette(87);
    Canvas->cd(1);
    TH2D *accumulator_xz = TMS_TrackFinder::GetFinder().AccumulatorToTH2D(false);
    accumulator_xz->Draw("colz");

    Canvas->cd(2);
    TH2D *accumulator_yz = TMS_TrackFinder::GetFinder().AccumulatorToTH2D(true);
    accumulator_yz->Draw("colz");
    Canvas->Print(CanvasName+".pdf");

    delete accumulator_xz;
    delete accumulator_yz;
  }

  nDraws++;
}

