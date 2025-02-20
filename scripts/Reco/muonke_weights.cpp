#include "paul_tol_colors.hpp"

// Helper to calculate Poisson
double CalcPoisson(TH1D *data, TH1D *mc) {
  double llh = 0;
  for (int i = 0; i < data->GetXaxis()->GetNbins(); ++i) {
    double data_c = data->GetBinContent(i+1);
    double mc_c = mc->GetBinContent(i+1);
    if (mc_c > 0 && data_c > 0) {
      llh += (mc_c - data_c + data_c * TMath::Log(data_c/mc_c));
    } else if (mc_c > 0 && data_c == 0) {
      llh += mc_c;
    }
  }
  return 2*llh;
}

void muonke_weights(std::string filename, 
    int nLinesCut = 1, 
    int nClustersCut = 0, 
    float OccupancyCut = 0.50, 
    bool CCmuOnly = false, 
    float ClusterEnergyCut = 0, 
    bool AtLeastOneLine = true, 
    bool AllDet = true) {
  /*
  int nLinesCut = 1; // How many lines can our events have?
  int nClustersCut = 0; // Only allow for this many clusters
  float OccupancyCut = 0.50; // Only allow for higher occupancy tracks than this
  bool CCmuOnly = false; // Only include events with a true CC muon (not necessarily selected as the track, but present in the event)
  float ClusterEnergyCut = 0; // How many MeV of energy in all clusters do we cut on (greater than this number gets excluded)
  bool AtLeastOneLine = true; // Do we require at least one line? (necessary for track length measurement)
  bool AllDet = true; // Muon starts in the whole detector? Or just thin region
  */

  TFile *f = new TFile(filename.c_str());
  TTree *truth = (TTree*)f->Get("Truth_Info");
  TTree *reco = (TTree*)f->Get("Line_Candidates");

  // Reco info
  const int MAX_LINES = 20;
  const int MAX_HITS = 1500;
  const int MAX_CLUSTERS = 25;
  int nLines;
  //float TotalTrackEnergy[MAX_LINES];
  float TrackLength[MAX_LINES];
  float Occupancy[MAX_LINES];
  float FirstHoughHit[MAX_LINES][2];
  float LastHoughHit[MAX_LINES][2];
  int EventNum_reco;

  //float RecoHitEnergy[MAX_HITS];
  float ClusterEnergy[MAX_CLUSTERS];
  int nClusters;

  reco->SetBranchStatus("*", false);

  reco->SetBranchStatus("nLines", true);
  reco->SetBranchAddress("nLines", &nLines);
  //reco->SetBranchStatus("TotalTrackEnergy", true);
  //reco->SetBranchAddress("TotalTrackEnergy", TotalTrackEnergy);
  reco->SetBranchStatus("TrackLength", true);
  reco->SetBranchAddress("TrackLength", TrackLength);
  reco->SetBranchStatus("Occupancy", true);
  reco->SetBranchAddress("Occupancy", Occupancy);
  reco->SetBranchStatus("FirstHoughHit", true);
  reco->SetBranchAddress("FirstHoughHit", FirstHoughHit);
  reco->SetBranchStatus("LastHoughHit", true);
  reco->SetBranchAddress("LastHoughHit", LastHoughHit);
  //reco->SetBranchStatus("RecoHitEnergy", true);
  //reco->SetBranchAddress("RecoHitEnergy", RecoHitEnergy);

  reco->SetBranchStatus("ClusterEnergy", true);
  reco->SetBranchAddress("ClusterEnergy", ClusterEnergy);
  reco->SetBranchStatus("nClusters", true);
  reco->SetBranchAddress("nClusters", &nClusters);

  reco->SetBranchStatus("EventNo", true);
  reco->SetBranchAddress("EventNo", &EventNum_reco);

  // Truth info
  float MuonP4[4];
  float Muon_Vertex[4];
  float Muon_Death[4];
  float Muon_TrueKE;
  int nParticles;
  int NeutrinoPDG;
  float NeutrinoP4[4];
  int EventNum_true;

  truth->SetBranchStatus("*", false);

  truth->SetBranchStatus("MuonP4", true);
  truth->SetBranchAddress("MuonP4", MuonP4);
  truth->SetBranchStatus("Muon_Vertex", true);
  truth->SetBranchAddress("Muon_Vertex", Muon_Vertex);
  truth->SetBranchStatus("Muon_Death", true);
  truth->SetBranchAddress("Muon_Death", Muon_Death);
  truth->SetBranchStatus("Muon_TrueKE", true);
  truth->SetBranchAddress("Muon_TrueKE", &Muon_TrueKE);
  //truth->SetBranchStatus("nParticles", true);
  //truth->SetBranchAddress("nParticles", &nParticles);
  //truth->SetBranchStatus("NeutrinoPDG", true);
  //truth->SetBranchAddress("NeutrinoPDG", &NeutrinoPDG);
  truth->SetBranchStatus("NeutrinoP4", true);
  truth->SetBranchAddress("NeutrinoP4", NeutrinoP4);
  truth->SetBranchStatus("EventNo", true);
  truth->SetBranchAddress("EventNo", &EventNum_true);

  // Set up the weights file
  TFile *weightfile = new TFile("tms_beamMonitoring_weights.root");
  const int nweights = 11;
  TH1D *weights[nweights];
  // All histograms are in true Enu, in GeV
  // Convert into TGraphs so can use interpolation?
  weights[0] = (TH1D*)weightfile->Get("h_xProjtargetDensity");
  weights[1] = (TH1D*)weightfile->Get("h_xProjBeamOffsetX;1");
  weights[2] = (TH1D*)weightfile->Get("h_xProjBeamTheta;1");
  weights[3] = (TH1D*)weightfile->Get("h_xProjBeamThetaPhi;1");
  weights[4] = (TH1D*)weightfile->Get("h_xProjHC;1");
  weights[5] = (TH1D*)weightfile->Get("h_xProjWL;1");
  weights[6] = (TH1D*)weightfile->Get("h_xProjDPR;1");
  weights[7] = (TH1D*)weightfile->Get("h_xProjHorn1_XShift;1");
  weights[8] = (TH1D*)weightfile->Get("h_xProjHorn1_YShift;1");
  weights[9] = (TH1D*)weightfile->Get("h_xProjHorn2_XShift;1");
  weights[10] = (TH1D*)weightfile->Get("h_xProjHorn2_YShift;1");

  TH2D *KE = new TH2D("KE", "KE;True muon KE (MeV);Track length of best track (g/cm^{2})", 100, 0, 5000, 50, 0, 2500);
  TH1D *h_Occupancy = new TH1D("Occ", "Occupancy; Occupancy of longest track; Number of events", 110, 0, 1.1);

  TH2D *KEest = new TH2D("KEest", "KE estimator; True muon KE (MeV); KE estimate", 50, 0, 5000, 50, 0, 5000);
  TH1D *TrueHad = new TH1D("TrueEHad", "True E_{had}; True E_{#nu}-E_{#mu} (GeV); N_{events}", 100, 0, 1);
  KE->GetYaxis()->SetTitleOffset(KE->GetYaxis()->GetTitleOffset()*1.4);
  KE->GetZaxis()->SetTitleOffset(KE->GetZaxis()->GetTitleOffset()*1.4);

  KEest->GetYaxis()->SetTitleOffset(KEest->GetYaxis()->GetTitleOffset()*1.4);
  KEest->GetZaxis()->SetTitleOffset(KEest->GetZaxis()->GetTitleOffset()*1.4);

  TH1D *MuonKEreco = new TH1D("MuonKEreco", "Muon KE reco;Muon reco KE (GeV); Number of events", 50, 0, 5);
  TH1D *MuonKEreco_w[nweights];
  for (int i = 0; i < nweights; ++i) {
    MuonKEreco_w[i] = new TH1D(Form("MuonKEreco_w_%i", i), Form("%s;Muon reco KE (GeV); Number of events, weighted", weights[i]->GetName()), 50, 0, 5);
    MuonKEreco_w[i]->SetLineColor(10000+i);
    MuonKEreco_w[i]->SetMarkerSize(0);
    MuonKEreco_w[i]->SetFillStyle(0);
    if (i > 5) MuonKEreco_w[i]->SetLineStyle(kDashed);
  }
  MuonKEreco->SetLineColor(kBlack);
  MuonKEreco->SetLineStyle(kDashed);
  MuonKEreco->SetMarkerSize(0);
  MuonKEreco->SetFillStyle(0);

  TH1D *hTrueEnu = new TH1D("hTrueEnu", "True E_{#nu};True E_{#nu} (GeV); Number of events", 50, 0, 10);
  TH1D *hTrueEnu_w[nweights];
  for (int i = 0; i < nweights; ++i) {
    hTrueEnu_w[i] = new TH1D(Form("hTrueEnu_w_%i", i), Form("%s;True E_{#nu} (GeV); Number of events, weighted", weights[i]->GetName()), 50, 0, 10);
    hTrueEnu_w[i]->SetLineColor(10000+i);
    hTrueEnu_w[i]->SetMarkerSize(0);
    hTrueEnu_w[i]->SetFillStyle(0);
    if (i > 5) hTrueEnu_w[i]->SetLineStyle(kDashed);
  }
  hTrueEnu->SetLineColor(kBlack);
  hTrueEnu->SetLineStyle(kDashed);
  hTrueEnu->SetMarkerSize(0);
  hTrueEnu->SetFillStyle(0);

  int ngood = 0;
  int nentries = truth->GetEntries();
  int trklen_counter = 0;
  std::cout << nentries << " events..." << std::endl;
  int true_entry = 0;
  int reco_entry = 0;
  for (int i = 0; i < nentries; ++i, ++true_entry, ++reco_entry) {
    truth->GetEntry(true_entry);
    reco->GetEntry(reco_entry);

    if (EventNum_reco != EventNum_true) {
      std::cout << "Event " << i << " " << EventNum_reco << " " << EventNum_true << std::endl;
      std::cout << true_entry << " " << reco_entry << std::endl;
      true_entry++;
      continue;
    }

    if (i % int(nentries/100.) == 0) std::cout << "Event " << i << std::endl;

    // The event has to have a muon
    if (Muon_TrueKE < 0) continue;

    // Ask for a true muon in the event if requested
    if (CCmuOnly && Muon_Vertex[2] < 0) continue;

    // Only include events with lines
    if (AtLeastOneLine && nLines < 1) continue;

    // Run the cuts
    if (nLines > nLinesCut) continue;
    if (nClusters > nClustersCut) continue;
    // Sum up the total cluster energy
    float cluster_en = 0;
    for (int j = 0; j < nClusters; ++j) {
      cluster_en += ClusterEnergy[j];
    }
    if (cluster_en > ClusterEnergyCut) continue;

    // Find the best track
    int besttrack = 0;
    for (int j = 0; j < nLines; ++j) {
      if (Occupancy[j] > Occupancy[besttrack]) besttrack = j;
    }

    // Also check the track with the longest track length
    int lon_trklen = 0;
    for (int j = 0; j < nLines; ++j) {
      if (TrackLength[j] > TrackLength[lon_trklen]) lon_trklen = j;
    }

    // And also check longest track
    float longest = 0;
    int longtrack = 0;
    for (int j = 0; j < nLines; ++j) {
      float xdist = LastHoughHit[j][0]-FirstHoughHit[j][0];
      float ydist = LastHoughHit[j][1]-FirstHoughHit[j][1];
      float dist = sqrt(xdist*xdist+ydist*ydist);
      if (dist > longest) longtrack = j;
    }

    if (longtrack != lon_trklen) {
      trklen_counter++;
    }
    longtrack = lon_trklen;

    // Look only at events with true muons that die inside the detector
    if (Muon_Death[1] > 1159 || Muon_Death[1] < -3864) continue;

    // Check that the longest track stops in the detector, and starts in the detector FV
    if (AllDet) {
      if (FirstHoughHit[longtrack][0] < 11185+55*2) continue;
    } else {
      if (FirstHoughHit[longtrack][0] < 11185+55*2 || FirstHoughHit[longtrack][0] > 14435) continue;  // Changed to new TMS_Thick_Start here
    }
    if (LastHoughHit[longtrack][0] > 18535-80*2) continue;
    //if (LastHoughHit[longtrack][0] > 14435) continue; // Same here

    // 20 cm inwards
    if (fabs(FirstHoughHit[longtrack][1]) > 3520-200) continue;
    if (fabs(LastHoughHit[longtrack][1]) > 3520-200) continue;

    h_Occupancy->Fill(Occupancy[longtrack]);

    // Ask for only small amount of other energy deposits
    if (Occupancy[longtrack] < OccupancyCut) continue;

    float best_tracklength = TrackLength[longtrack];
    KEest->Fill(Muon_TrueKE, 82+1.75*best_tracklength);
    KE->Fill(Muon_TrueKE, best_tracklength);

    MuonKEreco->Fill((82+1.75*best_tracklength)/1.E3);
    // Get the weighted
    // First true Enu, in GeV
    double TrueEnu = NeutrinoP4[3];
    for (int j = 0; j < nweights; ++j) {
      double weight = weights[j]->GetBinContent(weights[j]->FindBin(TrueEnu));
      MuonKEreco_w[j]->Fill((82+1.75*best_tracklength)/1.E3, weight);
      //MuonKEreco_w[j]->Fill((88.6+1.72*best_tracklength)/1.E3, weight);
      hTrueEnu_w[j]->Fill(TrueEnu, weight);
    }
    hTrueEnu->Fill(TrueEnu);

    TrueHad->Fill(NeutrinoP4[3] - MuonP4[3]/1000.);
    ngood++;
  }
  std::cout << ngood << "/" << nentries << std::endl;
  std::cout << trklen_counter << std::endl;

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->SetLeftMargin(canv->GetLeftMargin()*1.5);
  canv->SetRightMargin(canv->GetRightMargin()*1.4);

  // Fix filename

  while (filename.find("/") != std::string::npos) {
    filename = filename.substr(filename.find("/")+1, filename.size());
  }

  TString canvname = Form("MuonKE_%s", filename.c_str());
  canvname += Form("_nLinesCut%i_nClusterCut%i_OccupancyCut%.2f_CCmuOnly%o_ClusterEnCut%.2f_AtLeastOneLine%o_AllDet%o", nLinesCut, nClustersCut, OccupancyCut, CCmuOnly, ClusterEnergyCut, AtLeastOneLine, AllDet);
  canvname += ".pdf";
  canv->Print(canvname+"[");

  gStyle->SetPalette(55);
  KE->Draw("colz");
  canv->Print(canvname);
  h_Occupancy->Draw();
  canv->Print(canvname);

  TH1D *arith = new TH1D("arith", "arith", KE->GetXaxis()->GetNbins(), KE->GetXaxis()->GetBinLowEdge(1), KE->GetXaxis()->GetBinLowEdge(KE->GetXaxis()->GetNbins()+1));
  // Now make the muon KE
  gStyle->SetOptStat(1111);
  for (int i = 0; i < KE->GetXaxis()->GetNbins(); ++i) {
    double center = KE->GetXaxis()->GetBinCenter(i);
    TH1D *proj = KE->ProjectionY(Form("KE %.2f", center), i, i);
    double mean = proj->GetMean();
    double mode = proj->GetBinCenter(proj->GetMaximumBin());
    double rms = proj->GetRMS();
    arith->SetBinContent(i+1, mode);
    arith->SetBinError(i+1, rms);
    proj->Draw();
    proj->SetStats(1);
    canv->Print(canvname);
  }
  TF1 *fit = new TF1("fit", "[0]+[1]*x", KE->GetYaxis()->GetBinLowEdge(1), KE->GetYaxis()->GetBinLowEdge(KE->GetYaxis()->GetNbins()+1));
  gStyle->SetOptStat(0);
  KE->Draw("colz");
  arith->Draw("same");
  arith->Fit(fit, "S", "", 700, 2500);
  TLegend *leg = new TLegend(0.2, 0.5, 0.6, 0.9);
  leg->AddEntry(arith, "Artihmetic mean and RMS", "le");
  leg->AddEntry(fit, Form("y=%.2f x + %.2f", fit->GetParameter(1), fit->GetParameter(0)), "l");
  fit->SetLineColor(kRed);
  arith->GetXaxis()->SetTitle("Muon True KE");
  leg->Draw("same");

  canv->Print(canvname);

  KEest->Draw("colz");
  canv->Print(canvname);

  TrueHad->Draw("hist");
  canv->Print(canvname);

  MuonKEreco->Draw("hist");
  for (int i = 0; i < nweights; ++i) {
    double n2llh = CalcPoisson(MuonKEreco_w[i], MuonKEreco);
    MuonKEreco_w[i]->SetTitle(Form("%s #chi^{2}=%.2f", MuonKEreco_w[i]->GetTitle(), n2llh));
    MuonKEreco_w[i]->Draw("same,hist");
  }
  canv->BuildLegend();
  MuonKEreco->Draw("same,hist");
  canv->Print(canvname);

  TH1D *ratios[nweights];
  for (int i = 0; i < nweights; ++i) {
    ratios[i] = (TH1D*)MuonKEreco_w[i]->Clone(Form("%s_r", MuonKEreco_w[i]->GetName()));
    ratios[i]->Divide(MuonKEreco);
  }

  for (int i = 0; i < nweights; ++i) {
    if (i == 0) ratios[i]->Draw("hist");
    else ratios[i]->Draw("same,hist");
  }
  ratios[0]->GetYaxis()->SetRangeUser(0.8, 1.2);
  ratios[0]->GetYaxis()->SetTitle("Ratio to nom.");
  TLine *line = new TLine(ratios[0]->GetXaxis()->GetBinLowEdge(1), 1, ratios[0]->GetXaxis()->GetBinLowEdge(ratios[0]->GetXaxis()->GetNbins()+1), 1);
  line->SetLineWidth(2);
  line->SetLineColor(kRed);
  line->SetLineStyle(kDashed);
  canv->BuildLegend();
  line->Draw("same");
  canv->Print(canvname);

  // Write the enu histograms
  hTrueEnu->Draw("hist");
  for (int i = 0; i < nweights; ++i) {
    double n2llh = CalcPoisson(hTrueEnu_w[i], hTrueEnu);
    hTrueEnu_w[i]->SetTitle(Form("%s #chi^{2}=%.2f", hTrueEnu_w[i]->GetTitle(), n2llh));
    hTrueEnu_w[i]->Draw("same,hist");
  }
  canv->BuildLegend();
  hTrueEnu->Draw("same,hist");
  canv->Print(canvname);

  TH1D *ratioEnu[nweights];
  for (int i = 0; i < nweights; ++i) {
    ratioEnu[i] = (TH1D*)hTrueEnu_w[i]->Clone(Form("%s_r", hTrueEnu_w[i]->GetName()));
    ratioEnu[i]->Divide(hTrueEnu);
  }

  for (int i = 0; i < nweights; ++i) {
    if (i == 0) ratioEnu[i]->Draw("hist");
    else ratioEnu[i]->Draw("same,hist");
  }
  ratioEnu[0]->GetYaxis()->SetRangeUser(0.8, 1.2);
  ratioEnu[0]->GetYaxis()->SetTitle("Ratio to nom.");
  TLine *line2 = new TLine(ratioEnu[0]->GetXaxis()->GetBinLowEdge(1), 1, ratioEnu[0]->GetXaxis()->GetBinLowEdge(ratioEnu[0]->GetXaxis()->GetNbins()+1), 1);
  line2->SetLineWidth(2);
  line2->SetLineColor(kRed);
  line2->SetLineStyle(kDashed);
  canv->BuildLegend();
  line2->Draw("same");
  canv->Print(canvname);

  // Write all historams to a file too
  TString outname = canvname;
  outname = outname.ReplaceAll(".pdf", ".root");
  TFile *output = new TFile(outname, "recreate");

  KE->Write();
  h_Occupancy->Write();
  KEest->Write();
  TrueHad->Write();
  MuonKEreco->Write();
  hTrueEnu->Write();
  for (int i = 0; i < nweights; ++i) {
    MuonKEreco_w[i]->Write();
    ratios[i]->Write();
    hTrueEnu_w[i]->Write();
    ratioEnu[i]->Write();
  }
  output->Close();

  canv->Print(canvname+"]");
}
