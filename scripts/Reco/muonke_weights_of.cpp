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

void muonke_weights_of(std::string filename,
    int nLinesCut = 100,
    int nClustersCut = 100,
    float OccupancyCut = 0.50,
    bool CCmuOnly = false,
    float ClusterEnergyCut = 100,
    bool AtLeastOneLine = true,
    bool AllDet = true) {
  /*
  int nLinesCut = 100; // How many lines can our events have?
  int nClustersCut = 100; // Only allow for this many clusters
  float OccupancyCut = 0.50; // Only allow for higher occupancy tracks than this
  bool CCmuOnly = false; // Only include events with a true CC muon (not necessarily selected as the track, but present in the event)
  float ClusterEnergyCut = 1000; // How many MeV of energy in all clusters do we cut on (greater than this number gets excluded)
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
  //int FirstHoughPlane[MAX_LINES];
  //int LastHoughPlane[MAX_LINES];
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

  //reco->SetBranchStatus("FirstHoughPlane", true);
  //reco->SetBranchAddress("FirstHoughPlane", FirstHoughPlane);
  //reco->SetBranchStatus("LastHoughPlane", true);
  //reco->SetBranchAddress("LastHoughPlane", LastHoughPlane);

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

  TH1D *h_Occupancy = new TH1D("Occ", "Occupancy; Occupancy of longest track; Number of events", 110, 0, 1.1);
  TH1D *TrueEnu = new TH1D("TrueEnu", "True E_{#nu};True E_{#nu} (GeV)", 200, 0, 20);
  TH1D *TrueHad = new TH1D("TrueEHad", "True E_{had}; True E_{#nu}-E_{#mu} (GeV); N_{events}", 300, 0, 3);
  TH1D *TrueMuKE = new TH1D("TrueKEMu", "True Muon KE; True Muon KE (GeV); N_{events}", 200, 0, 20);

  TH1D *TrueMuKE_w[nweights];
  for (int i = 0; i < nweights; ++i) {
    TString name = weights[i]->GetName();
    name.ReplaceAll("h_xProj","");
    TrueMuKE_w[i] = new TH1D(Form("TrueMuKE_w_%i", i), Form("%s;Muon true KE (GeV); N_{events}, weighted", name.Data()), 200, 0, 20);
    std::cout << TrueMuKE_w[i]->GetTitle() << std::endl;
    TrueMuKE_w[i]->SetLineColor(10000+i);
    TrueMuKE_w[i]->SetMarkerSize(0);
    TrueMuKE_w[i]->SetFillStyle(0);
    if (i > 5) TrueMuKE_w[i]->SetLineStyle(kDashed);
  }

  TH1D *hClusterEnergy = new TH1D("ClusterEn", "Total cluster energy; Cluster energy (MeV); N_{events}", 100, 0, 100);

  TrueMuKE->SetLineColor(kBlack);
  TrueMuKE->SetLineStyle(kDashed);
  TrueMuKE->SetMarkerSize(0);
  TrueMuKE->SetFillStyle(0);

  int ngood = 0;
  int nentries = truth->GetEntries();
  //nentries = 5.5E6;
  int trklen_counter = 0;
  // Sometimes offset the entry
  int true_entry = 0;
  int reco_entry = 0;
  std::cout << nentries << " events..." << std::endl;
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
    hClusterEnergy->Fill(cluster_en);
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
    // Not needed for exiting
    //if (Muon_Death[1] > 1159 || Muon_Death[1] < -3864) continue;

    // Check that the longest track stops in the detector, and starts in the detector FV
    if (AllDet) {
      if (FirstHoughHit[longtrack][0] < 11362+55*2) continue;
    } else {
      if (FirstHoughHit[longtrack][0] < 11362+55*2 || FirstHoughHit[longtrack][0] > 13600) continue;
    }

    // Look for exiting tracks now
    if (LastHoughHit[longtrack][0] < 18300) continue;

    // 10 cm inwards
    if (fabs(FirstHoughHit[longtrack][1]) > 3520-100) continue;
    //if (fabs(LastHoughHit[longtrack][1]) > 3520-200) continue;

    h_Occupancy->Fill(Occupancy[longtrack]);

    // Ask for only small amount of other energy deposits
    if (Occupancy[longtrack] < OccupancyCut) continue;

    // Get the weighted
    // First true Enu, in GeV
    double TrueEnu_d = NeutrinoP4[3];
    TrueEnu->Fill(TrueEnu_d);
    double muonke = (MuonP4[3]-106.5)/1.E3;
    TrueMuKE->Fill(muonke);
    for (int j = 0; j < nweights; ++j) {
      double weight = weights[j]->GetBinContent(weights[j]->FindBin(TrueEnu_d));
      TrueMuKE_w[j]->Fill(muonke, weight);
    }

    double truehad = NeutrinoP4[3] - MuonP4[3]/1000.;
    TrueHad->Fill(truehad);
    ngood++;
    /*
    if (NeutrinoP4[3] < 2) {
      std::cout << "event " << i << std::endl;
      std::cout << "enu: " << NeutrinoP4[3] << std::endl;
      std::cout << "emu: " << MuonP4[3]/1E3 << std::endl;
      std::cout << "First last hough hit: " << FirstHoughHit[longtrack][0] << "-" << LastHoughHit[longtrack][0] << std::endl;
      std::cout << "tracklength: " << TrackLength[longtrack] << std::endl;
      std::cout << FirstHoughPlane[longtrack] << "-" << LastHoughPlane[longtrack] << std::endl;
      std::cout << EventNum_reco << " " << EventNum_true << std::endl;
    }
    */
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

  TString canvname = Form("MuonKE_OVERFLOW_SIDEINC_%s", filename.c_str());
  canvname += Form("_nLinesCut%i_nClusterCut%i_OccupancyCut%.2f_CCmuOnly%o_ClusterEnCut%.2f_AtLeastOneLine%o_AllDet%o", nLinesCut, nClustersCut, OccupancyCut, CCmuOnly, ClusterEnergyCut, AtLeastOneLine, AllDet);
  canvname += ".pdf";
  canv->Print(canvname+"[");

  gStyle->SetPalette(55);

  h_Occupancy->Draw("hist");
  canv->Print(canvname);

  hClusterEnergy->Draw("hist");
  canv->Print(canvname);

  TrueEnu->Draw("hist");
  canv->Print(canvname);

  TrueHad->Draw("hist");
  canv->Print(canvname);

  TrueMuKE->Draw("hist");
  for (int i = 0; i < nweights; ++i) {
    double n2llh = CalcPoisson(TrueMuKE_w[i], TrueMuKE);
    TrueMuKE_w[i]->SetTitle(Form("%s #chi^{2}=%.2f", TrueMuKE_w[i]->GetTitle(), n2llh));
    TrueMuKE_w[i]->Draw("same,hist");
  }
  canv->BuildLegend();
  TrueMuKE->Draw("same,hist");
  canv->Print(canvname);

  TH1D *ratios[nweights];
  for (int i = 0; i < nweights; ++i) {
    ratios[i] = (TH1D*)TrueMuKE_w[i]->Clone(Form("%s_r", TrueMuKE_w[i]->GetName()));
    ratios[i]->Divide(TrueMuKE);
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

  // Final single bin comparisons
  TH1D *onebin_nom = new TH1D("nom", "nom;;Number of events", 1, 0, 1);
  onebin_nom->SetBinContent(1, TrueMuKE->Integral());
  onebin_nom->SetLineColor(kBlack);
  onebin_nom->SetLineStyle(kDashed);
  onebin_nom->SetMarkerSize(0);
  onebin_nom->SetFillStyle(0);

  TH1D *onebin_w[nweights];
  for (int i = 0; i < nweights; ++i) {
    onebin_w[i] = new TH1D(Form("weight_%i", i), Form("%s;;Number of events", weights[i]->GetName()), 1, 0, 1);
    onebin_w[i]->SetLineColor(10000+i);
    onebin_w[i]->SetMarkerSize(0);
    onebin_w[i]->SetFillStyle(0);
    if (i > 5) onebin_w[i]->SetLineStyle(kDashed);
    onebin_w[i]->SetBinContent(1, TrueMuKE_w[i]->Integral());
    // Now make the llh
    double n2llh = CalcPoisson(onebin_w[i], onebin_nom);
    onebin_w[i]->SetTitle(Form("%s #chi^{2}=%.2f", onebin_w[i]->GetTitle(), n2llh));
  }

  // Draw them
  onebin_nom->Draw("hist");
  for (int i = 0; i < nweights; ++i) onebin_w[i]->Draw("same,hist");
  canv->BuildLegend();
  canv->Print(canvname);

  // Now make the single bin ratio
  TH1D *onebin_ratios[nweights];
  for (int i = 0; i < nweights; ++i) {
    onebin_ratios[i] = (TH1D*)onebin_w[i]->Clone(Form("%s_r", onebin_w[i]->GetName()));
    onebin_ratios[i]->Divide(onebin_nom);
  }
  onebin_ratios[0]->GetYaxis()->SetRangeUser(0.8, 1.2);
  onebin_ratios[0]->GetYaxis()->SetTitle("Ratio to nom.");
  TLine *line2 = new TLine(0, 1, 1, 1);
  line2->SetLineWidth(2);
  line2->SetLineColor(kRed);
  line2->SetLineStyle(kDashed);
  for (int i = 0; i < nweights; ++i) {
    if (i == 0) onebin_ratios[i]->Draw("hist");
    else onebin_ratios[i]->Draw("same,hist");
  }
  canv->BuildLegend();
  line2->Draw("same");
  canv->Print(canvname);

  // Write all historams to a file too
  TString outname = canvname;
  outname = outname.ReplaceAll(".pdf", ".root");
  TFile *output = new TFile(outname, "recreate");

  h_Occupancy->Write();
  TrueEnu->Write();
  TrueHad->Write();
  TrueMuKE->Write();
  for (int i = 0; i < nweights; ++i) {
    TrueMuKE_w[i]->Write();
    ratios[i]->Write();
  }
  hClusterEnergy->Write();

  onebin_nom->Write();
  for (int i = 0; i < nweights; ++i) {
    onebin_w[i]->Write();
    onebin_ratios[i]->Write();
  }
  output->Close();

  canv->Print(canvname+"]");
}
