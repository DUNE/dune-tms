void GetFitLimit(double &minfit, double &maxfit, int type) {

  if (type == 1) {
    // 0.3-1.5 KE range
    minfit = 0.1;
    maxfit = 1.0;
  } else if (type == 2) {
    // 0.3-2 KE range
    minfit = 1.0;
    maxfit = 2.8;

  } else if (type == 3) {
    // 0.3-3 KE range
    minfit = 0.5;
    maxfit = 5.5;

  } else if (type == 4) {
    // 0.3-4 KE range
    minfit = 1.2;
    maxfit = 2;

    /*
  } else if (type == 5) {
    // 0.5-1.5 KE range
    minfit=0.5;
    maxfit=1.5;

  } else if (type == 6) {
    // 0.5-2.5 KE range
    minfit = 0.5;
    maxfit = 2.5;

  } else if (type == 7) {
    // 0.5-4 KE range
    minfit=0.5;
    maxfit=4;

  } else if (type == 8) {
    // 0.0-4.0 KE range
    minfit = 0.5;
    maxfit=4.0;

  } else if (type == 9) {
    // 0.5-3.5 KE range
    minfit = 0.5;
    maxfit = 3.5;

  } else if (type == 10) {
    // 1-4 KE range
    minfit = 1;
    maxfit = 4;
    */
  }

}

void makeplotsenergy_unified_laronly(std::string filename, int scinttype) {
  gStyle->SetPalette(55);
  gErrorIgnoreLevel = kError;

  TFile *file = new TFile(filename.c_str());
  TTree *tree = (TTree*)file->Get("tree");

  int nEntries = tree->GetEntries();

  float vtx[3];
  int ievt;
  int muonReco;
  float Enu;
  float lepE;
  int lepPdg;
  float p3lep[3];
  float muonDeath[3];
  float muonBirth[3];
  float muonExitPt[3];
  float muScintLen;
  float muScintEnergy;
  float muonExitKE;
  float TMSKE;
  float muLArLen;
  float lepKE;

  tree->SetBranchStatus("*", false);
  tree->SetBranchStatus("vtx", true);
  tree->SetBranchAddress("vtx", vtx);
  tree->SetBranchStatus("muonDeath", true);
  tree->SetBranchAddress("muonDeath", muonDeath);
  tree->SetBranchStatus("muonBirth", true);
  tree->SetBranchAddress("muonBirth", muonBirth);
  tree->SetBranchStatus("muonExitPt", true);
  tree->SetBranchAddress("muonExitPt", muonExitPt);
  tree->SetBranchStatus("muonExitKE", true);
  tree->SetBranchAddress("muonExitKE", &muonExitKE);
  tree->SetBranchStatus("p3lep", true);
  tree->SetBranchAddress("p3lep", p3lep);
  tree->SetBranchStatus("ievt", true);
  tree->SetBranchAddress("ievt", &ievt);
  tree->SetBranchStatus("Ev", true);
  tree->SetBranchAddress("Ev", &Enu);
  tree->SetBranchStatus("lepE", true);
  tree->SetBranchAddress("lepE", &lepE);
  tree->SetBranchStatus("lepPdg", true);
  tree->SetBranchAddress("lepPdg", &lepPdg);
  tree->SetBranchStatus("muScintLen", true);
  tree->SetBranchAddress("muScintLen", &muScintLen);
  tree->SetBranchStatus("muScintEnergy", true);
  tree->SetBranchAddress("muScintEnergy", &muScintEnergy);
  tree->SetBranchStatus("rmmsKE", true);
  tree->SetBranchAddress("rmmsKE", &TMSKE);
  tree->SetBranchStatus("muLArLen_seg", true);
  tree->SetBranchAddress("muLArLen_seg", &muLArLen);
  tree->SetBranchStatus("lepKE", true);
  tree->SetBranchAddress("lepKE", &lepKE);

  std::vector<float> *zpt;
  tree->SetBranchStatus("zpt", true);
  tree->SetBranchAddress("zpt", &zpt);

  tree->SetBranchStatus("muonReco", true);
  tree->SetBranchAddress("muonReco", &muonReco);

  TVector3 beam_angle(0,0,1);
  beam_angle.RotateX(0.101);

  double MaxKE = 1.5;
  //double MaxScintLen = 2200+1050;
  double MaxScintLen = 1050;
  double MaxScintEn = 250;

  TH2D *ScintLenvsKE = new TH2D("scintlen_ke", "scintlen_ke; True Muon KE Vertex (GeV); Scintillation Length (g/cm^{2})", 50, 0, MaxKE, 50, 0, MaxScintLen);

  TH2D *ScintEnvsNplanes = new TH2D("scinten_planes", "scinten_planes; Number of planes; Scintillation Energy", 100, 0, 100, 50, 0, MaxScintEn);

  TH2D *ScintEnvsThick = new TH2D("scinten_thickplanes", "scinten_thickplanes; Number of planes; Scintillation Energy", 100, 0, 100, 50, 0, MaxScintEn);
  TH2D *ScintEnvsThin = new TH2D("scinten_thinplanes", "scinten_thinplanes; Number of planes; Scintillation Energy", 100, 0, 100, 50, 0, MaxScintEn);

  TH2D *KEExitvsKEEnt = new TH2D("keexit_keentr", "keexit_keentr; Muon KE vertex LAr (GeV); Muon KE (GeV) enter TMS", 50, 0, MaxKE, 50, 0, MaxKE);

  TH2D *KEvsCrudeKE =  new TH2D("ke_crudeke", "ke_crudeke; True vertex LAr KE (GeV); Muon KE crude (GeV)", 50, 0, MaxKE, 50, 0, MaxKE);

  int nPass = 0;
  int nTheta30 = 0;
  int nLAr = 0;
  int nInTMS = 0;
  int nNotCont = 0;
  int type = scinttype;

  for (int i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);

    TVector3 mu(p3lep[0], p3lep[1], p3lep[2]);
    double theta = mu.Angle(beam_angle)*180./3.1415;

    if (i % 100000 == 0) std::cout << "Event " << i << "/" << nEntries << " (" << i*100./nEntries << "%)" << std::endl;
    //if (abs(lepPdg) != 13) continue;
    if (lepPdg != 13) continue;
    if (theta < 30) nTheta30++;

    // If the event doesn't stop in the LAr
    nNotCont++;

    // Fiducial cut in LAr
    if (abs(vtx[0]) > 300. || abs(vtx[1]) > 100. || vtx[2] < 50. || vtx[2] > 350.) continue;
    nLAr++;
    // If the event has a stopping point in TMS
    if (muonReco == 2) nInTMS++;

    if (muonReco != 1) continue;

    // Add in the Argon track length
    muScintLen = muLArLen;

    muonExitKE=lepKE;
    KEExitvsKEEnt->Fill(muonExitKE/1E3, TMSKE/1E3);

    double crudeKE = (muScintLen/1.9);
    KEvsCrudeKE->Fill(muonExitKE/1E3, crudeKE/1E3);

    ScintLenvsKE->Fill(muonExitKE/1E3, muScintLen);

    nPass++;
  }

  std::cout << "nPass/total=" << nPass << "/" << nEntries << "=" << nPass*100./nEntries << "%" << std::endl;
  std::cout << "nTheta30/total=" << nTheta30 << "/" << nEntries << "=" << nTheta30*100./nEntries << "%" << std::endl;
  std::cout << "nPass/nLAr=" << nPass << "/" << nLAr << "=" << nPass*100./nLAr << "%" << std::endl;
  std::cout << "nPass/nTheta30=" << nPass << "/" << nTheta30 << "=" << nPass*100./nTheta30 << "%" << std::endl;
  std::cout << "nPass/nInTMS=" << nPass << "/" << nInTMS << "=" << nPass*100./nInTMS << "%" << std::endl;
  std::cout << "nPass/nNotCont=" << nPass << "/" << nNotCont << "=" << nPass*100./nNotCont << "%" << std::endl;

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  TString canvname = filename.c_str();
  canvname.ReplaceAll(".root", "muonKEest_xtratrackfix_uni_larlen_seg_laronly");
  canvname += Form("_type%i", type);

  std::cout << canvname << std::endl;
  canv->Print(canvname+".pdf[");

  // Draw scintillation length vs KE
  ScintLenvsKE->Draw("colz");
  canv->Print(canvname+".pdf");

  // Draw the crude estimate of KE
  KEvsCrudeKE->Draw("colz");
  canv->Print(canvname+".pdf");

  // Draw exit KE vs entrance KE
  KEExitvsKEEnt->Draw("colz");
  canv->Print(canvname+".pdf");

  // Draw the KE distribtuoin
  TH1D *kedist = ScintLenvsKE->ProjectionX();
  kedist->Draw();
  canv->Print(canvname+".pdf");

  // Make the projected slices too
  TH1D *GaussEst = new TH1D("GaussEst", "GaussEst;Muon KE; Gauss Profile Using Scint. Length", ScintLenvsKE->GetXaxis()->GetNbins(), ScintLenvsKE->GetXaxis()->GetBinLowEdge(1), ScintLenvsKE->GetXaxis()->GetBinLowEdge(ScintLenvsKE->GetXaxis()->GetNbins()+1)); 
  TH1D *ArithEst = new TH1D("ArithEst", "ArithEst;Muon KE; Arithmetic Profile Using Scint. Length", ScintLenvsKE->GetXaxis()->GetNbins(), ScintLenvsKE->GetXaxis()->GetBinLowEdge(1), ScintLenvsKE->GetXaxis()->GetBinLowEdge(ScintLenvsKE->GetXaxis()->GetNbins()+1)); 

  for (int i = 0; i < ScintLenvsKE->GetXaxis()->GetNbins(); ++i) {
    TH1D *proj = ScintLenvsKE->ProjectionY("_px", i, i);
    TFitResultPtr result = proj->Fit("gaus", "QS");
    if (result.Get() == NULL) continue;
    GaussEst->SetBinContent(i+1, result->GetParams()[1]);
    GaussEst->SetBinError(i+1, result->GetParams()[2]);
    double mean = proj->GetMean();
    double rms = proj->GetRMS();
    ArithEst->SetBinContent(i+1, mean);
    ArithEst->SetBinError(i+1, rms);
    //proj->Draw();
    //canv->Print(canvname+".pdf");
  }

  //GaussEst->GetYaxis()->SetRangeUser(0, ScintLenvsKE->GetXaxis()->GetBinLowEdge(ScintLenvsKE->GetXaxis()->GetNbins()+1));
  GaussEst->GetYaxis()->SetRangeUser(0, ScintLenvsKE->GetYaxis()->GetBinLowEdge(ScintLenvsKE->GetYaxis()->GetNbins()+1));
  GaussEst->SetMarkerColor(kBlue);
  GaussEst->Draw();
  ArithEst->SetMarkerColor(kRed);
  ArithEst->Draw("same");

  // Get the fit limits
  double minfit, maxfit;
  GetFitLimit(minfit, maxfit, type);

  // Define the fitting functions
  TF1 *fitting = new TF1("fitting", "[0]+[1]*x", minfit, maxfit);
  fitting->SetLineColor(GaussEst->GetMarkerColor());
  TF1 *fitting2 = new TF1("fitting2", "[0]+[1]*x", minfit, maxfit);
  fitting2->SetLineColor(ArithEst->GetMarkerColor());
  TF1 *fitting3 = new TF1("fitting2", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", minfit, maxfit);
  fitting3->SetLineColor(kBlack);

  GaussEst->Fit(fitting, "QS", "", minfit, maxfit);
  ArithEst->Fit(fitting2, "QS", "", minfit, maxfit);
  ArithEst->Fit(fitting3, "QS", "", minfit, maxfit);

  fitting->Draw("same");
  fitting2->Draw("same");
  fitting3->Draw("same");
  TLegend *leg = new TLegend(0.1, 0.5, 0.6, 0.9);
  // Now get the parameters
  leg->AddEntry(fitting, Form("Gauss Lin Est, a=%2.5f, m=%2.5f", fitting->GetParameter(0), fitting->GetParameter(1)), "l");
  leg->AddEntry(fitting2, Form("Arith Lin Est, a=%2.5f, m=%2.5f", fitting2->GetParameter(0), fitting2->GetParameter(1)), "l");
  leg->AddEntry(fitting3, Form("#splitline{Arith 3rd Est, a=%2.5f, b=%2.5f}{c=%2.5f, d=%2.5f}", fitting3->GetParameter(0), fitting3->GetParameter(1), fitting3->GetParameter(2), fitting3->GetParameter(3)), "l");
  leg->Draw("same");
  canv->Print(canvname+".pdf");

  // The fit parameters from the arithmetic projection
  double intercept = fitting2->GetParameter(0);
  double slope = fitting2->GetParameter(1);

  ScintEnvsNplanes->Draw("colz");
  canv->Print(canvname+".pdf");

  ScintEnvsNplanes->ProjectionX()->Draw();
  canv->Print(canvname+".pdf");

  ScintEnvsThick->Draw("colz");
  canv->Print(canvname+".pdf");

  ScintEnvsThick->ProjectionX()->Draw();
  canv->Print(canvname+".pdf");

  ScintEnvsThin->Draw("colz");
  canv->Print(canvname+".pdf");

  ScintEnvsThin->ProjectionX()->Draw();
  canv->Print(canvname+".pdf");

  /*
     for (int i = 0; i < 20; ++i) {
     KEestvsScintEnAr[i]->Draw("colz");
     KEestvsScintEnAr[i]->SetMaximum(KEestvsScintEn->GetMaximum());
     canv->Print(canvname+".gif+50");
     }
   */

  // Make the estimator plots
  TH2D *KEestvsKE = new TH2D("keest_ke", "keest_ke; True Muon KE Vertex LAr(GeV); Muon KE est from length", 50, 0, MaxKE, 50, 0, MaxKE);
  TH2D *KEestvsScintEn = new TH2D("keest_scint", "keest_scint; True KE - Muon KE est from length (GeV); Scint Energy", 50, -1., 1.0, 50, 0, MaxScintEn);
  TH2D *KEestvsScintEnAr[20];
  for (int i = 0; i < 20; ++i) {
    KEestvsScintEnAr[i] = new TH2D(Form("keest_scint_%i",i), Form("keest_scint_%i-%i; True KE - Muon KE est from length (GeV); Scint Energy", i*200,(i+1)*200), 50, -1., 1.0, 50, 0, MaxScintEn);
  }


  // Start the second loop where we used the coefficients to extract KE
  std::cout << "Running on second loop..." << std::endl;
  for (int i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    if (lepPdg != 13) continue;
    if (i % 100000 == 0) std::cout << "Event " << i << "/" << nEntries << " (" << i*100./nEntries << "%)" << std::endl;

    // Fiducial cut in LAr
    if (abs(vtx[0]) > 300. || abs(vtx[1]) > 100. || vtx[2] < 50. || vtx[2] > 350.) continue;
    nLAr++;
    // If the event has a stopping point in TMS
    if (muonReco == 2) nInTMS++;

    if (muonReco != 1) continue;

    muScintLen = muLArLen;

    double KEestimator = (muScintLen-intercept)/slope;
    muonExitKE = lepKE;

    KEestvsKE->Fill(muonExitKE/1E3, KEestimator);
    KEestvsScintEn->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);

    if (muonExitKE < 200) {
      KEestvsScintEnAr[0]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 400) {
      KEestvsScintEnAr[1]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 600) {
      KEestvsScintEnAr[2]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 800) {
      KEestvsScintEnAr[3]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 1000) {
      KEestvsScintEnAr[4]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 1200) {
      KEestvsScintEnAr[5]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 1400) {
      KEestvsScintEnAr[6]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 1600) {
      KEestvsScintEnAr[7]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 1800) {
      KEestvsScintEnAr[8]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 2000) {
      KEestvsScintEnAr[9]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 2200) {
      KEestvsScintEnAr[10]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 2400) {
      KEestvsScintEnAr[11]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 2600) {
      KEestvsScintEnAr[12]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 2800) {
      KEestvsScintEnAr[13]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 3000) {
      KEestvsScintEnAr[14]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 3200) {
      KEestvsScintEnAr[15]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 3400) {
      KEestvsScintEnAr[16]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 3600) {
      KEestvsScintEnAr[17]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else if (muonExitKE < 3800) {
      KEestvsScintEnAr[18]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    } else {
      KEestvsScintEnAr[19]->Fill(muonExitKE/1E3-KEestimator, muScintEnergy);
    }
  } // end event loop

  KEestvsKE->Draw("colz");
  canv->Print(canvname+".pdf");

  // Make the estimate plot
  TH1D *ke_est = new TH1D("ke_est", "ke_est;Muon KE LAr vertex est (GeV);Resolution (RMS/Mean) in %", KEestvsKE->GetXaxis()->GetNbins(), KEestvsKE->GetXaxis()->GetBinLowEdge(1), KEestvsKE->GetXaxis()->GetBinLowEdge(KEestvsKE->GetXaxis()->GetNbins()+1)); 

  TH1D *ke_bias = new TH1D("ke_bias", "ke_bias;Muon KE LAr vertex est (GeV);Bias [Estimate-True)/True] in %", KEestvsKE->GetXaxis()->GetNbins(), KEestvsKE->GetXaxis()->GetBinLowEdge(1), KEestvsKE->GetXaxis()->GetBinLowEdge(KEestvsKE->GetXaxis()->GetNbins()+1)); 
  // Get the estimates
  for (int i = 0; i < KEestvsKE->GetXaxis()->GetNbins(); ++i) {
    TH1D *proj2 = KEestvsKE->ProjectionY("_py", i, i);
    double mean2 = proj2->GetMean();
    double rms2 = proj2->GetRMS();
    if (mean2 != 0) ke_est->SetBinContent(i+1, 100*rms2/mean2);

    double trueval = ke_bias->GetBinCenter(i+1);
    if (trueval != 0) ke_bias->SetBinContent(i+1, 100*(mean2-trueval)/trueval);
  }
  ke_est->Draw();
  ke_est->GetYaxis()->SetRangeUser(0, 15);
  canv->Print(canvname+".pdf");

  ke_bias->Draw();
  ke_bias->GetYaxis()->SetRangeUser(-20, 20);
  canv->Print(canvname+".pdf");

  KEestvsScintEn->Draw("colz");
  canv->Print(canvname+".pdf");

  canv->Print(canvname+".pdf]");
}
