bool HasField;

void GetFitLimit(double &minfit, double &maxfit, int type) {

  if (type == 1) {
    // 0.3-1.5 KE range
    minfit = 0.3;
    maxfit = 1.5;
  } else if (type == 2) {
    // 0.3-2 KE range
    minfit = 0.3;
    maxfit = 2;

  } else if (type == 3) {
    // 0.3-3 KE range
    minfit = 0.3;
    maxfit = 3;

  } else if (type == 4) {
    // 0.3-4 KE range
    minfit = 0.3;
    maxfit = 4;

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
    // 0.5-4.5 KE range
    minfit = 0.5;
    maxfit = 4.0;

  } else if (type == 9) {
    // 0.5-3.5 KE range
    minfit = 0.5;
    maxfit = 3.5;

  } else if (type == 10) {
    // 1-4 KE range
    minfit = 1;
    maxfit = 4;
  }

}

double GetKEEstimate(float MuScintLen, int type) {
  double a = 0;
  double m = 0;

  // Magnetic field, reduced TMS region
  if (HasField) {

    if (type == 1) {
      // 0.3-1.5 KE range
      a = 50.53969;
      m = 606.48666;

    } else if (type == 2) {
      // 0.3-2 KE range
      a = 49.42859;
      m = 604.82862;

    } else if (type == 3) {
      // 0.3-3 KE range
      a = 45.18268;
      m = 599.36078;

    } else if (type == 4) {
      // 0.3-4 KE range
      a = 29.68949;
      m = 582.08401;

    } else if (type == 5) {
      // 0.5-1.5 KE range
      a = 39.08261;
      m = 595.25355;

    } else if (type == 6) {
      // 0.5-2.5 KE range
      a = 41.23895;
      m = 598.10703;

    } else if (type == 7) {
      // 0.5-4 KE range
      a = 19.71407;
      m = 577.29484;

    } else if (type == 8) {
      // 0.5-4.5 KE range
      a = 19.71407;
      m = 577.29484;

    } else if (type == 9) {
      // 0.5-3.5 KE range
      a = 36.26671;
      m = 592.97914;

    } else if (type == 10) {
      // 1-4 KE range
      a = 1.37366;
      m = 568.40946;
    }

  } else {

    if (type == 1) {
      // 0.3-1.5 KE range
      a = 51.89720;
      m = 611.16292;

    } else if (type == 2) {
      // 0.3-2 KE range
      a = 50.52472;
      m = 609.13586;

    } else if (type == 3) {
      // 0.3-3 KE range
      a = 45.19770;
      m = 602.32946;

    } else if (type == 4) {
      // 0.3-4 KE range
      a = 35.70390;
      m = 591.71539;

    } else if (type == 5) {
      // 0.5-1.5 KE range
      a = 41.76659;
      m = 601.39715;

    } else if (type == 6) {
      // 0.5-2.5 KE range
      a = 41.86355;
      m = 601.74399;

    } else if (type == 7) {
      // 0.5-4 KE range
      a = 27.13122;
      m = 587.34937;

    } else if (type == 8) {
      // 0.5-4.5 KE range
      a = 27.13122;
      m = 587.23937;

    } else if (type == 9) {
      // 0.5-3.5 KE range
      a = 35.48371;
      m = 595.23184;

    } else if (type == 10) {
      // 1-4 KE range
      a = 12.49478;
      m = 580.76844;
    }

  }

    // NO magnetic field notes:
    // KE range signficantly lower: from 0.1-3.5 GeV and then nothing
    // As a result, fitting to some kinematic ranges is rubbish!
    // e.g. 8 range
    // But the smear is much much larger
    //
    // 1 range (0.3-1.5)
    // -33.16954, 433.45964
    //
    // 2 range (0.3-2)
    // -26.84369, 420.22236
    //
    // 3 range (0.3-3)
    // -15.25318, 397.39066
    //
    // 4 range (0.3-4) (pretty rubbish--non linear)
    // 13.67681, 343.93573
    //
    // 5 range (0.5-1.5)
    // 14.01698, 380.30037
    //
    // 6 range (0.5-2.5)
    // 15.5296, 378.43798
    //
    // 7 range (0.5-4)
    // 70.43364, 310.90149 (pretty rubbish though)
    //
    // 8 range (0.5-4.5) RUBBISH
    // 109.67304, 264.51349
    //
    // 9 range (0.5-3.5)
    // 28.6099, 361.74179
    //
    // 10 range (1-4) rubbish!
    // 188.03465, 259.24341

  return (a + MuScintLen)/m;
}

void makeplotsenergy_tmske(std::string filename, int scinttype) {
  HasField = true;
  gStyle->SetPalette(55);
  gErrorIgnoreLevel = kError;

  TFile *file = new TFile(filename.c_str());
  TTree *tree = (TTree*)file->Get("tree");

  if (filename.find("0_field") != std::string::npos) HasField = false;
  else HasField = true;

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

  std::vector<float> *zpt;
  tree->SetBranchStatus("zpt", true);
  tree->SetBranchAddress("zpt", &zpt);

  tree->SetBranchStatus("muonReco", true);
  tree->SetBranchAddress("muonReco", &muonReco);

  TVector3 beam_angle(0,0,1);
  beam_angle.RotateX(0.101);

  double MaxKE = 4.0;
  double MaxScintLen = 2200;
  double MaxScintEn = 250;

  TH2D *ScintLenvsKE = new TH2D("scintlen_ke", "scintlen_ke; Muon KE (GeV); Scintillation Length (g/cm^{2})", 50, 0, MaxKE, 50, 0, MaxScintLen);

  TH2D *KEestvsKE = new TH2D("keest_ke", "keest_ke; Muon KE (GeV); Muon KE est from length", 50, 0, MaxKE, 50, 0, MaxKE);
  TH2D *KEestvsScintEn = new TH2D("keest_scint", "keest_scint; True KE - Muon KE est from length (GeV); Scint Energy", 50, -1., 1.0, 50, 0, MaxScintEn);
  TH2D *KEestvsScintEnAr[20];
  for (int i = 0; i < 20; ++i) {
    KEestvsScintEnAr[i] = new TH2D(Form("keest_scint_%i",i), Form("keest_scint_%i-%i; True KE - Muon KE est from length (GeV); Scint Energy", i*200,(i+1)*200), 50, -1., 1.0, 50, 0, MaxScintEn);
  }

  TH2D *KEvsCrudeKE =  new TH2D("ke_crudeke", "ke_crudeke; True KE (GeV); Muon KE crude (GeV)", 50, 0, MaxKE, 50, 0, MaxKE);

  TH2D *ScintEnvsNplanes = new TH2D("scinten_planes", "scinten_planes; Number of planes; Scintillation Energy", 100, 0, 100, 50, 0, MaxScintEn);

  TH2D *ScintEnvsThick = new TH2D("scinten_thickplanes", "scinten_thickplanes; Number of planes; Scintillation Energy", 100, 0, 100, 50, 0, MaxScintEn);
  TH2D *ScintEnvsThin = new TH2D("scinten_thinplanes", "scinten_thinplanes; Number of planes; Scintillation Energy", 100, 0, 100, 50, 0, MaxScintEn);

  TH2D *KEExitvsKEEnt = new TH2D("keexit_keentr", "keexit_keentr; Muon KE (GeV) exit LAr; Muon KE (GeV) enter TMS", 50, 0, MaxKE, 50, 0, MaxKE);

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

    //std::cout << "Event " << i << "/" << nEntries << " (" << i*100./nEntries << "%)" << std::endl;
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

    // Fiducial cut in TMS in z
    if (muonBirth[2] > 740. || muonDeath[2] > 1375. || 
        muonDeath[1] > 25.+50. || muonDeath[1] < -260.+50. ||
        abs(muonDeath[0]) > 300. ||
        muonExitKE <= 0 ||
        abs(muonDeath[0]) > 165. || abs(muonDeath[0]) < 10. ) continue;

    //abs(muonDeath[0]) < 190. ) continue;
    //if (muonReco == 0) std::cout << "AAAAAH" << std::endl;

    // subtract the extra track length added in
    muScintLen -= 24.35;
    double dedx_total = muScintEnergy/muScintLen;

    // Calculate the nplanes from birth and death point in z
    double distz = (muonDeath[2]-muonBirth[2]);
    // We know the planes change at z=949:
    // Between x and 949 we have thin layer: gap is 5.5 (1cm scint, 1.5cm steel, 3cm air)
    // between 949.3 and 949? There is a transition module with gap 9.5cm 
    // After 949.3 we have 8cm gap (1cm scint, 4cm steel, 3 cm air)
    // Check if the muon died in thick or thin regoin
    bool DiedInThick = false;
    if (muonDeath[2] > 949) DiedInThick = true;
    double distthin = 0;
    double distthick = 0;
    double nplanes = 0;
    double nplanes_thick = 0;
    double nplanes_thin = 0;
    if (DiedInThick) {
      distthin = 949-muonBirth[2];
      distthick = muonDeath[2]-949;
      // Add the transition plane
      nplanes = distthin/5.5. + distthick/8. + 2;
      nplanes_thick = distthick/8.;
      nplanes_thin = distthin/5.5;
      //std::cout << " nplanes thin: " << distthin/5.5 << std::endl;
      //std::cout << " nplanes thick: " << distthick/8.0 << std::endl;
    } else {
      distthin = muonDeath[2]-muonBirth[2];
      nplanes = distthin/5.5;
      nplanes_thin = distthin/5.5;
      //std::cout << "nplanes thin: " << distthin/5.5 << std::endl;
    }

    int nplanes_int = nplanes;
    ScintEnvsNplanes->Fill(nplanes, muScintEnergy);
    if (nplanes_thick > 0) ScintEnvsThick->Fill(nplanes_thick, muScintEnergy);
    ScintEnvsThin->Fill(nplanes_thin, muScintEnergy);

    KEExitvsKEEnt->Fill(muonExitKE/1E3, TMSKE/1E3);

    /*
       std::cout << "***" << std::endl;
       std::cout << "Event: " << ievt << std::endl;
       std::cout << "Birth: " << muonBirth[2] << " Death: " << muonDeath[2] << std::endl;
       std::cout << "Distance in thin: " << distthin << std::endl;
       std::cout << "Distance in thick: " << distthick << std::endl;
       std::cout << "Total distance: " << distz << std::endl;
       std::cout << "Died in thick? " << DiedInThick << std::endl;
       std::cout << "nplanes: " << nplanes << " int: " << nplanes_int << std::endl;
       */



    double crudeKE = (muScintLen/1.9)+(muScintEnergy-nplanes/2)*5; //*5; //(7.85/1.05);
    //double crudeKE = (muScintLen/1.9);
    /*
    int hits = (*zpt).size();
    double lastpoint = 0;
    for (int j = 0; j < hits; ++j) {
      //std::cout << (*zpt)[j] << std::endl;
      if (j == 0) lastpoint = (*zpt)[j];
      //std::cout << "dist z: " << (*zpt)[j]-lastpoint << std::endl;
      lastpoint = (*zpt)[j];
    }
    */
    //std::cout << muonExitKE/1E3 << " " << crudeKE << std::endl;
    muonExitKE=TMSKE;

    KEvsCrudeKE->Fill(muonExitKE/1E3, crudeKE/1E3);

    ScintLenvsKE->Fill(muonExitKE/1E3, muScintLen);
    KEestvsKE->Fill(muonExitKE/1E3, GetKEEstimate(muScintLen, type));
    KEestvsScintEn->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);

    if (muonExitKE < 200) {
      KEestvsScintEnAr[0]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 400) {
      KEestvsScintEnAr[1]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 600) {
      KEestvsScintEnAr[2]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 800) {
      KEestvsScintEnAr[3]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 1000) {
      KEestvsScintEnAr[4]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 1200) {
      KEestvsScintEnAr[5]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 1400) {
      KEestvsScintEnAr[6]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 1600) {
      KEestvsScintEnAr[7]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 1800) {
      KEestvsScintEnAr[8]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 2000) {
      KEestvsScintEnAr[9]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 2200) {
      KEestvsScintEnAr[10]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 2400) {
      KEestvsScintEnAr[11]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 2600) {
      KEestvsScintEnAr[12]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 2800) {
      KEestvsScintEnAr[13]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 3000) {
      KEestvsScintEnAr[14]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 3200) {
      KEestvsScintEnAr[15]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 3400) {
      KEestvsScintEnAr[16]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 3600) {
      KEestvsScintEnAr[17]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else if (muonExitKE < 3800) {
      KEestvsScintEnAr[18]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    } else {
      KEestvsScintEnAr[19]->Fill(muonExitKE/1E3-GetKEEstimate(muScintLen, type), muScintEnergy);
    }

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
  canvname.ReplaceAll(".root", "muon_TMS_KEest");
  canvname += Form("_type%i", type);
  canv->Print(canvname+".pdf[");
  ScintLenvsKE->Draw("colz");
  canv->Print(canvname+".pdf");

  KEvsCrudeKE->Draw("colz");
  canv->Print(canvname+".pdf");

  KEExitvsKEEnt->Draw("colz");
  canv->Print(canvname+".pdf");

  // Draw the KE distribtuoin
  TH1D *kedist = ScintLenvsKE->ProjectionX();
  kedist->Draw();
  canv->Print(canvname+".pdf");

  // Make the projected slices too
  TH1D *GaussEst = new TH1D("GaussEst", "GaussEst;Muon KE; Gauss Profile Using Scint. Length", ScintLenvsKE->GetXaxis()->GetNbins(), ScintLenvsKE->GetXaxis()->GetBinLowEdge(1), ScintLenvsKE->GetXaxis()->GetBinLowEdge(ScintLenvsKE->GetXaxis()->GetNbins()+1)); 
  TH1D *ArithEst = new TH1D("ArithEst", "ArithEst;Muon KE; Arithmetic Profile Using Scint. Length", ScintLenvsKE->GetXaxis()->GetNbins(), ScintLenvsKE->GetXaxis()->GetBinLowEdge(1), ScintLenvsKE->GetXaxis()->GetBinLowEdge(ScintLenvsKE->GetXaxis()->GetNbins()+1)); 
  //TH1D *GaussEst = new TH1D("GaussEst", "GaussEst;Gauss Profile Using Scint. Length;Muon KE", ScintLenvsKE->GetYaxis()->GetNbins(), ScintLenvsKE->GetYaxis()->GetBinLowEdge(1), ScintLenvsKE->GetYaxis()->GetBinLowEdge(ScintLenvsKE->GetYaxis()->GetNbins()+1)); 
  //TH1D *ArithEst = new TH1D("ArithEst", "ArithEst;Arithmetic Profile Using Scint. Length;Muon KE", ScintLenvsKE->GetYaxis()->GetNbins(), ScintLenvsKE->GetYaxis()->GetBinLowEdge(1), ScintLenvsKE->GetYaxis()->GetBinLowEdge(ScintLenvsKE->GetYaxis()->GetNbins()+1)); 

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

  double minfit, maxfit;
  GetFitLimit(minfit, maxfit, type);

  // Now get the parameters
  TF1 *fitting = new TF1("fitting", "[0]+[1]*x", minfit, maxfit);
  fitting->SetLineColor(GaussEst->GetMarkerColor());
  //fitting->SetLineStyle(kDashed);
  TF1 *fitting2 = new TF1("fitting2", "[0]+[1]*x", minfit, maxfit);
  fitting2->SetLineColor(ArithEst->GetMarkerColor());
  //fitting2->SetLineStyle(kDashed);
  TF1 *fitting3 = new TF1("fitting2", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", minfit, maxfit);
  fitting3->SetLineColor(kBlack);
  //fitting3->SetLineStyle(kDashed);

  GaussEst->Fit(fitting, "QS", "", minfit, maxfit);
  ArithEst->Fit(fitting2, "QS", "", minfit, maxfit);
  ArithEst->Fit(fitting3, "QS", "", minfit, maxfit);

  fitting->Draw("same");
  fitting2->Draw("same");
  fitting3->Draw("same");
  TLegend *leg = new TLegend(0.1, 0.5, 0.6, 0.9);
  leg->AddEntry(fitting, Form("Gauss Lin Est, a=%2.5f, m=%2.5f", fitting->GetParameter(0), fitting->GetParameter(1)), "l");
  leg->AddEntry(fitting2, Form("Arith Lin Est, a=%2.5f, m=%2.5f", fitting2->GetParameter(0), fitting2->GetParameter(1)), "l");
  leg->AddEntry(fitting3, Form("#splitline{Arith 3rd Est, a=%2.5f, b=%2.5f}{c=%2.5f, d=%2.5f}", fitting3->GetParameter(0), fitting3->GetParameter(1), fitting3->GetParameter(2), fitting3->GetParameter(3)), "l");
  leg->Draw("same");
  canv->Print(canvname+".pdf");

  double intercept = fitting->GetParameter(0);
  double slope = fitting->GetParameter(1);

  KEestvsKE->Draw("colz");
  canv->Print(canvname+".pdf");

  // Make the estimate plot
  TH1D *ke_est = new TH1D("ke_est", "ke_est;Muon KE est (GeV);Resolution (RMS/Mean) in %", KEestvsKE->GetXaxis()->GetNbins(), KEestvsKE->GetXaxis()->GetBinLowEdge(1), KEestvsKE->GetXaxis()->GetBinLowEdge(KEestvsKE->GetXaxis()->GetNbins()+1)); 

  TH1D *ke_bias = new TH1D("ke_bias", "ke_bias;Muon KE est (GeV);Bias [Estimate-True)/True] in %", KEestvsKE->GetXaxis()->GetNbins(), KEestvsKE->GetXaxis()->GetBinLowEdge(1), KEestvsKE->GetXaxis()->GetBinLowEdge(KEestvsKE->GetXaxis()->GetNbins()+1)); 
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

  canv->Print(canvname+".pdf]");
}
