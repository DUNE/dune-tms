void GetFitLimit(double &minfit, double &maxfit, int type) {

  if (type == 1) {
    // 0.3-1.5 KE range
    minfit = 0.3;
    maxfit = 1.5;
  } else if (type == 2) {
    // 0.3-2 KE range
    minfit = 0.2;
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
    // 0.0-4.0 KE range
    minfit = 0.0;
    maxfit=4.0;

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

void makeplotsenergy_tmske_unified(std::string filename, int scinttype) {
  gStyle->SetPalette(55);
  gErrorIgnoreLevel = kError;

  TFile *file = new TFile(filename.c_str());
  //TTree *tree = (TTree*)file->Get("tree");
  TChain *tree = new TChain("tree");
  TString fileset = filename.c_str();
  fileset = fileset(0,fileset.Last('/')+1);
  tree->Add(fileset+"*FlatTree.root");

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
  std::vector<float> *xpt = NULL;
  std::vector<float> *ypt = NULL;
  std::vector<float> *zpt = NULL;

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

  tree->SetBranchStatus("xpt", true);
  tree->SetBranchAddress("xpt", &xpt);
  tree->SetBranchStatus("ypt", true);
  tree->SetBranchAddress("ypt", &ypt);
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

  TH2D *ScintEnvsNplanes = new TH2D("scinten_planes", "scinten_planes; Number of planes; Scintillation Energy", 100, 0, 100, 50, 0, MaxScintEn);

  TH2D *ScintEnvsThick = new TH2D("scinten_thickplanes", "scinten_thickplanes; Number of planes; Scintillation Energy", 100, 0, 100, 50, 0, MaxScintEn);
  TH2D *ScintEnvsThin = new TH2D("scinten_thinplanes", "scinten_thinplanes; Number of planes; Scintillation Energy", 100, 0, 100, 50, 0, MaxScintEn);

  TH2D *KEExitvsKEEnt = new TH2D("keexit_keentr", "keexit_keentr; Muon KE (GeV) exit LAr; Muon KE (GeV) enter TMS", 50, 0, MaxKE, 50, 0, MaxKE);

  TH2D *KEvsCrudeKE =  new TH2D("ke_crudeke", "ke_crudeke; True KE (GeV); Muon KE crude (GeV)", 50, 0, MaxKE, 50, 0, MaxKE);


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

    // Fiducial cut in TMS in z
    if (muonBirth[2] > 735. || muonDeath[2] > 1365. || 
        muonDeath[1] > 87-50. || muonDeath[1] < -234+50. ||
        abs(muonDeath[0]) > 160. || abs(muonDeath[0]) < 10 ||
        muonExitKE <= 0 ) continue;
    /*
    if (muonBirth[2] > 733. || muonDeath[2] > 1375. || 
        muonDeath[1] > 25.+50. || muonDeath[1] < -260.+50. ||
        abs(muonDeath[0]) > 300. ||
        muonExitKE <= 0 ||
        abs(muonDeath[0]) > 165. || abs(muonDeath[0]) < 10. ) continue;
        */

    int nhits = (*xpt).size();
    bool bad = false;
    for (int j = 0; j < nhits; ++j) {
      double x = (*xpt)[j];
      if (fabs(x) > 160) {
        bad=true;
        break;
      }
    }
    if (bad) continue;

    //abs(muonDeath[0]) < 190. ) continue;
    //if (muonReco == 0) std::cout << "AAAAAH" << std::endl;

    // Recalculate the scintillation length
    //muScintLen -= 24.35;
    // Calculate costheta by taking LAr vertex
    double xdist = muonExitPt[0]-vtx[0];
    double ydist = muonExitPt[1]-vtx[1];
    double zdist = muonExitPt[2]-vtx[2];

    double dist = sqrt(xdist*xdist+ydist*ydist+zdist*zdist);
    double costheta = fabs(zdist/dist);

    // Or calculate from exit point in LAr vs birth point in TMS (seems better?)
    double xdist2 = muonBirth[0]-muonExitPt[0];
    double ydist2 = muonBirth[1]-muonExitPt[1];
    double zdist2 = muonBirth[2]-muonExitPt[2];
    double dist2 = sqrt(xdist2*xdist2+ydist2*ydist2+zdist2*zdist2);
    double costheta2 = fabs(zdist2/dist2);

    // this has already been added in, so now need to remove it
    muScintLen -= 24.35/costheta2;
    // Then also add in the first layer of steel that is missing
    // Also account for angle
    muScintLen -= 1.5*7.85/costheta2;

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
      nplanes = distthin/5.5 + distthick/8. + 2;
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

    muonExitKE=TMSKE;

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
  canvname.ReplaceAll(".root", "muonKEest_notrackfixinscript_uni_tmsKE");
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
  TH2D *KEestvsKE = new TH2D("keest_ke", "keest_ke; Muon KE (GeV); Muon KE est from length", 50, 0, MaxKE, 50, 0, MaxKE);
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

    // Fiducial cut in TMS in z
    /*
    if (muonBirth[2] > 733. || muonDeath[2] > 1375. || 
        muonDeath[1] > 25.+50. || muonDeath[1] < -260.+50. ||
        abs(muonDeath[0]) > 300. ||
        muonExitKE <= 0 ||
        abs(muonDeath[0]) > 165. || abs(muonDeath[0]) < 10. ) continue;
    */
    if (muonBirth[2] > 735. || muonDeath[2] > 1365. || 
        muonDeath[1] > 87-50. || muonDeath[1] < -234+50. ||
        abs(muonDeath[0]) > 160. || abs(muonDeath[0]) < 10. ||
        muonExitKE <= 0 ) continue;

    int nhits = (*xpt).size();
    bool bad = false;
    for (int j = 0; j < nhits; ++j) {
      double x = (*xpt)[j];
      if (fabs(x) > 160) {
        bad=true;
        break;
      }
    }
    if (bad) continue;
    // Need to shift the scintillator length again
    
    //muScintLen -= 24.35;
    // Or calculate from exit point in LAr vs birth point in TMS (seems better?)
    double xdist2 = muonBirth[0]-muonExitPt[0];
    double ydist2 = muonBirth[1]-muonExitPt[1];
    double zdist2 = muonBirth[2]-muonExitPt[2];
    double dist2 = sqrt(xdist2*xdist2+ydist2*ydist2+zdist2*zdist2);
    double costheta2 = fabs(zdist2/dist2);

    muScintLen -= 24.35/costheta2;
    // Then also add in the first layer of steel that is missing
    // Also account for angle
    muScintLen -= 1.5*7.85/costheta2;

    muonExitKE=TMSKE;

    double KEestimator = (muScintLen-intercept)/slope;

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

  canv->Print(canvname+".pdf]");

  TFile *output = new TFile(canvname+".root", "recreate");
  output->cd();
  ScintLenvsKE->Write();
  KEvsCrudeKE->Write();
  KEExitvsKEEnt->Write();
  kedist->Write();
  GaussEst->Write();
  ArithEst->Write();
  fitting->Write();
  fitting2->Write();
  fitting3->Write();

  KEestvsKE->Write();
  ke_est->Write();
  ke_bias->Write();
  KEestvsScintEn->Write();
  output->Close();
}
