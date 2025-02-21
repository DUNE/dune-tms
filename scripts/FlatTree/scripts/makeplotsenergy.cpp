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

double GetKEEstimate(float MuScintLen, int type) {
  double a = 0;
  double m = 0;

  // Magnetic field, reduced TMS region
  if (HasField) {

    // Numbers latest updated correcting for extra track length costheta, and extra 1.5cm steel layer and costheta
    if (type == 1) {
      // 0.3-1.5 KE range
      //a = 70.08876;
      //m = 588.16938;
      a = 55.01993;
      m = 586.75049;

    } else if (type == 2) {
      // 0.3-2 KE range
      //a = 67.72533;
      //m = 584.53105;
      a = 52.73001;
      m = 583.19945;

    } else if (type == 3) {
      // 0.3-3 KE range
      //a = 63.37710;
      //m = 578.75542;
      a = 48.63773;
      m = 577.72723;

    } else if (type == 4) {
      // 0.3-4 KE range
      //a = 53.30983;
      //m = 567.10362;
      a = 41.73991;
      m = 569.67691;

    } else if (type == 5) {
      // 0.5-1.5 KE range
      //a = 56.78289;
      //m = 575.05039;
      a = 42.57383;
      m = 574.46354;

    } else if (type == 6) {
      // 0.5-2.5 KE range
      //a = 57.89044;
      //m = 576.50358;
      a = 43.44416;
      m = 575.60326;

    } else if (type == 7) {
      // 0.5-4 KE range
      //a = 43.34080;
      //m = 561.98283;
      a = 33.31205;
      m = 565.43834;

    } else if (type == 8) {
      // 0.5-4.5 KE range
      //a = 47.21737;
      //m = 563.18167;
      a = 46.65463;
      m = 572.25792;

    } else if (type == 9) {
      // 0.5-3.5 KE range
      //a = 52.98120;
      //m = 571.36314;
      a = 39.18512;
      m = 571.12740;

    } else if (type == 10) {
      // 1-4 KE range
      //a = 26.55431;
      //m = 554.28554;
      a = 22.98453;
      m = 560.70198;
    }

  } else {

    if (type == 1) {
      // 0.3-1.5 KE range
      a = 72.14711;
      m = 592.34453;

    } else if (type == 2) {
      // 0.3-2 KE range
      a = 59.55622;
      m = 588.36491;

    } else if (type == 3) {
      // 0.3-3 KE range
      a = 64.18097;
      m = 581.26033;

    } else if (type == 4) {
      // 0.3-4 KE range
      a = 57.29929;
      m = 573.25095;

    } else if (type == 5) {
      // 0.5-1.5 KE range
      a = 60.70602;
      m = 581.25857;

    } else if (type == 6) {
      // 0.5-2.5 KE range
      a = 59.44676;
      m = 579.80508;

    } else if (type == 7) {
      // 0.5-4 KE range
      a = 48.75966;
      m = 568.98775;

    } else if (type == 8) {
      // 0.5-4.5 KE range
      a = 46.54097;
      m = 566.75401;

    } else if (type == 9) {
      // 0.5-3.5 KE range
      a = 53.77334;
      m = 573.84698;

    } else if (type == 10) {
      // 1-4 KE range
      a = 37.97820;
      m = 564.07516;
    }

  }

  return (a + MuScintLen)/m;
}

void makeplotsenergy(std::string filename, int scinttype) {
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
    if (muonBirth[2] > 733. || muonDeath[2] > 1375. || 
        muonDeath[1] > 25.+50. || muonDeath[1] < -260.+50. ||
        abs(muonDeath[0]) > 300. ||
        muonExitKE <= 0 ||
        abs(muonDeath[0]) > 165. || abs(muonDeath[0]) < 10. ) continue;
    //abs(muonDeath[0]) < 190. ) continue;
    //if (muonReco == 0) std::cout << "AAAAAH" << std::endl;

    // Recalculate the scintillation length
    muScintLen -= 24.35;
    // Calculate costheta by taking LAr vertex
    double xdist = muonExitPt[0]-vtx[0];
    double ydist = muonExitPt[1]-vtx[1];
    double zdist = muonExitPt[2]-vtx[2];

    double dist = sqrt(xdist*xdist+ydist*ydist+zdist*zdist);
    double costheta = fabs(zdist/dist);

    double xdist2 = muonBirth[0]-muonExitPt[0];
    double ydist2 = muonBirth[1]-muonExitPt[1];
    double zdist2 = muonBirth[2]-muonExitPt[2];
    double dist2 = sqrt(xdist2*xdist2+ydist2*ydist2+zdist2*zdist2);
    double costheta2 = fabs(zdist2/dist2);
    //std::cout << costheta << " " << costheta2 << std::endl;

    //std::cout << "Adding 24.35/" << costheta << "=" << 24.35/costheta << " to length" << std::endl;
    muScintLen += 24.35/costheta2;
    // Then also add in the first layer of steel that is missing
    // Also account for angle
    muScintLen += 1.5*7.85/costheta2;

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
    int hits = (*zpt).size();
    double lastpoint = 0;
    for (int j = 0; j < hits; ++j) {
      //std::cout << (*zpt)[j] << std::endl;
      if (j == 0) lastpoint = (*zpt)[j];
      //std::cout << "dist z: " << (*zpt)[j]-lastpoint << std::endl;
      lastpoint = (*zpt)[j];
    }
    //std::cout << muonExitKE/1E3 << " " << crudeKE << std::endl;

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
  canvname.ReplaceAll(".root", "muonKEest_xtratrackfix_");
  canvname += Form("_type%i", type);
  std::cout << canvname << std::endl;
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
