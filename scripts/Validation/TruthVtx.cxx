{
  REGISTER_AXIS(GenericHadE,
                std::make_tuple("Hadronic Energy (MeV)", 51, 0, 100));
  REGISTER_AXIS(GenericVisE,
                std::make_tuple("Visible Energy (MeV)", 51, 0, 100));
  REGISTER_AXIS(ENu, std::make_tuple("E_{#nu} (GeV)", 40, 0, 20));
  REGISTER_AXIS(PNu, std::make_tuple("P_{#nu} (GeV)", 40, 0, 20));
  REGISTER_AXIS(PNuZ, std::make_tuple("P_{#nu Z} (GeV)", 40, 0, 20));
  if (on_new_spill) {
    GetHist("basic__truth_vtx__TrueVtxN", "TrueVtxN", "n0-120")
        ->Fill(truth.TrueVtxN);
    for (int ivtx = 0; ivtx < truth.TrueVtxN; ivtx++) {
      // Position info
      GetHist("basic__truth_vtx__TrueVtxX", "TrueVtxX", "X")
          ->Fill(truth.TrueVtxX[ivtx] * CM);
      GetHist("basic__truth_vtx__TrueVtxY", "TrueVtxY", "Y_full")
          ->Fill(truth.TrueVtxY[ivtx] * CM);
      GetHist("basic__truth_vtx__TrueVtxZ", "TrueVtxZ", "Z_full")
          ->Fill(truth.TrueVtxZ[ivtx] * CM);
      GetHist("basic__truth_vtx__TrueVtxT", "TrueVtxT", "T")
          ->Fill(truth.TrueVtxT[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxXY", "TrueVtx X vs Y", "X", "Y_full")
          ->Fill(truth.TrueVtxX[ivtx] * CM, truth.TrueVtxY[ivtx] * CM);
      GetHist("basic__truth_vtx__TrueVtxYZ", "TrueVtx Y vs Z", "Z_full",
              "Y_full")
          ->Fill(truth.TrueVtxZ[ivtx] * CM, truth.TrueVtxY[ivtx] * CM);

      GetHist("basic__truth_vtx__TrueVtxPz", "TrueVtxPz", "PNuZ")
          ->Fill(truth.TrueVtxPz[ivtx] * GEV);
      double P = TVector3(truth.TrueVtxPx[ivtx], truth.TrueVtxPy[ivtx],
                          truth.TrueVtxPz[ivtx])
                     .Mag() *
                 GEV;
      GetHist("basic__truth_vtx__TrueVtxP", "TrueVtxP", "PNu")->Fill(P);
      GetHist("basic__truth_vtx__TrueVtxE", "TrueVtxE", "ENu")
          ->Fill(truth.TrueVtxE[ivtx] * GEV);

      GetHist("basic__truth_vtx__TrueVtxHadronicELarShell",
              "TrueVtxHadronicELarShell", "GenericHadE")
          ->Fill(truth.TrueVtxHadronicELarShell[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxHadronicELarShell",
              "TrueVtxHadronicELAr", "GenericHadE")
          ->Fill(truth.TrueVtxHadronicELAr[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxHadronicELAr", "TrueVtxHadronicELAr",
              "GenericHadE")
          ->Fill(truth.TrueVtxHadronicELAr[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxHadronicETMS", "TrueVtxHadronicETMS",
              "GenericHadE")
          ->Fill(truth.TrueVtxHadronicETMS[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxHadronicE", "TrueVtxHadronicE",
              "GenericHadE")
          ->Fill(truth.TrueVtxHadronicE[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxVisibleETMS", "TrueVtxVisibleETMS",
              "GenericVisE")
          ->Fill(truth.TrueVtxVisibleETMS[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxVisibleELAr", "TrueVtxVisibleELAr",
              "GenericVisE")
          ->Fill(truth.TrueVtxVisibleELAr[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxVisibleE", "TrueVtxVisibleE",
              "GenericVisE")
          ->Fill(truth.TrueVtxVisibleE[ivtx]);

      if (truth.TrueVtxFiducialCut[ivtx]) {
        GetHist("basic__truth_vtx__FidCut__TrueVtxX", "TrueVtxX", "X")
            ->Fill(truth.TrueVtxX[ivtx] * CM);
        GetHist("basic__truth_vtx__FidCut__TrueVtxY", "TrueVtxY", "Y_full")
            ->Fill(truth.TrueVtxY[ivtx] * CM);
        GetHist("basic__truth_vtx__FidCut__TrueVtxZ", "TrueVtxZ", "Z_full")
            ->Fill(truth.TrueVtxZ[ivtx] * CM);
        GetHist("basic__truth_vtx__FidCut__TrueVtxXY", "TrueVtx X vs Y", "X",
                "Y_full")
            ->Fill(truth.TrueVtxX[ivtx] * CM, truth.TrueVtxY[ivtx] * CM);
      }
      if (truth.TrueVtxShellEnergyCut[ivtx]) {
        GetHist("basic__truth_vtx__ShellCut__TrueVtxX", "TrueVtxX", "X")
            ->Fill(truth.TrueVtxX[ivtx] * CM);
        GetHist("basic__truth_vtx__ShellCut__TrueVtxY", "TrueVtxY", "Y_full")
            ->Fill(truth.TrueVtxY[ivtx] * CM);
        GetHist("basic__truth_vtx__ShellCut__TrueVtxZ", "TrueVtxZ", "Z_full")
            ->Fill(truth.TrueVtxZ[ivtx] * CM);
        GetHist("basic__truth_vtx__ShellCut__TrueVtxXY", "TrueVtx X vs Y", "X",
                "Y_full")
            ->Fill(truth.TrueVtxX[ivtx] * CM, truth.TrueVtxY[ivtx] * CM);

        GetHist("basic__truth_vtx__ShellCut__TrueVtxHadronicELarShell",
                "TrueVtxHadronicELarShell", "GenericHadE")
            ->Fill(truth.TrueVtxHadronicELarShell[ivtx]);
      }
      if (truth.TrueVtxNDPhysicsCut[ivtx]) {
        GetHist("basic__truth_vtx__NDPhysicsCut__TrueVtxX", "TrueVtxX", "X")
            ->Fill(truth.TrueVtxX[ivtx] * CM);
        GetHist("basic__truth_vtx__NDPhysicsCut__TrueVtxY", "TrueVtxY",
                "Y_full")
            ->Fill(truth.TrueVtxY[ivtx] * CM);
        GetHist("basic__truth_vtx__NDPhysicsCut__TrueVtxZ", "TrueVtxZ",
                "Z_full")
            ->Fill(truth.TrueVtxZ[ivtx] * CM);
        GetHist("basic__truth_vtx__NDPhysicsCut__TrueVtxXY", "TrueVtx X vs Y",
                "X", "Y_full")
            ->Fill(truth.TrueVtxX[ivtx] * CM, truth.TrueVtxY[ivtx] * CM);
      }

      // Pdg code
      GetHist("basic__truth_vtx__Pdg", "TrueVtxPDG", "nu_pdg")
          ->Fill(NuPDGtoIndex(truth.TrueVtxPDG[ivtx]));

      // Cut info
      GetHist("basic__truth_vtx__TrueVtxFiducialCut", "TrueVtxFiducialCut",
              "falsetrue")
          ->Fill(truth.TrueVtxFiducialCut[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxShellEnergyCut",
              "TrueVtxShellEnergyCut", "falsetrue")
          ->Fill(truth.TrueVtxShellEnergyCut[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxNDPhysicsCut", "TrueVtxNDPhysicsCut",
              "falsetrue")
          ->Fill(truth.TrueVtxNDPhysicsCut[ivtx]);
    }
  }
}
