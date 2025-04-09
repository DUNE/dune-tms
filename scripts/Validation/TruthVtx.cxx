{
  REGISTER_AXIS(GenericHadE,
                std::make_tuple("Hadronic Energy (MeV)", 51, 0, 100));
  REGISTER_AXIS(GenericVisE,
                std::make_tuple("Visible Energy (MeV)", 51, 0, 100));
  REGISTER_AXIS(ENu, std::make_tuple("E_{#nu} (GeV)", 40, 0, 20));
  REGISTER_AXIS(PNu, std::make_tuple("P_{#nu} (GeV)", 40, 0, 20));
  REGISTER_AXIS(PNuZ, std::make_tuple("P_{#nu Z} (GeV)", 40, 0, 20));
  if (on_new_spill) {
    GetHist("basic__truth_vtx__TrueVtxN", "TrueVtxN", "n0-120", "#N Spills")
        ->Fill(truth.TrueVtxN);
    for (int ivtx = 0; ivtx < truth.TrueVtxN; ivtx++) {
      // Position info
      GetHist("basic__truth_vtx__TrueVtxX", "TrueVtxX", "X", "#N Vertices")
          ->Fill(truth.TrueVtxX[ivtx] * CM);
      GetHist("basic__truth_vtx__TrueVtxY", "TrueVtxY", "Y_full", "#N Vertices")
          ->Fill(truth.TrueVtxY[ivtx] * CM);
      GetHist("basic__truth_vtx__TrueVtxZ", "TrueVtxZ", "Z_full", "#N Vertices")
          ->Fill(truth.TrueVtxZ[ivtx] * CM);
      GetHist("basic__truth_vtx__TrueVtxT", "TrueVtxT", "T", "#N Vertices")
          ->Fill(truth.TrueVtxT[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxXY", "TrueVtx X vs Y", "X", "Y_full")
          ->Fill(truth.TrueVtxX[ivtx] * CM, truth.TrueVtxY[ivtx] * CM);
      GetHist("basic__truth_vtx__TrueVtxYZ", "TrueVtx Y vs Z", "Z_full",
              "Y_full")
          ->Fill(truth.TrueVtxZ[ivtx] * CM, truth.TrueVtxY[ivtx] * CM);

      GetHist("basic__truth_vtx__TrueVtxPz", "TrueVtxPz", "PNuZ", "#N Vertices")
          ->Fill(truth.TrueVtxPz[ivtx] * GEV);
      double P = TVector3(truth.TrueVtxPx[ivtx], truth.TrueVtxPy[ivtx],
                          truth.TrueVtxPz[ivtx])
                     .Mag() *
                 GEV;
      GetHist("basic__truth_vtx__TrueVtxP", "TrueVtxP", "PNu", "#N Vertices")->Fill(P);
      GetHist("basic__truth_vtx__TrueVtxE", "TrueVtxE", "ENu", "#N Vertices")
          ->Fill(truth.TrueVtxE[ivtx] * GEV);

      GetHist("basic__truth_vtx__TrueVtxHadronicELarShell",
              "TrueVtxHadronicELarShell", "GenericHadE", "#N Vertices")
          ->Fill(truth.TrueVtxHadronicELarShell[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxHadronicELarShell",
              "TrueVtxHadronicELAr", "GenericHadE", "#N Vertices")
          ->Fill(truth.TrueVtxHadronicELAr[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxHadronicELAr", "TrueVtxHadronicELAr",
              "GenericHadE", "#N Vertices")
          ->Fill(truth.TrueVtxHadronicELAr[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxHadronicETMS", "TrueVtxHadronicETMS",
              "GenericHadE", "#N Vertices")
          ->Fill(truth.TrueVtxHadronicETMS[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxHadronicE", "TrueVtxHadronicE",
              "GenericHadE", "#N Vertices")
          ->Fill(truth.TrueVtxHadronicE[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxVisibleETMS", "TrueVtxVisibleETMS",
              "GenericVisE", "#N Vertices")
          ->Fill(truth.TrueVtxVisibleETMS[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxVisibleELAr", "TrueVtxVisibleELAr",
              "GenericVisE", "#N Vertices")
          ->Fill(truth.TrueVtxVisibleELAr[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxVisibleE", "TrueVtxVisibleE",
              "GenericVisE", "#N Vertices")
          ->Fill(truth.TrueVtxVisibleE[ivtx]);

      if (truth.TrueVtxFiducialCut[ivtx]) {
        GetHist("basic__truth_vtx__FidCut__TrueVtxX", "TrueVtxX", "X", "#N Vertices")
            ->Fill(truth.TrueVtxX[ivtx] * CM);
        GetHist("basic__truth_vtx__FidCut__TrueVtxY", "TrueVtxY", "Y_full", "#N Vertices")
            ->Fill(truth.TrueVtxY[ivtx] * CM);
        GetHist("basic__truth_vtx__FidCut__TrueVtxZ", "TrueVtxZ", "Z_full", "#N Vertices")
            ->Fill(truth.TrueVtxZ[ivtx] * CM);
        GetHist("basic__truth_vtx__FidCut__TrueVtxXY", "TrueVtx X vs Y", "X",
                "Y_full")
            ->Fill(truth.TrueVtxX[ivtx] * CM, truth.TrueVtxY[ivtx] * CM);
      }
      if (truth.TrueVtxShellEnergyCut[ivtx]) {
        GetHist("basic__truth_vtx__ShellCut__TrueVtxX", "TrueVtxX", "X", "#N Vertices")
            ->Fill(truth.TrueVtxX[ivtx] * CM);
        GetHist("basic__truth_vtx__ShellCut__TrueVtxY", "TrueVtxY", "Y_full", "#N Vertices")
            ->Fill(truth.TrueVtxY[ivtx] * CM);
        GetHist("basic__truth_vtx__ShellCut__TrueVtxZ", "TrueVtxZ", "Z_full", "#N Vertices")
            ->Fill(truth.TrueVtxZ[ivtx] * CM);
        GetHist("basic__truth_vtx__ShellCut__TrueVtxXY", "TrueVtx X vs Y", "X",
                "Y_full")
            ->Fill(truth.TrueVtxX[ivtx] * CM, truth.TrueVtxY[ivtx] * CM);

        GetHist("basic__truth_vtx__ShellCut__TrueVtxHadronicELarShell",
                "TrueVtxHadronicELarShell", "GenericHadE", "#N Vertices")
            ->Fill(truth.TrueVtxHadronicELarShell[ivtx]);
      }
      if (truth.TrueVtxNDPhysicsCut[ivtx]) {
        GetHist("basic__truth_vtx__NDPhysicsCut__TrueVtxX", "TrueVtxX", "X", "#N Vertices")
            ->Fill(truth.TrueVtxX[ivtx] * CM);
        GetHist("basic__truth_vtx__NDPhysicsCut__TrueVtxY", "TrueVtxY",
                "Y_full", "#N Vertices")
            ->Fill(truth.TrueVtxY[ivtx] * CM);
        GetHist("basic__truth_vtx__NDPhysicsCut__TrueVtxZ", "TrueVtxZ",
                "Z_full", "#N Vertices")
            ->Fill(truth.TrueVtxZ[ivtx] * CM);
        GetHist("basic__truth_vtx__NDPhysicsCut__TrueVtxXY", "TrueVtx X vs Y",
                "X", "Y_full", "#N Vertices")
            ->Fill(truth.TrueVtxX[ivtx] * CM, truth.TrueVtxY[ivtx] * CM);
      }

      // Pdg code
      GetHist("basic__truth_vtx__Pdg", "TrueVtxPDG", "nu_pdg", "#N Vertices")
          ->Fill(NuPDGtoIndex(truth.TrueVtxPDG[ivtx]));

      // Cut info
      GetHist("basic__truth_vtx__TrueVtxFiducialCut", "TrueVtxFiducialCut",
              "falsetrue", "#N Vertices")
          ->Fill(truth.TrueVtxFiducialCut[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxShellEnergyCut",
              "TrueVtxShellEnergyCut", "falsetrue", "#N Vertices")
          ->Fill(truth.TrueVtxShellEnergyCut[ivtx]);
      GetHist("basic__truth_vtx__TrueVtxNDPhysicsCut", "TrueVtxNDPhysicsCut",
              "falsetrue", "#N Vertices")
          ->Fill(truth.TrueVtxNDPhysicsCut[ivtx]);
    }
  }
}
