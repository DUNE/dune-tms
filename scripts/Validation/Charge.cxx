{
  for (int it = 0; it < reco.nTracks; it++) {
    // Only consider those that have enough visible energy in the TMS to count
    bool tms_touch = truth.RecoTrackPrimaryParticleTrueVisibleEnergy[it] >=
                     MINIMUM_VISIBLE_ENERGY;
    bool nd_physics_muon = NDPhysicsMuon(truth, reco, it);
    bool ismuon = truth.RecoTrackPrimaryParticlePDG[it] == 13;
    bool isantimuon = truth.RecoTrackPrimaryParticlePDG[it] == -13;
    bool ispion = abs(truth.RecoTrackPrimaryParticlePDG[it]) == 211;
    int true_charge = (truth.RecoTrackPrimaryParticlePDG[it] < 0) ? -1 : 1;
    int charge = reco.Charge[it];
    bool correct_charge_id = true_charge * charge > 0;
    double particle_starting_ke =
        truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][3] * GEV;
    if (ismuon || isantimuon) {
      // Correct-charge fraction in the plane of true KE and true x-z angle at
      // the TMS front face (the ND PDR muon-spectrometer figure). Numerator and
      // denominator are kept as separate 2D histograms so files can be merged
      // with hadd and divided afterwards. Filled for two selections: "ndphys"
      // (ND-physics muons) and "touch" (any muon with enough TMS energy).
      {
        const auto &p = truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it];
        double angle_xz_deg = atan2(p[0], p[2]) * 180.0 / M_PI;
        // Component 3 of this branch is the muon's TOTAL energy (sqrt(p^2 + m^2)), as are the
        // KE values the 1D histograms in this file and Reco_Eff.cxx / EnergyResolution.cxx fill
        // (so those axes sit 105.7 MeV above the true kinetic energy). These maps use the real
        // kinetic energy.
        const double kMuonMassMeV = 105.6584;
        const double ke_mev = p[3] - kMuonMassMeV;
        const char *species = ismuon ? "muon" : "antimuon";
        auto fill2d = [&](const char *selection, bool selected) {
          if (!selected)
            return;
          for (const char *kind : {"den", "num"}) {
            if (kind[0] == 'n' && !correct_charge_id)
              continue;
            std::string name = std::string("charge_id__angle_ke__") + selection +
                               "_" + species + "_" + kind;
            if (mapForGetHist.find(name) == mapForGetHist.end()) {
              mapForGetHist[name] = new TH2D(
                  name.c_str(),
                  ";True KE entering TMS (MeV);True x-z angle entering TMS (deg)",
                  50, 0, 5000, 60, -60, 60);
            }
            ((TH2D *)mapForGetHist[name])->Fill(ke_mev, angle_xz_deg);
          }
        };
        fill2d("ndphys", nd_physics_muon);
        fill2d("touch", tms_touch);
      }
      // incorrect charge id
      if (nd_physics_muon && ismuon && correct_charge_id)
        GetSpecialHist("special__charge_id_muon_numerator",
                       "charge_id__true_muon_numerator",
                       "#mu^{-} Charge ID Efficiency vs True KE_{Enter}",
                       "ke_tms_enter", "#Efficiency")
            ->Fill(particle_starting_ke);
      if (nd_physics_muon && isantimuon && correct_charge_id)
        GetSpecialHist("special__charge_id_antimuon_numerator",
                       "charge_id__true_antimuon_numerator",
                       "#mu^{+} Charge ID Efficiency vs True KE_{Enter}",
                       "ke_tms_enter", "#Efficiency")
            ->Fill(particle_starting_ke);
      if (nd_physics_muon && ismuon)
        GetSpecialHist("special__charge_id_muon_denominator",
                       "charge_id__true_muon_denominator",
                       "#mu^{-} Charge ID Efficiency vs True KE_{Enter}",
                       "ke_tms_enter", "#Efficiency")
            ->Fill(particle_starting_ke);
      if (nd_physics_muon && isantimuon)
        GetSpecialHist("special__charge_id_antimuon_denominator",
                       "charge_id__true_antimuon_denominator",
                       "#mu^{+} Charge ID Efficiency vs True KE_{Enter}",
                       "ke_tms_enter", "#Efficiency")
            ->Fill(particle_starting_ke);

      if (tms_touch && ismuon && correct_charge_id)
        GetHist("charge_id__all_true_muon_numerator",
                "#mu^{-} Charge ID Efficiency vs True KE_{Enter}",
                "ke_tms_enter", "#Efficiency")
            ->Fill(particle_starting_ke);
      if (tms_touch && isantimuon && correct_charge_id)
        GetHist("charge_id__all_true_antimuon_numerator",
                "#mu^{+} Charge ID Efficiency vs True KE_{Enter}",
                "ke_tms_enter", "#Efficiency")
            ->Fill(particle_starting_ke);
      if (tms_touch && ismuon)
        GetHist("charge_id__all_true_muon_denominator",
                "#mu^{-} Charge ID Efficiency vs True KE_{Enter}",
                "ke_tms_enter", "#Efficiency")
            ->Fill(particle_starting_ke);
      if (tms_touch && isantimuon)
        GetHist("charge_id__all_true_antimuon_denominator",
                "#mu^{+} Charge ID Efficiency vs True KE_{Enter}",
                "ke_tms_enter", "#Efficiency")
            ->Fill(particle_starting_ke);
      if (nd_physics_muon && !correct_charge_id) {
        double sz = reco.StartPos[it][2] * CM;
        double sx = reco.StartPos[it][0] * CM;
        double ez = reco.EndPos[it][2] * CM;
        double ex = reco.EndPos[it][0] * CM;
        if (isantimuon)
          GetHist("charge_id__validation__startpoint_antimuons",
                  "Startpoint of Incorrectedly Identified Antimuons", "Z", "X")
              ->Fill(sz, sx);
        if (ismuon)
          GetHist("charge_id__validation__startpoint_muon",
                  "Startpoint of Incorrectedly Identified Muons", "Z", "X")
              ->Fill(sz, sx);
        if (isantimuon)
          GetHist("charge_id__validation__endpoint_antimuons",
                  "Endpoint of Incorrectedly Identified Antimuons", "Z", "X")
              ->Fill(ez, ex);
        if (ismuon)
          GetHist("charge_id__validation__endpoint_muon",
                  "Endpoint of Incorrectedly Identified Muons", "Z", "X")
              ->Fill(ez, ex);
      }
    } // if muon || antimuon
    if (ispion && tms_touch) {
      if (correct_charge_id)
        GetHist("charge_id__any_pion_numerator",
                "#pi Charge ID Efficiency vs True KE_{Enter}",
                "pion_ke_tms_enter", "#Efficiency")
            ->Fill(particle_starting_ke);
      GetHist("charge_id__any_pion_denominator",
              "#pi Charge ID Efficiency vs True KE_{Enter}",
              "pion_ke_tms_enter", "#Efficiency")
          ->Fill(particle_starting_ke);
    }
  }
}
