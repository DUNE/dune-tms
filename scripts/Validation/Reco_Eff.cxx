// Add scope to avoid cross talk with other scripts
{
  REGISTER_AXIS(ke_tms, std::make_tuple("True Muon KE (GeV)", 28, 0.0, 7.0));
  REGISTER_AXIS(
      pion_ke_tms_enter,
      std::make_tuple("True Pion KE Entering TMS (GeV)", 20, 0.0, 5.0));
  REGISTER_AXIS(muon_endpoint_z, std::make_tuple("True Muon Endpoint Z (cm)",
                                                 20, 1136.0, 1831.0));
  REGISTER_AXIS(muon_endpoint_y,
                std::make_tuple("True Muon Endpoint Y (cm)", 20, -216.0, 16.0));
  REGISTER_AXIS(muon_endpoint_x, std::make_tuple("True Muon Endpoint X (cm)",
                                                 20, -330.0, 330.0));
  REGISTER_AXIS(
      muon_startpoint_z,
      std::make_tuple("True Muon Startpoint Z (cm)", 20, 1036.0, 1831.0));
  REGISTER_AXIS(
      muon_startpoint_y,
      std::make_tuple("True Muon Startpoint Y (cm)", 20, -216.0, 16.0));
  REGISTER_AXIS(
      muon_startpoint_x,
      std::make_tuple("True Muon Startpoint X (cm)", 20, -330.0, 330.0));
  if (on_new_spill) {
    for (int ip = 0; ip < truth.nTrueParticles; ip++) {
      // Only consider those that have enough visible energy in the TMS to count
      bool tms_touch = truth.TrueVisibleEnergy[ip] >= MINIMUM_VISIBLE_ENERGY;
      if (tms_touch) {
        bool ismuon = abs(truth.PDG[ip]) == 13;
        bool ispion = abs(truth.PDG[ip]) == 211;
        bool lar_start = truth.LArFiducialStart[ip];
        bool tms_end = truth.TMSFiducialEnd[ip];
        double particle_starting_ke = truth.MomentumTMSStart[ip][3] * GEV;
        double particle_ke = truth.BirthMomentum[ip][3] * GEV;
        if (ismuon && tms_touch) {
          GetHist(
              "reco_eff__no_lar_tms_cuts__all_muon_ke_tms_enter_denominator",
              "Reconstruction Efficiency: Denominator", "ke_tms_enter")
              ->Fill(particle_starting_ke);
          GetHist("reco_eff__all_muon_ke_tms_enter_denominator",
                  "Reco Efficiency vs True KE, All Muons: Denominator",
                  "ke_tms_enter")
              ->Fill(particle_starting_ke);
          GetHist(
              "reco_eff__all_muon_ke_tms_enter_including_doubles_denominator",
              "Reco Efficiency w/2x vs True KE, All Muons: Denominator",
              "ke_tms_enter")
              ->Fill(particle_starting_ke);
          GetHist("reco_eff__all_muon_ke_denominator",
                  "Reco Efficiency vs True KE, All Muons: Denominator",
                  "ke_tms")
              ->Fill(particle_ke);

          GetHist("reco_eff__endpoint__muon_endpoint_x_denominator",
                  "Reconstruction Efficiency: Denominator", "muon_endpoint_x")
              ->Fill(truth.PositionTMSEnd[ip][0] * CM);
          GetHist("reco_eff__endpoint__muon_endpoint_y_denominator",
                  "Reconstruction Efficiency: Denominator", "muon_endpoint_y")
              ->Fill(truth.PositionTMSEnd[ip][1] * CM);
          GetHist("reco_eff__endpoint__muon_endpoint_z_denominator",
                  "Reconstruction Efficiency: Denominator", "muon_endpoint_z")
              ->Fill(truth.PositionTMSEnd[ip][2] * CM);

          GetHist("reco_eff__startpoint__muon_startpoint_x_denominator",
                  "Reconstruction Efficiency: Denominator", "muon_startpoint_x")
              ->Fill(truth.PositionTMSStart[ip][0] * CM);
          GetHist("reco_eff__startpoint__muon_startpoint_y_denominator",
                  "Reconstruction Efficiency: Denominator", "muon_startpoint_y")
              ->Fill(truth.PositionTMSStart[ip][1] * CM);
          GetHist("reco_eff__startpoint__muon_startpoint_z_denominator",
                  "Reconstruction Efficiency: Denominator", "muon_startpoint_z")
              ->Fill(truth.PositionTMSStart[ip][2] * CM);
        }
        if (ispion)
          GetHist(
              "reco_eff__no_lar_tms_cuts__all_pion_ke_tms_enter_denominator",
              "Reconstruction Efficiency: Denominator", "pion_ke_tms_enter")
              ->Fill(particle_starting_ke);
        if (ismuon && lar_start && tms_end) {
          GetSpecialHist("special__reco_eff_muon_ke_tms_enter_denominator", 
                  "reco_eff__muon_ke_tms_enter_denominator",
                  "Reconstruction Efficiency: Denominator", "ke_tms_enter")
              ->Fill(particle_starting_ke);
          GetHist("reco_eff__good_reco_muon_ke_tms_enter_denominator",
                  "Reconstruction Efficiency, with Good Reco: Denominator", "ke_tms_enter")
              ->Fill(particle_starting_ke);
          GetHist("reco_eff__muon_ke_denominator",
                  "Reco Efficiency vs True KE: Denominator", "ke_tms")
              ->Fill(particle_ke);
          double particle_starting_angle =
              std::atan2(truth.MomentumTMSStart[ip][0],
                         truth.MomentumTMSStart[ip][2]) *
              DEG;
          GetHist("reco_eff__angle_tms_enter_denominator",
                  "Reconstruction Efficiency: Denominator", "angle_tms_enter")
              ->Fill(particle_starting_angle);
        }
      } // end if (tms_touch)
    }
  }

  // Fill numerators of reco_eff plots here
  for (int it = 0; it < reco.nTracks; it++) {
    // Only consider those that have enough visible energy in the TMS to count
    bool tms_touch = truth.RecoTrackPrimaryParticleTrueVisibleEnergy[it] >=
                     MINIMUM_VISIBLE_ENERGY;
    if (tms_touch) {
      bool ismuon = abs(truth.RecoTrackPrimaryParticlePDG[it]) == 13;
      bool ispion = abs(truth.RecoTrackPrimaryParticlePDG[it]) == 211;
      bool not_double_reco = true;
      bool lar_start = truth.RecoTrackPrimaryParticleLArFiducialStart[it];
      bool tms_end = truth.RecoTrackPrimaryParticleTMSFiducialEnd[it];
      double particle_starting_ke =
          truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][3] * GEV;
      double particle_ke = -1;
      // truth.RecoTrackPrimaryParticleTrueMomentum[it][3] * GEV;
      int ip = truth.RecoTrackPrimaryParticleIndex[it];
      if (ip >= truth.nTrueParticles || ip < 0) {
        std::cout << "Warning: RecoTrackPrimaryParticleIndex is broken: Found "
                     "no true particle in reco eff with ip="
                  << ip << std::endl;
        ip = -1;
      }
      if (ip >= 0) {
        particle_ke = truth.BirthMomentum[ip][3] * GEV;
        particle_indices_reconstructed[ip]++;
        if (particle_indices_reconstructed[ip] > 1) {
          // Found this particle already at least once
          not_double_reco = false;
          /*if (ismuon)
            std::cout
                << "Warning: Found muon reconstructed multiple times in spill: "
                << particle_starting_ke << " GeV" << std::endl;
          else
            std::cout << "Warning: Found non-muon reconstructed multiple times "
                         "in spill: "
                      << particle_starting_ke << " GeV" << std::endl;*/
          if (ismuon) {
            GetHist("reco_eff__multi_reco__muons",
                    "Muons which were reconstructed more than once",
                    "ke_tms_enter")
                ->Fill(particle_starting_ke);
            GetHist("reco_eff__multi_reco__probability_multi_reco_numerator",
                    "Chance of Getting Reco'd more than Once", "ke_tms_enter")
                ->Fill(particle_starting_ke);
          } else
            GetHist("reco_eff__multi_reco__nonmuon",
                    "Non-muons which were reconstructed more than once",
                    "ke_tms")
                ->Fill(particle_starting_ke);
          if (ismuon)
            DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                      "reco_eff/multi_reco_muon",
                      TString::Format("Particle %d reco'd %dx", ip,
                                      particle_indices_reconstructed[ip])
                          .Data(),
                      reco, lc, truth, DrawSliceN::many);
          else
            DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                      "reco_eff/multi_reco_nonmuon",
                      TString::Format("Particle %d reco'd %dx", ip,
                                      particle_indices_reconstructed[ip])
                          .Data(),
                      reco, lc, truth, DrawSliceN::many);
        }
      }
      if (ismuon) {
        GetHist("reco_eff__all_muon_ke_tms_enter_including_doubles_numerator",
                "Reco Efficiency w/2x vs True TMS-Entering KE, All Muons: "
                "Numerator",
                "ke_tms_enter")
            ->Fill(particle_starting_ke);
      }
      if (ismuon && not_double_reco) {
        GetHist("reco_eff__all_muon_ke_tms_enter_numerator",
                "Reco Efficiency vs True TMS-Entering KE, All Muons: Numerator",
                "ke_tms_enter")
            ->Fill(particle_starting_ke);
        GetHist("reco_eff__all_muon_ke_numerator",
                "Reco Efficiency vs True KE, All Muons: Numerator", "ke_tms")
            ->Fill(particle_ke);
        GetHist("reco_eff__multi_reco__probability_multi_reco_denominator",
                "Chance of Getting Reco'd more than Once", "ke_tms_enter")
            ->Fill(particle_starting_ke);

        GetHist("reco_eff__endpoint__muon_endpoint_x_numerator",
                "Reconstruction Efficiency: Numerator", "muon_endpoint_x")
            ->Fill(truth.PositionTMSEnd[ip][0] * CM);
        GetHist("reco_eff__endpoint__muon_endpoint_y_numerator",
                "Reconstruction Efficiency: Numerator", "muon_endpoint_y")
            ->Fill(truth.PositionTMSEnd[ip][1] * CM);
        GetHist("reco_eff__endpoint__muon_endpoint_z_numerator",
                "Reconstruction Efficiency: Numerator", "muon_endpoint_z")
            ->Fill(truth.PositionTMSEnd[ip][2] * CM);

        GetHist("reco_eff__startpoint__muon_startpoint_x_numerator",
                "Reconstruction Efficiency: Numerator", "muon_startpoint_x")
            ->Fill(truth.PositionTMSStart[ip][0] * CM);
        GetHist("reco_eff__startpoint__muon_startpoint_y_numerator",
                "Reconstruction Efficiency: Numerator", "muon_startpoint_y")
            ->Fill(truth.PositionTMSStart[ip][1] * CM);
        GetHist("reco_eff__startpoint__muon_startpoint_z_numerator",
                "Reconstruction Efficiency: Numerator", "muon_startpoint_z")
            ->Fill(truth.PositionTMSStart[ip][2] * CM);
            
        if (particle_starting_ke < 1)
            DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                      "reco_eff/reco_low_e_muon",
                      TString::Format("Low E muon: %.1f GeV",
                                      particle_starting_ke)
                          .Data(),
                      reco, lc, truth, DrawSliceN::few);
        if (particle_starting_ke < 0.25)
            DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                      "reco_eff/reco_very_low_e_muon",
                      TString::Format("Low E muon: %.2f GeV",
                                      particle_starting_ke)
                          .Data(),
                      reco, lc, truth, DrawSliceN::few);
      }
      if (ispion && not_double_reco) {
        GetHist("reco_eff__no_lar_tms_cuts__all_pion_ke_tms_enter_numerator",
                "Reco Efficiency, All Pions: Numerator", "pion_ke_tms_enter")
            ->Fill(particle_starting_ke);
        DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                  "reco_eff/reco_pion",
                  TString::Format("Pion with %.1f GeV",
                                  particle_starting_ke)
                      .Data(),
                  reco, lc, truth, DrawSliceN::few);
      }
      if (ismuon && lar_start && tms_end && not_double_reco) {
        GetSpecialHist("special__reco_eff_muon_ke_tms_enter_numerator",
                "reco_eff__muon_ke_tms_enter_numerator",
                "Reco Efficiency vs True TMS-Entering KE: Numerator",
                "ke_tms_enter")
            ->Fill(particle_starting_ke);
        int charge = reco.Charge[it];
        int true_charge = (truth.RecoTrackPrimaryParticlePDG[it] < 0) ? -1 : 1;
        bool correct_charge_id = true_charge * charge > 0;
        /*double true_muon_starting_ke =
          truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][3] * 1e-3;
        double reco_ke = length_to_energy(reco.Length[it]);
        bool correct_energy = std::abs((reco_ke - true_muon_starting_ke) / true_muon_starting_ke) < 0.05;
        
        GetHist("reco_eff__debugging__energy_reco",
                "Reco Energy",
                "ke_tms_enter")
            ->Fill(reco_ke);
        GetHist("reco_eff__debugging__energy_true",
                "Reco Energy",
                "ke_tms_enter")
            ->Fill(true_muon_starting_ke);
        REGISTER_AXIS(fractional_e_resolution, std::make_tuple("Fractional E", 41, -1.0, 1.0));
        GetHist("reco_eff__debugging__energy_resolution_fraction",
                "Fractional E Resolution",
                "fractional_e_resolution")
            ->Fill((reco_ke - true_muon_starting_ke) / true_muon_starting_ke);*/
        double true_z = truth.RecoTrackPrimaryParticleTruePositionEnd[it][2] * CM;
        double reco_z = reco.EndPos[it][2] * CM;
        double dz = reco_z - true_z;
        bool correct_endpoint = -30 < dz && dz < 0;
        bool is_good_reco = correct_charge_id && correct_endpoint;
        if (is_good_reco)
          GetHist("reco_eff__good_reco_muon_ke_tms_enter_numerator",
                  "Reco Efficiency vs True TMS-Entering KE: Numerator",
                  "ke_tms_enter")
              ->Fill(particle_starting_ke);
        GetHist("reco_eff__muon_ke_numerator",
                "Reco Efficiency vs True KE, TMS ending: Numerator", "ke_tms")
            ->Fill(particle_ke);
        double particle_starting_angle =
            std::atan2(
                truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][0],
                truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][2]) *
            DEG;
        GetHist("reco_eff__angle_tms_enter_numerator",
                "Reco Efficiency vs True TMS-Entering Angle: Numerator",
                "angle_tms_enter")
            ->Fill(particle_starting_angle);
      }
    } // end if (tms_touch)
  }
}
