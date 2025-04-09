// Add scope to avoid cross talk with other scripts
{
  // True muon that starts in LAr and ends in TMS
  bool passes_truth_cuts = true;
  // Reco starts in front of TMS, and end is inside containment zone far from
  // edges of TMS
  bool passes_reco_cuts = true;

  // Event truly has a muon
  bool has_muon = false;
  // Event truly has a muon that starts in LAr and ends in TMS
  bool has_good_muon = false;
  // Event truly has a muon that ends in TMS
  bool has_contained_muon = false;
  double highest_true_muon_starting_ke = -9e-12;
  double highest_reco_starting_muon_ke = -9e-12;
  double true_areal_density_of_highest = -9e-12;
  double reco_areal_density_of_highest = -9e-12;
  double highest_reco_starting_muon_ke_uncorrected = -9e-12;
  double highest_reco_starting_muon_ke_without_fit = -9e-12;
  double true_muon_ke_lar = -9e-12;
  double true_muon_areal_density_lar = -9e12; // starting from lar to end
  double true_muon_areal_density_lar_only =
      -9e12; // only lar (and window) component

  if (on_new_spill) {

    for (int ip = 0; ip < truth.nTrueParticles; ip++) {
      // Only consider those that have enough visible energy in the TMS to count
      /*bool tms_touch = truth.TrueVisibleEnergy[ip] >= MINIMUM_VISIBLE_ENERGY;
      if (tms_touch) {
        bool ismuon = abs(truth.PDG[ip]) == 13;
        bool tms_end = truth.TMSFiducialEnd[ip];
        double particle_starting_ke = truth.MomentumTMSStart[ip][3] * GEV;
        //double particle_ke = truth.BirthMomentum[ip][3] * GEV;
        GetHist("energy_resolution__selection_eff__selection_eff_ke_tms_enter_denominator",
                "Selection Efficiency vs True TMS-Entering KE: Denominator",
                "ke_tms_enter")->Fill(particle_starting_ke);
      }*/

      // TVector3 birth_pos(truth.BirthPosition[ip][0],
      // truth.BirthPosition[ip][1], truth.BirthPosition[ip][2]); bool start_lar
      // = LArFiducialCut(birth_pos); TVector3
      // end_pos(truth.DeathPosition[ip][0],
      //                  truth.DeathPosition[ip][1],
      //                  truth.DeathPosition[ip][2]);
      bool start_lar = truth.LArFiducialStart[ip];
      bool end_tms = true; // isTMSContained(end_pos) == 0;
      if (start_lar && end_tms &&
          abs(truth.RecoTrackPrimaryParticlePDG[ip]) == 13) {
        double particle_starting_ke = truth.MomentumTMSStart[ip][3] * GEV;
        GetHist("energy_resolution__selection_eff__selection_eff_ke_tms_enter_"
                "denominator",
                "Selection Efficiency vs True TMS-Entering KE: Denominator",
                "ke_tms_enter")
            ->Fill(particle_starting_ke);
        GetHist("energy_resolution__selection_eff__selection_eff_contained_ke_"
                "tms_enter_denominator",
                "Selection Efficiency vs True TMS-Entering KE: Denominator",
                "ke_tms_enter")
            ->Fill(particle_starting_ke);
      }
    }
  }

  for (int it = 0; it < reco.nTracks; it++) {
    bool ismuon = abs(truth.RecoTrackPrimaryParticlePDG[it]) == 13;
    if (ismuon) {
      // Event has a muon
      has_muon = true;

      double true_muon_starting_ke =
          truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][3] * 1e-3;
      if (highest_true_muon_starting_ke < true_muon_starting_ke) {
        highest_true_muon_starting_ke = true_muon_starting_ke;
        true_areal_density_of_highest =
            truth.RecoTrackPrimaryParticleTrueTrackLengthInTMS[it];

        int particle_index = truth.RecoTrackPrimaryParticleIndex[it];
        if (particle_index >= 0 && particle_index < truth.nTrueParticles) {
          true_muon_ke_lar = truth.BirthMomentum[particle_index][3] * 1e-3;
        } else {
          std::cout << "Warning: Didn't get a valid particle index: "
                    << particle_index << std::endl;
          true_muon_ke_lar =
              truth.RecoTrackPrimaryParticleTrueMomentum[it][3] * 1e-3;
        }
        true_muon_areal_density_lar =
            truth.RecoTrackPrimaryParticleTrueTrackLength[it];
        true_muon_areal_density_lar_only =
            true_muon_areal_density_lar - true_areal_density_of_highest;
      }
      double length_to_use = reco.Length[it];
      TVector3 direction(reco.EndDirection[it][0], 0, reco.EndDirection[it][2]);
      if (direction.Mag() == 0)
        direction.SetZ(1); // Normalize to 1 in z since we don't know direction
      else
        direction.SetMag(1); // Normalize to 1
      #ifdef add_half_steel
      // add half thickness of steel
      // plus the avg we're off because we're failing to reco the track end
      // note that length is in g/cm^2
      double density_to_use = -99999.0;
      if (reco.EndPos[it][2] < 13500) {        // thin section
        length_to_use += 5.89 * direction.Z(); // g/cm^2
        density_to_use = 4.43;                 // g/cm^2, 1/2 steel, 1/2 plastic
      } else {                                 // thick section
        length_to_use += 11.8 * direction.Z(); // g/cm^2
        density_to_use = 5.57; // g/cm^2, 2/3rds steel, 1/3rd plastic
      }
      // On avg we're off by about this much in cm
      // direction.SetMag(30);
      length_to_use += density_to_use * direction.Z();
      #endif // ifdef add_half_steel
      double reco_muon_starting_ke = length_to_energy(length_to_use);
      if (highest_reco_starting_muon_ke < reco_muon_starting_ke) {
        highest_reco_starting_muon_ke = reco_muon_starting_ke;
        highest_reco_starting_muon_ke_uncorrected =
            length_to_energy(reco.Length[it]);
        highest_reco_starting_muon_ke_without_fit =
            default_length_to_energy(length_to_use);
        reco_areal_density_of_highest = length_to_use;
      }

      // Now check if it started in LAr and ended in TMS
      /*TVector3
      birth_pos(truth.RecoTrackPrimaryParticleTruePositionStart[it][0],
          truth.RecoTrackPrimaryParticleTruePositionStart[it][1],
          truth.RecoTrackPrimaryParticleTruePositionTrackStart[it][2]);
      bool start_lar = LArFiducialCut(birth_pos); */
      /*TVector3 end_pos(truth.RecoTrackPrimaryParticleTruePositionEnd[it][0],
                       truth.RecoTrackPrimaryParticleTruePositionEnd[it][1],
                       truth.RecoTrackPrimaryParticleTruePositionEnd[it][2]);
      bool end_tms = isTMSContained(end_pos) == 0;*/
      bool start_lar = truth.RecoTrackPrimaryParticleLArFiducialStart[it];
      bool end_tms = true; // truth.RecoTrackPrimaryParticleTMSFiducialEnd[it];
      if (!start_lar || !end_tms) {
        passes_truth_cuts = false;
      } else
        has_good_muon = true;

      if (truth.RecoTrackPrimaryParticleTMSFiducialEnd[it])
        has_contained_muon = true;

      // Now check reco cuts
      TVector3 birth_pos_reco(reco.StartPos[it][0], reco.StartPos[it][1],
                              reco.StartPos[it][2]);
      bool start_tms_reco =
          birth_pos_reco.Z() < 11362 + 30; // Check that z is in front of TMS
      TVector3 end_pos_reco(reco.EndPos[it][0], reco.EndPos[it][1],
                            reco.EndPos[it][2]);
      bool end_tms_reco = TMSEndpointUncontained(end_pos_reco) == 0;
      if (!start_tms_reco || !end_tms_reco) {
        passes_reco_cuts = false;
      }
      GetHist("energy_resolution__cut_validation__passes_reco_cuts",
              "passes_reco_cuts", "yesno")
          ->Fill(passes_reco_cuts ? 0 : 1);
      GetHist("energy_resolution__cut_validation__passes_truth_cuts",
              "passes_truth_cuts", "yesno")
          ->Fill(passes_truth_cuts ? 0 : 1);
      GetHist("energy_resolution__cut_validation__start_tms_reco",
              "start_tms_reco", "yesno")
          ->Fill(start_tms_reco ? 0 : 1);
      GetHist("energy_resolution__cut_validation__end_tms_reco", "end_tms_reco",
              "yesno")
          ->Fill(end_tms_reco ? 0 : 1);
      GetHist("energy_resolution__cut_validation__truth_start_lar", "start_lar",
              "yesno")
          ->Fill(start_lar ? 0 : 1);
      GetHist("energy_resolution__cut_validation__truth_end_tms", "end_tms",
              "yesno")
          ->Fill(end_tms ? 0 : 1);
      if (start_tms_reco && end_tms_reco) {
        GetHist("energy_resolution__cuts__reco_end_pos_z_passes",
                "Reco End Pos Z: Passes Reco Cuts", "Z")
            ->Fill(end_pos_reco.Z() * CM);
        GetHist("energy_resolution__cuts__reco_start_pos_z_passes",
                "Reco Start Pos Z: Passes Reco Cuts", "Z")
            ->Fill(birth_pos_reco.Z() * CM);
      }
      if (end_tms_reco)
        GetHist("energy_resolution__cuts__reco_end_pos_z_passes_end",
                "Reco End Pos Z: Passes Reco Cuts", "Z")
            ->Fill(end_pos_reco.Z() * CM);
      if (start_tms_reco)
        GetHist("energy_resolution__cuts__reco_start_pos_z_passes_start",
                "Reco Start Pos Z: Passes Reco Cuts", "Z")
            ->Fill(birth_pos_reco.Z() * CM);
      GetHist("energy_resolution__cuts__reco_end_pos_z_all",
              "Reco End Pos Z: All Muons", "Z")
          ->Fill(end_pos_reco.Z() * CM);
      GetHist("energy_resolution__cuts__reco_start_pos_z_all",
              "Reco Start Pos Z: All Muons", "Z")
          ->Fill(birth_pos_reco.Z() * CM);
    }
  }

  // Only want events with a muon
  if (!has_muon)
    passes_truth_cuts = false;

  double resolution =
      highest_reco_starting_muon_ke - highest_true_muon_starting_ke;
  double fractional_resolution = resolution / highest_true_muon_starting_ke;
  double fractional_resolution_uncorrected =
      (highest_reco_starting_muon_ke_uncorrected -
       highest_true_muon_starting_ke) /
      highest_true_muon_starting_ke;
  double fractional_resolution_without_fit =
      (highest_reco_starting_muon_ke_without_fit -
       highest_true_muon_starting_ke) /
      highest_true_muon_starting_ke;

  // double lar_component_reco_ke_estimate =
  // lar_length_to_energy(true_muon_areal_density_lar_only);
  double lar_component_true = true_muon_ke_lar - highest_true_muon_starting_ke;
  double lar_component_reco_ke_estimate = lar_component_true;
  double lar_reco_ke =
      highest_reco_starting_muon_ke + lar_component_reco_ke_estimate;
  double lar_fractional_resolution =
      (lar_reco_ke - true_muon_ke_lar) / true_muon_ke_lar;
  double lar_component_fractional_resolution =
      (lar_component_reco_ke_estimate - lar_component_true) /
      lar_component_true;

  // Plot some areal density info first, since it doesn't rely on our formula
  REGISTER_AXIS(
      true_areal_density,
      std::make_tuple("True Areal Density (g/cm^2)", 50, 0.0, 3000.0));
  REGISTER_AXIS(
      reco_areal_density,
      std::make_tuple("Reco Areal Density (g/cm^2)", 50, 0.0, 3000.0));
  REGISTER_AXIS(residual_ke,
                std::make_tuple("Residual TMS-Entering KE (reco - true) (GeV)",
                                40, -1.0, 1.0));
  REGISTER_AXIS(residual_ke_slice,
                std::make_tuple("Residual TMS-Entering KE (reco - true) (GeV)",
                                21, -1.0, 1.0));
  REGISTER_AXIS(
      energy_resolution_slice,
      std::make_tuple("KE Resolution (reco - true) / true", 31, -0.4, 0.4));
  REGISTER_AXIS(basic_true_ke_enter,
                std::make_tuple("True TMS-Entering KE (GeV)", 50, 0.0, 5.0));
  REGISTER_AXIS(basic_reco_ke_enter,
                std::make_tuple("Reco TMS-Entering KE (GeV)", 50, 0.0, 5.0));
  REGISTER_AXIS(basic_true_ke,
                std::make_tuple("True Muon KE (GeV)", 30, 0.0, 5.0));
  REGISTER_AXIS(basic_reco_ke,
                std::make_tuple("Reco Muon KE (GeV)", 30, 0.0, 5.0));
  if (has_muon) {
    // All muons
    GetHist(
        "energy_resolution__areal_density__all_muons__all_areal_density_true",
        "Areal Density: True", "true_areal_density")
        ->Fill(true_areal_density_of_highest);
    GetHist(
        "energy_resolution__areal_density__all_muons__all_areal_density_reco",
        "Areal Density: Reco", "reco_areal_density")
        ->Fill(reco_areal_density_of_highest);
    GetHist("energy_resolution__areal_density__all_muons__all_areal_density_"
            "comparison",
            "Areal Density", "true_areal_density", "reco_areal_density")
        ->Fill(true_areal_density_of_highest, reco_areal_density_of_highest);
    GetHist("energy_resolution__areal_density__all_muons__all_areal_density_"
            "comparison_column_normalize",
            "Areal Density", "true_areal_density", "reco_areal_density")
        ->Fill(true_areal_density_of_highest, reco_areal_density_of_highest);
  }
  if (passes_truth_cuts && passes_reco_cuts) {
    // Selected sample
    GetHist("energy_resolution__areal_density__areal_density_true",
            "Areal Density: True", "true_areal_density")
        ->Fill(true_areal_density_of_highest);
    GetHist("energy_resolution__areal_density__areal_density_reco",
            "Areal Density: Reco", "reco_areal_density")
        ->Fill(reco_areal_density_of_highest);
    GetHist("energy_resolution__areal_density__areal_density_comparison",
            "Areal Density", "true_areal_density", "reco_areal_density")
        ->Fill(true_areal_density_of_highest, reco_areal_density_of_highest);
    GetHist("energy_resolution__areal_density__areal_density_comparison_column_"
            "normalize",
            "Areal Density", "true_areal_density", "reco_areal_density")
        ->Fill(true_areal_density_of_highest, reco_areal_density_of_highest);

    //
    GetHist("energy_resolution__areal_density_and_energy__vs_energy_comparison_"
            "true",
            "True Areal Density Compared to Energy", "true_areal_density",
            "basic_true_ke_enter")
        ->Fill(true_areal_density_of_highest, highest_true_muon_starting_ke);
    GetHist("energy_resolution__areal_density_and_energy__vs_energy_comparison_"
            "column_normalize_true",
            "True Areal Density Compared to Energy, Col Norm.",
            "true_areal_density", "basic_true_ke_enter")
        ->Fill(true_areal_density_of_highest, highest_true_muon_starting_ke);

    REGISTER_AXIS(
        true_areal_density_zoom,
        std::make_tuple("True Areal Density (g/cm^2)", 70, 0.0, 2100.0));
    REGISTER_AXIS(fractional_error_ke_over_density,
                  std::make_tuple("True KE / Areal Density (MeV / g/cm^2)", 30,
                                  1.4, 2.4));
    REGISTER_AXIS(fractional_error_density_over_ke,
                  std::make_tuple("True Areal Density / True KE (g/cm^2 / MeV)",
                                  30, 0.3, 0.8));
    std::vector<double> de_to_try = {0.5, 0.25, 0.1, 0.01, 0.001, 0.0001};
    int de_index = 0;
    for (auto de : de_to_try) {
      std::vector<double> energy_slices = {0.5, 1, 1.5, 2, 2.5, 3};
      std::vector<double> energy_slices_high;
      std::vector<double> energy_slices_low;
      for (auto slice_e : energy_slices) {
        energy_slices_high.push_back(slice_e + de);
        energy_slices_low.push_back(slice_e - de);
      }
      for (size_t ie = 0; ie < energy_slices_low.size(); ie++) {
        double e_low = energy_slices_low.at(ie);
        double e_high = energy_slices_high.at(ie);
        if (e_low < highest_true_muon_starting_ke &&
            highest_true_muon_starting_ke <= e_high) {
          int precision = 5;
          if (de >= 0.0001)
            precision = 4;
          if (de >= 0.001)
            precision = 3;
          if (de >= 0.01)
            precision = 2;
          if (de >= 0.1)
            precision = 1;
          if (de == 0.25)
            precision = 2;
          if (de >= 1)
            precision = 0;
          GetHist(
              TString::Format(
                  "energy_resolution__areal_density_and_energy__energy_slices__"
                  "energy_slices_dev%d_normalize_to_one_nostack_%ld",
                  de_index, ie)
                  .Data(),
              TString::Format("True Areal Density in Slices of TMS-Entering "
                              "KE: %.*lf < E / GeV <= %.*lf",
                              precision, e_low, precision, e_high)
                  .Data(),
              "true_areal_density_zoom", "Normalized Frequency")
              ->Fill(true_areal_density_of_highest);
          double fractional_ke_over_density =
              true_areal_density_of_highest /
              (highest_true_muon_starting_ke * 1000);
          GetHist(TString::Format(
                      "energy_resolution__areal_density_and_energy__energy_"
                      "slices__energy_fractional_ke_over_density_slices_dev%d_"
                      "normalize_to_probability_nostack_%ld",
                      de_index, ie)
                      .Data(),
                  TString::Format("True Areal Density / KE in Slices of "
                                  "TMS-Entering KE: %.*lf < E / GeV <= %.*lf",
                                  precision, e_low, precision, e_high)
                      .Data(),
                  "fractional_error_density_over_ke", "Probability")
              ->Fill(fractional_ke_over_density);
        }
      }
      de_index += 1;
    }

    REGISTER_AXIS(true_kinetic_energy_zoom,
                  std::make_tuple("True TMS-Entering KE (GeV)", 60, 0.0, 4.5));
    int da_index = 0;
    std::vector<double> density_slices;
    double density_change = 2400 / 6;
    for (int id = 0; id < 6; id++) {
      double density = (id + 0.5) * density_change;
      density_slices.push_back(density);
    }
    std::vector<double> da_to_try = {density_change, density_change * 0.05, 50,
                                     20, 10};
    for (auto dd : da_to_try) {
      std::vector<double> density_slices_high;
      std::vector<double> density_slices_low;
      for (auto slice_d : density_slices) {
        double offset_to_use = 0;
        if (da_index == 0)
          offset_to_use = -0.5 * density_change;
        if (da_index == 1)
          offset_to_use = 0;
        if (da_index == 2)
          offset_to_use = 0;
        if (da_index == 3)
          offset_to_use = 0;
        if (da_index == 4)
          offset_to_use = 0;
        density_slices_high.push_back(slice_d + dd + offset_to_use);
        density_slices_low.push_back(slice_d + offset_to_use);
      }
      for (size_t id = 0; id < density_slices_low.size(); id++) {
        double d_low = density_slices_low.at(id);
        double d_high = density_slices_high.at(id);
        if (d_low < true_areal_density_of_highest &&
            true_areal_density_of_highest <= d_high) {
          int precision = 0;
          if (dd >= 1)
            precision = 0;
          GetHist(
              TString::Format(
                  "energy_resolution__areal_density_and_energy__density_slices_"
                  "_density_slices_dev%d_normalize_to_one_nostack_%ld",
                  da_index, id)
                  .Data(),
              TString::Format("True TMS-Entering KE in Slices of True Areal "
                              "Density: %.*lf < #rho / g/cm^2 <= %.*lf",
                              precision, d_low, precision, d_high)
                  .Data(),
              "true_kinetic_energy_zoom", "Normalized Frequency")
              ->Fill(highest_true_muon_starting_ke);
          double fractional_density_over_ke =
              true_areal_density_of_highest /
              (highest_true_muon_starting_ke * 1000);
          double fractional_ke_over_density = highest_true_muon_starting_ke *
                                              1000 /
                                              true_areal_density_of_highest;
          GetHist(
              TString::Format(
                  "energy_resolution__areal_density_and_energy__density_slices_"
                  "_density_fractional_ke_over_density_slices_dev%d_normalize_"
                  "to_probability_nostack_%ld",
                  da_index, id)
                  .Data(),
              TString::Format("True KE / Areal Density in Slices of True Areal "
                              "Density: %.*lf < #rho / g/cm^2 <= %.*lf",
                              precision, d_low, precision, d_high)
                  .Data(),
              "fractional_error_ke_over_density", "Probability")
              ->Fill(fractional_ke_over_density);
          GetHist(
              TString::Format(
                  "energy_resolution__areal_density_and_energy__density_slices_"
                  "_density_fractional_density_over_ke_slices_dev%d_normalize_"
                  "to_probability_nostack_%ld",
                  da_index, id)
                  .Data(),
              TString::Format("True Areal Density / KE in Slices of True Areal "
                              "Density: %.*lf < #rho / g/cm^2 <= %.*lf",
                              precision, d_low, precision, d_high)
                  .Data(),
              "fractional_error_density_over_ke", "Probability")
              ->Fill(fractional_density_over_ke);
        }
      }
      da_index += 1;
    }

    //
    double target_e = (1.75 * true_areal_density_of_highest) * 1e-3 + 0.5;
    if (target_e < highest_true_muon_starting_ke) {
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                "true_areal_density_not_correlated_with_ke",
                TString::Format("n tracks = %d", reco.nTracks).Data(), reco, lc,
                truth, DrawSliceN::many);
      GetHist("energy_resolution__areal_density_and_energy__vs_energy_"
              "comparison_drawn_slices",
              "True Areal Density Compared to Energy", "true_areal_density",
              "basic_true_ke_enter")
          ->Fill(true_areal_density_of_highest, highest_true_muon_starting_ke);
    }
  }

  // Now plot some resolution information
  if (has_muon) {

    GetHist("energy_resolution__resolution__all_muons__muon_starting_ke_"
            "resolution_raw",
            "Muon Resolution: All muons", "ke_tms_enter_true",
            "ke_tms_enter_reco")
        ->Fill(highest_true_muon_starting_ke, highest_reco_starting_muon_ke);
    GetHist("energy_resolution__resolution__all_muons__muon_starting_ke_"
            "resolution_raw_column_normalized",
            "Muon Resolution: All muons", "ke_tms_enter_true",
            "ke_tms_enter_reco")
        ->Fill(highest_true_muon_starting_ke, highest_reco_starting_muon_ke);
    GetHist("energy_resolution__resolution__all_muons__muon_starting_ke_"
            "fractional_resolution_raw",
            "Muon Resolution: All muons", "energy_resolution")
        ->Fill(fractional_resolution);

    GetHist("energy_resolution__resolution__all_muon_starting_ke_fractional_"
            "resolution_nostack_0_raw",
            "Muon Resolution: All muons", "energy_resolution")
        ->Fill(fractional_resolution);
    GetHist("energy_resolution__resolution__only_two_muon_starting_ke_"
            "fractional_resolution_nostack_0_raw",
            "Muon Resolution: All muons", "energy_resolution")
        ->Fill(fractional_resolution);
    GetHist("energy_resolution__resolution__uncontained_comparison_muon_"
            "starting_ke_fractional_resolution_nostack_0_raw",
            "Muon Resolution: All muons", "energy_resolution")
        ->Fill(fractional_resolution);
  }
  if (has_good_muon) {
    GetHist("energy_resolution__resolution__all_muon_starting_ke_fractional_"
            "resolution_nostack_2_good",
            "Muon Resolution: Start LAr, End TMS", "energy_resolution")
        ->Fill(fractional_resolution);
  }

  if (has_contained_muon) {
    GetHist("energy_resolution__resolution__all_muon_starting_ke_fractional_"
            "resolution_nostack_1_contained",
            "Muon Resolution: End TMS", "energy_resolution")
        ->Fill(fractional_resolution);
    GetHist("energy_resolution__resolution__uncontained_comparison_muon_"
            "starting_ke_fractional_resolution_nostack_1_contained",
            "Muon Resolution: Contained Muons", "energy_resolution")
        ->Fill(fractional_resolution);

    if (passes_truth_cuts)
      GetHist("energy_resolution__selection_eff__selection_eff_contained_ke_"
              "tms_enter_numerator",
              "Selection Efficiency vs True TMS-Entering KE: Numerator",
              "ke_tms_enter")
          ->Fill(highest_true_muon_starting_ke);
  } else if (has_muon) {
    // Uncontained muon
    GetHist("energy_resolution__resolution__uncontained_comparison_muon_"
            "starting_ke_fractional_resolution_nostack_2_uncontained",
            "Muon Resolution: Uncontained Muons", "energy_resolution")
        ->Fill(fractional_resolution);
  }

  if (passes_truth_cuts && passes_reco_cuts) {
    GetHist("energy_resolution__resolution__muon_starting_ke_resolution",
            "Muon Resolution: Selected Muons", "ke_tms_enter_true",
            "ke_tms_enter_reco")
        ->Fill(highest_true_muon_starting_ke, highest_reco_starting_muon_ke);
    GetHist("energy_resolution__resolution__muon_starting_ke_resolution_column_"
            "normalized",
            "Muon Resolution: Column Normalized", "ke_tms_enter_true",
            "ke_tms_enter_reco")
        ->Fill(highest_true_muon_starting_ke, highest_reco_starting_muon_ke);
    GetHist("energy_resolution__resolution__muon_starting_ke_resolution_row_"
            "normalized",
            "Muon Resolution: Row Normalized", "ke_tms_enter_true",
            "ke_tms_enter_reco")
        ->Fill(highest_true_muon_starting_ke, highest_reco_starting_muon_ke);
    GetHist(
        "energy_resolution__resolution__muon_starting_ke_fractional_resolution",
        "Muon Resolution: Selected Muons", "energy_resolution")
        ->Fill(fractional_resolution);
    GetHist("energy_resolution__resolution__muon_starting_ke_resolution_"
            "residual_column_maximize",
            "Residual Muon KE, Column Peak-Normalized", "basic_true_ke_enter",
            "residual_ke")
        ->Fill(highest_true_muon_starting_ke, resolution);

    GetHist("energy_resolution__selection_eff__selection_eff_ke_tms_enter_"
            "numerator",
            "Selection Efficiency vs True TMS-Entering KE: Numerator",
            "ke_tms_enter")
        ->Fill(highest_true_muon_starting_ke);

    GetHist("energy_resolution__resolution__all_muon_starting_ke_fractional_"
            "resolution_nostack_3_nd_physics_sample",
            "Muon Resolution: Selected Muons", "energy_resolution")
        ->Fill(fractional_resolution);
    GetHist("energy_resolution__resolution__only_two_muon_starting_ke_"
            "fractional_resolution_nostack_3_nd_physics_sample",
            "Muon Resolution: Selected Muons", "energy_resolution")
        ->Fill(fractional_resolution);

    GetHist("energy_resolution__resolution__corrected_starting_ke_fractional_"
            "resolution_nostack_1_corrected",
            "Muon Resolution: Corrected", "energy_resolution")
        ->Fill(fractional_resolution);
    GetHist("energy_resolution__resolution__corrected_starting_ke_fractional_"
            "resolution_nostack_3_uncorrected",
            "Muon Resolution: Uncorrected", "energy_resolution")
        ->Fill(fractional_resolution_uncorrected);
    GetHist("energy_resolution__resolution__fit_starting_ke_fractional_"
            "resolution_nostack_1_with_fit",
            "Muon Resolution: With Fit", "energy_resolution")
        ->Fill(fractional_resolution);
    GetHist("energy_resolution__resolution__fit_starting_ke_fractional_"
            "resolution_nostack_3_no_fit",
            "Muon Resolution: Without Fit", "energy_resolution")
        ->Fill(fractional_resolution_without_fit);

    GetHist("energy_resolution__resolution__contained_starting_ke_fractional_"
            "resolution_nostack_0_all",
            "Muon Resolution: All Muons", "energy_resolution")
        ->Fill(fractional_resolution);
    if (has_contained_muon) {
      GetHist("energy_resolution__resolution__contained_starting_ke_fractional_"
              "resolution_nostack_1_contained",
              "Muon Resolution: Contained", "energy_resolution")
          ->Fill(fractional_resolution);
    } else {
      GetHist("energy_resolution__resolution__contained_starting_ke_fractional_"
              "resolution_nostack_2_uncontained",
              "Muon Resolution: Uncontained", "energy_resolution")
          ->Fill(fractional_resolution);
    }

    GetHist("energy_resolution__lar_resolution__lar_muon_fractional_resolution",
            "LAr-starting Muon Resolution", "energy_resolution")
        ->Fill(lar_fractional_resolution);
    GetHist("energy_resolution__lar_resolution__lar_muon_resolution_column_"
            "normalized",
            "LAr-starting Muon Resolution, Column Normalized", "basic_true_ke",
            "basic_reco_ke")
        ->Fill(true_muon_ke_lar, lar_reco_ke);
    GetHist("energy_resolution__lar_resolution__lar_muon_resolution",
            "LAr-starting Muon Resolution", "basic_true_ke", "basic_reco_ke")
        ->Fill(true_muon_ke_lar, lar_reco_ke);
    REGISTER_AXIS(basic_lar_component_ke_true,
                  std::make_tuple("True Muon KE (GeV)", 30, 0.0, 5.0));
    REGISTER_AXIS(basic_lar_component_ke_reco,
                  std::make_tuple("Reco Muon KE (GeV)", 30, 0.0, 5.0));
    REGISTER_AXIS(
        basic_lar_component_areal_density,
        std::make_tuple("True Muon Areal Density (g/cm^2)", 30, 0.0, 3000.0));
    GetHist(
        "energy_resolution__lar_resolution__lar_component_estimate_comparison",
        "LAr-component Muon Resolution", "basic_lar_component_ke_true",
        "basic_lar_component_ke_reco")
        ->Fill(lar_component_true, lar_component_reco_ke_estimate);
    GetHist("energy_resolution__lar_resolution__lar_component_areal_density",
            "LAr Muon True KE vs areal density", "basic_lar_component_ke_true",
            "basic_lar_component_areal_density")
        ->Fill(true_muon_ke_lar - highest_true_muon_starting_ke,
               true_muon_areal_density_lar_only);
    GetHist("energy_resolution__lar_resolution__lar_component_muon_fractional_"
            "resolution",
            "LAr Component Muon Resolution", "energy_resolution")
        ->Fill(lar_component_fractional_resolution);

    if (highest_reco_starting_muon_ke + 1 < highest_true_muon_starting_ke)
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                "poor_reco_starting_muon_ke",
                TString::Format("n tracks = %d", reco.nTracks).Data(), reco, lc,
                truth, DrawSliceN::many);

    std::vector<double> cuts = {0.0,  0.25, 0.5,  0.75, 1.0,  1.25, 1.5,
                                1.75, 2.0,  2.25, 2.5,  2.75, 3.0,  3.25,
                                3.5,  3.75, 4.0,  4.25, 4.5,  4.75, 5.0};
    for (size_t ic = 0; ic < cuts.size() - 1; ic++) {
      double min = cuts[ic];
      double max = cuts[ic + 1];
      if (min <= highest_true_muon_starting_ke &&
          highest_true_muon_starting_ke < max) {
        GetHist(Form("energy_resolution__resolution__slices__muon_%ld", ic),
                Form("Muon Resolution: %.2f < True KE / GeV< %.2f", min, max),
                "energy_resolution_slice")
            ->Fill(fractional_resolution);
        GetHist(Form("energy_resolution__resolution__residual_slices__muon_%ld",
                     ic),
                Form("Muon Resolution: %.2f < True KE / GeV < %.2f", min, max),
                "residual_ke_slice")
            ->Fill(resolution);
      }
    }
  }
}
