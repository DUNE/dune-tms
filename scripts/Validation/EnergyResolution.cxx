// Add scope to avoid cross talk with other scripts
{
  // True muon that starts in LAr and ends in TMS
  bool passes_truth_cuts = true;
  // Reco starts in front of TMS, and end is inside containment zone far from edges of TMS
  bool passes_reco_cuts  = true;

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
      
      TVector3 birth_pos(truth.BirthPosition[ip][0], truth.BirthPosition[ip][1], truth.BirthPosition[ip][2]);
      bool start_lar = LArFiducialCut(birth_pos);
      TVector3 end_pos(truth.DeathPosition[ip][0],
                       truth.DeathPosition[ip][1],
                       truth.DeathPosition[ip][2]);
      bool end_tms = isTMSContained(end_pos) == 0;
      if (start_lar && end_tms) {
        double particle_starting_ke = truth.MomentumTMSStart[ip][3] * GEV;
        GetHist("energy_resolution__selection_eff__selection_eff_ke_tms_enter_denominator",
                "Selection Efficiency vs True TMS-Entering KE: Denominator",
                "ke_tms_enter")->Fill(particle_starting_ke);
        GetHist("energy_resolution__selection_eff__selection_eff_contained_ke_tms_enter_denominator",
                "Selection Efficiency vs True TMS-Entering KE: Denominator",
                "ke_tms_enter")->Fill(particle_starting_ke);
      }
    }   
  }

  for (int it = 0; it < reco.nTracks; it++) {
    bool ismuon = abs(truth.RecoTrackPrimaryParticlePDG[it]) == 13;
    if (ismuon) {
      // Event has a muon
      has_muon = true;

      double true_muon_starting_ke = truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][3] * 1e-3;
      if (highest_true_muon_starting_ke < true_muon_starting_ke) { 
        highest_true_muon_starting_ke = true_muon_starting_ke;
        true_areal_density_of_highest = truth.RecoTrackPrimaryParticleTrueTrackLengthInTMS[it];
      }
      double reco_muon_starting_ke = (82+1.75*reco.Length[it])*1e-3;
      if (highest_reco_starting_muon_ke < reco_muon_starting_ke) {  
        highest_reco_starting_muon_ke = reco_muon_starting_ke;
        reco_areal_density_of_highest = reco.Length[it];
      }

      // Now check if it started in LAr and ended in TMS
      /*TVector3 birth_pos(truth.RecoTrackPrimaryParticleTruePositionStart[it][0], 
          truth.RecoTrackPrimaryParticleTruePositionStart[it][1], 
          truth.RecoTrackPrimaryParticleTruePositionTrackStart[it][2]);
      bool start_lar = LArFiducialCut(birth_pos);
      TVector3 end_pos(truth.RecoTrackPrimaryParticleTruePositionEnd[it][0],
                       truth.RecoTrackPrimaryParticleTruePositionEnd[it][1],
                       truth.RecoTrackPrimaryParticleTruePositionEnd[it][2]);
      bool end_tms = isTMSContained(end_pos) == 0;*/
      bool start_lar = truth.RecoTrackPrimaryParticleLArFiducialStart[it];
      bool end_tms = truth.RecoTrackPrimaryParticleTMSFiducialEnd[it];
      if (!start_lar || !end_tms) { 
        passes_truth_cuts = false;
      }
      else has_good_muon = true;
      
      if (end_tms) has_contained_muon = true;

      // Now check reco cuts
      TVector3 birth_pos_reco(reco.StartPos[it][0], reco.StartPos[it][1], reco.StartPos[it][2]);
      bool start_tms_reco = birth_pos_reco.Z() < 11362 + 30; // Check that z is in front of TMS
      TVector3 end_pos_reco(reco.EndPos[it][0],
                            reco.EndPos[it][1],
                            reco.EndPos[it][2]);
      bool end_tms_reco = isTMSContained(end_pos_reco) == 0;
      if (!start_tms_reco || !end_tms_reco) { 
        passes_reco_cuts = false;
      }
      GetHist("energy_resolution__cut_validation__passes_reco_cuts", "passes_reco_cuts", "yesno")->Fill(passes_reco_cuts ? 0 : 1);
      GetHist("energy_resolution__cut_validation__passes_truth_cuts", "passes_truth_cuts", "yesno")->Fill(passes_truth_cuts ? 0 : 1);
      GetHist("energy_resolution__cut_validation__start_tms_reco", "start_tms_reco", "yesno")->Fill(start_tms_reco ? 0 : 1);
      GetHist("energy_resolution__cut_validation__end_tms_reco", "end_tms_reco", "yesno")->Fill(end_tms_reco ? 0 : 1);
      GetHist("energy_resolution__cut_validation__truth_start_lar", "start_lar", "yesno")->Fill(start_lar ? 0 : 1);
      GetHist("energy_resolution__cut_validation__truth_end_tms", "end_tms", "yesno")->Fill(end_tms ? 0 : 1);
      if (start_tms_reco && end_tms_reco) {
        GetHist("energy_resolution__cuts__reco_end_pos_z_passes", "Reco End Pos Z: Passes Reco Cuts", "Z")->Fill(end_pos_reco.Z() * CM);
        GetHist("energy_resolution__cuts__reco_start_pos_z_passes", "Reco Start Pos Z: Passes Reco Cuts", "Z")->Fill(birth_pos_reco.Z() * CM);
      }
      if (end_tms_reco)
        GetHist("energy_resolution__cuts__reco_end_pos_z_passes_end", "Reco End Pos Z: Passes Reco Cuts", "Z")->Fill(end_pos_reco.Z() * CM);
      if (start_tms_reco)
        GetHist("energy_resolution__cuts__reco_start_pos_z_passes_start", "Reco Start Pos Z: Passes Reco Cuts", "Z")->Fill(birth_pos_reco.Z() * CM);
      GetHist("energy_resolution__cuts__reco_end_pos_z_all", "Reco End Pos Z: All Muons", "Z")->Fill(end_pos_reco.Z() * CM);
      GetHist("energy_resolution__cuts__reco_start_pos_z_all", "Reco Start Pos Z: All Muons", "Z")->Fill(birth_pos_reco.Z() * CM);
    }
  }

  // Only want events with a muon
  if (!has_muon) passes_truth_cuts = false;
  
  double resolution = highest_reco_starting_muon_ke - highest_true_muon_starting_ke;
  double fractional_resolution = resolution / highest_true_muon_starting_ke;
  double residual_areal_density = reco_areal_density_of_highest - true_areal_density_of_highest;
  
  // Plot some areal density info first, since it doesn't rely on our formula
  REGISTER_AXIS(true_areal_density, std::make_tuple("True Areal Density (g/cm^2)", 50, 0.0, 3000.0));
  REGISTER_AXIS(reco_areal_density, std::make_tuple("Reco Areal Density (g/cm^2)", 50, 0.0, 3000.0));
  REGISTER_AXIS(residual_areal_density, std::make_tuple("Residual Areal Density (reco - true) (g/cm^2)", 40, -200.0, 200.0));
  REGISTER_AXIS(residual_ke, std::make_tuple("Residual TMS-Entering KE (reco - true) (GeV)", 40, -1.0, 1.0));
  REGISTER_AXIS(basic_true_ke_enter, std::make_tuple("True TMS-Entering KE (GeV)", 50, 0.0, 5.0));
  if (has_muon) {
    GetHist("energy_resolution__areal_density__all_areal_density_true",
            "Areal Density: True",
            "true_areal_density")->Fill(true_areal_density_of_highest);
    GetHist("energy_resolution__areal_density__all_areal_density_reco",
            "Areal Density: Reco",
            "reco_areal_density")->Fill(reco_areal_density_of_highest);
    GetHist("energy_resolution__areal_density__all_areal_density_comparison",
            "Areal Density", "true_areal_density",
            "reco_areal_density")->Fill(true_areal_density_of_highest, reco_areal_density_of_highest);
  }
  if (has_contained_muon) {
    GetHist("energy_resolution__areal_density__contained_areal_density_true",
            "Areal Density: True",
            "true_areal_density")->Fill(true_areal_density_of_highest);
    GetHist("energy_resolution__areal_density__contained_areal_density_reco",
            "Areal Density: Reco",
            "reco_areal_density")->Fill(reco_areal_density_of_highest);
    GetHist("energy_resolution__areal_density__contained_areal_density_comparison",
            "Areal Density", "true_areal_density",
            "reco_areal_density")->Fill(true_areal_density_of_highest, reco_areal_density_of_highest);
    GetHist("energy_resolution__areal_density__contained_areal_density_comparison_column_normalize",
            "Areal Density", "true_areal_density",
            "reco_areal_density")->Fill(true_areal_density_of_highest, reco_areal_density_of_highest);
    GetHist("energy_resolution__areal_density__contained_areal_density_residual",
            "Residual Areal Density", "true_areal_density",
            "residual_areal_density")->Fill(true_areal_density_of_highest, residual_areal_density);
    GetHist("energy_resolution__areal_density__contained_areal_density_residual_column_maximize",
            "Residual Areal Density, Column Peak-Normalized", "true_areal_density",
            "residual_areal_density")->Fill(true_areal_density_of_highest, residual_areal_density);
  }
  if (passes_truth_cuts && passes_reco_cuts) {
    GetHist("energy_resolution__areal_density__areal_density_true",
            "Areal Density: True",
            "true_areal_density")->Fill(true_areal_density_of_highest);
    GetHist("energy_resolution__areal_density__areal_density_reco",
            "Areal Density: Reco",
            "reco_areal_density")->Fill(reco_areal_density_of_highest);
    GetHist("energy_resolution__areal_density__areal_density_comparison",
            "Areal Density", "true_areal_density",
            "reco_areal_density")->Fill(true_areal_density_of_highest, reco_areal_density_of_highest);
    GetHist("energy_resolution__areal_density__areal_density_comparison_column_normalize",
            "Areal Density", "true_areal_density",
            "reco_areal_density")->Fill(true_areal_density_of_highest, reco_areal_density_of_highest);
  }
  

  // Now plot some resolution information
  if (has_muon) {
    GetHist("energy_resolution__resolution__muon_starting_ke_resolution_raw",
            "Muon Resolution: All muons", "ke_tms_enter_true",
            "ke_tms_enter_reco")->Fill(highest_true_muon_starting_ke, highest_reco_starting_muon_ke);
    GetHist("energy_resolution__resolution__muon_starting_ke_fractional_resolution_raw",
            "Muon Resolution: All muons", "energy_resolution")->Fill(fractional_resolution);
  }
  if (has_good_muon) GetHist("energy_resolution__resolution__muon_starting_ke_resolution_fid",
                              "Muon Resolution: All True Fiducial Muons", "ke_tms_enter_true",
                             "ke_tms_enter_reco")->Fill(highest_true_muon_starting_ke, highest_reco_starting_muon_ke);
                             
                             
  if (has_contained_muon) { 
    GetHist("energy_resolution__resolution__muon_starting_ke_resolution_contained",
            "Muon Resolution: True Contained Muon", "ke_tms_enter_true",
            "ke_tms_enter_reco")->Fill(highest_true_muon_starting_ke, highest_reco_starting_muon_ke);
    GetHist("energy_resolution__resolution__muon_starting_ke_resolution_residual_column_maximize_contained",
            "Residual Muon KE, Column Peak-Normalized", "basic_true_ke_enter",
            "residual_ke")->Fill(highest_true_muon_starting_ke, resolution);
    GetHist("energy_resolution__resolution__muon_starting_ke_fractional_resolution_contained",
            "Muon Resolution: True Contained Muon", "energy_resolution")->Fill(fractional_resolution);
            
    if (passes_truth_cuts)
      GetHist("energy_resolution__selection_eff__selection_eff_contained_ke_tms_enter_numerator",
              "Selection Efficiency vs True TMS-Entering KE: Numerator",
              "ke_tms_enter")->Fill(highest_true_muon_starting_ke);
  }

  if (passes_truth_cuts && passes_reco_cuts) {
    GetHist("energy_resolution__resolution__muon_starting_ke_resolution",
            "Muon Resolution: ND Physics Sample", "ke_tms_enter_true",
            "ke_tms_enter_reco")->Fill(highest_true_muon_starting_ke, highest_reco_starting_muon_ke);
    GetHist("energy_resolution__resolution__muon_starting_ke_resolution_column_normalized",
            "Muon Resolution: Column Normalized", "ke_tms_enter_true",
            "ke_tms_enter_reco")->Fill(highest_true_muon_starting_ke, highest_reco_starting_muon_ke);
    GetHist("energy_resolution__resolution__muon_starting_ke_resolution_row_normalized",
            "Muon Resolution: Row Normalized", "ke_tms_enter_true",
            "ke_tms_enter_reco")->Fill(highest_true_muon_starting_ke, highest_reco_starting_muon_ke);
    GetHist("energy_resolution__resolution__muon_starting_ke_fractional_resolution",
            "Muon Resolution: ND Physics Sample", "energy_resolution")->Fill(fractional_resolution);
    GetHist("energy_resolution__resolution__muon_starting_ke_resolution_residual_column_maximize",
            "Residual Muon KE, Column Peak-Normalized", "basic_true_ke_enter",
            "residual_ke")->Fill(highest_true_muon_starting_ke, resolution);
            
    GetHist("energy_resolution__selection_eff__selection_eff_ke_tms_enter_numerator",
            "Selection Efficiency vs True TMS-Entering KE: Numerator",
            "ke_tms_enter")->Fill(highest_true_muon_starting_ke);
  }
}
