// Add scope to avoid cross talk with other scripts
{
  bool passes_truth_cuts = true;
  bool passes_reco_cuts  = true;

  bool has_muon = false;
  bool has_good_muon = false;
  bool has_contained_muon = false;
  double highest_true_muon_starting_ke = -9e-12;
  double highest_reco_starting_muon_ke = -9e-12;

  for (int it = 0; it < reco.nTracks; it++) {
    bool ismuon = abs(truth.RecoTrackPrimaryParticlePDG[it]) == 13;
    if (ismuon) {
      // Event has a muon
      has_muon = true;

      double true_muon_starting_ke = truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][3] * 1e-3;
      if (highest_true_muon_starting_ke < true_muon_starting_ke) highest_true_muon_starting_ke = true_muon_starting_ke;
      double reco_muon_starting_ke = (82+1.75*reco.Length[it])*1e-3;
      if (highest_reco_starting_muon_ke < reco_muon_starting_ke) highest_reco_starting_muon_ke = reco_muon_starting_ke;

      // Now check if it started in LAr and ended in TMS
      TVector3 birth_pos(truth.LeptonX4[0], truth.LeptonX4[1], truth.LeptonX4[2]);
      bool start_lar = LArFiducialCut(birth_pos);
      TVector3 end_pos(truth.RecoTrackPrimaryParticleTruePositionEnd[it][0],
                       truth.RecoTrackPrimaryParticleTruePositionEnd[it][1],
                       truth.RecoTrackPrimaryParticleTruePositionEnd[it][2]);
      bool end_tms = isTMSContained(end_pos) == 0;
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
    GetHist("energy_resolution__resolution__muon_starting_ke_fractional_resolution_contained",
            "Muon Resolution: True Contained Muon", "energy_resolution")->Fill(fractional_resolution);
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
  }
}
