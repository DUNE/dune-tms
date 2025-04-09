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
    double particle_starting_ke = truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][3] * GEV;
    if (ismuon || isantimuon) {
      // incorrect charge id
      if (nd_physics_muon && ismuon && correct_charge_id)
        GetSpecialHist("nd_physics_sample__charge_id_muon_numerator", "charge_id__true_muon_numerator",
                "#mu^{-} Charge ID Efficiency vs True KE_{Enter}",
                "ke_tms_enter", "#Efficiency")
            ->Fill(particle_starting_ke);
      if (nd_physics_muon && isantimuon && correct_charge_id)
        GetSpecialHist("nd_physics_sample__charge_id_antimuon_numerator", "charge_id__true_antimuon_numerator",
                "#mu^{+} Charge ID Efficiency vs True KE_{Enter}",
                "ke_tms_enter", "#Efficiency")
            ->Fill(particle_starting_ke);
      if (nd_physics_muon && ismuon)
        GetSpecialHist("nd_physics_sample__charge_id_muon_denominator", "charge_id__true_muon_denominator",
                "#mu^{-} Charge ID Efficiency vs True KE_{Enter}",
                "ke_tms_enter", "#Efficiency")
            ->Fill(particle_starting_ke);
      if (nd_physics_muon && isantimuon)
        GetSpecialHist("nd_physics_sample__charge_id_antimuon_denominator", "charge_id__true_antimuon_denominator",
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
                  "Startpoint of Incorrectedly Identified Antimuons",
                  "Z", "X")
              ->Fill(sz, sx);
        if (ismuon)
          GetHist("charge_id__validation__startpoint_muon",
                  "Startpoint of Incorrectedly Identified Muons",
                  "Z", "X")
              ->Fill(sz, sx);
        if (isantimuon)
          GetHist("charge_id__validation__endpoint_antimuons",
                  "Endpoint of Incorrectedly Identified Antimuons",
                  "Z", "X")
              ->Fill(ez, ex);
        if (ismuon)
          GetHist("charge_id__validation__endpoint_muon",
                  "Endpoint of Incorrectedly Identified Muons",
                  "Z", "X")
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
