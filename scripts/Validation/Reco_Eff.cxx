// Add scope to avoid cross talk with other scripts
{
  REGISTER_AXIS(ke_tms, std::make_tuple("True Muon KE (GeV)", 28, 0.0, 7.0));
  REGISTER_AXIS(pion_ke_tms_enter, std::make_tuple("True Pion KE Entering TMS (GeV)", 20, 0.0, 5.0));
  REGISTER_AXIS(muon_endpoint_z, std::make_tuple("True Muon Endpoint Z (cm)", 20, 1136.0, 1831.0));
  REGISTER_AXIS(muon_endpoint_y, std::make_tuple("True Muon Endpoint Y (cm)", 20, -216.0, 16.0));
  REGISTER_AXIS(muon_endpoint_x, std::make_tuple("True Muon Endpoint X (cm)", 20, -330.0, 330.0));
  REGISTER_AXIS(muon_startpoint_z, std::make_tuple("True Muon Startpoint Z (cm)", 20, 1136.0, 1831.0));
  REGISTER_AXIS(muon_startpoint_y, std::make_tuple("True Muon Startpoint Y (cm)", 20, -216.0, 16.0));
  REGISTER_AXIS(muon_startpoint_x, std::make_tuple("True Muon Startpoint X (cm)", 20, -330.0, 330.0));
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
          GetHist("reco_eff__no_lar_tms_cuts__all_muon_ke_tms_enter_denominator", "Reconstruction Efficiency: Denominator",
            "ke_tms_enter")->Fill(particle_starting_ke);
          GetHist("reco_eff__all_muon_ke_denominator", "Reco Efficiency vs True KE, All Muons: Denominator", 
            "ke_tms")->Fill(particle_ke);
            
          GetHist("reco_eff__endpoint__muon_endpoint_x_denominator", "Reconstruction Efficiency: Denominator",
            "muon_endpoint_x")->Fill(truth.PositionTMSEnd[ip][0] * CM);
          GetHist("reco_eff__endpoint__muon_endpoint_y_denominator", "Reconstruction Efficiency: Denominator",
            "muon_endpoint_y")->Fill(truth.PositionTMSEnd[ip][1] * CM);
          GetHist("reco_eff__endpoint__muon_endpoint_z_denominator", "Reconstruction Efficiency: Denominator",
            "muon_endpoint_z")->Fill(truth.PositionTMSEnd[ip][2] * CM);
            
          GetHist("reco_eff__startpoint__muon_startpoint_x_denominator", "Reconstruction Efficiency: Denominator",
            "muon_startpoint_x")->Fill(truth.PositionTMSStart[ip][0] * CM);
          GetHist("reco_eff__startpoint__muon_startpoint_y_denominator", "Reconstruction Efficiency: Denominator",
            "muon_startpoint_y")->Fill(truth.PositionTMSStart[ip][1] * CM);
          GetHist("reco_eff__startpoint__muon_startpoint_z_denominator", "Reconstruction Efficiency: Denominator",
            "muon_startpoint_z")->Fill(truth.PositionTMSStart[ip][2] * CM);
            
        }
        if (ispion)
          GetHist("reco_eff__no_lar_tms_cuts__all_pion_ke_tms_enter_denominator", "Reconstruction Efficiency: Denominator",
            "pion_ke_tms_enter")->Fill(particle_starting_ke);
        if (ismuon && lar_start && tms_end) {
          GetHist("reco_eff__muon_ke_tms_enter_denominator", "Reconstruction Efficiency: Denominator",
            "ke_tms_enter")->Fill(particle_starting_ke);
          GetHist("reco_eff__muon_ke_denominator", "Reco Efficiency vs True KE: Denominator", 
            "ke_tms")->Fill(particle_ke);
          double particle_starting_angle = std::atan2(truth.MomentumTMSStart[ip][0], truth.MomentumTMSStart[ip][2]) * DEG;
          GetHist("reco_eff__angle_tms_enter_denominator", "Reconstruction Efficiency: Denominator",
            "angle_tms_enter")->Fill(particle_starting_angle);
          
        }
      } // end if (tms_touch)
    } 
  }
      
  // Fill numerators of reco_eff plots here
  for (int it = 0; it < reco.nTracks; it++) {
    // Only consider those that have enough visible energy in the TMS to count
    bool tms_touch = truth.RecoTrackPrimaryParticleTrueVisibleEnergy[it] >= MINIMUM_VISIBLE_ENERGY;
    if (tms_touch) {
      bool ismuon = abs(truth.RecoTrackPrimaryParticlePDG[it]) == 13;
      bool ispion = abs(truth.RecoTrackPrimaryParticlePDG[it]) == 211;
      bool lar_start = truth.RecoTrackPrimaryParticleLArFiducialStart[it];
      bool tms_end = truth.RecoTrackPrimaryParticleTMSFiducialEnd[it];
      double particle_starting_ke = truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][3] * GEV;
      double particle_ke = -1;
      //truth.RecoTrackPrimaryParticleTrueMomentum[it][3] * GEV;
      int ip = truth.RecoTrackPrimaryParticleIndex[it];
      if (ip >= truth.nTrueParticles || ip < 0) {
        std::cout<<"Warning: RecoTrackPrimaryParticleIndex is broken: Found no true particle in reco eff with ip="<<ip<<std::endl;
        ip = -1;
      }
      if (ip >= 0) {
        particle_ke = truth.BirthMomentum[ip][3] * GEV;
      }
      if (ismuon) {
        GetHist("reco_eff__all_muon_ke_tms_enter_numerator", "Reco Efficiency vs True TMS-Entering KE, All Muons: Numerator", 
          "ke_tms_enter")->Fill(particle_starting_ke);
        GetHist("reco_eff__all_muon_ke_numerator", "Reco Efficiency vs True KE, All Muons: Numerator", 
          "ke_tms")->Fill(particle_ke);
            
        GetHist("reco_eff__endpoint__muon_endpoint_x_numerator", "Reconstruction Efficiency: Numerator",
          "muon_endpoint_x")->Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[it][0] * CM);
        GetHist("reco_eff__endpoint__muon_endpoint_y_numerator", "Reconstruction Efficiency: Numerator",
          "muon_endpoint_y")->Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[it][1] * CM);
        GetHist("reco_eff__endpoint__muon_endpoint_z_numerator", "Reconstruction Efficiency: Numerator",
          "muon_endpoint_z")->Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[it][2] * CM);
            
        GetHist("reco_eff__startpoint__muon_startpoint_x_numerator", "Reconstruction Efficiency: Numerator",
          "muon_startpoint_x")->Fill(truth.RecoTrackPrimaryParticleTruePositionStart[it][0] * CM);
        GetHist("reco_eff__startpoint__muon_startpoint_y_numerator", "Reconstruction Efficiency: Numerator",
          "muon_startpoint_y")->Fill(truth.RecoTrackPrimaryParticleTruePositionStart[it][1] * CM);
        GetHist("reco_eff__startpoint__muon_startpoint_z_numerator", "Reconstruction Efficiency: Numerator",
          "muon_startpoint_z")->Fill(truth.RecoTrackPrimaryParticleTruePositionStart[it][2] * CM);
      }
      if (ispion)
        GetHist("reco_eff__no_lar_tms_cuts__all_pion_ke_tms_enter_numerator", "Reco Efficiency, All Pions: Numerator", 
          "pion_ke_tms_enter")->Fill(particle_starting_ke);
      if (ismuon && lar_start && tms_end) {
        GetHist("reco_eff__muon_ke_tms_enter_numerator", "Reco Efficiency vs True TMS-Entering KE: Numerator", 
          "ke_tms_enter")->Fill(particle_starting_ke);
        GetHist("reco_eff__muon_ke_numerator", "Reco Efficiency vs True KE, TMS ending: Numerator", 
          "ke_tms")->Fill(particle_ke);
        double particle_starting_angle = std::atan2(truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][0],
                                                truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][2]) * DEG;
        GetHist("reco_eff__angle_tms_enter_numerator", "Reco Efficiency vs True TMS-Entering Angle: Numerator", 
          "angle_tms_enter")->Fill(particle_starting_angle);
      }
    } // end if (tms_touch)
  }
}
