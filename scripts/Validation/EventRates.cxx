// Add scope to avoid cross talk with other scripts
{
  // Same but with LAr shell cut
  bool lar_shell_cut = truth.LArOuterShellEnergyFromVertex < 30; // MeV
  if (true) { // on_new_spill || lar_shell_cut // shell energy plot needs every
              // event
    std::map<int, double> vertices_with_visible_e;
    std::map<int, int> vertices_seen;
    std::map<int, double> vertices_with_visible_e_above_minimum;
    for (int ip = 0; ip < truth.nTrueParticles; ip++) {
      bool tms_energy = truth.TrueVisibleEnergy[ip] >= MINIMUM_VISIBLE_ENERGY;
      bool ismuon = abs(truth.PDG[ip]) == 13;
      //bool ispion = abs(truth.PDG[ip]) == 211;
      bool lar_end = truth.LArFiducialEnd[ip];
      bool lar_start = truth.LArFiducialStart[ip];
      TVector3 birth_position(truth.BirthPosition[ip][0] * CM,
                              truth.BirthPosition[ip][1] * CM,
                              truth.BirthPosition[ip][2] * CM);
      bool lar_fiducial_start = IsInLAr(birth_position, LArFiducial);
      bool tms_touch = truth.TMSFiducialTouch[ip];
      bool tms_end = truth.TMSFiducialEnd[ip];
      bool is_rock_muon = false;
      TVector3 birth(truth.BirthPosition[ip][0], truth.BirthPosition[ip][1],
                     truth.BirthPosition[ip][2]);
      bool start_in_lar_full = IsInLAr(birth * CM, LArFull);
      bool start_in_tms_full = IsInTMSFull(birth);
      is_rock_muon = ismuon && !start_in_lar_full && !start_in_tms_full;
      // Check that this particle truly came from interaction added up in shell
      // cut
      bool particle_from_shell_cut_interaction =
          truth.VertexIdOfMostEnergyInEvent == truth.VertexID[ip];

      if (on_new_spill) {
        int vertex_id = truth.VertexID[ip];
        if (truth.TrueVisibleEnergy[ip] > 0) {
          vertices_with_visible_e[vertex_id] += truth.TrueVisibleEnergy[ip];
          if (tms_energy)
            vertices_with_visible_e_above_minimum[vertex_id] +=
                truth.TrueVisibleEnergy[ip];
        }
        vertices_seen[vertex_id] += 1;
      }

      // if (ismuon) std::cout<<TString::Format("Rock muon status: %d, %d, %d,
      // %d", is_rock_muon, ismuon, start_in_lar_full,
      // start_in_tms_full)<<std::endl;
      if (is_rock_muon && on_new_spill) {
        GetHist("event_rates__rock_muons__rock_muons", "Rock Muon Event Rates",
                "rock_muon_event_rate")
            ->Fill(0);
        if (lar_end)
          GetHist("event_rates__rock_muons__rock_muons",
                  "Rock Muon Event Rates", "rock_muon_event_rate")
              ->Fill(1);
        else if (tms_end)
          GetHist("event_rates__rock_muons__rock_muons",
                  "Rock Muon Event Rates", "rock_muon_event_rate")
              ->Fill(2);
        else
          GetHist("event_rates__rock_muons__rock_muons",
                  "Rock Muon Event Rates", "rock_muon_event_rate")
              ->Fill(3);

        GetHist("event_rates__rock_muons__validation__rock_muon_x",
                "Rock Muon X Position", "X")
            ->Fill(birth.X() * CM);
        GetHist("event_rates__rock_muons__validation__rock_muon_y",
                "Rock Muon Y Position", "Y_full")
            ->Fill(birth.Y() * CM);
        GetHist("event_rates__rock_muons__validation__rock_muon_z",
                "Rock Muon Z Position", "Z_full")
            ->Fill(birth.Z() * CM);
        GetHist("event_rates__rock_muons__validation__rock_muon_xy",
                "Rock Muon XY Position", "X", "Y_full")
            ->Fill(birth.X() * CM, birth.Y() * CM);
        GetHist("event_rates__rock_muons__validation__rock_muon_xz",
                "Rock Muon XZ Position", "Z_full", "X")
            ->Fill(birth.Z() * CM, birth.X() * CM);
        GetHist("event_rates__rock_muons__validation__rock_muon_yz",
                "Rock Muon YZ Position", "Z_full", "Y_full")
            ->Fill(birth.Z() * CM, birth.Y() * CM);
      }

      bool is_nd_physics_muon = lar_fiducial_start && ismuon;
      if (is_nd_physics_muon && on_new_spill) {
        GetHist("event_rates__tms__tms_event_rates", "TMS Event Rates",
                "tms_event_rate")
            ->Fill(0);
        if (tms_touch)
          GetHist("event_rates__tms__tms_event_rates", "TMS Event Rates",
                  "tms_event_rate")
              ->Fill(1);
        if (tms_end)
          GetHist("event_rates__tms__tms_event_rates", "TMS Event Rates",
                  "tms_event_rate")
              ->Fill(2);
        if (tms_touch && !tms_end)
          GetHist("event_rates__tms__tms_event_rates", "TMS Event Rates",
                  "tms_event_rate")
              ->Fill(3);
      }

      if (is_nd_physics_muon && lar_shell_cut &&
          particle_from_shell_cut_interaction) {
        // This is the full ND physics sample
        GetHist("event_rates__tms__tms_event_rates_with_shell_cut",
                "TMS Event Rates", "tms_event_rate")
            ->Fill(0);
        if (tms_touch)
          GetHist("event_rates__tms__tms_event_rates_with_shell_cut",
                  "TMS Event Rates", "tms_event_rate")
              ->Fill(1);
        if (tms_end)
          GetHist("event_rates__tms__tms_event_rates_with_shell_cut",
                  "TMS Event Rates", "tms_event_rate")
              ->Fill(2);
        if (tms_touch && !tms_end)
          GetHist("event_rates__tms__tms_event_rates_with_shell_cut",
                  "TMS Event Rates", "tms_event_rate")
              ->Fill(3);
        GetHist("event_rates__tms__validation__shell_energy_cut",
                "LAr Visible Hadronic Energy from Primary Vertex in 30cm Shell",
                "LArOuterShellEnergyFromVertex")
            ->Fill(truth.LArOuterShellEnergyFromVertex);
      }
      if (is_nd_physics_muon && particle_from_shell_cut_interaction) {
        GetHist("event_rates__tms__validation__shell_energy",
                "LAr Visible Hadronic Energy from Primary Vertex in 30cm Shell",
                "LArOuterShellEnergyFromVertex")
            ->Fill(truth.LArOuterShellEnergyFromVertex);
      }

      if (lar_start && ismuon && on_new_spill) {
        GetHist("event_rates__tms__tms_event_rates_lar_active",
                "TMS Event Rates", "tms_event_rate")
            ->Fill(0);
        if (tms_touch)
          GetHist("event_rates__tms__tms_event_rates_lar_active",
                  "TMS Event Rates", "tms_event_rate")
              ->Fill(1);
        if (tms_end)
          GetHist("event_rates__tms__tms_event_rates_lar_active",
                  "TMS Event Rates", "tms_event_rate")
              ->Fill(2);
        if (tms_touch && !tms_end)
          GetHist("event_rates__tms__tms_event_rates_lar_active",
                  "TMS Event Rates", "tms_event_rate")
              ->Fill(3);
      }
    } // end for loop over particles

    // Now loop over vertices once per spill
    if (on_new_spill) {
      REGISTER_AXIS(true_vertex_energy,
                    std::make_tuple("True Visible E (GeV)", 51, 0, 2));
      for (auto& vi : vertices_seen) {
        GetHist("event_rates__tms__tms_event_rates_by_vertex",
                "TMS Event Rates by Vertex", "tms_event_rate_by_vertex")
            ->Fill(0);
        (void)vi; // unused
      }
      for (auto &vi : vertices_with_visible_e) {
        GetHist("event_rates__tms__validation__true_vertex_energy",
                "True Visible E in TMS from Single Vtx", "true_vertex_energy")
            ->Fill(vi.second * GEV);
        GetHist("event_rates__tms__tms_event_rates_by_vertex",
                "TMS Event Rates by Vertex", "tms_event_rate_by_vertex")
            ->Fill(1);
      }
      for (auto &vi : vertices_with_visible_e_above_minimum) {
        GetHist(
            "event_rates__tms__validation__true_vertex_energy_with_minimum_cut",
            "True Vis. E in TMS from Single Vtx with Min. Part. E",
            "true_vertex_energy")
            ->Fill(vi.second * GEV);
        GetHist("event_rates__tms__tms_event_rates_by_vertex",
                "TMS Event Rates by Vertex", "tms_event_rate_by_vertex")
            ->Fill(2);
      }
    }
  } // end if on new spill or lar shell cut

  for (int it = 0; it < truth.RecoTrackN; it++) {
    bool ismuon = abs(truth.RecoTrackPrimaryParticlePDG[it]) == 13;
    bool lar_start = truth.RecoTrackPrimaryParticleLArFiducialStart[it];
    TVector3 birth_position(
        truth.RecoTrackPrimaryParticleTruePositionStart[it][0] * CM,
        truth.RecoTrackPrimaryParticleTruePositionStart[it][1] * CM,
        truth.RecoTrackPrimaryParticleTruePositionStart[it][2] * CM);
    bool lar_fiducial_start = IsInLAr(birth_position, LArFiducial);
    bool tms_touch = truth.RecoTrackPrimaryParticleTMSFiducialTouch[it];
    bool tms_end = truth.RecoTrackPrimaryParticleTMSFiducialEnd[it];
    bool is_nd_physics_muon = lar_fiducial_start && ismuon;
    bool lar_shell_cut = truth.LArOuterShellEnergyFromVertex < 30; // MeV
    if (is_nd_physics_muon) {
      GetHist("event_rates__tms__tms_muon_reco_rates",
              "TMS N Muons Reconstructed", "tms_muon_reco_rates")
          ->Fill(0);
      if (tms_touch)
        GetHist("event_rates__tms__tms_muon_reco_rates",
                "TMS N Muons Reconstructed", "tms_muon_reco_rates")
            ->Fill(1);
      if (tms_end)
        GetHist("event_rates__tms__tms_muon_reco_rates",
                "TMS N Muons Reconstructed", "tms_muon_reco_rates")
            ->Fill(2);
      if (tms_touch && !tms_end)
        GetHist("event_rates__tms__tms_muon_reco_rates",
                "TMS N Muons Reconstructed", "tms_muon_reco_rates")
            ->Fill(3);
    }
    if (is_nd_physics_muon && lar_shell_cut) {
      GetHist("event_rates__tms__tms_muon_reco_rates_with_shell_cut",
              "TMS N Muons Reconstructed", "tms_muon_reco_rates")
          ->Fill(0);
      if (tms_touch)
        GetHist("event_rates__tms__tms_muon_reco_rates_with_shell_cut",
                "TMS N Muons Reconstructed", "tms_muon_reco_rates")
            ->Fill(1);
      if (tms_end)
        GetHist("event_rates__tms__tms_muon_reco_rates_with_shell_cut",
                "TMS N Muons Reconstructed", "tms_muon_reco_rates")
            ->Fill(2);
      if (tms_touch && !tms_end)
        GetHist("event_rates__tms__tms_muon_reco_rates_with_shell_cut",
                "TMS N Muons Reconstructed", "tms_muon_reco_rates")
            ->Fill(3);
    }
    if (lar_start && ismuon) {
      GetHist("event_rates__tms__tms_muon_reco_rates_lar_active_region",
              "TMS N Muons Reconstructed", "tms_muon_reco_rates")
          ->Fill(0);
      if (tms_touch)
        GetHist("event_rates__tms__tms_muon_reco_rates_lar_active_region",
                "TMS N Muons Reconstructed", "tms_muon_reco_rates")
            ->Fill(1);
      if (tms_end)
        GetHist("event_rates__tms__tms_muon_reco_rates_lar_active_region",
                "TMS N Muons Reconstructed", "tms_muon_reco_rates")
            ->Fill(2);
      if (tms_touch && !tms_end)
        GetHist("event_rates__tms__tms_muon_reco_rates_lar_active_region",
                "TMS N Muons Reconstructed", "tms_muon_reco_rates")
            ->Fill(3);
    }
    if (lar_fiducial_start) {
      GetHist("event_rates__reco_pdg__pdg", "Primary PDG Reconstructed", "pdg")
          ->Fill(PDGtoIndex(truth.RecoTrackPrimaryParticlePDG[it]));
      GetHist(
          "event_rates__reco_pdg__fiducial_comparison_pdg_nostack_1_fiducial",
          "Primary PDG Reco: Fiducial", "pdg")
          ->Fill(PDGtoIndex(truth.RecoTrackPrimaryParticlePDG[it]));
    }
    if (lar_fiducial_start && lar_shell_cut) {
      GetHist("event_rates__reco_pdg__pdg_with_shell_cut",
              "Primary PDG Reconstructed", "pdg")
          ->Fill(PDGtoIndex(truth.RecoTrackPrimaryParticlePDG[it]));
      GetHist(
          "event_rates__reco_pdg__fiducial_comparison_pdg_nostack_2_nd_physics",
          "Primary PDG Reco: ND Physics Muons", "pdg")
          ->Fill(PDGtoIndex(truth.RecoTrackPrimaryParticlePDG[it]));
      GetHist("event_rates__tms__validation__shell_energy_reco_muons_cut",
              "LAr Visible Hadronic Energy from Primary Vertex in 30cm Shell",
              "LArOuterShellEnergyFromVertex")
          ->Fill(truth.LArOuterShellEnergyFromVertex);
    }
    if (lar_fiducial_start) {
      GetHist("event_rates__tms__validation__shell_energy_reco_muons",
              "LAr Visible Hadronic Energy from Primary Vertex in 30cm Shell",
              "LArOuterShellEnergyFromVertex")
          ->Fill(truth.LArOuterShellEnergyFromVertex);
    }
    if (lar_start) {
      GetHist("event_rates__reco_pdg__pdg_lar_active_region",
              "Primary PDG Reco., LAr Active Region", "pdg")
          ->Fill(PDGtoIndex(truth.RecoTrackPrimaryParticlePDG[it]));
      GetHist("event_rates__reco_pdg__fiducial_comparison_pdg_nostack_0_all",
              "Primary PDG Reco: Active Region", "pdg")
          ->Fill(PDGtoIndex(truth.RecoTrackPrimaryParticlePDG[it]));
    }
  }
}
