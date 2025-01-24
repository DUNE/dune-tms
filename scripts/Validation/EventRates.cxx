// Add scope to avoid cross talk with other scripts
{
  if (on_new_spill) {
    std::map<int, double> vertices_with_visible_e;
    std::map<int, int> vertices_seen;
    std::map<int, double> vertices_with_visible_e_above_minimum;
    for (int ip = 0; ip < truth.nTrueParticles; ip++) {
      bool tms_energy = truth.TrueVisibleEnergy[ip] >= MINIMUM_VISIBLE_ENERGY;
      bool ismuon = abs(truth.PDG[ip]) == 13;
      bool ispion = abs(truth.PDG[ip]) == 211;
      bool lar_end = truth.LArFiducialEnd[ip];
      bool lar_start = truth.LArFiducialStart[ip];
      TVector3 birth_position(truth.BirthPosition[ip][0]*CM,
                              truth.BirthPosition[ip][1]*CM,
                              truth.BirthPosition[ip][2]*CM);
      bool lar_fiducial_start = IsInLAr(birth_position, LArFiducial);
      bool tms_touch = truth.TMSFiducialTouch[ip];
      bool tms_end = truth.TMSFiducialEnd[ip];
      bool is_rock_muon = false;
      TVector3 birth(truth.BirthPosition[ip][0], truth.BirthPosition[ip][1], truth.BirthPosition[ip][2]);
      bool start_in_lar_full = IsInLAr(birth*CM, LArFull);
      bool start_in_tms_full = IsInTMSFull(birth);
      is_rock_muon = ismuon && !start_in_lar_full && !start_in_tms_full;
      
      int vertex_id = truth.VertexID[ip];
      if (truth.TrueVisibleEnergy[ip] > 0) {
        vertices_with_visible_e[vertex_id] += truth.TrueVisibleEnergy[ip];
        if (tms_energy) vertices_with_visible_e_above_minimum[vertex_id] += truth.TrueVisibleEnergy[ip];
      }
      vertices_seen[vertex_id] += 1;
      
      //if (ismuon) std::cout<<TString::Format("Rock muon status: %d, %d, %d, %d", is_rock_muon, ismuon, start_in_lar_full, start_in_tms_full)<<std::endl;
      if (is_rock_muon) {
        GetHist("event_rates__rock_muons__rock_muons", "Rock Muon Event Rates",
          "rock_muon_event_rate")->Fill(0);
        if (lar_end) GetHist("event_rates__rock_muons__rock_muons", "Rock Muon Event Rates",
          "rock_muon_event_rate")->Fill(1);
        else if (tms_end) GetHist("event_rates__rock_muons__rock_muons", "Rock Muon Event Rates",
          "rock_muon_event_rate")->Fill(2);
        else GetHist("event_rates__rock_muons__rock_muons", "Rock Muon Event Rates",
          "rock_muon_event_rate")->Fill(3);
          
        GetHist("event_rates__rock_muons__validation__rock_muon_x", "Rock Muon X Position",
          "X")->Fill(birth.X() * CM);
        GetHist("event_rates__rock_muons__validation__rock_muon_y", "Rock Muon Y Position",
          "Y_full")->Fill(birth.Y() * CM);
        GetHist("event_rates__rock_muons__validation__rock_muon_z", "Rock Muon Z Position",
          "Z_full")->Fill(birth.Z() * CM);
        GetHist("event_rates__rock_muons__validation__rock_muon_xy", "Rock Muon XY Position",
          "X", "Y_full")->Fill(birth.X() * CM, birth.Y() * CM);
        GetHist("event_rates__rock_muons__validation__rock_muon_xz", "Rock Muon XZ Position",
          "Z_full", "X")->Fill(birth.Z() * CM, birth.X() * CM);
        GetHist("event_rates__rock_muons__validation__rock_muon_yz", "Rock Muon YZ Position",
          "Z_full", "Y_full")->Fill(birth.Z() * CM, birth.Y() * CM);
      }
      
      bool is_nd_physics_muon = lar_fiducial_start && ismuon; 
      if (is_nd_physics_muon) {
        GetHist("event_rates__tms__tms_event_rates", "TMS Event Rates",
          "tms_event_rate")->Fill(0);
        if (tms_touch) GetHist("event_rates__tms__tms_event_rates", "TMS Event Rates",
          "tms_event_rate")->Fill(1);
        if (tms_end) GetHist("event_rates__tms__tms_event_rates", "TMS Event Rates",
          "tms_event_rate")->Fill(2);
        if (tms_touch && !tms_end) GetHist("event_rates__tms__tms_event_rates", "TMS Event Rates",
          "tms_event_rate")->Fill(3);
      }
      if (lar_start && ismuon) {
        GetHist("event_rates__tms__tms_event_rates_lar_active", "TMS Event Rates",
          "tms_event_rate")->Fill(0);
        if (tms_touch) GetHist("event_rates__tms__tms_event_rates_lar_active", "TMS Event Rates",
          "tms_event_rate")->Fill(1);
        if (tms_end) GetHist("event_rates__tms__tms_event_rates_lar_active", "TMS Event Rates",
          "tms_event_rate")->Fill(2);
        if (tms_touch && !tms_end) GetHist("event_rates__tms__tms_event_rates_lar_active", "TMS Event Rates",
          "tms_event_rate")->Fill(3);
      }
    } // end for loop over particles
    
    // Now loop over vertices
    REGISTER_AXIS(true_vertex_energy, std::make_tuple("True Visible E (GeV)", 51, 0, 2));
    for (auto& vi : vertices_seen) {
      GetHist("event_rates__tms__tms_event_rates_by_vertex", "TMS Event Rates by Vertex",
        "tms_event_rate_by_vertex")->Fill(0);
    }
    for (auto& vi : vertices_with_visible_e) {
      GetHist("event_rates__tms__validation__true_vertex_energy", "True Visible E in TMS from Single Vtx",
        "true_vertex_energy")->Fill(vi.second * GEV);
      GetHist("event_rates__tms__tms_event_rates_by_vertex", "TMS Event Rates by Vertex",
        "tms_event_rate_by_vertex")->Fill(1);
    }
    for (auto& vi : vertices_with_visible_e_above_minimum) {
      GetHist("event_rates__tms__validation__true_vertex_energy_with_minimum_cut", "True Vis. E in TMS from Single Vtx with Min. Part. E",
        "true_vertex_energy")->Fill(vi.second * GEV);
      GetHist("event_rates__tms__tms_event_rates_by_vertex", "TMS Event Rates by Vertex",
        "tms_event_rate_by_vertex")->Fill(2);
    }
  } // end if on new spill
  
  for (int it = 0; it < truth.RecoTrackN; it++) {
    bool ismuon = abs(truth.RecoTrackPrimaryParticlePDG[it]) == 13;
    bool lar_start = truth.RecoTrackPrimaryParticleLArFiducialStart[it];
    TVector3 birth_position(truth.RecoTrackPrimaryParticleTruePositionStart[it][0]*CM,
                            truth.RecoTrackPrimaryParticleTruePositionStart[it][1]*CM,
                            truth.RecoTrackPrimaryParticleTruePositionStart[it][2]*CM);
    bool lar_fiducial_start = IsInLAr(birth_position, LArFiducial);
    bool tms_touch = truth.RecoTrackPrimaryParticleTMSFiducialTouch[it];
    bool tms_end = truth.RecoTrackPrimaryParticleTMSFiducialEnd[it];
    bool is_nd_physics_muon = lar_fiducial_start && ismuon; 
    if (is_nd_physics_muon) {
      GetHist("event_rates__tms__tms_muon_reco_rates", "TMS N Muons Reconstructed",
        "tms_muon_reco_rates")->Fill(0);
      if (tms_touch) GetHist("event_rates__tms__tms_muon_reco_rates", "TMS N Muons Reconstructed",
        "tms_muon_reco_rates")->Fill(1);
      if (tms_end) GetHist("event_rates__tms__tms_muon_reco_rates", "TMS N Muons Reconstructed",
        "tms_muon_reco_rates")->Fill(2);
      if (tms_touch && !tms_end) GetHist("event_rates__tms__tms_muon_reco_rates", "TMS N Muons Reconstructed",
        "tms_muon_reco_rates")->Fill(3);
    }
    if (lar_start && ismuon) {
      GetHist("event_rates__tms__tms_muon_reco_rates_lar_active_region", "TMS N Muons Reconstructed",
        "tms_muon_reco_rates")->Fill(0);
      if (tms_touch) GetHist("event_rates__tms__tms_muon_reco_rates_lar_active_region", "TMS N Muons Reconstructed",
        "tms_muon_reco_rates")->Fill(1);
      if (tms_end) GetHist("event_rates__tms__tms_muon_reco_rates_lar_active_region", "TMS N Muons Reconstructed",
        "tms_muon_reco_rates")->Fill(2);
      if (tms_touch && !tms_end) GetHist("event_rates__tms__tms_muon_reco_rates_lar_active_region", "TMS N Muons Reconstructed",
        "tms_muon_reco_rates")->Fill(3);
    }
    if (lar_fiducial_start) {
      GetHist("event_rates__reco_pdg__pdg", "Primary PDG Reconstructed",
        "pdg")->Fill(PDGtoIndex(truth.RecoTrackPrimaryParticlePDG[it]));
      GetHist("event_rates__reco_pdg__fiducial_comparison_pdg_nostack_1_fiducial", "Primary PDG Reco: Fiducial",
        "pdg")->Fill(PDGtoIndex(truth.RecoTrackPrimaryParticlePDG[it]));
    }
    if (lar_start) {
      GetHist("event_rates__reco_pdg__pdg_lar_active_region", "Primary PDG Reco., LAr Active Region",
        "pdg")->Fill(PDGtoIndex(truth.RecoTrackPrimaryParticlePDG[it]));
      GetHist("event_rates__reco_pdg__fiducial_comparison_pdg_nostack_0_all", "Primary PDG Reco: Active Region",
        "pdg")->Fill(PDGtoIndex(truth.RecoTrackPrimaryParticlePDG[it]));
    }
  }
}
