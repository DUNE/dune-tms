{
  // raw variables
  // Fill some basic "raw" variables from the reco tree
  GetHist("basic__raw__EventNo", "EventNo", "EventNo", "#N Slices")
      ->Fill(reco.EventNo);
  GetHist("basic__raw__SliceNo", "SliceNo", "SliceNo", "#N Slices")
      ->Fill(reco.SliceNo);
  GetHist("basic__raw__SpillNo", "SpillNo", "SpillNo", "#N Slices")
      ->Fill(reco.SpillNo);

  GetSpecialHist("special__track_n", "basic__raw__ntracks", "N Reco Tracks", "ntracks", "#N Slices")
      ->Fill(reco.nTracks);
  for (int it = 0; it < reco.nTracks; it++) {
    GetSpecialHist("special__track_nhits", "basic__raw__nHits", "N Hits per Track", "n0-120", "#N Tracks")->Fill(reco.nHits[it]);
    GetHist("basic__raw__nKalmanNodes", "nKalmanNodes", "n0-120")
        ->Fill(reco.nKalmanNodes[it]);
    for (int ih = 0; ih < reco.nHits[it]; ih++) {
      GetHist("basic__raw__TrackHitPos_X", "TrackHitPos X", "X", "#N Hits")
          ->Fill(reco.TrackHitPos[it][ih][0] * CM);
      GetHist("basic__raw__TrackHitPos_Y", "TrackHitPos Y", "Y", "#N Hits")
          ->Fill(reco.TrackHitPos[it][ih][1] * CM);
      GetHist("basic__raw__TrackHitPos_Z", "TrackHitPos Z", "Z", "#N Hits")
          ->Fill(reco.TrackHitPos[it][ih][2] * CM);
    }
    for (int ih = 0; ih < reco.nKalmanNodes[it]; ih++) {
      GetHist("basic__raw__KalmanPos_X", "KalmanPos X", "X", "#N Nodes")
          ->Fill(reco.KalmanPos[it][ih][0] * CM);
      GetHist("basic__raw__KalmanPos_Y", "KalmanPos Y", "Y", "#N Nodes")
          ->Fill(reco.KalmanPos[it][ih][1] * CM);
      GetHist("basic__raw__KalmanPos_Z", "KalmanPos Z", "Z", "#N Nodes")
          ->Fill(reco.KalmanPos[it][ih][2] * CM);
      GetHist("basic__raw__KalmanTruePos_X", "KalmanTruePos X", "X", "#N Nodes")
          ->Fill(reco.KalmanTruePos[it][ih][0] * CM);
      GetHist("basic__raw__KalmanTruePos_Y", "KalmanTruePos Y", "Y", "#N Nodes")
          ->Fill(reco.KalmanTruePos[it][ih][1] * CM);
      GetHist("basic__raw__KalmanTruePos_Z", "KalmanTruePos Z", "Z", "#N Nodes")
          ->Fill(reco.KalmanTruePos[it][ih][2] * CM);
    }
    GetHist("basic__raw__StartDirection_X", "StartDirection X", "dx",
            "#N Tracks")
        ->Fill(reco.StartDirection[it][0]);
    GetHist("basic__raw__StartDirection_Y", "StartDirection Y", "dy",
            "#N Tracks")
        ->Fill(reco.StartDirection[it][1]);
    GetHist("basic__raw__StartDirection_Z", "StartDirection Z", "dz",
            "#N Tracks")
        ->Fill(reco.StartDirection[it][2]);
    GetHist("basic__raw__EndDirection_X", "EndDirection X", "dx", "#N Tracks")
        ->Fill(reco.EndDirection[it][0]);
    GetHist("basic__raw__EndDirection_Y", "EndDirection Y", "dy", "#N Tracks")
        ->Fill(reco.EndDirection[it][1]);
    GetHist("basic__raw__EndDirection_Z", "EndDirection Z", "dz", "#N Tracks")
        ->Fill(reco.EndDirection[it][2]);

    REGISTER_AXIS(DirectionSanityCheck,
                  std::make_tuple("Direction Mag", 51, 0, 10));
    double start_direction_mag =
        std::sqrt(reco.StartDirection[it][0] * reco.StartDirection[it][0] +
                  reco.StartDirection[it][1] * reco.StartDirection[it][1] +
                  reco.StartDirection[it][2] * reco.StartDirection[it][2]);
    double end_direction_mag =
        std::sqrt(reco.EndDirection[it][0] * reco.EndDirection[it][0] +
                  reco.EndDirection[it][1] * reco.EndDirection[it][1] +
                  reco.EndDirection[it][2] * reco.EndDirection[it][2]);
    // const double epsilon = 1e-6;
    //  if (std::abs(start_direction_mag - 1) > epsilon)
    //    std::cout<<"Warning: Found start_direction_mag outside of 1.
    //    "<<start_direction_mag<<std::endl;
    //  if (std::abs(end_direction_mag - 1) > epsilon)
    //    std::cout<<"Warning: Found end_direction_mag outside of 1.
    //    "<<end_direction_mag<<std::endl;
    GetHist("basic__sanity__StartDirectionMag", "StartDirection Mag",
            "DirectionSanityCheck", "#N Tracks")
        ->Fill(start_direction_mag);
    GetHist("basic__sanity__EndDirectionMag", "EndDirection Mag",
            "DirectionSanityCheck", "#N Tracks")
        ->Fill(end_direction_mag);
    /*if (std::abs(reco.StartDirection[it][1]) > 1)
      std::cout << "big y dir: " << reco.StartDirection[it][0] << ","
                << reco.StartDirection[it][1] << ","
                << reco.StartDirection[it][2] << "\t" << start_direction_mag
                << std::endl;*/

    const double big_epsilon = 1e-1;
    bool has_y_start_direction_zero =
        (std::abs(reco.StartDirection[it][1]) < big_epsilon);
    if (has_y_start_direction_zero)
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                "issues/y_start_direction_zero",
                TString::Format("n tracks = %d", reco.nTracks).Data(), reco, lc,
                truth, DrawSliceN::many);
    bool has_y_end_direction_zero =
        (std::abs(reco.EndDirection[it][1]) < big_epsilon);
    if (has_y_end_direction_zero)
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                "issues/y_end_direction_zero",
                TString::Format("n tracks = %d", reco.nTracks).Data(), reco, lc,
                truth, DrawSliceN::many);

    GetHist("basic__raw__StartPos_X", "StartPos X", "X", "#N Tracks")
        ->Fill(reco.StartPos[it][0] * CM);
    GetHist("basic__raw__StartPos_Y", "StartPos Y", "Y", "#N Tracks")
        ->Fill(reco.StartPos[it][1] * CM);
    GetHist("basic__raw__StartPos_Z", "StartPos Z", "Z", "#N Tracks")
        ->Fill(reco.StartPos[it][2] * CM);
    GetHist("basic__raw__EndPos_X", "EndPos X", "X", "#N Tracks")
        ->Fill(reco.EndPos[it][0] * CM);
    GetHist("basic__raw__EndPos_Y", "EndPos Y", "Y", "#N Tracks")
        ->Fill(reco.EndPos[it][1] * CM);
    GetHist("basic__raw__EndPos_Z", "EndPos Z", "Z", "#N Tracks")
        ->Fill(reco.EndPos[it][2] * CM);

    GetHist("basic__raw__Charge", "Charge", "charge", "#N Tracks")
        ->Fill(reco.Charge[it]);
    GetHist("basic__raw__Length", "Areal Density, AKA Length", "areal_density",
            "#N Tracks")
        ->Fill(reco.Length[it]);
    GetHist("basic__raw__Momentum", "Momentum", "momentum", "#N Tracks")
        ->Fill(reco.Momentum[it]);
    GetHist("basic__raw__EnergyRange", "EnergyRange", "energy_range",
            "#N Tracks")
        ->Fill(reco.EnergyRange[it]);
    GetHist("basic__raw__EnergyDeposit", "EnergyDeposit", "energy_deposit",
            "#N Tracks")
        ->Fill(reco.EnergyDeposit[it]);
  }

  // Adjust to kalman if needed
  if (!has_kalman) {
    for (int itrack = 0; itrack < reco.nTracks; itrack++) {
      for (int ihit = 0; ihit < reco.nHits[itrack]; ihit++) {
        for (int i = 0; i < 3; i++) {
          reco.KalmanPos[itrack][ihit][i] = reco.TrackHitPos[itrack][ihit][i];
          reco.KalmanTruePos[itrack][ihit][i] =
              truth.RecoTrackTrueHitPosition[itrack][ihit][i];
        }
      }
      reco.nKalmanNodes[itrack] = reco.nHits[itrack];
    }
  }
  // Fix n kalman nodes which is always zero.
  for (int itrack = 0; itrack < reco.nTracks; itrack++) {
    int n_nonzero = 0;
    for (int ihit = 0; ihit < reco.nHits[itrack]; ihit++) {
      // Since kalman pos need to be above z = 11000, we can use this to find
      // zero points safely
      if (reco.KalmanPos[itrack][ihit][2] > 1000)
        n_nonzero++;
    }
    reco.nKalmanNodes[itrack] = n_nonzero;

    /*if (reco.nHits[itrack] > 0) {
      const int LOOKBACK_WINDOW = 10;
      int n_nodes = n_nonzero > LOOKBACK_WINDOW ? LOOKBACK_WINDOW : n_nonzero;
      n_nodes -= 1;
      reco.StartDirection[itrack][0] = reco.KalmanPos[itrack][0][0] -
    reco.KalmanPos[itrack][n_nodes][0]; reco.StartDirection[itrack][1] =
    reco.KalmanPos[itrack][0][1] - reco.KalmanPos[itrack][n_nodes][1];
      reco.StartDirection[itrack][2] = reco.KalmanPos[itrack][0][2] -
    reco.KalmanPos[itrack][n_nodes][2];

      int last_index = reco.nKalmanNodes[itrack] - 1;
      reco.EndDirection[itrack][0] = reco.KalmanPos[itrack][last_index][0] -
    reco.KalmanPos[itrack][last_index - n_nodes][0];
      reco.EndDirection[itrack][1] = reco.KalmanPos[itrack][last_index][1] -
    reco.KalmanPos[itrack][last_index - n_nodes][1];
      reco.EndDirection[itrack][2] = reco.KalmanPos[itrack][last_index][2] -
    reco.KalmanPos[itrack][last_index - n_nodes][2];
    }*/
  }

  // "Fixed" variables
  for (int it = 0; it < reco.nTracks; it++) {
    GetHist("basic__fixed__nKalmanNodes", "nKalmanNodes", "n0-120", "#N Tracks")
        ->Fill(reco.nKalmanNodes[it]);
    for (int ih = 0; ih < reco.nKalmanNodes[it]; ih++) {
      GetHist("basic__fixed__KalmanPos_X", "KalmanPos X", "X", "#N Nodes")
          ->Fill(reco.KalmanPos[it][ih][0] * CM);
      GetHist("basic__fixed__KalmanPos_Y", "KalmanPos Y", "Y", "#N Nodes")
          ->Fill(reco.KalmanPos[it][ih][1] * CM);
      GetHist("basic__fixed__KalmanPos_Z", "KalmanPos Z", "Z", "#N Nodes")
          ->Fill(reco.KalmanPos[it][ih][2] * CM);
      GetHist("basic__fixed__KalmanTruePos_X", "KalmanTruePos X", "X",
              "#N Nodes")
          ->Fill(reco.KalmanTruePos[it][ih][0] * CM);
      GetHist("basic__fixed__KalmanTruePos_Y", "KalmanTruePos Y", "Y",
              "#N Nodes")
          ->Fill(reco.KalmanTruePos[it][ih][1] * CM);
      GetHist("basic__fixed__KalmanTruePos_Z", "KalmanTruePos Z", "Z",
              "#N Nodes")
          ->Fill(reco.KalmanTruePos[it][ih][2] * CM);
    }
    GetHist("basic__fixed__StartDirection_X", "StartDirection X", "dx",
            "#N Tracks")
        ->Fill(reco.StartDirection[it][0]);
    GetHist("basic__fixed__StartDirection_Y", "StartDirection Y", "dy",
            "#N Tracks")
        ->Fill(reco.StartDirection[it][1]);
    GetHist("basic__fixed__StartDirection_Z", "StartDirection Z", "dz",
            "#N Tracks")
        ->Fill(reco.StartDirection[it][2]);
    GetHist("basic__fixed__EndDirection_X", "EndDirection X", "dx", "#N Tracks")
        ->Fill(reco.EndDirection[it][0]);
    GetHist("basic__fixed__EndDirection_Y", "EndDirection Y", "dy", "#N Tracks")
        ->Fill(reco.EndDirection[it][1]);
    GetHist("basic__fixed__EndDirection_Z", "EndDirection Z", "dz", "#N Tracks")
        ->Fill(reco.EndDirection[it][2]);

    // Check for reco hits outside the TMS scint volume
    for (int ih = 0; ih < reco.nKalmanNodes[it]; ih++) {
      int hit_flag_x = -1;
      int hit_flag_y = -1;
      int hit_flag_z = -1;
      if (std::abs(reco.KalmanPos[it][ih][2]) < 0.01)
        hit_flag_z = 0;
      if (reco.KalmanPos[it][ih][2] < -990)
        hit_flag_z = 1;
      if (reco.KalmanPos[it][ih][2] < 11362)
        hit_flag_z = 2;
      if (reco.KalmanPos[it][ih][2] > 18313)
        hit_flag_z = 3;
      // if (hit_flag_z == 3) std::cout<<"reco.KalmanPos[it][ih][2]:
      // "<<reco.KalmanPos[it][ih][2]<<std::endl;
      if (reco.KalmanPos[it][ih][1] < -2950)
        hit_flag_y = 2;
      if (reco.KalmanPos[it][ih][1] > 240)
        hit_flag_y = 3;
      if (reco.KalmanPos[it][ih][0] < -3350)
        hit_flag_x = 2;
      if (reco.KalmanPos[it][ih][0] > 3350)
        hit_flag_x = 3;
      // Only check x and y == 0 if z is strange
      if (hit_flag_z != -1) {
        if (std::abs(reco.KalmanPos[it][ih][0]) < 0.01)
          hit_flag_x = 0;
        if (std::abs(reco.KalmanPos[it][ih][1]) < 0.01)
          hit_flag_y = 0;
      }
      if (hit_flag_z != -1)
        GetHist("basic__sanity__Unusual_Reco_Hit_Locations_Z",
                "Unusual Reco Hit Locations Z", "unusual_hit_locations",
                "#N Nodes")
            ->Fill(hit_flag_z);
      if (hit_flag_x != -1)
        GetHist("basic__sanity__Unusual_Reco_Hit_Locations_X",
                "Unusual Reco Hit Locations X", "unusual_hit_locations",
                "#N Nodes")
            ->Fill(hit_flag_x);
      if (hit_flag_y != -1)
        GetHist("basic__sanity__Unusual_Reco_Hit_Locations_Y",
                "Unusual Reco Hit Locations Y", "unusual_hit_locations",
                "#N Nodes")
            ->Fill(hit_flag_y);
    }
  }
  
  // Truth variables
  REGISTER_AXIS(HitEnergy, std::make_tuple("Hit Energy (MeV)", 51, 0, 100));
  REGISTER_AXIS(
      TrueVisibleEnergy,
      std::make_tuple("Particle True Visible Energy (MeV)", 51, 0, 1000));
  REGISTER_AXIS(HitHadronicEnergy,
                std::make_tuple("Hit Hadronic Energy (MeV)", 51, 0, 100));
  REGISTER_AXIS(HitEnergy_zoom, std::make_tuple("Hit Energy (MeV)", 20, 0, 1));
  REGISTER_AXIS(Hitdedx, std::make_tuple("Hit dEdx (MeV / cm)", 51, 0, 100));
  REGISTER_AXIS(TotalEnergy,
                std::make_tuple("Total Hit Energy (GeV)", 51, 0, 100));
  REGISTER_AXIS(TotalHadEnergy,
                std::make_tuple("Total Hit Hadronic Energy (GeV)", 51, 0, 100));
  REGISTER_AXIS(TotalRatio, std::make_tuple("Ratio E_{Had}/E", 51, 0, 1.0));
  REGISTER_AXIS(
      TotalEnergyPerVertex,
      std::make_tuple("Total Hit Energy Per Vertex (GeV)", 51, 0, 20));
  REGISTER_AXIS(
      TotalHadronicEnergyPerVertex,
      std::make_tuple("Total Hit Hadronic Energy Per Vertex (GeV)", 51, 0, 20));
  REGISTER_AXIS(ShellHadronicEnergyPerVertex,
                std::make_tuple("Hit Hadronic Energy Per Vertex In Shell (GeV)",
                                51, 0, 20));
  REGISTER_AXIS(ShellHadronicEnergyPerVertexZoom,
                std::make_tuple("Hit Hadronic Energy Per Vertex In Shell (MeV)",
                                20, 0, 200));
  if (on_new_spill) {
    GetHist("basic__truth__nTrueParticles", "nTrueParticles", "n0-500",
            "#N Spills")
        ->Fill(truth.nTrueParticles);
    GetHist("basic__truth__nTruePrimaryParticles", "nTruePrimaryParticles",
            "n0-500", "#N Spills")
        ->Fill(truth.nTruePrimaryParticles);
    GetHist("basic__truth__nTrueForgottenParticles", "nTrueForgottenParticles",
            "n0-120", "#N Spills")
        ->Fill(truth.nTrueForgottenParticles);

    std::map<int, bool> vertex_id_to_fiducial;
    for (int ip = 0; ip < truth.nTrueParticles; ip++) {
      int pdg = truth.PDG[ip];
      if (std::abs(pdg) != 13) {
        int vid = truth.VertexID[ip];
        TVector3 position(truth.BirthPosition[ip][0] * CM,
                          truth.BirthPosition[ip][1] * CM,
                          truth.BirthPosition[ip][2] * CM);
        bool result = IsInLAr(position, LArFiducial);
        vertex_id_to_fiducial[vid] = result;
      }
    }

    for (int ip = 0; ip < truth.nTrueParticles; ip++) {
      GetHist("basic__truth__PDG", "PDG", "pdg", "#N Particles")
          ->Fill(PDGtoIndex(truth.PDG[ip]));
      if (truth.IsPrimary[ip])
        GetHist("basic__truth__PDG_Primary", "PDG Primary Particles", "pdg",
                "#N Particles")
            ->Fill(PDGtoIndex(truth.PDG[ip]));
      if (!truth.IsPrimary[ip])
        GetHist("basic__truth__PDG_Secondary", "PDG Secondary Particles", "pdg",
                "#N Particles")
            ->Fill(PDGtoIndex(truth.PDG[ip]));
      GetHist("basic__truth__TrueVisibleEnergy", "TrueVisibleEnergy",
              "TrueVisibleEnergy", "#N Particles")
          ->Fill(truth.TrueVisibleEnergy[ip]);
      // TODO these are looping outside n reco tracks, right?
      GetHist("basic__truth__RecoTrackPrimaryParticleTruePositionStart_X",
              "RecoTrackPrimaryParticleTruePositionStart X", "X",
              "#N Particles")
          ->Fill(truth.RecoTrackPrimaryParticleTruePositionStart[ip][0] * CM);
      GetHist("basic__truth__RecoTrackPrimaryParticleTruePositionStart_Y",
              "RecoTrackPrimaryParticleTruePositionStart Y", "Y_full",
              "#N Particles")
          ->Fill(truth.RecoTrackPrimaryParticleTruePositionStart[ip][1] * CM);
      GetHist("basic__truth__RecoTrackPrimaryParticleTruePositionStart_Z",
              "RecoTrackPrimaryParticleTruePositionStart Z", "Z_full",
              "#N Particles")
          ->Fill(truth.RecoTrackPrimaryParticleTruePositionStart[ip][2] * CM);
      GetHist("basic__truth__RecoTrackPrimaryParticleTruePositionEnd_X",
              "RecoTrackPrimaryParticleTruePositionEnd X", "X", "#N Particles")
          ->Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[ip][0] * CM);
      GetHist("basic__truth__RecoTrackPrimaryParticleTruePositionEnd_Y",
              "RecoTrackPrimaryParticleTruePositionEnd Y", "Y_full",
              "#N Particles")
          ->Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[ip][1] * CM);
      GetHist("basic__truth__RecoTrackPrimaryParticleTruePositionEnd_Z",
              "RecoTrackPrimaryParticleTruePositionEnd Z", "Z_full",
              "#N Particles")
          ->Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[ip][2] * CM);

      GetHist("basic__truth__RecoTrackPrimaryParticleLArFiducialStart",
              "RecoTrackPrimaryParticleLArFiducialStart", "yesno",
              "#N Particles")
          ->Fill(truth.RecoTrackPrimaryParticleLArFiducialStart[ip] ? 0 : 1);
      GetHist("basic__truth__RecoTrackPrimaryParticleTMSFiducialEnd",
              "RecoTrackPrimaryParticleTMSFiducialEnd", "yesno", "#N Particles")
          ->Fill(truth.RecoTrackPrimaryParticleTMSFiducialEnd[ip] ? 0 : 1);
    }
  }
}
