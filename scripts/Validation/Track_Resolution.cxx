// Add scope to avoid cross talk with other scripts
{
  REGISTER_AXIS(endpoint_dz,
                std::make_tuple("dZ (reco - true) (cm)", 31, -100.0, 50.0));
  REGISTER_AXIS(endpoint_dx,
                std::make_tuple("dX (reco - true) (cm)", 31, -50.0, 50.0));
  REGISTER_AXIS(endpoint_dy,
                std::make_tuple("dY (reco - true) (cm)", 31, -50.0, 50.0));
  REGISTER_AXIS(
      endpoint_dboth,
      std::make_tuple("2D XZ resolution (reco - true) (cm)", 31, 0, 100.0));
  for (int it = 0; it < reco.nTracks; it++) {
    bool has_contained_muon = false;
    bool end_tms = truth.RecoTrackPrimaryParticleTMSFiducialEnd[it];
    bool is_muon = abs(truth.RecoTrackPrimaryParticlePDG[it]) == 13;
    if (end_tms && is_muon)
      has_contained_muon = true;

    if (has_contained_muon) {
      double true_endpoint_x =
          truth.RecoTrackPrimaryParticleTruePositionEnd[it][0];
      double reco_endpoint_x = reco.EndPos[it][0];
      double true_endpoint_y =
          truth.RecoTrackPrimaryParticleTruePositionEnd[it][1];
      double reco_endpoint_y = reco.EndPos[it][1];
      double true_endpoint_z =
          truth.RecoTrackPrimaryParticleTruePositionEnd[it][2];
      double reco_endpoint_z = reco.EndPos[it][2];
      double endpoint_dz = (reco_endpoint_z - true_endpoint_z) * CM;
      double endpoint_dy = (reco_endpoint_y - true_endpoint_y) * CM;
      double endpoint_dx = (reco_endpoint_x - true_endpoint_x) * CM;
      double endpoint_dxz =
          std::sqrt(endpoint_dz * endpoint_dz + endpoint_dx * endpoint_dx);
      GetSpecialHist("special__track_resolution_endpoint_z", "reco_track__track_resolution__endpoint_z",
              "Endpoint Resolution dz", "endpoint_dz", "#N Tracks")
          ->Fill(endpoint_dz);
      GetSpecialHist("special__track_resolution_endpoint_x", "reco_track__track_resolution__endpoint_x",
              "Endpoint Resolution dx", "endpoint_dx", "#N Tracks")
          ->Fill(endpoint_dx);
      GetSpecialHist("special__track_resolution_endpoint_y", "reco_track__track_resolution__endpoint_y",
              "Endpoint Resolution dy", "endpoint_dy", "#N Tracks")
          ->Fill(endpoint_dy);
      GetHist("reco_track__track_resolution__endpoint_xz",
              "Endpoint Resolution 2d", "endpoint_dboth")
          ->Fill(endpoint_dxz);

      TVector3 direction(reco.EndDirection[it][0], 0, reco.EndDirection[it][2]);
      // TVector3 direction(0, 0, 1);
      if (direction.Mag() > 0)
        direction.SetMag(30);
      else
        direction.SetZ(30);
      GetHist("reco_track__track_resolution__corrected_endpoint_z",
              "Endpoint Resolution dz", "endpoint_dz")
          ->Fill(endpoint_dz + direction.Z());
    }
  }
}
