// Add scope to avoid cross talk with other scripts
{
  REGISTER_AXIS(endpoint_dz,
                std::make_tuple("dZ (reco - true) (cm)", 31, -100.0, 50.0));
  REGISTER_AXIS(endpoint_dz_wide,
                std::make_tuple("dZ (reco - true) (cm)", 100, -400.0, 100.0));
  REGISTER_AXIS(endpoint_dx,
                std::make_tuple("dX (reco - true) (cm)", 31, -50.0, 50.0));
  REGISTER_AXIS(endpoint_dy,
                std::make_tuple("dY (reco - true) (cm)", 31, -50.0, 50.0));
  REGISTER_AXIS(
      endpoint_dboth,
      std::make_tuple("2D XZ resolution (reco - true) (cm)", 31, 0, 100.0));
  REGISTER_AXIS(Z_reco, std::make_tuple("Reco Z (cm)", 100, 1100, 1900));
  REGISTER_AXIS(Z_true, std::make_tuple("True Z (cm)", 100, 1100, 1900));
  for (int it = 0; it < reco.nTracks; it++) {
    bool has_contained_muon = false;
    bool end_tms = truth.RecoTrackPrimaryParticleTMSFiducialEnd[it];
    bool is_muon = abs(truth.RecoTrackPrimaryParticlePDG[it]) == 13;
    if (end_tms && is_muon)
      has_contained_muon = true;

    if (has_contained_muon) {
      double reco_startpoint_x = reco.StartPos[it][0];
      double reco_startpoint_y = reco.StartPos[it][1];
      double reco_startpoint_z = reco.StartPos[it][2];
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
      GetSpecialHist("special__track_resolution_endpoint_z",
                     "reco_track__track_resolution__endpoint_z",
                     "Endpoint Resolution dz", "endpoint_dz", "#N Tracks")
          ->Fill(endpoint_dz);
      GetSpecialHist("special__track_resolution_endpoint_x",
                     "reco_track__track_resolution__endpoint_x",
                     "Endpoint Resolution dx", "endpoint_dx", "#N Tracks")
          ->Fill(endpoint_dx);
      GetSpecialHist("special__track_resolution_endpoint_y",
                     "reco_track__track_resolution__endpoint_y",
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

      GetHist("reco_track__track_positions__track_endpoint_xz",
              "Track Endpoint XZ", "Z", "X")
          ->Fill(reco_endpoint_z * CM, reco_endpoint_x * CM);
      GetHist("reco_track__track_positions__track_endpoint_yz",
              "Track Endpoint YZ", "Z", "Y")
          ->Fill(reco_endpoint_z * CM, reco_endpoint_y * CM);
      GetHist("reco_track__track_positions__track_endpoint_yz",
              "Track Endpoint XY", "X", "Y")
          ->Fill(reco_endpoint_x * CM, reco_endpoint_y * CM);

      // First because end is start of kalman track
      int plane = reco.RecoTrackKalmanFirstPlaneBarView[it][0];
      int bar = reco.RecoTrackKalmanFirstPlaneBarView[it][1];
      int view = reco.RecoTrackKalmanFirstPlaneBarView[it][2];
      if (view == 0)
        GetHist("reco_track__track_positions__track_endpoint_plane_bar_view_x",
                "Track Endpoint Plane Bar X View", "plane", "bar")
            ->Fill(plane, bar);
      if (view == 1)
        GetHist("reco_track__track_positions__track_endpoint_plane_bar_view_y",
                "Track Endpoint Plane Bar Y View", "plane", "bar")
            ->Fill(plane, bar);
      if (view == 2)
        GetHist("reco_track__track_positions__track_endpoint_plane_bar_view_u",
                "Track Endpoint Plane Bar U View", "plane", "bar")
            ->Fill(plane, bar);
      if (view == 3)
        GetHist("reco_track__track_positions__track_endpoint_plane_bar_view_v",
                "Track Endpoint Plane Bar V View", "plane", "bar")
            ->Fill(plane, bar);

      GetHist("reco_track__track_positions__track_endpoint_z_vs_plane",
              "Track Endpoint Plane vs Z", "Z", "plane")
          ->Fill(reco_endpoint_z * CM, plane);
      int true_plane = reco.RecoTrackKalmanFirstPlaneBarViewTrue[it][0];
      GetHist("reco_track__track_positions__track_endpoint_plane_true_vs_reco",
              "Track Endpoint True vs Reco Plane", "true_plane", "plane")
          ->Fill(true_plane, plane);
      GetHist(
          "reco_track__track_positions__track_endpoint_plane_vs_z_resolution",
          "Track Endpoint True vs Reco Plane", "plane", "endpoint_dz_wide")
          ->Fill(plane, endpoint_dz);

      GetHist("reco_track__track_positions__track_endpoint_z_reco_vs_true",
              "Track Endpoint Z, Reco vs True", "Z_true", "Z_reco")
          ->Fill(true_endpoint_z * CM, reco_endpoint_z * CM);
      if (view == 0)
        GetHist(
            "reco_track__track_positions__track_endpoint_z_reco_vs_true_view_x",
            "Track Endpoint Z, Reco vs True, X View Endpoint", "Z_true",
            "Z_reco")
            ->Fill(true_endpoint_z * CM, reco_endpoint_z * CM);
      if (view == 1)
        GetHist(
            "reco_track__track_positions__track_endpoint_z_reco_vs_true_view_y",
            "Track Endpoint Z, Reco vs True, Y View Endpoint", "Z_true",
            "Z_reco")
            ->Fill(true_endpoint_z * CM, reco_endpoint_z * CM);
      if (view == 2)
        GetHist(
            "reco_track__track_positions__track_endpoint_z_reco_vs_true_view_u",
            "Track Endpoint Z, Reco vs True, U View Endpoint", "Z_true",
            "Z_reco")
            ->Fill(true_endpoint_z * CM, reco_endpoint_z * CM);
      if (view == 3)
        GetHist(
            "reco_track__track_positions__track_endpoint_z_reco_vs_true_view_v",
            "Track Endpoint Z, Reco vs True, V View Endpoint", "Z_true",
            "Z_reco")
            ->Fill(true_endpoint_z * CM, reco_endpoint_z * CM);

      GetHist("reco_track__track_positions__track_endpoint_z_vs_z_resolution",
              "Track Endpoint Z vs Track Z Resolution", "Z", "endpoint_dz_wide")
          ->Fill(reco_endpoint_z * CM, endpoint_dz);
      GetHist("reco_track__track_positions__track_endpoint_z_vs_xz_resolution",
              "Track Endpoint Z vs XZ Resolution", "Z", "endpoint_dboth")
          ->Fill(reco_endpoint_z * CM, endpoint_dxz);

      GetHist("reco_track__track_positions__track_startpoint_xz",
              "Track Startpoint XZ", "Z", "X")
          ->Fill(reco_startpoint_z * CM, reco_startpoint_x * CM);
      GetHist("reco_track__track_positions__track_startpoint_yz",
              "Track Startpoint YZ", "Z", "Y")
          ->Fill(reco_startpoint_z * CM, reco_startpoint_y * CM);
      GetHist("reco_track__track_positions__track_startpoint_yz",
              "Track Startpoint XY", "X", "Y")
          ->Fill(reco_startpoint_x * CM, reco_startpoint_y * CM);

      if (endpoint_dz > 10)
        DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                  "resolution/track_endpoint/z_resolution",
                  TString::Format("endpoint_dz = %.1fcm", endpoint_dz).Data(),
                  reco, lc, truth, DrawSliceN::few);
      if (endpoint_dz > 75)
        DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                  "resolution/track_endpoint/poor_z_resolution",
                  TString::Format("endpoint_dz = %.1fcm", endpoint_dz).Data(),
                  reco, lc, truth, DrawSliceN::few);
      if (endpoint_dz > 400)
        DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                  "resolution/track_endpoint/terrible_z_resolution",
                  TString::Format("endpoint_dz = %.1fcm", endpoint_dz).Data(),
                  reco, lc, truth, DrawSliceN::few);
      if (endpoint_dz > 800)
        DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                  "resolution/track_endpoint/horrific_z_resolution",
                  TString::Format("endpoint_dz = %.1fcm", endpoint_dz).Data(),
                  reco, lc, truth, DrawSliceN::few);
    }
  }
}
