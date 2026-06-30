// Add scope to avoid cross talk with other scripts
{
  REGISTER_AXIS(ke_tms, std::make_tuple("True Muon KE (GeV)", 28, 0.0, 7.0));
  REGISTER_AXIS(
      pion_ke_tms_enter,
      std::make_tuple("True Pion KE Entering TMS (GeV)", 20, 0.0, 5.0));
  REGISTER_AXIS(muon_endpoint_z, std::make_tuple("True Muon Endpoint Z (cm)",
                                                 20, 1136.0, 1831.0));
  REGISTER_AXIS(muon_endpoint_y,
                std::make_tuple("True Muon Endpoint Y (cm)", 20, -216.0, 16.0));
  REGISTER_AXIS(muon_endpoint_x, std::make_tuple("True Muon Endpoint X (cm)",
                                                 20, -330.0, 330.0));
  REGISTER_AXIS(
      muon_startpoint_z,
      std::make_tuple("True Muon Startpoint Z (cm)", 20, 1036.0, 1831.0));
  REGISTER_AXIS(
      muon_startpoint_y,
      std::make_tuple("True Muon Startpoint Y (cm)", 20, -216.0, 16.0));
  REGISTER_AXIS(
      muon_startpoint_x,
      std::make_tuple("True Muon Startpoint X (cm)", 20, -330.0, 330.0));
  REGISTER_AXIS(
      reco_eff_visible_energy,
      std::make_tuple("True Visible Energy in TMS (MeV)", 100, -50.0, 2000.0));
  REGISTER_AXIS(reco_eff_track_count_delta,
                std::make_tuple("Truth_Info RecoTrackN - Reco_Tree nTracks",
                                41, -20.5, 20.5));
  REGISTER_AXIS(reco_eff_spill_count,
                std::make_tuple("Count per spill", 101, -0.5, 100.5));
  REGISTER_AXIS(reco_eff_denominator_cut,
                std::make_tuple("Denominator cut", 9, 0.5, 9.5));
  REGISTER_AXIS(reco_eff_denominator_fail,
                std::make_tuple("Denominator fail reason", 7, 0.5, 7.5));
  REGISTER_AXIS(reco_eff_numerator_cut,
                std::make_tuple("Numerator cut", 9, 0.5, 9.5));
  REGISTER_AXIS(reco_eff_numerator_fail,
                std::make_tuple("Numerator fail reason", 9, 0.5, 9.5));

  static const std::vector<std::string> denominator_cut_labels = {
      "all truth",      "muon PDG",       "TMS visible",
      "KE in bounds",   "LAr branch",     "LAr position",
      "TMS end",        "branch strict",  "position strict"};
  static const std::vector<std::string> denominator_fail_labels = {
      "non-muon PDG",   "no TMS visible", "KE out of bounds",
      "no LAr branch",  "no LAr position", "no TMS end",
      "LAr disagreement"};
  static const std::vector<std::string> numerator_cut_labels = {
      "all reco",       "muon PDG",       "TMS visible",
      "KE in bounds",   "not duplicate",  "LAr branch",
      "LAr position",   "TMS end",        "position strict"};
  static const std::vector<std::string> numerator_fail_labels = {
      "RecoTrackN != nTracks", "non-muon PDG",   "no TMS visible",
      "KE out of bounds",     "duplicate",      "no LAr branch",
      "no LAr position",      "no TMS end",     "LAr disagreement"};

  auto configure_labeled_axis = [](TH1 *hist,
                                   const std::vector<std::string> &labels) {
    hist->SetStats(0);
    hist->SetNdivisions(labels.size());
    for (size_t i = 0; i < labels.size(); ++i)
      hist->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    hist->LabelsOption("v", "X");
  };
  auto fill_labeled_bin = [&](const std::string &name, const std::string &title,
                              const std::string &axis,
                              const std::vector<std::string> &labels,
                              int bin, const std::string &ytitle) {
    TH1 *hist = GetHist(name, title, axis, ytitle);
    configure_labeled_axis(hist, labels);
    hist->Fill(bin);
  };
  auto fill_bool = [](const std::string &name, const std::string &title,
                      bool value) {
    GetHist(name, title, "falsetrue", "# entries")->Fill(value ? 1.0 : 0.0);
  };
  auto is_lar_start_position = [](const Float_t position[4]) {
    const double lar_start[] = {-3478.48, -2166.71, 4179.24};
    const double lar_end[] = {3478.48, 829.282, 9135.88};
    return position[0] >= lar_start[0] && position[0] <= lar_end[0] &&
           position[1] >= lar_start[1] && position[1] <= lar_end[1] &&
           position[2] >= lar_start[2] && position[2] <= lar_end[2];
  };
  auto quantize_reco_truth = [](float value) -> long long {
    return static_cast<long long>(std::llround(value * 1000.0));
  };
  auto reco_truth_fingerprint = [&](int track) -> std::string {
    std::string key = std::to_string(truth.RecoTrackPrimaryParticlePDG[track]);
    for (int dim = 0; dim < 4; ++dim)
      key += ":" + std::to_string(quantize_reco_truth(
                       truth.RecoTrackPrimaryParticleTrueMomentum[track][dim]));
    for (int dim = 0; dim < 4; ++dim)
      key += ":" + std::to_string(quantize_reco_truth(
                       truth.RecoTrackPrimaryParticleTruePositionStart[track][dim]));
    for (int dim = 0; dim < 4; ++dim)
      key += ":" + std::to_string(quantize_reco_truth(
                       truth.RecoTrackPrimaryParticleTruePositionEnd[track][dim]));
    return key;
  };

  // Per-spill debug counters for reco-eff fills. These are intentionally
  // separate from the histograms so we can print the actual tree/spill/entry
  // range when numerator and denominator bookkeeping disagree.
  static bool have_previous_spill = false;
  static long spill_counter = 0;
  static int previous_tree = -1;
  static int previous_spill = -1;
  static Long64_t previous_first_entry = -1;
  static Long64_t previous_last_entry = -1;
  static int cnt_spill_rows = 0;
  static int cnt_mu_den = 0, cnt_mu_num = 0;
  static int cnt_pi_den = 0, cnt_pi_num = 0;

  if (on_new_spill) {
    if (have_previous_spill) {
      GetHist("reco_eff__diagnostics__per_spill_mu_denominator",
              "Reco-eff muon denominators per spill", "reco_eff_spill_count")
          ->Fill(cnt_mu_den);
      GetHist("reco_eff__diagnostics__per_spill_mu_numerator",
              "Reco-eff muon numerators per spill", "reco_eff_spill_count")
          ->Fill(cnt_mu_num);
      GetHist("reco_eff__diagnostics__per_spill_rows",
              "Tracking_Validation rows per spill", "reco_eff_spill_count")
          ->Fill(cnt_spill_rows);
      bool mu_inconsistent =
          (cnt_mu_num > cnt_mu_den) || (cnt_mu_den == 0 && cnt_mu_num > 0);
      bool pi_inconsistent =
          (cnt_pi_num > cnt_pi_den) || (cnt_pi_den == 0 && cnt_pi_num > 0);
      if (mu_inconsistent || pi_inconsistent) {
        std::cout << "[Reco_Eff] Potential inconsistency in tree "
                  << previous_tree << ", SpillNo " << previous_spill
                  << " (spill index " << spill_counter
                  << ", entries " << previous_first_entry << "-"
                  << previous_last_entry << ", rows=" << cnt_spill_rows
                  << "): mu_den=" << cnt_mu_den
                  << ", mu_num=" << cnt_mu_num
                  << ", pi_den=" << cnt_pi_den
                  << ", pi_num=" << cnt_pi_num << std::endl;
      }
    }

    spill_counter++;
    have_previous_spill = true;
    previous_tree = current_tree;
    previous_spill = current_spill;
    previous_first_entry = entry_number;
    previous_last_entry = entry_number;
    cnt_spill_rows = 0;
    cnt_mu_den = cnt_mu_num = cnt_pi_den = cnt_pi_num = 0;

    const int n_true_particles = std::max(
        0, std::min(truth.nTrueParticles, Truth_Info::kMaxTrueParticles));
    for (int ip = 0; ip < n_true_particles; ip++) {
      fill_labeled_bin("reco_eff__diagnostics__denominator_cutflow",
                       "Reco-eff denominator cutflow",
                       "reco_eff_denominator_cut", denominator_cut_labels, 1,
                       "# truth particles");
      GetHist("reco_eff__diagnostics__den_pdg",
              "Reco-eff denominator truth-particle PDG", "pdg",
              "# truth particles")
          ->Fill(PDGtoIndex(truth.PDG[ip]));
      const bool ismuon = abs(truth.PDG[ip]) == 13;
      const bool ispion = abs(truth.PDG[ip]) == 211;
      if (ismuon)
        fill_labeled_bin("reco_eff__diagnostics__denominator_cutflow",
                         "Reco-eff denominator cutflow",
                         "reco_eff_denominator_cut", denominator_cut_labels, 2,
                         "# truth particles");
      else
        fill_labeled_bin("reco_eff__diagnostics__denominator_fail_reason",
                         "Reco-eff denominator failed cuts",
                         "reco_eff_denominator_fail", denominator_fail_labels,
                         1, "# truth particles");

      const bool tms_touch = truth.TrueVisibleEnergy[ip] >= MINIMUM_VISIBLE_ENERGY;
      if (!tms_touch) {
        if (ismuon || ispion)
          fill_labeled_bin("reco_eff__diagnostics__denominator_fail_reason",
                           "Reco-eff denominator failed cuts",
                           "reco_eff_denominator_fail", denominator_fail_labels,
                           2, "# truth particles");
        continue;
      }
      if (ismuon)
        fill_labeled_bin("reco_eff__diagnostics__denominator_cutflow",
                         "Reco-eff denominator cutflow",
                         "reco_eff_denominator_cut", denominator_cut_labels, 3,
                         "# truth particles");

      const bool lar_start_branch = truth.LArFiducialStart[ip];
      const bool lar_start = is_lar_start_position(truth.BirthPosition[ip]);
      const bool tms_end = truth.TMSFiducialEnd[ip];
      const double particle_starting_ke = truth.MomentumTMSStart[ip][3] * GEV;
      const double particle_ke = truth.BirthMomentum[ip][3] * GEV;
      const bool ke_in_bounds = particle_starting_ke >= muon_ke_bins[0] &&
                                particle_starting_ke <
                                    muon_ke_bins[n_muon_ke_bins];
      bool should_include = true;
      if (only_nd_physics)
        should_include = tms_end && lar_start;

      if (ismuon) {
        GetHist("reco_eff__diagnostics__den_muon_visible_energy",
                "Denominator muon visible energy", "reco_eff_visible_energy",
                "# muons")
            ->Fill(truth.TrueVisibleEnergy[ip]);
        GetHist("reco_eff__diagnostics__den_muon_ke_tms_enter",
                "Denominator muon KE entering TMS", "ke_tms_enter")
            ->Fill(particle_starting_ke);
        fill_bool("reco_eff__diagnostics__den_lar_start_branch",
                  "Denominator LArFiducialStart branch", lar_start_branch);
        fill_bool("reco_eff__diagnostics__den_lar_start_position",
                  "Denominator true start in LAr box", lar_start);
        fill_bool("reco_eff__diagnostics__den_tms_end",
                  "Denominator TMSFiducialEnd branch", tms_end);
        fill_bool("reco_eff__diagnostics__den_lar_branch_matches_position",
                  "Denominator LAr branch agrees with position",
                  lar_start_branch == lar_start);
        if (ke_in_bounds)
          fill_labeled_bin("reco_eff__diagnostics__denominator_cutflow",
                           "Reco-eff denominator cutflow",
                           "reco_eff_denominator_cut", denominator_cut_labels,
                           4, "# truth particles");
        else
          fill_labeled_bin("reco_eff__diagnostics__denominator_fail_reason",
                           "Reco-eff denominator failed cuts",
                           "reco_eff_denominator_fail", denominator_fail_labels,
                           3, "# truth particles");
        if (lar_start_branch)
          fill_labeled_bin("reco_eff__diagnostics__denominator_cutflow",
                           "Reco-eff denominator cutflow",
                           "reco_eff_denominator_cut", denominator_cut_labels,
                           5, "# truth particles");
        else
          fill_labeled_bin("reco_eff__diagnostics__denominator_fail_reason",
                           "Reco-eff denominator failed cuts",
                           "reco_eff_denominator_fail", denominator_fail_labels,
                           4, "# truth particles");
        if (lar_start)
          fill_labeled_bin("reco_eff__diagnostics__denominator_cutflow",
                           "Reco-eff denominator cutflow",
                           "reco_eff_denominator_cut", denominator_cut_labels,
                           6, "# truth particles");
        else
          fill_labeled_bin("reco_eff__diagnostics__denominator_fail_reason",
                           "Reco-eff denominator failed cuts",
                           "reco_eff_denominator_fail", denominator_fail_labels,
                           5, "# truth particles");
        if (tms_end)
          fill_labeled_bin("reco_eff__diagnostics__denominator_cutflow",
                           "Reco-eff denominator cutflow",
                           "reco_eff_denominator_cut", denominator_cut_labels,
                           7, "# truth particles");
        else
          fill_labeled_bin("reco_eff__diagnostics__denominator_fail_reason",
                           "Reco-eff denominator failed cuts",
                           "reco_eff_denominator_fail", denominator_fail_labels,
                           6, "# truth particles");
        if (lar_start_branch != lar_start)
          fill_labeled_bin("reco_eff__diagnostics__denominator_fail_reason",
                           "Reco-eff denominator failed cuts",
                           "reco_eff_denominator_fail", denominator_fail_labels,
                           7, "# truth particles");
        if (lar_start_branch && tms_end)
          fill_labeled_bin("reco_eff__diagnostics__denominator_cutflow",
                           "Reco-eff denominator cutflow",
                           "reco_eff_denominator_cut", denominator_cut_labels,
                           8, "# truth particles");
        if (lar_start && tms_end)
          fill_labeled_bin("reco_eff__diagnostics__denominator_cutflow",
                           "Reco-eff denominator cutflow",
                           "reco_eff_denominator_cut", denominator_cut_labels,
                           9, "# truth particles");
      }

      if (ismuon && should_include) {
        GetHist("reco_eff__no_lar_tms_cuts__all_muon_ke_tms_enter_denominator",
                "Reconstruction Efficiency: Denominator", "ke_tms_enter")
            ->Fill(particle_starting_ke);
        GetHist("reco_eff__all_muon_ke_tms_enter_denominator",
                "Reco Efficiency vs True KE, All Muons: Denominator",
                "ke_tms_enter")
            ->Fill(particle_starting_ke);
        cnt_mu_den++;
        GetHist(
            "reco_eff__all_muon_ke_tms_enter_including_doubles_denominator",
            "Reco Efficiency w/2x vs True KE, All Muons: Denominator",
            "ke_tms_enter")
            ->Fill(particle_starting_ke);
        GetHist("reco_eff__all_muon_ke_denominator",
                "Reco Efficiency vs True KE, All Muons: Denominator", "ke_tms")
            ->Fill(particle_ke);

        GetHist("reco_eff__endpoint__muon_endpoint_x_denominator",
                "Reconstruction Efficiency: Denominator", "muon_endpoint_x")
            ->Fill(truth.PositionTMSEnd[ip][0] * CM);
        GetHist("reco_eff__endpoint__muon_endpoint_y_denominator",
                "Reconstruction Efficiency: Denominator", "muon_endpoint_y")
            ->Fill(truth.PositionTMSEnd[ip][1] * CM);
        GetHist("reco_eff__endpoint__muon_endpoint_z_denominator",
                "Reconstruction Efficiency: Denominator", "muon_endpoint_z")
            ->Fill(truth.PositionTMSEnd[ip][2] * CM);

        GetHist("reco_eff__startpoint__muon_startpoint_x_denominator",
                "Reconstruction Efficiency: Denominator", "muon_startpoint_x")
            ->Fill(truth.PositionTMSStart[ip][0] * CM);
        GetHist("reco_eff__startpoint__muon_startpoint_y_denominator",
                "Reconstruction Efficiency: Denominator", "muon_startpoint_y")
            ->Fill(truth.PositionTMSStart[ip][1] * CM);
        GetHist("reco_eff__startpoint__muon_startpoint_z_denominator",
                "Reconstruction Efficiency: Denominator", "muon_startpoint_z")
            ->Fill(truth.PositionTMSStart[ip][2] * CM);
      }
      if (ispion && should_include) {
        GetHist("reco_eff__no_lar_tms_cuts__all_pion_ke_tms_enter_denominator",
                "Reconstruction Efficiency: Denominator", "pion_ke_tms_enter")
            ->Fill(particle_starting_ke);
        cnt_pi_den++;
      }
      if (ismuon && lar_start && tms_end && should_include) {
        GetSpecialHist("special__reco_eff_muon_ke_tms_enter_denominator",
                       "reco_eff__muon_ke_tms_enter_denominator",
                       "Reconstruction Efficiency: Denominator",
                       "ke_tms_enter")
            ->Fill(particle_starting_ke);
        GetHist("reco_eff__good_reco_muon_ke_tms_enter_denominator",
                "Reconstruction Efficiency, with Good Reco: Denominator",
                "ke_tms_enter")
            ->Fill(particle_starting_ke);
        GetHist("reco_eff__muon_ke_denominator",
                "Reco Efficiency vs True KE: Denominator", "ke_tms")
            ->Fill(particle_ke);
        double particle_starting_angle =
            std::atan2(truth.MomentumTMSStart[ip][0],
                       truth.MomentumTMSStart[ip][2]) *
            DEG;
        GetHist("reco_eff__angle_tms_enter_denominator",
                "Reconstruction Efficiency: Denominator", "angle_tms_enter")
            ->Fill(particle_starting_angle);
      }
    }
  }

  if (have_previous_spill) {
    cnt_spill_rows++;
    previous_last_entry = entry_number;
  }

  const int reco_track_delta = truth.RecoTrackN - reco.nTracks;
  GetHist("reco_eff__diagnostics__recotrackn_minus_ntracks",
          "Truth_Info RecoTrackN - Reco_Tree nTracks",
          "reco_eff_track_count_delta", "# entries")
      ->Fill(reco_track_delta);
  if (truth.RecoTrackN != reco.nTracks)
    fill_labeled_bin("reco_eff__diagnostics__numerator_fail_reason",
                     "Reco-eff numerator failed cuts",
                     "reco_eff_numerator_fail", numerator_fail_labels, 1,
                     "# reco tracks");

  int n_reco_eff_tracks = truth.RecoTrackN;
  if (n_reco_eff_tracks > reco.nTracks)
    n_reco_eff_tracks = reco.nTracks;
  if (n_reco_eff_tracks > Truth_Info::kMaxRecoTracks)
    n_reco_eff_tracks = Truth_Info::kMaxRecoTracks;
  if (n_reco_eff_tracks < 0)
    n_reco_eff_tracks = 0;

  for (int it = 0; it < n_reco_eff_tracks; it++) {
    fill_labeled_bin("reco_eff__diagnostics__numerator_cutflow",
                     "Reco-eff numerator cutflow", "reco_eff_numerator_cut",
                     numerator_cut_labels, 1, "# reco tracks");
    GetHist("reco_eff__diagnostics__num_pdg",
            "Reco-eff numerator primary-particle PDG", "pdg", "# reco tracks")
        ->Fill(PDGtoIndex(truth.RecoTrackPrimaryParticlePDG[it]));
    const bool tms_touch = truth.RecoTrackPrimaryParticleTrueVisibleEnergy[it] >=
                           MINIMUM_VISIBLE_ENERGY;
    const bool ismuon = abs(truth.RecoTrackPrimaryParticlePDG[it]) == 13;
    const bool ispion = abs(truth.RecoTrackPrimaryParticlePDG[it]) == 211;
    bool not_double_reco = true;
    const bool lar_start_branch =
        truth.RecoTrackPrimaryParticleLArFiducialStart[it];
    const bool lar_start =
        is_lar_start_position(truth.RecoTrackPrimaryParticleTruePositionStart[it]);
    const bool tms_end = truth.RecoTrackPrimaryParticleTMSFiducialEnd[it];
    const double particle_starting_ke =
        truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][3] * GEV;
    const double particle_ke =
        truth.RecoTrackPrimaryParticleTrueMomentum[it][3] * GEV;
    const bool ke_in_bounds = particle_starting_ke >= muon_ke_bins[0] &&
                              particle_starting_ke <
                                  muon_ke_bins[n_muon_ke_bins];

    if (ismuon)
      fill_labeled_bin("reco_eff__diagnostics__numerator_cutflow",
                       "Reco-eff numerator cutflow", "reco_eff_numerator_cut",
                       numerator_cut_labels, 2, "# reco tracks");
    else
      fill_labeled_bin("reco_eff__diagnostics__numerator_fail_reason",
                       "Reco-eff numerator failed cuts",
                       "reco_eff_numerator_fail", numerator_fail_labels, 2,
                       "# reco tracks");

    if (ismuon) {
      GetHist("reco_eff__diagnostics__num_muon_visible_energy",
              "Numerator muon visible energy", "reco_eff_visible_energy",
              "# muons")
          ->Fill(truth.RecoTrackPrimaryParticleTrueVisibleEnergy[it]);
      GetHist("reco_eff__diagnostics__num_muon_ke_tms_enter",
              "Numerator muon KE entering TMS", "ke_tms_enter")
          ->Fill(particle_starting_ke);
      fill_bool("reco_eff__diagnostics__num_lar_start_branch",
                "Numerator LArFiducialStart branch", lar_start_branch);
      fill_bool("reco_eff__diagnostics__num_lar_start_position",
                "Numerator true start in LAr box", lar_start);
      fill_bool("reco_eff__diagnostics__num_tms_end",
                "Numerator TMSFiducialEnd branch", tms_end);
      fill_bool("reco_eff__diagnostics__num_lar_branch_matches_position",
                "Numerator LAr branch agrees with position",
                lar_start_branch == lar_start);
      if (ke_in_bounds)
        fill_labeled_bin("reco_eff__diagnostics__numerator_cutflow",
                         "Reco-eff numerator cutflow",
                         "reco_eff_numerator_cut", numerator_cut_labels, 4,
                         "# reco tracks");
      else
        fill_labeled_bin("reco_eff__diagnostics__numerator_fail_reason",
                         "Reco-eff numerator failed cuts",
                         "reco_eff_numerator_fail", numerator_fail_labels, 4,
                         "# reco tracks");
      if (lar_start_branch)
        fill_labeled_bin("reco_eff__diagnostics__numerator_cutflow",
                         "Reco-eff numerator cutflow",
                         "reco_eff_numerator_cut", numerator_cut_labels, 6,
                         "# reco tracks");
      else
        fill_labeled_bin("reco_eff__diagnostics__numerator_fail_reason",
                         "Reco-eff numerator failed cuts",
                         "reco_eff_numerator_fail", numerator_fail_labels, 6,
                         "# reco tracks");
      if (lar_start)
        fill_labeled_bin("reco_eff__diagnostics__numerator_cutflow",
                         "Reco-eff numerator cutflow",
                         "reco_eff_numerator_cut", numerator_cut_labels, 7,
                         "# reco tracks");
      else
        fill_labeled_bin("reco_eff__diagnostics__numerator_fail_reason",
                         "Reco-eff numerator failed cuts",
                         "reco_eff_numerator_fail", numerator_fail_labels, 7,
                         "# reco tracks");
      if (tms_end)
        fill_labeled_bin("reco_eff__diagnostics__numerator_cutflow",
                         "Reco-eff numerator cutflow",
                         "reco_eff_numerator_cut", numerator_cut_labels, 8,
                         "# reco tracks");
      else
        fill_labeled_bin("reco_eff__diagnostics__numerator_fail_reason",
                         "Reco-eff numerator failed cuts",
                         "reco_eff_numerator_fail", numerator_fail_labels, 8,
                         "# reco tracks");
      if (lar_start_branch != lar_start)
        fill_labeled_bin("reco_eff__diagnostics__numerator_fail_reason",
                         "Reco-eff numerator failed cuts",
                         "reco_eff_numerator_fail", numerator_fail_labels, 9,
                         "# reco tracks");
    }

    if (!tms_touch) {
      if (ismuon || ispion)
        fill_labeled_bin("reco_eff__diagnostics__numerator_fail_reason",
                         "Reco-eff numerator failed cuts",
                         "reco_eff_numerator_fail", numerator_fail_labels, 3,
                         "# reco tracks");
      continue;
    }
    if (ismuon)
      fill_labeled_bin("reco_eff__diagnostics__numerator_cutflow",
                       "Reco-eff numerator cutflow", "reco_eff_numerator_cut",
                       numerator_cut_labels, 3, "# reco tracks");

    const std::string particle_fingerprint = reco_truth_fingerprint(it);
    particle_fingerprints_reconstructed[particle_fingerprint]++;
    if (particle_fingerprints_reconstructed[particle_fingerprint] > 1) {
      // Found this apparent truth particle already at least once
      not_double_reco = false;
      fill_labeled_bin("reco_eff__diagnostics__numerator_fail_reason",
                       "Reco-eff numerator failed cuts",
                       "reco_eff_numerator_fail", numerator_fail_labels, 5,
                       "# reco tracks");
      if (ismuon) {
        GetHist("reco_eff__multi_reco__muons",
                "Muons which were reconstructed more than once",
                "ke_tms_enter")
            ->Fill(particle_starting_ke);
        GetHist("reco_eff__multi_reco__probability_multi_reco_numerator",
                "Chance of Getting Reco'd more than Once", "ke_tms_enter")
            ->Fill(particle_starting_ke);
      } else
        GetHist("reco_eff__multi_reco__nonmuon",
                "Non-muons which were reconstructed more than once", "ke_tms")
            ->Fill(particle_starting_ke);
      if (ismuon)
        DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                  "reco_eff/multi_reco_muon",
                  TString::Format("Particle fingerprint reco'd %dx",
                                  particle_fingerprints_reconstructed
                                      [particle_fingerprint])
                      .Data(),
                  reco, lc, truth, DrawSliceN::many);
      else
        DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                  "reco_eff/multi_reco_nonmuon",
                  TString::Format("Particle fingerprint reco'd %dx",
                                  particle_fingerprints_reconstructed
                                      [particle_fingerprint])
                      .Data(),
                  reco, lc, truth, DrawSliceN::many);
    }
    if (not_double_reco)
      fill_labeled_bin("reco_eff__diagnostics__numerator_cutflow",
                       "Reco-eff numerator cutflow", "reco_eff_numerator_cut",
                       numerator_cut_labels, 5, "# reco tracks");

    if (ismuon) {
      GetHist("reco_eff__all_muon_ke_tms_enter_including_doubles_numerator",
              "Reco Efficiency w/2x vs True TMS-Entering KE, All Muons: "
              "Numerator",
              "ke_tms_enter")
          ->Fill(particle_starting_ke);
    }
    if (ismuon && not_double_reco) {
      GetHist("reco_eff__all_muon_ke_tms_enter_numerator",
              "Reco Efficiency vs True TMS-Entering KE, All Muons: Numerator",
              "ke_tms_enter")
          ->Fill(particle_starting_ke);
      cnt_mu_num++;
      GetHist("reco_eff__all_muon_ke_numerator",
              "Reco Efficiency vs True KE, All Muons: Numerator", "ke_tms")
          ->Fill(particle_ke);
      GetHist("reco_eff__multi_reco__probability_multi_reco_denominator",
              "Chance of Getting Reco'd more than Once", "ke_tms_enter")
          ->Fill(particle_starting_ke);

      GetHist("reco_eff__endpoint__muon_endpoint_x_numerator",
              "Reconstruction Efficiency: Numerator", "muon_endpoint_x")
          ->Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[it][0] * CM);
      GetHist("reco_eff__endpoint__muon_endpoint_y_numerator",
              "Reconstruction Efficiency: Numerator", "muon_endpoint_y")
          ->Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[it][1] * CM);
      GetHist("reco_eff__endpoint__muon_endpoint_z_numerator",
              "Reconstruction Efficiency: Numerator", "muon_endpoint_z")
          ->Fill(truth.RecoTrackPrimaryParticleTruePositionEnd[it][2] * CM);

      GetHist("reco_eff__startpoint__muon_startpoint_x_numerator",
              "Reconstruction Efficiency: Numerator", "muon_startpoint_x")
          ->Fill(truth.RecoTrackPrimaryParticleTruePositionEnteringTMS[it][0] *
                 CM);
      GetHist("reco_eff__startpoint__muon_startpoint_y_numerator",
              "Reconstruction Efficiency: Numerator", "muon_startpoint_y")
          ->Fill(truth.RecoTrackPrimaryParticleTruePositionEnteringTMS[it][1] *
                 CM);
      GetHist("reco_eff__startpoint__muon_startpoint_z_numerator",
              "Reconstruction Efficiency: Numerator", "muon_startpoint_z")
          ->Fill(truth.RecoTrackPrimaryParticleTruePositionEnteringTMS[it][2] *
                 CM);

      if (particle_starting_ke < 1)
        DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                  "reco_eff/reco_low_e_muon",
                  TString::Format("Low E muon: %.1f GeV", particle_starting_ke)
                      .Data(),
                  reco, lc, truth, DrawSliceN::few);
      if (particle_starting_ke < 0.25)
        DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                  "reco_eff/reco_very_low_e_muon",
                  TString::Format("Low E muon: %.2f GeV", particle_starting_ke)
                      .Data(),
                  reco, lc, truth, DrawSliceN::few);
    }
    if (ispion && not_double_reco) {
      GetHist("reco_eff__no_lar_tms_cuts__all_pion_ke_tms_enter_numerator",
              "Reco Efficiency, All Pions: Numerator", "pion_ke_tms_enter")
          ->Fill(particle_starting_ke);
      cnt_pi_num++;
      DrawSlice(TString::Format("entry_%lld", entry_number).Data(),
                "reco_eff/reco_pion",
                TString::Format("Pion with %.1f GeV", particle_starting_ke)
                    .Data(),
                reco, lc, truth, DrawSliceN::few);
    }
    if (ismuon && lar_start && tms_end && not_double_reco) {
      fill_labeled_bin("reco_eff__diagnostics__numerator_cutflow",
                       "Reco-eff numerator cutflow", "reco_eff_numerator_cut",
                       numerator_cut_labels, 9, "# reco tracks");
      GetSpecialHist("special__reco_eff_muon_ke_tms_enter_numerator",
                     "reco_eff__muon_ke_tms_enter_numerator",
                     "Reco Efficiency vs True TMS-Entering KE: Numerator",
                     "ke_tms_enter")
          ->Fill(particle_starting_ke);
      int charge = reco.Charge[it];
      int true_charge = (truth.RecoTrackPrimaryParticlePDG[it] < 0) ? -1 : 1;
      bool correct_charge_id = true_charge * charge > 0;
      /*double true_muon_starting_ke =
        truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][3] * 1e-3;
      double reco_ke = length_to_energy(reco.Length[it]);
      bool correct_energy = std::abs((reco_ke - true_muon_starting_ke) /
      true_muon_starting_ke) < 0.05;

      GetHist("reco_eff__debugging__energy_reco", "Reco Energy",
              "ke_tms_enter")
          ->Fill(reco_ke);
      GetHist("reco_eff__debugging__energy_true", "Reco Energy",
              "ke_tms_enter")
          ->Fill(true_muon_starting_ke);
      REGISTER_AXIS(fractional_e_resolution,
                    std::make_tuple("Fractional E", 41, -1.0, 1.0));
      GetHist("reco_eff__debugging__energy_resolution_fraction",
              "Fractional E Resolution", "fractional_e_resolution")
          ->Fill((reco_ke - true_muon_starting_ke) / true_muon_starting_ke);*/
      double true_z =
          truth.RecoTrackPrimaryParticleTruePositionEnd[it][2] * CM;
      double reco_z = reco.EndPos[it][2] * CM;
      double dz = reco_z - true_z;
      bool correct_endpoint = -30 < dz && dz < 0;
      bool is_good_reco = correct_charge_id && correct_endpoint;
      if (is_good_reco)
        GetHist("reco_eff__good_reco_muon_ke_tms_enter_numerator",
                "Reco Efficiency vs True TMS-Entering KE: Numerator",
                "ke_tms_enter")
            ->Fill(particle_starting_ke);
      GetHist("reco_eff__muon_ke_numerator",
              "Reco Efficiency vs True KE, TMS ending: Numerator", "ke_tms")
          ->Fill(particle_ke);
      double particle_starting_angle =
          std::atan2(
              truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][0],
              truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][2]) *
          DEG;
      GetHist("reco_eff__angle_tms_enter_numerator",
              "Reco Efficiency vs True TMS-Entering Angle: Numerator",
              "angle_tms_enter")
          ->Fill(particle_starting_angle);
    }
  }
}
