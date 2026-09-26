// Fitted momentum vs truth, and track completeness, for muons that start in
// the ND-LAr fiducial volume (added 2026-09-26).
//
// For each such true muon, only its most complete track in the slice is used
// (the one holding most of the muon's visible energy), so a muon split into
// fragments contributes once, through its best piece. Momentum is the track's
// fitted momentum (reco.Momentum) against the magnitude of the true momentum
// entering the TMS, split into muons that stop in the TMS (contained) and
// ones that leave it (exiting). Completeness is the muon's visible energy on
// the track over its visible energy in the whole spill -- the counterpart of
// reco_track__cleanliness_energy (purity).
// Add scope to avoid cross talk with other scripts
{
  REGISTER_AXIS(momentum_resolution,
                std::make_tuple("(Fit - True) / True Momentum", 40, -1.0, 1.0));
  REGISTER_AXIS(p_true_tms_enter,
                std::make_tuple("True Muon Momentum Entering TMS (GeV/c)", 20, 0.0, 5.0));
  REGISTER_AXIS(p_fit, std::make_tuple("Fitted Momentum (GeV/c)", 25, 0.0, 5.0));
  REGISTER_AXIS(track_completeness,
                std::make_tuple("Muon Visible Energy on Track / in Spill", 22, 0.0, 1.1));

  // Best (most complete) track per true muon in this slice.
  std::map<int, int> best_track_of_particle;
  for (int it = 0; it < reco.nTracks; it++) {
    if (std::abs(truth.RecoTrackPrimaryParticlePDG[it]) != 13) continue;
    if (!truth.RecoTrackPrimaryParticleLArFiducialStart[it]) continue;
    const int ip = truth.RecoTrackPrimaryParticleIndex[it];
    if (ip < 0 || ip >= truth_spill.nTrueParticles) continue;
    auto found = best_track_of_particle.find(ip);
    if (found == best_track_of_particle.end() ||
        truth.RecoTrackPrimaryParticleTrueVisibleEnergy[it] >
            truth.RecoTrackPrimaryParticleTrueVisibleEnergy[found->second])
      best_track_of_particle[ip] = it;
  }

  for (const auto &particle_and_track : best_track_of_particle) {
    const int ip = particle_and_track.first;
    const int it = particle_and_track.second;
    const bool contained = truth.RecoTrackPrimaryParticleTMSFiducialEnd[it];
    const std::string sample = contained ? "contained" : "exiting";

    // Completeness (contained and exiting alike).
    if (truth_spill.TrueVisibleEnergy[ip] > 0)
      GetHist("momentum__completeness__" + sample, "Muon Track Completeness, " + sample,
              "track_completeness", "#N Muons")
          ->Fill(truth.RecoTrackPrimaryParticleTrueVisibleEnergy[it] / truth_spill.TrueVisibleEnergy[ip]);

    // Fitted momentum.
    const double px = truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][0];
    const double py = truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][1];
    const double pz = truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][2];
    const double p_true = std::sqrt(px * px + py * py + pz * pz);  // MeV
    const double p_fit = reco.Momentum[it];                         // MeV
    if (!(p_true > 0) || !(p_fit > 0) || !std::isfinite(p_fit)) continue;
    const double resolution = (p_fit - p_true) / p_true;
    GetHist("momentum__resolution__" + sample, "Fitted Momentum Resolution, " + sample,
            "momentum_resolution", "#N Muons")
        ->Fill(resolution);
    GetHist("momentum__resolution__" + sample + "_vs_p", "Fitted Momentum Resolution vs True, " + sample,
            "p_true_tms_enter", "momentum_resolution")
        ->Fill(p_true * GEV, resolution);
    GetHist("momentum__fit_vs_true__" + sample, "Fitted vs True Momentum, " + sample,
            "p_true_tms_enter", "p_fit")
        ->Fill(p_true * GEV, p_fit * GEV);
  }
}
