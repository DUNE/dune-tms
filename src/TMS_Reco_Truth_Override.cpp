#include "TMS_Reco_Truth_Override.h"

#include "TMS_Hit.h"
#include "TMS_Kalman.h"
#include "TMS_Reco.h"
#include "TMS_Track.h"

std::vector<std::pair<int, int>>
TMS_Reco_Truth_Override::GetAllUniqueParticles(std::vector<TMS_Hit> &hits) {
  // Finds all the unique vertex/primary ids and returns them as a vector
  // They represent the unique particle indices
  std::map<std::pair<int, int>, double> counter;
  for (auto &hit : hits) {
    auto th = hit.GetTrueHit();
    for (size_t i = 0; i < th.GetNTrueParticles(); i++) {
      int vertex_id = th.GetVertexIds(i);
      int particle_id = th.GetPrimaryIds(i);
      double energy = th.GetEnergyShare(i);
      std::pair<int, int> key(vertex_id, particle_id);
      counter[key] += energy;
    }
  }

  // Now make a vector with unique pairs
  std::vector<std::pair<int, int>> out;
  for (auto it : counter) {
    out.push_back(it.first);
  }
  return out;
}

std::vector<std::vector<TMS_Hit>>
TMS_Reco_Truth_Override::GetSeparateHitsByParticle(
    std::vector<TMS_Hit> &hits, std::vector<std::pair<int, int>> particles) {
  std::map<std::pair<int, int>, std::vector<TMS_Hit>> map;

  // Finds all the unique vertex/primary ids and returns them as a vector
  // They represent the unique particle indices
  for (auto &hit : hits) {
    auto th = hit.GetTrueHit();
    for (size_t i = 0; i < th.GetNTrueParticles(); i++) {
      int vertex_id = th.GetVertexIds(i);
      int particle_id = th.GetPrimaryIds(i);
      std::pair<int, int> key(vertex_id, particle_id);
      // Only add hit if it's not there yet
      auto &vec = map[key];
      if (std::find(vec.begin(), vec.end(), hit) == vec.end())
        vec.push_back(hit);
    }
  }
  // Now make a vector with unique pairs
  std::vector<std::vector<TMS_Hit>> out;
  for (auto key : particles) {
    out.push_back(map[key]);
  }
  return out;
}

bool TMS_Reco_Truth_Override::CanKalmanProcess(
    std::vector<TMS_Hit> &hits_copy) {
  bool found_u = false;
  bool found_v = false;
  bool found_x = false;
  bool found_y = false;
  for (size_t i = 0; i < hits_copy.size(); i++) {
    if (hits_copy.at(i).GetBar().GetBarType() == TMS_Bar::kVBar)
      found_v = true;
    if (hits_copy.at(i).GetBar().GetBarType() == TMS_Bar::kUBar)
      found_u = true;
    if (hits_copy.at(i).GetBar().GetBarType() == TMS_Bar::kXBar)
      found_x = true;
    if (hits_copy.at(i).GetBar().GetBarType() == TMS_Bar::kYBar)
      found_y = true;
  }
  int n_found = 0;
  if (found_u)
    n_found += 1;
  if (found_v)
    n_found += 1;
  if (found_x)
    n_found += 1;
  if (found_y)
    n_found += 1;
  bool verbose = false;
  if (verbose)
    std::cout << "Found " << n_found << " bar types" << std::endl;
  if (verbose && found_u)
    std::cout << "Found u bar" << std::endl;
  if (verbose && found_v)
    std::cout << "Found v bar" << std::endl;
  if (verbose && found_x)
    std::cout << "Found x bar" << std::endl;
  if (verbose && found_y)
    std::cout << "Found y bar" << std::endl;
  bool can_process = n_found > 1 && hits_copy.size() > 3;
  if (verbose && can_process)
    std::cout << "Decided Kalman can process" << std::endl;
  else if (verbose)
    std::cout << "Decided Kalman can't process" << std::endl;
  return can_process;
}

TMS_Kalman TMS_Reco_Truth_Override::RunKalmanFilter(std::vector<TMS_Hit> hits,
                                                    double charge) {
  // Need separate copies for each since kalman filter doesn't keep hits cleanly
  // separate
  std::vector<TMS_Hit> hits_copy;
  for (auto hit : hits) {
    // std::cout<<"x: "<<hit.GetX()<<", reco x: "<<hit.GetRecoX()<<", reco x
    // uncertainty: "<<hit.GetRecoXUncertainty()<<std::endl;
    hits_copy.push_back(hit);
  }
  TMS_Kalman kalman(hits_copy, charge);
  return kalman;
}

TMS_Track TMS_Reco_Truth_Override::RunTrackReco(std::vector<TMS_Hit> hits) {
  TMS_Track trk;
  #define USE_THREE
  #ifdef USE_THREE
  TMS_Kalman kalman_neutral = RunKalmanFilter(hits, 0);
  const double charge_scale = 1;
  auto kalman_positive = RunKalmanFilter(hits, charge_scale);
  auto kalman_negative = RunKalmanFilter(hits, -charge_scale);

  double chi2_positive = kalman_positive.GetTrackChi2();
  double chi2_negative = kalman_negative.GetTrackChi2();

  TMS_Kalman &kalman_best_candidate = kalman_neutral;
  double chi2_best = kalman_best_candidate.GetTrackChi2();
  double charge = 0;
  if (chi2_positive < chi2_best) {
    kalman_best_candidate = kalman_positive;
    chi2_best = kalman_best_candidate.GetTrackChi2();
    charge = 13;
  }
  if (chi2_negative < chi2_best) {
    kalman_best_candidate = kalman_negative;
    chi2_best = kalman_best_candidate.GetTrackChi2();
    charge = -13;
  }

  // Set chi2 values
  trk.SetChi2(chi2_best);
  trk.SetChi2_plus(chi2_positive);
  trk.SetChi2_minus(chi2_negative);
  // Fill the KalmanNodes of the TMS_Track
  trk.KalmanNodes = kalman_best_candidate.GetKalmanNodes();
  trk.nKalmanNodes = trk.KalmanNodes.size();
  trk.SetMomentum(kalman_best_candidate.GetMomentum());
  trk.KalmanNodes_plus = kalman_positive.GetKalmanNodes();
  trk.KalmanNodes_minus = kalman_negative.GetKalmanNodes();
  #else
  const double MAX_CHARGE = 1.5;
  double current_charge = -MAX_CHARGE;
  double dc = 0.1;
  TMS_Kalman temp = RunKalmanFilter(hits, 0);
  TMS_Kalman temp_pos = RunKalmanFilter(hits, 1);
  TMS_Kalman temp_neg = RunKalmanFilter(hits, -1);
  TMS_Kalman &kalman_best_candidate = temp;
  TMS_Kalman &kalman_positive = temp_pos;
  TMS_Kalman &kalman_negative = temp_neg;
  
  double chi2_best = 1e12;
  double chi2_best_plus = 1e12;
  double chi2_best_minus = 1e12;
  double charge = -99999;
  do {
    auto kalman = RunKalmanFilter(hits, current_charge);
    double chi2 = kalman.GetTrackChi2();
    if (chi2 < chi2_best) {
      kalman_best_candidate = kalman;
      chi2_best = chi2;
      charge = current_charge * 10;
    }
    if (current_charge > 0 && chi2 < chi2_best_plus) { 
      kalman_positive = kalman;
      chi2_best_plus = chi2;
    }
    if (current_charge < 0 && chi2 < chi2_best_minus) {
      kalman_negative = kalman;
      chi2_best_minus = chi2;
    }
    current_charge += dc;
  } while (current_charge < MAX_CHARGE);
  // Set chi2 values
  trk.SetChi2(chi2_best);
  trk.SetChi2_plus(chi2_best_plus);
  trk.SetChi2_minus(chi2_best_minus);
  // Fill the KalmanNodes of the TMS_Track
  trk.KalmanNodes = kalman_best_candidate.GetKalmanNodes();
  trk.nKalmanNodes = trk.KalmanNodes.size();
  trk.SetMomentum(kalman_best_candidate.GetMomentum());
  trk.KalmanNodes_plus = kalman_positive.GetKalmanNodes();
  trk.KalmanNodes_minus = kalman_negative.GetKalmanNodes();
  #endif
  // Set track charge based on chi2
  trk.Charge_Kalman = charge;
  trk.Charge = charge;

  trk.Hits = hits;
  trk.nHits = trk.Hits.size();


  // TODO this is technically the back but we'll remain consistent with the
  // convention for now
  trk.SetStartPosition(kalman_best_candidate.Start[0],
                       kalman_best_candidate.Start[1],
                       kalman_best_candidate.Start[2]);
  // And this is the front of the track, the part closest to LAr
  trk.SetEndPosition(kalman_best_candidate.End[0], kalman_best_candidate.End[1],
                     kalman_best_candidate.End[2]);

  // Set directions
  trk.SetStartDirection(kalman_best_candidate.StartDirection[0],
                        kalman_best_candidate.StartDirection[1],
                        kalman_best_candidate.StartDirection[2]);
  trk.SetEndDirection(kalman_best_candidate.EndDirection[0],
                      kalman_best_candidate.EndDirection[1],
                      kalman_best_candidate.EndDirection[2]);

  trk.Length = TMS_TrackFinder::GetFinder().CalculateTrackLengthKalman(trk);
  //std::cout<<"trk.nKalmanNodes: "<<trk.nKalmanNodes<<std::endl;
  //std::cout<<"trk.Length: "<<trk.Length<<std::endl;
  // This is rough estimate for relationship
  trk.EnergyRange = 1.75 * trk.Length; 
  trk.EnergyDeposit = GetEnergyForOccupancy(hits);
  // Track time is min time
  trk.Time = 1e12;
  for (auto &hit : hits) {
    if (trk.Time > hit.GetT())
      trk.Time = hit.GetT();
  }

  return trk;
}

double
TMS_Reco_Truth_Override::GetEnergyForOccupancy(std::vector<TMS_Hit> &hits) {
  double total_energy = 0;
  // TODO do I want E or EVis?
  for (auto &hit : hits)
    total_energy += hit.GetEVis();
  return total_energy;
}

void TMS_Reco_Truth_Override::AdjustHitInformation(std::vector<TMS_Hit> &hits) {
  for (auto &hit : hits) {
    //  Use truth info for now
    hit.SetRecoX(hit.GetTrueHit().GetX());
    hit.SetRecoY(hit.GetTrueHit().GetY());
    double uncertainty_x = 50; // mm, 5cm for U,U,Y planes
    double uncertainty_y = 300; // mm, 30cm for U,U,Y planes, 5cm for X plane
    if (hit.GetBar().GetBarType() == TMS_Bar::kXBar) {
      uncertainty_y = 50;
      uncertainty_x = 300;
    }
    hit.SetRecoXUncertainty(uncertainty_x);
    hit.SetRecoYUncertainty(uncertainty_y);
  }
}

bool TMS_Reco_Truth_Override::FindTracksUsingTruth(
    std::vector<TMS_Hit> &hits, TMS_Event &event,
    std::vector<TMS_Track> &HoughTracks3D) {
  if (hits.size() > 100000) {
    std::cout << "Fatal: Expected reasonable hits. Got this many: "
              << hits.size() << std::endl;
    exit(0);
  }

  bool verbose = false;

  // Moves all the reco positions
  // Required for kalman
  AdjustHitInformation(hits);

  if (verbose)
    std::cout << "FindTracksUsingTruth: Processing this many hits: "
              << hits.size() << std::endl;
  // First get all the possible indices to look for
  auto indices = GetAllUniqueParticles(hits);
  // Now look up their truth info in the event (optional?)

  // Now split the hits by particle (how do I deal with overlap?)
  auto hits_by_particle = GetSeparateHitsByParticle(hits, indices);
  if (verbose)
    std::cout << "FindTracksUsingTruth: Separated hits into this many groups "
                 "by particle: "
              << hits_by_particle.size() << std::endl;

  double total_energy = GetEnergyForOccupancy(hits);
  if (verbose)
    std::cout << "Total energy:\t" << total_energy << "\t for this many hits:\t"
              << hits.size() << std::endl;

  // Now run each through the kalman filter
  // Put all the info in HoughTracks3D
  int i = 0;
  for (auto &particle_hits : hits_by_particle) {
    i += 1;
    if (particle_hits.size() > hits.size()) {
      std::cout << "Fatal: Found more hits for one particle than there are "
                   "hits in total. "
                << particle_hits.size() << " vs " << hits.size() << std::endl;
      exit(1);
    }
    //  TODO currently the kalman filter seems to have some crash conditions
    //  Like it's happy to process 4 hit u-only vectors but not 4 hit v-only
    //  vectors
    bool should_run_reco = CanKalmanProcess(particle_hits);
    if (should_run_reco) {
      auto track = RunTrackReco(particle_hits);
      track.Occupancy = track.EnergyDeposit / total_energy;
      if (verbose) {
        std::cout<<"Track number "<<i<<std::endl;
        track.Print();
        if (track.Length < 0) {
          std::cout<<"Fatal: Found length < 0: "<<track.Length<<std::endl;
          exit(0);
        }
      }
      HoughTracks3D.push_back(track);
    }
  }

  return true;
}
