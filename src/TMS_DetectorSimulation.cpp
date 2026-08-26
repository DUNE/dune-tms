#include "TMS_DetectorSimulation.h"
#include "TMS_Readout_Manager.h"

void TMS_DetectorSimulation::SimulateOpticalModel(TMS_Event &event, std::default_random_engine &generator) {
  // Steps:
  // Loop over hits
  // Convert hit E -> PE
  // Do a poisson throw to get the number of PE
  // Apply some effect of PE capture into the fiber (assumed to be part of E -> PE conversion)
  // Split PE in two randomly for short path vs the long way.
  std::vector<TMS_Hit> &TMS_Hits = event.GetHitsRawRef();

  // TODO add second exponential term using fast decay length
  const double birks_constant = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_BirksConstant(); // mm / MeV
  const bool should_simulate_poisson_throws = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_ShouldSimulatePoisson();
  const bool should_simulate_fiber_lengths = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_ShouldSimulateFiberLengths();

  const double wsf_attenuation_length = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_WSFAttenuationLength(); // m
  // In reality, light bounces so there's a length multiplier
  const double wsf_length_multiplier = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_WSFLengthMultiplier();
  const double wsf_decay_constant = 1/wsf_attenuation_length;
  const double wsf_fiber_reflection_eff = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_WSFEndReflectionEff(); // How much light will reflect at the end
  const double fiber_coupling_eff = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_AdditionalFiberCouplingEff();
  const double optic_fiber_attenuation_length = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_AdditionalFiberAttenuationLength();
  const double optic_fiber_decay_constant = 1/optic_fiber_attenuation_length;
  const double optic_fiber_length = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_AdditionalFiberLength();

  const double readout_coupling_eff = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_ReadoutCouplingEff();

  for (auto& hit : TMS_Hits) {
    double pe = hit.GetPE();

    // Applies birk's suppression
    double de = hit.GetTrueHit().GetE();
    double dx = hit.GetTrueHit().GetdX();
    double dedx = 0;
    if (dx > 1e-8) dedx = de / dx;
    else dedx = de / 1.0;
    pe *= 1.0 / (1.0 + birks_constant * dedx);

    double pe_short = pe;
    double pe_long = 0;
    if (should_simulate_poisson_throws) {
      // Do a poisson throw to get the number of PE
      std::poisson_distribution<int> poisson(pe);
      pe = poisson(generator);
      // Now split the photons into the long and short paths with 50% chance of each
      std::binomial_distribution<int> binomial(pe, 0.5);
      pe_short = binomial(generator);
      pe_long = pe - pe_short;
    }

    if (should_simulate_fiber_lengths) {

      // Calculate the long and short path lengths
#ifdef USE_OLD_CODE
      double true_y = hit.GetTrueHit().GetY() / 1000.0; // m
      // In case of orthogonal (X) layers change to GetX()
      if (hit.GetBar().GetBarType() == TMS_Bar::kXBar) true_y = hit.GetTrueHit().GetX() / 1000.0;
      // assuming 0 is center, and assume we're reading out from top, then top would be biased negative and bottom positive, so -true_y.
      // TODO manually found this center. Make function in geom tools that returns values about scint
      // TODO fix math
      double distance_from_middle = TMS_Manager::GetInstance().Get_Geometry_YMIDDLE() - true_y;  // -1.54799
      double distance_from_end = distance_from_middle + 2;
      double long_way_distance_from_end = 4 + (4 - distance_from_end);
#else
      double distance_from_end = hit.GetTrueDistanceFromReadout() * 1e-3; // m
      double long_way_distance_from_end = hit.GetTrueLongDistanceFromReadout() * 1e-3; // m
#endif
      // In reality, light bounces so there's a multiplier
      // TODO it may be more realistic to make this non-linear
      distance_from_end *= wsf_length_multiplier;
      long_way_distance_from_end *= wsf_length_multiplier;

      // Now do exponential decay
      pe_short = pe_short * std::exp(-wsf_decay_constant * distance_from_end);
      pe_long = pe_long * std::exp(-wsf_decay_constant * long_way_distance_from_end) * wsf_fiber_reflection_eff;

      // Now possibly couple to a regular optical fiber
      if (optic_fiber_length > 0) {
        pe_short = fiber_coupling_eff * pe_short * std::exp(-optic_fiber_decay_constant * optic_fiber_length);
        pe_long = fiber_coupling_eff * pe_long * std::exp(-optic_fiber_decay_constant * optic_fiber_length);
      }
    }

    // Now couple between the fibers and the readout
    pe_long = pe_long * readout_coupling_eff;
    pe_short = pe_short * readout_coupling_eff;

    // Now save this information
    pe = pe_long + pe_short;

    // We want to save info right after fibers but before electronic conversion noise
    // This is particularly useful for timing information which cares about the first photon to be detected
    hit.GetAdjustableTrueHit().SetPEAfterFibers(pe);
    hit.GetAdjustableTrueHit().SetPEAfterFibersLongPath(pe_long);
    hit.GetAdjustableTrueHit().SetPEAfterFibersShortPath(pe_short);

    // Now save the reconstructed information
    hit.SetPE(pe);
    // Need to convert from PE to MeV. Could use 1/LY but have to account for additional effects.
    // so get constant to remove the effect of poisson, birks, fiber length, and readout noise above
    // The largest effect is fiber length
    double calibration_constant = TMS_Manager::GetInstance().Get_RECO_CALIBRATION_EnergyCalibration();
    double reco_e = pe * calibration_constant;
    hit.SetE(reco_e);
    hit.SetEVis(reco_e);
  }
}

void TMS_DetectorSimulation::SimulateDarkCount(TMS_Event &event) {
  (void)event;
  // TODO Add noise hits. They can fake readout
  // One issue is that there's no truth info about the particles to save.
}

void TMS_DetectorSimulation::SimulateTimingModel(TMS_Event &event, std::default_random_engine &generator) {
  // List of timing effects to simulate:
  // Random electronic timing noise
  // Deadtime
  // Time skew from first PE to hit sensor
  // Optical fiber length delays (corrected to strip center)
  // Timing effects from random noise, cross talk, afterpulsing
  // TODO check constants or put in config
  std::vector<TMS_Hit> &TMS_Hits = event.GetHitsRawRef();

  std::normal_distribution<double> noise_distribution(0.0, 1); // Mean of 0.0 and standard deviation of 1ns
  std::uniform_int_distribution<int> coin_flip(0, 1); // 0 or 1 depending on if you went long or short
  double scintillator_decay_time = 3.0; // ns
  double wsf_decay_time = 20.0; // ns
  std::exponential_distribution<double> exp_scint(1 / scintillator_decay_time); // Decay time = 3ns for scintillator
  std::exponential_distribution<double> exp_wsf(1 / wsf_decay_time); // 20ns for wavelength shifting fiber
  const double SPEED_OF_LIGHT =  0.2998; // m/ns
  const double FIBER_N = 1.5; //
  const double SPEED_OF_LIGHT_IN_FIBER = SPEED_OF_LIGHT / FIBER_N;

  const double wsf_length_multiplier = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_WSFLengthMultiplier();

  //double avg = 0;
  //double maxy = -1e9;
  //double miny = 1e9;
  //int n = 0;
  for (auto& hit : TMS_Hits) {
    double t = 0;
    // Random electronic timing noise (~1ns or less)
    t += noise_distribution(generator);
    // Optical fiber length delay (corrected to strip center)
    // (up to 13.4ns assuming 4m from edge, but correlated with y position. If delta y = 1m spread, than relative error is only 3.3ns)
#ifdef USE_OLD_CODE
    double true_y = hit.GetTrueHit().GetY() / 1000.0; // m
    // Making sure this gets changed for orthogonal (X) layers
    if (hit.GetBar().GetBarType() == TMS_Bar::kXBar) true_y = hit.GetTrueHit().GetX() / 1000.0;
    //miny = std::min(miny, true_y);
    //maxy = std::max(maxy, true_y);
    // assuming 0 is center, and assume we're reading out from top, then top would be biased negative and bottom positive, so -true_y.
    // TODO manually found this center. Want a better way in case things change
    double distance_from_middle = TMS_Manager::GetInstance().Get_Geometry_YMIDDLE() - true_y;  //-1.54799
    double long_way_distance = distance_from_middle + 8;
#else
    double distance_from_middle = hit.GetTrueDistanceFromMiddle() * 1e-3; // m
    double long_way_distance = hit.GetTrueLongDistanceFromMiddle() * 1e-3; // m
#endif
    // In reality, light bounces so there's a multiplier to the distance
    // todo, it may be more realistic to make this non-linear
    distance_from_middle *= wsf_length_multiplier;
    long_way_distance *= wsf_length_multiplier;

    // Find the time correction
    double time_correction = distance_from_middle / SPEED_OF_LIGHT_IN_FIBER;
    // This is the time correction if you go the long way instead
    double time_correction_long_way = long_way_distance / SPEED_OF_LIGHT_IN_FIBER;

    // Time slew (up to 30ns for 1pe hits, 9ns for 5pe, ~2ns 22pe. Typically 22pe mips assuming 45 pe mips with half going the long way)
    double pe_short_path = hit.GetTrueHit().GetPEAfterFibersShortPath();
    double pe_long_path = hit.GetTrueHit().GetPEAfterFibersLongPath();
    double minimum_time_offset = 1e100;

    #define USE_GAMMA_DISTRIBUTION
    #ifdef USE_GAMMA_DISTRIBUTION
    // TODO test this version with gamma

    // Gamma distribution does throws without having to do the throws
    // So it answer the question, what's the lowest of a random exponential assuming N throws.
    std::gamma_distribution<double> gamma_scint_short_path(pe_short_path, 1.0 / scintillator_decay_time);
    std::gamma_distribution<double> gamma_wsf_short_path(pe_short_path, 1.0 / wsf_decay_time);
    double minimum_time_gamma_scint_short_path = gamma_scint_short_path(generator);
    double minimum_time_gamma_wsf_short_path = gamma_wsf_short_path(generator);
    double minimum_time_offset_short_path = minimum_time_gamma_scint_short_path + minimum_time_gamma_wsf_short_path + time_correction;
    minimum_time_offset = std::min(minimum_time_offset_short_path, minimum_time_offset);

    std::gamma_distribution<double> gamma_scint_long_path(pe_long_path, 1.0 / scintillator_decay_time);
    std::gamma_distribution<double> gamma_wsf_long_path(pe_long_path, 1.0 / wsf_decay_time);
    double minimum_time_gamma_scint_long_path = gamma_scint_long_path(generator);
    double minimum_time_gamma_wsf_long_path = gamma_wsf_long_path(generator);
    double minimum_time_offset_long_path = minimum_time_gamma_scint_long_path + minimum_time_gamma_wsf_long_path + time_correction_long_way;
    minimum_time_offset = std::min(minimum_time_offset_long_path, minimum_time_offset);

    #else
    // We don't have to do 1000s of throws. The time will be very close to zero.
    // Assuming 1k PE, the mean time is ~0.02ns vs ~0.06ns for 300 PE.
    const double MAX_PE_THROWS = 300;
    if (pe_short_path > MAX_PE_THROWS) {
      pe_short_path = MAX_PE_THROWS;
    }
    while (pe_short_path > 0) {
      double time_offset = time_correction;
      time_offset += exp_scint(generator);
      time_offset += exp_wsf(generator);
      minimum_time_offset = std::min(time_offset, minimum_time_offset);
      pe_short_path -= 1;
    }
    if (pe_long_path > MAX_PE_THROWS) {
      pe_long_path = MAX_PE_THROWS;
    }
    while (pe_long_path > 0) {
      double time_offset = time_correction_long_way;
      time_offset += exp_scint(generator);
      time_offset += exp_wsf(generator);
      minimum_time_offset = std::min(time_offset, minimum_time_offset);
      pe_long_path -= 1;
    }
    #endif
    t += minimum_time_offset;
    double hit_time = hit.GetT();
    //std::cout<<"dt: "<<t<<", hit t: "<<hit_time<<", reco t: "<<hit_time + t<<", min t offset: "<<minimum_time_offset<<", t corr: "<<time_correction<<", dist from middle: "<<distance_from_middle<<", long way t corr: "<<time_correction_long_way<<", long way dist: "<<long_way_distance<<", hit pe: "<<hit.GetPE()<<std::endl;
    //std::cout<<"Hit time: "<<hit_time<<std::endl;
    //std::cout<<"Adjusted hit time: "<<hit_time + t<<std::endl;
    // Finally set the time
    hit.SetT(hit_time + t);
  }
  //avg /= n;
  //std::cout<<"Avg middle: "<<avg<<std::endl;
  /*std::cout<<"Max y: "<<maxy<<std::endl;
  std::cout<<"Min y: "<<miny<<std::endl;
  std::cout<<"Center y: "<<0.5*(maxy - miny)<<std::endl;*/
}

void TMS_DetectorSimulation::SimulateDeadtime(TMS_Event &event) {
  // Simulates readout windows, deadtime and zombie time.
  // |----  readout -----|-------deadtime-------{zombie time}]
  // The actual merging of hits in a readout window happens in TMS_SignalProcessing::MergeCoincidentHits.
  std::vector<TMS_Hit> &TMS_Hits = event.GetHitsRawRef();

  // How long a channel can read out
  double readout_time = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_ReadoutTime();;// TMS_Const::TMS_TimeThreshold; // ns
  // How long a channel is dead before it can read out again
  double deadtime = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_Deadtime();; // ns
  // Imagine a hit before the end of deadtime, but the system is close enough to resetting that you
  // channel starts recording energy. So when readout is ready, it looks like a hit hit at t=readout_ready_time.
  // That's zombie time.
  double zombie_time = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_ZombieTime();;

  const bool deadtime_verbose = false;


  if (deadtime > 0) {
    // Want sorted hits by T so we can find the first hit in a channel to be the start of readout windows and deadtime.
    std::sort(TMS_Hits.begin(), TMS_Hits.end(), TMS_Hit::SortByT);

    int n_dead_hits = 0;
    int n_zombie_hits = 0;

    // These store the end time for each window by the per-channel id.
    // If there's no matching id, then we haven't seen that id yet and this hit can be the start of a readout window.
    std::map<int, double> readout_map;
    std::map<int, double> deadtime_map;
    std::map<int, double> zombie_map;
    std::map<int, bool> has_zombie_map;
    std::map<int, double> x_map;
    std::map<int, double> z_map;
    std::map<int, double> t_map;
    //for (auto& hit : TMS_Hits) {
    for (size_t i = 0; i < TMS_Hits.size(); ++i) {
      auto& hit = TMS_Hits[i];
      double t = hit.GetT();
      // For a per-channel deadtime, this is a unique id for a channel.
      // But some detectors have deadtime for a whole board, in which case this should
      // return a single id for the whole board.
      const int id = hit.GetNotZ() + 100000 * hit.GetZ(); // TODO make sure it's unique
      auto it_read = readout_map.find(id);
      auto it_dead = deadtime_map.find(id);
      auto it_zombie = zombie_map.find(id);
      bool should_reset_times = false;
      if (it_read == readout_map.end()) {
        // This id hasn't been seen before. We can read with no issue
        should_reset_times = true;
        if (deadtime_verbose) std::cout<<"New channel -> Do regular read"<<std::endl;
      }
      else {
        // We have seen this channel. So next we need to see if we're in the readout time or deadtime
        double t_read = it_read->second;
        double t_dead = it_dead->second;
        double t_zombie = it_zombie->second;
        if (deadtime_verbose) std::cout<<"Found channel we found already with Notz: "<<hit.GetNotZ()<<", z: "<<hit.GetZ()<<", id: "<<id<<"\n";
        if (deadtime_verbose) std::cout<<"Compare with previous channel Notz: "<<x_map[id]<<", z: "<<z_map[id]<<", t: "<<t_map[id]<<"\n";
        if (x_map[id] != hit.GetNotZ() || z_map[id] != hit.GetZ()) std::cout<<"\n** Found mismatch in Notz,z **\n"<<std::endl;
        if (deadtime_verbose) std::cout<<"i="<<i<<", t="<<t<<", t_read="<<t_read<<", t_dead="<<t_dead<<", t_zombie="<<t_zombie<<", dt="<<(t-t_read+readout_time)<<", dt_map: "<<(t-t_map[id])<<std::endl;
        if (t < t_read) {
          // We can do a regular read
          if (deadtime_verbose) std::cout<<"t < t_read -> Do regular read"<<std::endl;
        }
        else if (zombie_time > 0 && t < t_zombie) {
          // Zombie time sets the time to the start of the next readout which is the end of deadtime
          hit.SetT(t_dead);
          n_zombie_hits += 1;
          has_zombie_map[id] = true;
          if (deadtime_verbose) std::cout<<"t < t_zombie -> Treat as zombie"<<std::endl;
        }
        else if (t < t_dead) {
          // Suppress channel since it's dead
          hit.SetPedSup(true);
          n_dead_hits += 1;
          if (deadtime_verbose) std::cout<<"t < t_dead -> Kill channel"<<std::endl;
        }
        else {
          // This channel is past the deadtime so reset our windows
          if (deadtime_verbose) std::cout<<"t >= t_dead -> Regular read and reset windows"<<std::endl;
          should_reset_times = true;
        }
      }

      // Calculate the times that this id can be read or is dead or is zombie
      if (should_reset_times) {
        // If there's a zombie time, that hit should be the start of the next readout
        // And its time will be the end of the current deadtime
        // But we don't need to do this check if this hit is outside the deadtime window of the next hit
        // So calculate that window first
        // But also this hit now needs to be checked against the zombie times's windows, so redo this hit
        // Only relevant if we've seen this channel id before -- for a brand new channel, it_dead is
        // deadtime_map.end() (there's no prior zombie state to re-check against), so skip entirely
        // rather than dereferencing end().
        if (it_dead != deadtime_map.end() && has_zombie_map[id] == true) {
          double deadtime_window_starting_from_end_of_deadtime = it_dead->second + deadtime;
          if (t < deadtime_window_starting_from_end_of_deadtime) {
            t = deadtime_map[id];
            has_zombie_map[id] = false;
            // Need to redo this hit to check that it isn't in the deadtime of the zombie hit
            i--;
          }
        }
        // Since this is the first hit for an id, the end of the readout window is t + readout_time
        double t_read = t + readout_time;
        readout_map[id] = t_read;
        // Deadtime starts at the end of the readout period
        double t_dead = t_read + deadtime;
        deadtime_map[id] = t_dead;
        // Zombie time is the time before the end of deadtime
        // So 20ns of zombie time in a 500ns deadtime would start at 480ns.
        double t_zombie = t_dead - zombie_time;
        zombie_map[id] = t_zombie;

        x_map[id] = hit.GetNotZ();
        z_map[id] = hit.GetZ();
        t_map[id] = t;

        auto position = std::make_pair(hit.GetNotZ(), hit.GetZ());
        auto deadtime_range = std::make_pair(t_read, t_dead);
        auto readout_range = std::make_pair(t, t_read);
        event.AddDeadtimeChannelRecord(position, deadtime_range, readout_range);

        #ifdef RECORD_HIT_DEADTIME
        // In theory merging hits should capture correct deadtimes based on the first deadtime recorded here
        // But if doing things by sipm or some larger group, then we're not recording all the info
        hit.SetDeadtimeStart(t_read);
        hit.SetDeadtimeStop(t_dead);
        #endif
      }

      // TODO remove
      //if (i > 100) exit(0);

    }
    double n_dead_hits_as_percent = 100.0 * n_dead_hits / (float)TMS_Hits.size();
    std::cout<<"N dead hits: "<<n_dead_hits<<" out of "<<TMS_Hits.size()<<" hits. That's "<<n_dead_hits_as_percent<<"%"<<std::endl;
    if (zombie_time > 0) std::cout<<"N zombie hits: "<<n_zombie_hits<<std::endl;
  }
}

void TMS_DetectorSimulation::SimulateReadoutNoise(TMS_Event &event, std::default_random_engine &generator) {
  // Only want to simulate the little bit of electronic noise from reading out after merging hits
  // Otherwise we're adding together a bunch of random numbers centered around zero, leading to an average of zero
  std::vector<TMS_Hit> &TMS_Hits = event.GetHitsRawRef();

  const double readout_noise = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_ReadoutNoise();
  if (readout_noise > 0) {
    for (auto& hit : TMS_Hits) {
      double pe = hit.GetPE();
      // Skip 0 pe hits, they'll be removed in ped sup step
      if (pe > 0) {
        double E = hit.GetE();
        // Now take into account electronic readout noise
        std::normal_distribution<double> normal(pe, readout_noise);
        double newpe = normal(generator);
        double newE = E * newpe / pe;
        hit.SetPE(newpe);
        hit.SetE(newE);
      }
    }
  }
}
