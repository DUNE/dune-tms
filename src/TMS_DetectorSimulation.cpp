#include "TMS_DetectorSimulation.h"
#include "TMS_Readout_Manager.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <tuple>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
constexpr double SPEED_OF_LIGHT = 0.2998; // m/ns

// Relocated from TMS_Hit.cpp (Phase III) -- these need the true hit position, which is no
// longer embedded in TMS_Hit. Only ever called from this file, so kept file-local rather than
// added to the TMS_DetectorSimulation public interface.

// Distance along the bar from the point (x, y) to the bar's readout end. Position-based so the
// response-element path can evaluate it per optical deposit; the true-hit wrappers below are
// the original functions.
double DistanceFromReadout(const TMS_Hit &hit, double x, double y) {
  const double barLength = hit.GetBar().GetBarLength();
  const double barCenter = hit.GetBar().GetAxisReadoutCenter();
  // Note that you want to do always do more positive - less positive, or else you get a sign error
  if (hit.GetBar().GetBarType() == TMS_Bar::kXBar) {
    // Readout from sides
    if (x < 0) return x - TMS_Geom::GetInstance().XBarNegReadoutLocation(barCenter, barLength);
    else return TMS_Geom::GetInstance().XBarPosReadoutLocation(barCenter, barLength) - x;
  }
  else {
    // Readout from top. Assuming U ~ V ~ Y for now
    return TMS_Geom::GetInstance().YBarReadoutLocation(barCenter, barLength) - y;
  }
}

double LongDistanceFromReadout(const TMS_Hit &hit, double x, double y) {
  const double barLength = hit.GetBar().GetBarLength();
  double additional_length;
  if (hit.GetBar().GetBarType() == TMS_Bar::kXBar) {
    // Readout from sides
    additional_length = 2 * TMS_Geom::GetInstance().XBarLength(barLength);
  }
  else {
    // Readout from top. Assuming U ~ V ~ Y for now
    additional_length = 2 * TMS_Geom::GetInstance().YBarLength(barLength);
  }
  return additional_length - DistanceFromReadout(hit, x, y);
}

// Distance from a bar's own readout end to its own geometric center is half its length,
// for both bar types -- X-bars and Y-bars are both single-ended readout with a reflecting
// far end (X-bars split into two mirror-image halves at the detector's central gap, Y-bars
// not split), so both are defined symmetrically about their own center (readout locations
// are barCenter +/- 0.5*barLength in XBarPosReadoutLocation/XBarNegReadoutLocation/
// YBarReadoutLocation). Written directly rather than routed through XBarLength()/
// YBarLength() since those two are now identical (both just return barLength).
double DistanceFromMiddle(const TMS_Hit &hit, double x, double y) {
  return DistanceFromReadout(hit, x, y) - 0.5 * hit.GetBar().GetBarLength();
}

// Same bar-type-independent center offset as DistanceFromMiddle() above.
double LongDistanceFromMiddle(const TMS_Hit &hit, double x, double y) {
  return LongDistanceFromReadout(hit, x, y) - 0.5 * hit.GetBar().GetBarLength();
}

double GetTrueDistanceFromReadout(const TMS_Hit &hit, const TMS_TrueHit &true_hit) {
  return DistanceFromReadout(hit, true_hit.GetX(), true_hit.GetY());
}

double GetTrueLongDistanceFromReadout(const TMS_Hit &hit, const TMS_TrueHit &true_hit) {
  return LongDistanceFromReadout(hit, true_hit.GetX(), true_hit.GetY());
}

double GetTrueDistanceFromMiddle(const TMS_Hit &hit, const TMS_TrueHit &true_hit) {
  return DistanceFromMiddle(hit, true_hit.GetX(), true_hit.GetY());
}

double GetTrueLongDistanceFromMiddle(const TMS_Hit &hit, const TMS_TrueHit &true_hit) {
  return LongDistanceFromMiddle(hit, true_hit.GetX(), true_hit.GetY());
}
} // namespace

void TMS_DetectorSimulation::SimulateOpticalModel(TMS_Event &event, std::default_random_engine &generator) {
  // Steps:
  // Loop over hits
  // Convert hit E -> PE
  // Apply some effect of PE capture into the fiber (assumed to be part of E -> PE conversion)
  // Default: do a poisson throw to get the number of PE, split it in two randomly for the
  // short path vs the long way, then attenuate each path.
  // Sim.Optical.PoissonAfterAttenuation: attenuate the expected PE of each path first, then
  // do a poisson throw per path.
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
  const bool poisson_after_attenuation = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_PoissonAfterAttenuation();
  const bool use_response_elements = TMS_Readout_Manager::GetInstance().Get_Sim_DetSim_UseResponseElements();

  // Path efficiencies for light produced at a given (length-multiplied) distance from the
  // readout: WLS attenuation (and end reflection for the long way), optional additional optical
  // fiber, then fiber-to-readout coupling. Applied as a sequence of multiplications so the
  // default (draw-then-attenuate) path reproduces the previous arithmetic exactly.
  auto attenuate_short = [&](double x, double distance) {
    if (should_simulate_fiber_lengths) {
      // Now do exponential decay
      x = x * std::exp(-wsf_decay_constant * distance);
      // Now possibly couple to a regular optical fiber
      if (optic_fiber_length > 0) x = fiber_coupling_eff * x * std::exp(-optic_fiber_decay_constant * optic_fiber_length);
    }
    // Now couple between the fibers and the readout
    return x * readout_coupling_eff;
  };
  auto attenuate_long = [&](double x, double distance) {
    if (should_simulate_fiber_lengths) {
      x = x * std::exp(-wsf_decay_constant * distance) * wsf_fiber_reflection_eff;
      if (optic_fiber_length > 0) x = fiber_coupling_eff * x * std::exp(-optic_fiber_decay_constant * optic_fiber_length);
    }
    return x * readout_coupling_eff;
  };
  // std::poisson_distribution requires a positive mean
  auto poisson_draw = [&](double mean) {
    if (mean <= 0) return 0.0;
    std::poisson_distribution<int> poisson(mean);
    return static_cast<double>(poisson(generator));
  };

  // Response-element path: light is simulated per optical deposit and each detected photon
  // gets an explicit sensor arrival time (used by SimulateTimingModel()).
  const double light_yield = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_LightYield();
  const double deposit_bin_length = TMS_Readout_Manager::GetInstance().Get_Sim_DetSim_DepositBinLength();
  const double passage_max_gap = TMS_Readout_Manager::GetInstance().Get_Sim_DetSim_PassageMaxGap();
  const double speed_of_light_in_fiber = SPEED_OF_LIGHT / TMS_Readout_Manager::GetInstance().Get_Sim_Timing_FiberRefractiveIndex();
  std::exponential_distribution<double> exp_scint(1 / TMS_Readout_Manager::GetInstance().Get_Sim_Timing_ScintillatorDecayTime());
  std::exponential_distribution<double> exp_wsf(1 / TMS_Readout_Manager::GetInstance().Get_Sim_Timing_WLSDecayTime());
  std::uniform_real_distribution<double> uniform_share(0.0, 1.0);

  for (auto& hit : TMS_Hits) {
    if (use_response_elements) {
      const std::vector<TMS_Passage::Segment>* segments = event.GetResponseSegments(hit.GetHitId());
      if (segments == nullptr) throw std::runtime_error("Fatal: SimulateOpticalModel() found a hit with no recorded steps while Sim.DetSim.UseResponseElements is on");
      double pe_short = 0;
      double pe_long = 0;
      // Vertex of each contributing trajectory, for per-photon provenance
      std::map<int, long long> vertex_of_trajectory;
      if (const TMS_TrueHit* truth = event.GetTrueHit(hit.GetHitId())) {
        for (size_t i = 0; i < truth->GetNTrueParticles(); ++i) vertex_of_trajectory[truth->GetPrimaryIds(i)] = truth->GetVertexGlobalIds(i);
      }
      for (const auto& passage_indices : TMS_Passage::BuildPassages(*segments, passage_max_gap)) {
        std::vector<TMS_Passage::Segment> passage;
        for (size_t i : passage_indices) passage.push_back((*segments)[i]);
        for (const auto& deposit : TMS_Passage::Resegment(passage, deposit_bin_length)) {
          // Local Birks suppression from the deposit's own dE/dx
          const double dedx = (deposit.dx > 1e-8) ? deposit.energy / deposit.dx : deposit.energy / 1.0;
          const double pe_produced = deposit.energy * light_yield / (1.0 + birks_constant * dedx);
          const double x = deposit.position[0];
          const double y = deposit.position[1];
          double distance = 0, long_distance = 0, middle = 0, long_middle = 0;
          if (should_simulate_fiber_lengths) {
            distance = DistanceFromReadout(hit, x, y) * 1e-3 * wsf_length_multiplier; // m
            long_distance = LongDistanceFromReadout(hit, x, y) * 1e-3 * wsf_length_multiplier;
          }
          middle = DistanceFromMiddle(hit, x, y) * 1e-3 * wsf_length_multiplier;
          long_middle = LongDistanceFromMiddle(hit, x, y) * 1e-3 * wsf_length_multiplier;
          // Half the light goes each way; detected photons per path are Poisson with the fully
          // attenuated mean (as Sim.Optical.PoissonAfterAttenuation)
          const double mean_short = attenuate_short(0.5 * pe_produced, distance);
          const double mean_long = attenuate_long(0.5 * pe_produced, long_distance);
          const double n_short = should_simulate_poisson_throws ? poisson_draw(mean_short) : mean_short;
          const double n_long = should_simulate_poisson_throws ? poisson_draw(mean_long) : mean_long;
          pe_short += n_short;
          pe_long += n_long;
          // Particle that produced each photon: drawn in proportion to the deposit's energy shares
          // (no draw when there is a single contributor)
          auto photon_source = [&]() {
            int trajectory = deposit.shares.empty() ? -999 : deposit.shares.front().first;
            if (deposit.shares.size() > 1) {
              double pick = uniform_share(generator) * deposit.energy;
              for (const auto& share : deposit.shares) {
                trajectory = share.first;
                pick -= share.second;
                if (pick <= 0) break;
              }
            }
            const auto it = vertex_of_trajectory.find(trajectory);
            return std::make_pair(trajectory, it == vertex_of_trajectory.end() ? -999LL : it->second);
          };
          // Arrival time at the sensor, corrected to the strip center as in SimulateTimingModel()
          const int photons_short = static_cast<int>(std::ceil(n_short));
          const int photons_long = static_cast<int>(std::ceil(n_long));
          for (int i = 0; i < photons_short; ++i) {
            const double t = deposit.t + middle / speed_of_light_in_fiber + exp_scint(generator) + exp_wsf(generator);
            const auto source = photon_source();
            event.AddPhotonArrival(hit.GetHitId(), t, hit.GetHitId(), false, source.first, source.second);
          }
          for (int i = 0; i < photons_long; ++i) {
            const double t = deposit.t + long_middle / speed_of_light_in_fiber + exp_scint(generator) + exp_wsf(generator);
            const auto source = photon_source();
            event.AddPhotonArrival(hit.GetHitId(), t, hit.GetHitId(), true, source.first, source.second);
          }
        }
      }
      event.SortPhotonArrivals(hit.GetHitId());
      const double pe = pe_short + pe_long;
      TMS_TrueHit* adjustable_true_hit = event.GetAdjustableTrueHit(hit.GetHitId());
      if (adjustable_true_hit == nullptr) throw std::runtime_error("Fatal: SimulateOpticalModel() found a hit with no truth -- this stage is MC-only");
      adjustable_true_hit->SetPEAfterFibers(pe);
      adjustable_true_hit->SetPEAfterFibersLongPath(pe_long);
      adjustable_true_hit->SetPEAfterFibersShortPath(pe_short);
      hit.SetPE(pe);
      double reco_e = pe * TMS_Manager::GetInstance().Get_RECO_CALIBRATION_EnergyCalibration();
      hit.SetE(reco_e);
      hit.SetEVis(reco_e);
      continue;
    }

    double pe = hit.GetPE();

    // Applies birk's suppression
    // Expected present: this stage only ever runs on MC input, where truth was just
    // constructed for every hit in TMS_Event::ProcessTG4Event(). Fail fast rather than
    // segfault if that invariant is ever violated (e.g. this stage invoked on a truthless event).
    const TMS_TrueHit* true_hit = event.GetTrueHit(hit.GetHitId());
    if (true_hit == nullptr) throw std::runtime_error("Fatal: SimulateOpticalModel() found a hit with no truth -- this stage is MC-only");
    double de = true_hit->GetE();
    double dx = true_hit->GetdX();
    double dedx = 0;
    if (dx > 1e-8) dedx = de / dx;
    else dedx = de / 1.0;
    pe *= 1.0 / (1.0 + birks_constant * dedx);

    // Path efficiencies: WLS attenuation (and end reflection for the long way), optional
    // additional optical fiber, then fiber-to-readout coupling. Applied as a sequence of
    // multiplications so the default (draw-then-attenuate) path reproduces the previous
    // arithmetic exactly.
    double distance_from_end = 0;
    double long_way_distance_from_end = 0;
    if (should_simulate_fiber_lengths) {
      // Calculate the long and short path lengths
#ifdef USE_OLD_CODE
      double true_y = true_hit->GetY() / 1000.0; // m
      // In case of orthogonal (X) layers change to GetX()
      if (hit.GetBar().GetBarType() == TMS_Bar::kXBar) true_y = true_hit->GetX() / 1000.0;
      // assuming 0 is center, and assume we're reading out from top, then top would be biased negative and bottom positive, so -true_y.
      // TODO manually found this center. Make function in geom tools that returns values about scint
      // TODO fix math
      double distance_from_middle = TMS_Manager::GetInstance().Get_Geometry_YMIDDLE() - true_y;  // -1.54799
      distance_from_end = distance_from_middle + 2;
      long_way_distance_from_end = 4 + (4 - distance_from_end);
#else
      distance_from_end = GetTrueDistanceFromReadout(hit, *true_hit) * 1e-3; // m
      long_way_distance_from_end = GetTrueLongDistanceFromReadout(hit, *true_hit) * 1e-3; // m
#endif
      // In reality, light bounces so there's a multiplier
      // TODO it may be more realistic to make this non-linear
      distance_from_end *= wsf_length_multiplier;
      long_way_distance_from_end *= wsf_length_multiplier;
    }
    double pe_short = pe;
    double pe_long = 0;
    if (should_simulate_poisson_throws && poisson_after_attenuation) {
      // Photons are produced and each one independently survives to the sensor, so the
      // number detected on each path is Poisson with the fully attenuated expected mean
      // (half the produced light goes each way). Integer PE, same mean as below.
      pe_short = poisson_draw(attenuate_short(0.5 * pe, distance_from_end));
      pe_long = poisson_draw(attenuate_long(0.5 * pe, long_way_distance_from_end));
    } else {
      if (should_simulate_poisson_throws) {
        // Do a poisson throw to get the number of PE
        std::poisson_distribution<int> poisson(pe);
        pe = poisson(generator);
        // Now split the photons into the long and short paths with 50% chance of each
        std::binomial_distribution<int> binomial(pe, 0.5);
        pe_short = binomial(generator);
        pe_long = pe - pe_short;
      }
      pe_short = attenuate_short(pe_short, distance_from_end);
      pe_long = attenuate_long(pe_long, long_way_distance_from_end);
    }

    // Now save this information
    pe = pe_long + pe_short;

    // We want to save info right after fibers but before electronic conversion noise
    // This is particularly useful for timing information which cares about the first photon to be detected
    TMS_TrueHit* adjustable_true_hit = event.GetAdjustableTrueHit(hit.GetHitId());
    if (adjustable_true_hit == nullptr) throw std::runtime_error("Fatal: SimulateOpticalModel() found a hit with no truth -- this stage is MC-only");
    adjustable_true_hit->SetPEAfterFibers(pe);
    adjustable_true_hit->SetPEAfterFibersLongPath(pe_long);
    adjustable_true_hit->SetPEAfterFibersShortPath(pe_short);

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

  std::normal_distribution<double> noise_distribution(0.0, TMS_Readout_Manager::GetInstance().Get_Sim_Timing_ElectronicTimeNoise()); // ns
  double scintillator_decay_time = TMS_Readout_Manager::GetInstance().Get_Sim_Timing_ScintillatorDecayTime(); // ns
  double wsf_decay_time = TMS_Readout_Manager::GetInstance().Get_Sim_Timing_WLSDecayTime(); // ns
  std::exponential_distribution<double> exp_scint(1 / scintillator_decay_time);
  std::exponential_distribution<double> exp_wsf(1 / wsf_decay_time); // wavelength shifting fiber
  const double FIBER_N = TMS_Readout_Manager::GetInstance().Get_Sim_Timing_FiberRefractiveIndex();
  const double SPEED_OF_LIGHT_IN_FIBER = SPEED_OF_LIGHT / FIBER_N;

  if (TMS_Readout_Manager::GetInstance().Get_Sim_DetSim_UseResponseElements()) {
    // Response-element path: SimulateOpticalModel() already generated every detected photon
    // with its sensor arrival time; the hit time is the first arrival. The electronic time noise
    // is added once per readout in SimulateChannelReadout().
    for (auto& hit : TMS_Hits) {
      const std::vector<TMS_PhotonArrival>* arrivals = event.GetPhotonArrivals(hit.GetHitId());
      // No detected photons: keep the true time, as the legacy path does
      if (arrivals != nullptr && !arrivals->empty()) hit.SetT(arrivals->front().Time);
    }
    return;
  }

  const double wsf_length_multiplier = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_WSFLengthMultiplier();

  //double avg = 0;
  //double maxy = -1e9;
  //double miny = 1e9;
  //int n = 0;
  for (auto& hit : TMS_Hits) {
    double t = 0;
    // Expected present: this stage only ever runs on MC input, right after
    // SimulateOpticalModel() populated PEAfterFibers* for every hit. Fail fast rather than
    // segfault if that invariant is ever violated (e.g. this stage invoked on a truthless event).
    const TMS_TrueHit* true_hit = event.GetTrueHit(hit.GetHitId());
    if (true_hit == nullptr) throw std::runtime_error("Fatal: SimulateTimingModel() found a hit with no truth -- this stage is MC-only");
    // Random electronic timing noise (~1ns or less)
    t += noise_distribution(generator);
    // Optical fiber length delay (corrected to strip center)
    // (up to 13.4ns assuming 4m from edge, but correlated with y position. If delta y = 1m spread, than relative error is only 3.3ns)
#ifdef USE_OLD_CODE
    double true_y = true_hit->GetY() / 1000.0; // m
    // Making sure this gets changed for orthogonal (X) layers
    if (hit.GetBar().GetBarType() == TMS_Bar::kXBar) true_y = true_hit->GetX() / 1000.0;
    //miny = std::min(miny, true_y);
    //maxy = std::max(maxy, true_y);
    // assuming 0 is center, and assume we're reading out from top, then top would be biased negative and bottom positive, so -true_y.
    // TODO manually found this center. Want a better way in case things change
    double distance_from_middle = TMS_Manager::GetInstance().Get_Geometry_YMIDDLE() - true_y;  //-1.54799
    double long_way_distance = distance_from_middle + 8;
#else
    double distance_from_middle = GetTrueDistanceFromMiddle(hit, *true_hit) * 1e-3; // m
    double long_way_distance = GetTrueLongDistanceFromMiddle(hit, *true_hit) * 1e-3; // m
#endif
    // In reality, light bounces so there's a multiplier to the distance
    // todo, it may be more realistic to make this non-linear
    distance_from_middle *= wsf_length_multiplier;
    long_way_distance *= wsf_length_multiplier;

    // Find the time correction
    double time_correction = distance_from_middle / SPEED_OF_LIGHT_IN_FIBER;
    // This is the time correction if you go the long way instead
    double time_correction_long_way = long_way_distance / SPEED_OF_LIGHT_IN_FIBER;

    // Simulate every timing photon and use the first arrival as this hit's
    // representative time. A previous std::gamma_distribution shortcut was
    // wrong: Gamma(shape=N) describes a sum of N exponential draws, whereas
    // this model needs the minimum of N independent scintillator+WLS delays.
    // Keep the established 300-photon timing cap and ceil(mean PE) behavior;
    // this is intentionally a targeted timing correction, not a readout-model
    // redesign. PhotonArrivals records every sampled sensor arrival for a
    // future threshold/window model.
    double pe_short_path = true_hit->GetPEAfterFibersShortPath();
    double pe_long_path = true_hit->GetPEAfterFibersLongPath();
    double minimum_time_offset = 1e100;
    const double MAX_PE_THROWS = TMS_Readout_Manager::GetInstance().Get_Sim_Timing_MaxTimingPhotons();
    const int n_short_photons = std::min(static_cast<int>(std::ceil(pe_short_path)),
                                         static_cast<int>(MAX_PE_THROWS));
    const int n_long_photons = std::min(static_cast<int>(std::ceil(pe_long_path)),
                                        static_cast<int>(MAX_PE_THROWS));
    const double hit_time = hit.GetT();
    for (int i = 0; i < n_short_photons; ++i) {
      double time_offset = time_correction;
      time_offset += exp_scint(generator);
      time_offset += exp_wsf(generator);
      minimum_time_offset = std::min(time_offset, minimum_time_offset);
      event.AddPhotonArrival(hit.GetHitId(), hit_time + time_offset, hit.GetHitId(), false);
    }
    for (int i = 0; i < n_long_photons; ++i) {
      double time_offset = time_correction_long_way;
      time_offset += exp_scint(generator);
      time_offset += exp_wsf(generator);
      minimum_time_offset = std::min(time_offset, minimum_time_offset);
      event.AddPhotonArrival(hit.GetHitId(), hit_time + time_offset, hit.GetHitId(), true);
    }

    event.SortPhotonArrivals(hit.GetHitId());
    // Both paths had 0 PE: retain the existing no-slew fallback rather than
    // propagating the sentinel into the reconstructed hit time.
    if (minimum_time_offset == 1e100) minimum_time_offset = 0;
    t += minimum_time_offset;
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
    // Channel = (plane, bar, view). (z, NotZ) alone does not separate the two halves of a split
    // X-bar, which are read out at opposite ends.
    using ChannelKey = std::tuple<int, int, int>;
    std::map<ChannelKey, double> readout_map;
    std::map<ChannelKey, double> deadtime_map;
    std::map<ChannelKey, double> zombie_map;
    std::map<ChannelKey, bool> has_zombie_map;
    std::map<ChannelKey, double> x_map;
    std::map<ChannelKey, double> z_map;
    std::map<ChannelKey, double> t_map;
    //for (auto& hit : TMS_Hits) {
    for (size_t i = 0; i < TMS_Hits.size(); ++i) {
      auto& hit = TMS_Hits[i];
      double t = hit.GetT();
      // For a per-channel deadtime, this is a unique id for a channel.
      // But some detectors have deadtime for a whole board, in which case this should
      // return a single id for the whole board.
      const ChannelKey id(hit.GetPlaneNumber(), hit.GetBarNumber(), hit.GetBar().GetBarTypeNumber());
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
        if (deadtime_verbose) std::cout<<"Found channel we found already with Notz: "<<hit.GetNotZ()<<", z: "<<hit.GetZ()<<", plane: "<<std::get<0>(id)<<", bar: "<<std::get<1>(id)<<", view: "<<std::get<2>(id)<<"\n";
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
        auto it_has_zombie = has_zombie_map.find(id);
        if (it_dead != deadtime_map.end() && it_has_zombie != has_zombie_map.end() && it_has_zombie->second == true) {
          double deadtime_window_starting_from_end_of_deadtime = it_dead->second + deadtime;
          if (t < deadtime_window_starting_from_end_of_deadtime) {
            t = it_dead->second;
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
    double n_dead_hits_as_percent = TMS_Hits.empty() ? 0.0 : 100.0 * n_dead_hits / (float)TMS_Hits.size();
    std::cout<<"N dead hits: "<<n_dead_hits<<" out of "<<TMS_Hits.size()<<" hits. That's "<<n_dead_hits_as_percent<<"%"<<std::endl;
    if (zombie_time > 0) std::cout<<"N zombie hits: "<<n_zombie_hits<<std::endl;
  }
}

void TMS_DetectorSimulation::SimulateChannelReadout(TMS_Event &event, std::default_random_engine &generator) {
  // Same readout model as SimulateDeadtime() + TMS_SignalProcessing::MergeCoincidentHits(), done
  // once per channel:  |----  readout -----|------- deadtime ----{zombie time}]
  // The first hit in a channel opens a readout window of Sim.Readout.ReadoutTime; later hits in
  // the window are merged into it. The channel is then dead for Sim.Readout.Deadtime (if > 0):
  // hits in the deadtime are lost, except that hits in the final Sim.Readout.ZombieTime (if > 0)
  // are read at the end of the deadtime as the start of the next readout.
  std::vector<TMS_Hit> &TMS_Hits = event.GetHitsRawRef();
  const double readout_time = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_ReadoutTime();
  const double deadtime = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_Deadtime();
  const double zombie_time = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_ZombieTime();
  std::normal_distribution<double> time_noise(0.0, TMS_Readout_Manager::GetInstance().Get_Sim_Timing_ElectronicTimeNoise());

  std::map<TMS_ChannelId, std::vector<size_t>> hits_by_channel;
  for (size_t i = 0; i < TMS_Hits.size(); ++i) {
    if (TMS_Hits[i].GetPedSup()) continue;
    hits_by_channel[TMS_Hits[i].GetChannelId()].push_back(i);
  }

  std::vector<bool> remove(TMS_Hits.size(), false);
  std::vector<size_t> readouts;
  for (auto& channel : hits_by_channel) {
    std::vector<size_t>& idx = channel.second;
    std::stable_sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return TMS_Hits[a].GetT() < TMS_Hits[b].GetT(); });
    bool open = false;
    size_t current = 0;
    double t_read = 0, t_dead = 0;
    auto start_readout = [&](size_t i, double t) {
      current = i;
      open = true;
      readouts.push_back(i);
      t_read = t + readout_time;
      t_dead = t_read + std::max(deadtime, 0.0);
      if (deadtime > 0) {
        const TMS_Hit& hit = TMS_Hits[i];
        event.AddDeadtimeChannelRecord(std::make_pair(hit.GetNotZ(), hit.GetZ()), std::make_pair(t_read, t_dead), std::make_pair(t, t_read));
      }
    };
    for (size_t i : idx) {
      TMS_Hit& hit = TMS_Hits[i];
      const double t = hit.GetT();
      if (!open || t >= t_dead) {
        start_readout(i, t);
      } else if (t < t_read) {
        TMS_Hits[current].MergeWith(hit);
        event.MergeTrueHit(TMS_Hits[current].GetHitId(), hit.GetHitId());
        event.MergePhotonArrivals(TMS_Hits[current].GetHitId(), hit.GetHitId());
        event.MergeResponseSegments(TMS_Hits[current].GetHitId(), hit.GetHitId());
        remove[i] = true;
      } else if (zombie_time > 0 && t >= t_dead - zombie_time) {
        // Read at the end of the deadtime, as the start of the next readout
        hit.SetT(t_dead);
        start_readout(i, t_dead);
      } else {
        // Lost in the deadtime
        remove[i] = true;
      }
    }
  }

  // One electronic time measurement per readout
  for (size_t i : readouts) TMS_Hits[i].SetT(TMS_Hits[i].GetT() + time_noise(generator));

  std::vector<TMS_Hit> remaining_hits;
  for (size_t i = 0; i < TMS_Hits.size(); ++i) {
    if (!remove[i]) {
      remaining_hits.push_back(TMS_Hits[i]);
    } else {
      event.EraseTrueHit(TMS_Hits[i].GetHitId());
      event.ErasePhotonArrivals(TMS_Hits[i].GetHitId());
      event.EraseResponseSegments(TMS_Hits[i].GetHitId());
    }
  }
  // Same final ordering as MergeCoincidentHits()
  std::sort(remaining_hits.begin(), remaining_hits.end(), TMS_Hit::SortByZThenT);
  TMS_Hits = std::move(remaining_hits);
}

void TMS_DetectorSimulation::SimulateFrontEndTimingMode(TMS_Event &event, std::default_random_engine &generator) {
  std::vector<TMS_Hit> &TMS_Hits = event.GetHitsRawRef();
  TMS_Readout_Manager& config = TMS_Readout_Manager::GetInstance();
  const double peaking_time = config.Get_Sim_FrontEnd_FastShaperPeakingTime();
  const int order = config.Get_Sim_FrontEnd_FastShaperOrder();
  const double threshold = config.Get_Sim_FrontEnd_DiscriminatorThreshold();
  const double tdc_step = config.Get_Sim_FrontEnd_TDCStep();
  const double crosstalk = config.Get_Sim_FrontEnd_SiPMCrosstalkProbability();
  const double pixels = config.Get_Sim_FrontEnd_SiPMPixels();
  const bool energy_from_tot = config.Get_Sim_FrontEnd_EnergyFromToT();
  const double calibration = TMS_Manager::GetInstance().Get_RECO_CALIBRATION_EnergyCalibration();
  std::normal_distribution<double> jitter(0.0, config.Get_Sim_FrontEnd_TimingJitter());
  std::uniform_real_distribution<double> uniform(0.0, 1.0);

  // Single-photoelectron fast-shaper pulse, unit height at t = peaking_time
  auto pulse = [&](double dt) {
    if (dt <= 0) return 0.0;
    const double x = dt / peaking_time;
    return std::pow(x, order) * std::exp(order * (1.0 - x));
  };
  // Beyond this a pulse contributes < 1e-4 pe
  double tail = peaking_time;
  while (pulse(tail) > 1e-4) tail += 0.1 * peaking_time;
  const double grid = 0.1; // ns, waveform sampling
  auto floor_to_tdc = [&](double t) { return std::floor(t / tdc_step) * tdc_step; };

  // For EnergyFromToT: time over threshold of the mean waveform of A photoelectrons, tabulated once.
  // The mean single-photoelectron waveform is the fast-shaper pulse convolved with the photon
  // arrival-time spread (scintillator and WLS decay times, as in SimulateOpticalModel()), since
  // that spread lengthens the time over threshold; a real-data ToT calibration contains it too.
  // Depends only on the configuration, which is fixed for the job: build it once.
  static std::vector<std::pair<double, double>> tot_of_charge;
  if (energy_from_tot && tot_of_charge.empty()) {
    const double tau_scint = config.Get_Sim_Timing_ScintillatorDecayTime();
    const double tau_wls = config.Get_Sim_Timing_WLSDecayTime();
    auto arrival_pdf = [&](double s) {
      if (s <= 0) return 0.0;
      if (std::fabs(tau_wls - tau_scint) < 1e-9) return s * std::exp(-s / tau_wls) / (tau_wls * tau_wls);
      return (std::exp(-s / tau_wls) - std::exp(-s / tau_scint)) / (tau_wls - tau_scint);
    };
    const double step = 0.05; // ns
    const double span = tail + 15 * std::max(tau_scint, tau_wls);
    std::vector<double> mean_waveform;
    for (double t = 0; t < span; t += step) {
      double v = 0;
      for (double u = 0; u <= t; u += step) v += pulse(t - u) * arrival_pdf(u) * step;
      mean_waveform.push_back(v);
    }
    for (double a = threshold; a < 1e5; a *= 1.02) {
      int i0 = -1, i1 = -1;
      for (size_t i = 0; i < mean_waveform.size(); ++i) {
        const bool above = a * mean_waveform[i] >= threshold;
        if (above && i0 < 0) i0 = static_cast<int>(i);
        if (!above && i0 >= 0) { i1 = static_cast<int>(i); break; }
      }
      if (i0 >= 0 && i1 > i0) tot_of_charge.push_back({(i1 - i0) * step, a});
    }
  }
  auto charge_from_tot = [&](double tot) {
    if (tot_of_charge.empty() || tot <= tot_of_charge.front().first) return tot_of_charge.empty() ? threshold : tot_of_charge.front().second;
    for (size_t k = 1; k < tot_of_charge.size(); ++k) {
      if (tot <= tot_of_charge[k].first) {
        const auto& lo = tot_of_charge[k - 1];
        const auto& hi = tot_of_charge[k];
        return lo.second + (hi.second - lo.second) * (tot - lo.first) / (hi.first - lo.first);
      }
    }
    return tot_of_charge.back().second;
  };

  std::map<TMS_ChannelId, std::vector<size_t>> hits_by_channel;
  for (size_t i = 0; i < TMS_Hits.size(); ++i) {
    if (TMS_Hits[i].GetPedSup()) continue;
    hits_by_channel[TMS_Hits[i].GetChannelId()].push_back(i);
  }

  // Hits merged into another hit's discriminator pulse are removed; hits whose photons never
  // crossed the threshold stay, flagged as suppressed like the default pedestal threshold does.
  std::vector<bool> merged_away(TMS_Hits.size(), false);
  std::vector<bool> fired(TMS_Hits.size(), false);
  for (auto& channel : hits_by_channel) {
    // Every photoelectron in the channel: (arrival time, charge in pe including crosstalk, hit)
    struct Photoelectron { double t; double charge; size_t hit; };
    std::vector<Photoelectron> pes;
    for (size_t i : channel.second) {
      const std::vector<TMS_PhotonArrival>* arrivals = event.GetPhotonArrivals(TMS_Hits[i].GetHitId());
      if (arrivals == nullptr) continue;
      for (const auto& a : *arrivals) {
        double charge = 1;
        while (uniform(generator) < crosstalk) charge += 1;
        pes.push_back({a.Time, charge, i});
      }
    }
    if (pes.empty()) continue;
    std::sort(pes.begin(), pes.end(), [](const Photoelectron& a, const Photoelectron& b) { return a.t < b.t; });

    // Discriminator pulses: [t_up, t_down) where the summed waveform is at or above threshold
    struct DiscriminatorPulse { double t_up; double t_down; };
    std::vector<DiscriminatorPulse> pulses;
    size_t first = 0; // earliest photoelectron still within `tail` of the sampling time
    size_t next = 0;  // next photoelectron not yet started
    bool above = false;
    double t_up = 0;
    double t = pes.front().t;
    while (true) {
      while (next < pes.size() && pes[next].t <= t) ++next;
      while (first < next && t - pes[first].t > tail) ++first;
      if (first == next && !above) {
        // Nothing contributing: jump to the next photoelectron
        if (next >= pes.size()) break;
        t = pes[next].t;
        continue;
      }
      double v = 0;
      for (size_t k = first; k < next; ++k) v += pes[k].charge * pulse(t - pes[k].t);
      if (!above && v >= threshold) { above = true; t_up = t; }
      else if (above && v < threshold) { above = false; pulses.push_back({t_up, t}); }
      t += grid;
    }

    // Photoelectrons belong to the pulse they arrive before the end of (after the previous one)
    std::vector<std::map<size_t, double>> charge_by_hit(pulses.size());
    std::vector<double> pulse_charge(pulses.size(), 0);
    size_t p = 0;
    for (const auto& pe : pes) {
      while (p < pulses.size() && pe.t >= pulses[p].t_down) ++p;
      if (p == pulses.size()) break;
      charge_by_hit[p][pe.hit] += pe.charge;
      pulse_charge[p] += pe.charge;
    }
    // Each original hit goes (with its truth) to the pulse holding most of its charge; a pulse
    // no hit is assigned to (e.g. a late retrigger on one hit's tail) is dropped.
    std::map<size_t, size_t> pulse_of_hit;
    std::map<size_t, double> best_charge;
    for (size_t k = 0; k < pulses.size(); ++k) {
      for (const auto& hc : charge_by_hit[k]) {
        if (hc.second > best_charge[hc.first]) { best_charge[hc.first] = hc.second; pulse_of_hit[hc.first] = k; }
      }
    }
    for (size_t k = 0; k < pulses.size(); ++k) {
      size_t survivor = TMS_Hits.size();
      double survivor_charge = -1;
      for (const auto& hp : pulse_of_hit) {
        if (hp.second == k && charge_by_hit[k][hp.first] > survivor_charge) { survivor = hp.first; survivor_charge = charge_by_hit[k][hp.first]; }
      }
      if (survivor == TMS_Hits.size()) continue;
      TMS_Hit& hit = TMS_Hits[survivor];
      for (const auto& hp : pulse_of_hit) {
        if (hp.second != k || hp.first == survivor) continue;
        merged_away[hp.first] = true;
        event.MergeTrueHit(hit.GetHitId(), TMS_Hits[hp.first].GetHitId());
        event.MergePhotonArrivals(hit.GetHitId(), TMS_Hits[hp.first].GetHitId());
        event.MergeResponseSegments(hit.GetHitId(), TMS_Hits[hp.first].GetHitId());
      }
      const double t_stamp = floor_to_tdc(pulses[k].t_up + jitter(generator));
      const double tot = floor_to_tdc(pulses[k].t_down) - floor_to_tdc(pulses[k].t_up);
      // SiPM pixel saturation
      const double pe = pixels * (1.0 - std::exp(-pulse_charge[k] / pixels));
      hit.SetT(t_stamp);
      hit.SetToT(tot);
      hit.SetPE(pe);
      const double energy = (energy_from_tot ? charge_from_tot(tot) : pe) * calibration;
      hit.SetE(energy);
      hit.SetEVis(energy);
      fired[survivor] = true;
    }
  }

  std::vector<TMS_Hit> remaining_hits;
  for (size_t i = 0; i < TMS_Hits.size(); ++i) {
    if (merged_away[i]) {
      event.EraseTrueHit(TMS_Hits[i].GetHitId());
      event.ErasePhotonArrivals(TMS_Hits[i].GetHitId());
      event.EraseResponseSegments(TMS_Hits[i].GetHitId());
      continue;
    }
    if (!fired[i]) TMS_Hits[i].SetPedSup(true);
    remaining_hits.push_back(TMS_Hits[i]);
  }
  std::sort(remaining_hits.begin(), remaining_hits.end(), TMS_Hit::SortByZThenT);
  TMS_Hits = std::move(remaining_hits);
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
