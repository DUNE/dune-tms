#ifndef _TMS_DETECTORSIMULATION_H_SEEN_
#define _TMS_DETECTORSIMULATION_H_SEEN_

#include <random>

#include "TMS_Event.h"

// Sim-only detector-response steps: methods that only make sense when simulating a
// detector response from truth information, not applicable to real DAQ readout.
// Counterpart to TMS_SignalProcessing, which handles the real-or-simulated steps.
class TMS_DetectorSimulation {

  public:

    static TMS_DetectorSimulation& GetInstance() {
      static TMS_DetectorSimulation Instance;
      return Instance;
    }

    void SimulateOpticalModel(TMS_Event &event, std::default_random_engine &generator);
    void SimulateDarkCount(TMS_Event &event);
    void SimulateTimingModel(TMS_Event &event, std::default_random_engine &generator);
    void SimulateDeadtime(TMS_Event &event);
    // Response-element pipeline (Sim.DetSim.UseResponseElements): readout windows, deadtime and
    // zombie time per TMS_ChannelId, replacing SimulateDeadtime() + the post-simulation
    // MergeCoincidentHits(). Adds the electronic time noise once per readout.
    void SimulateChannelReadout(TMS_Event &event, std::default_random_engine &generator);
    // A5202 timing mode (Sim.DetSim.FrontEndTimingMode, with UseResponseElements): per channel,
    // sums the photons' fast-shaper pulses; each discriminator pulse becomes a hit with the
    // threshold-crossing time and time over threshold in TDC steps. Replaces
    // SimulateChannelReadout(), SimulateReadoutNoise() and the pedestal threshold.
    void SimulateFrontEndTimingMode(TMS_Event &event, std::default_random_engine &generator);
    void SimulateReadoutNoise(TMS_Event &event, std::default_random_engine &generator);

  private:
    TMS_DetectorSimulation() {};
    TMS_DetectorSimulation(TMS_DetectorSimulation const &) = delete;
    void operator=(TMS_DetectorSimulation const &) = delete;
    ~TMS_DetectorSimulation() {};

};

#endif
