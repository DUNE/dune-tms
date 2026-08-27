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
    void SimulateReadoutNoise(TMS_Event &event, std::default_random_engine &generator);

  private:
    TMS_DetectorSimulation() {};
    TMS_DetectorSimulation(TMS_DetectorSimulation const &) = delete;
    void operator=(TMS_DetectorSimulation const &) = delete;
    ~TMS_DetectorSimulation() {};

};

#endif
