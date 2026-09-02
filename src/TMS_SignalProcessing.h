#ifndef _TMS_SIGNALPROCESSING_H_SEEN_
#define _TMS_SIGNALPROCESSING_H_SEEN_

#include "TMS_Event.h"

// Real-or-simulated signal-processing steps: methods a real DAQ readout would also
// need to run on digitized channel hits, regardless of whether those hits came from
// TMS_DetectorSimulation or a real detector. This is the seam where real data can join.
class TMS_SignalProcessing {

  public:

    static TMS_SignalProcessing& GetInstance() {
      static TMS_SignalProcessing Instance;
      return Instance;
    }

    void MergeCoincidentHits(TMS_Event &event);
    void SimulatePedestalSubtraction(TMS_Event &event);

  private:
    TMS_SignalProcessing() {};
    TMS_SignalProcessing(TMS_SignalProcessing const &) = delete;
    void operator=(TMS_SignalProcessing const &) = delete;
    ~TMS_SignalProcessing() {};

};

#endif
