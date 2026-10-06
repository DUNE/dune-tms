#ifndef _TMS_DECODER_H_SEEN_
#define _TMS_DECODER_H_SEEN_

#include <vector>

#include "TMS_Hit.h"
#include "TMS_RawReadout.h"

// Turns raw readouts into the TMS_Hit objects that reconstruction uses. This is where calibration
// and the pedestal threshold belong, so that simulation and data go through the same step; the
// simulation's decoder is here, a decoder for real data will sit beside it once the raw format is
// defined.
class TMS_MCDecoder {
  public:
    // One hit per raw readout, in the same order, carrying the readout's HitId. Hits whose
    // channel did not trigger, and (in the charge mode) hits below Sim.Readout.PedestalSubtractionThreshold,
    // are kept and flagged as pedestal suppressed, as the simulation has always done. The
    // A5202 timing mode has no separate threshold: its discriminator is the threshold.
    static std::vector<TMS_Hit> Decode(const std::vector<TMS_RawReadout>& readouts);
};

#endif
