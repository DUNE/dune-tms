#include "TMS_Decoder.h"
#include "TMS_Readout_Manager.h"

std::vector<TMS_Hit> TMS_MCDecoder::Decode(const std::vector<TMS_RawReadout>& readouts) {
  TMS_Readout_Manager& config = TMS_Readout_Manager::GetInstance();
  const bool timing_mode = config.Get_Sim_DetSim_FrontEndTimingMode();
  const double threshold = config.Get_Sim_Readout_PedestalSubtractionThreshold();

  std::vector<TMS_Hit> hits;
  hits.reserve(readouts.size());
  for (const TMS_RawReadout& readout : readouts) {
    // No calibration yet: the readout's energies are the simulation's
    TMS_Hit hit(readout.Bar, readout.Energy, readout.Time, readout.Charge);
    hit.SetHitId(readout.HitId);
    hit.SetEVis(readout.EnergyVisible);
    hit.SetToT(readout.ToT);
    hit.SetPedSup(!readout.Triggered || (!timing_mode && readout.Charge < threshold));
    hits.push_back(std::move(hit));
  }
  return hits;
}
