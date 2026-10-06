#ifndef _TMS_RAWREADOUT_H_SEEN_
#define _TMS_RAWREADOUT_H_SEEN_

#include "TMS_Hit.h"

// What one electronics channel reports for one readout: the output of the response-element
// detector simulation (Sim.DetSim.UseResponseElements), and the stand-in for a record of the raw
// data, whose format is not settled. Calibration and the pedestal threshold are applied when a
// decoder turns it into a TMS_Hit (TMS_MCDecoder for the simulation; a data decoder will do the
// same job for real data), not by the electronics simulation.
//
// For now it carries both the charge (PE) and the time over threshold of the A5202 timing mode,
// plus the energies the simulation gives, which stand in for a calibration; the raw record will
// hold counts, and the decoder the conversions, once the DAQ format exists.
struct TMS_RawReadout {
  // The electronics channel: plane, bar (or bar half) and view
  TMS_ChannelId Channel;
  // Where the channel is. A placeholder for the channel map that a decoder will use to find a
  // bar from a channel id; copied from the simulated hit's bar.
  TMS_Bar Bar;
  // Key of this readout's truth (TMS_TrueHit, photon arrivals, response segments) in TMS_Event.
  // Simulation only: real data has none.
  int HitId;
  // Time stamp of the readout, ns
  double Time;
  // Charge, in photoelectrons
  double Charge;
  // Time over threshold, ns; -999 where not measured (charge mode)
  double ToT;
  // Energy, MeV, as the electronics simulation gives it (the decoder's calibration, for now)
  double Energy;
  double EnergyVisible;
  // False where the channel did not produce a readout above its hardware threshold (the A5202
  // discriminator never fired). The record is kept only so decoding reproduces the hits the
  // simulation has always written, flagged as pedestal suppressed.
  bool Triggered;

  TMS_RawReadout(const TMS_ChannelId& channel, const TMS_Bar& bar, int hit_id, double time, double charge,
                 double tot, double energy, double energy_visible, bool triggered) :
    Channel(channel), Bar(bar), HitId(hit_id), Time(time), Charge(charge), ToT(tot), Energy(energy),
    EnergyVisible(energy_visible), Triggered(triggered) {}
};

#endif
