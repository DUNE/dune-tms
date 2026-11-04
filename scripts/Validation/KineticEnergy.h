#ifndef TMS_VALIDATION_KINETIC_ENERGY_H
#define TMS_VALIDATION_KINETIC_ENERGY_H

#include <cmath>
#include <cstdlib>

// Component [3] of the true four-momentum branches (MomentumTMSStart, BirthMomentum,
// RecoTrackPrimaryParticleTrueMomentum, ...EnteringTMS, ...) is the particle's TOTAL energy,
// sqrt(p^2 + m^2). The validation plots are labeled "KE", and the reconstructed energy they are
// compared with (a range-based estimate) is a kinetic energy, so the true value has to be
// E - m. Using E directly puts every true "KE" axis one particle mass too high (105.7 MeV for a
// muon) and shows up as a spurious offset in the energy resolution.

// Mass in MeV of the particles the validation selections use, keyed by |PDG code|. Returns a
// negative number for a particle that is not listed, so the caller can tell.
inline double ParticleMassMeV(int pdg) {
  switch (std::abs(pdg)) {
  case 11:   return 0.510999;  // electron / positron
  case 13:   return 105.6584;  // muon
  case 22:   return 0.0;       // photon
  case 111:  return 134.9768;  // pi0
  case 211:  return 139.5704;  // charged pion
  case 321:  return 493.677;   // charged kaon
  case 2112: return 939.5654;  // neutron
  case 2212: return 938.2721;  // proton
  default:   return -1.0;
  }
}

// True kinetic energy in MeV from the total energy (MeV) of a particle with the given PDG code.
// A value that is not a physical total energy (the -999 default the trees use when a branch was
// not filled) is returned unchanged so it still fails the same validity checks as before, and so
// is the energy of a particle whose mass is not in the table above.
inline double KineticEnergyMeV(int pdg, double total_energy_mev) {
  if (!(total_energy_mev > 0.0)) return total_energy_mev;
  const double mass = ParticleMassMeV(pdg);
  if (mass < 0.0) return total_energy_mev;
  return total_energy_mev - mass;
}

#endif // TMS_VALIDATION_KINETIC_ENERGY_H
