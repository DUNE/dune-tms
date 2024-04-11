#ifndef __TMS_UTILS_H__
#define __TMS_UTILS_H__

#include <iostream>
#include <unordered_map>
#include <vector>

// Can only run this if we have linkage to CAF format
#ifdef DUNEANAOBJ_ENABLED
// Include standard record format for TMS
#include "duneanaobj/StandardRecord/SRTMS.h"
#endif

class TMS_Hit;

namespace TMS_Utils {
  struct ParticleInfo {
    double total_energy = 0;
    std::vector<int> indices;
    std::vector<double> energies;
  };

  // Get the sign of a number
  template <typename T> inline int sgn(T value) {
    return (T(0)<value - (value<T(0)));
  }
  
  ParticleInfo GetPrimaryIdsByEnergy(const std::vector<TMS_Hit>& hits);
  ParticleInfo GetSumAndHighest(const std::unordered_map<int, double>& map);

// Can only run this if we have linkage to CAF format
#ifdef DUNEANAOBJ_ENABLED
  caf::SRTMS ConvertEvent();
#endif

}

#endif
