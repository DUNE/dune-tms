#ifndef __TMS_UTILS_H__
#define __TMS_UTILS_H__

#include <iostream>
// could switch back to an unordered map if you define a hash function for pairs
#include <map>
#include <vector>

// Can only run this if we have linkage to CAF format
#ifdef DUNEANAOBJ_ENABLED
// Include standard record format for TMS
#include "duneanaobj/StandardRecord/SRTMS.h"
#endif

class TMS_Hit;

namespace TMS_Utils {
/*
  struct PairID {
    int vertexid;
    int trackid;
    
    PairID(v, t) : vertexid(v), trackid(t) {}
  };*/
  
  struct ParticleInfo {
    double total_energy = 0;
    std::vector<int> trackids;
    std::vector<int> vertexids;
    std::vector<double> energies;
    
    void Print() {
      std::cout<<"ParticleInfo.total_energy = "<<total_energy<<std::endl;
      std::cout<<"n trackids, vertexids and energies = "<<trackids.size()<<std::endl;
      for (size_t i = 0; i < trackids.size(); i++) {
        std::cout<<trackids.at(i)<<": "<<energies.at(i)<<std::endl;
        std::cout<<vertexids.at(i)<<": "<<energies.at(i)<<std::endl;
      }
    }
  };

  // Get the sign of a number
  template <typename T> inline int sgn(T value) {
    return (T(0)<value - (value<T(0)));
  }
  
  ParticleInfo GetPrimaryIdsByEnergy(const std::vector<TMS_Hit>& hits);
  ParticleInfo GetSumAndHighest(const std::map<std::pair<int, int>, double>& map);

// Can only run this if we have linkage to CAF format
#ifdef DUNEANAOBJ_ENABLED
  caf::SRTMS ConvertEvent();
#endif

}

#endif
