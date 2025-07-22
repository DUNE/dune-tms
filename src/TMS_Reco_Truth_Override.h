#ifndef _TMS_RECO_TRUTH_OVERRIDE_H_SEEN_
#define _TMS_RECO_TRUTH_OVERRIDE_H_SEEN_

#include <vector>

class TMS_Hit;
class TMS_Event;
class TMS_Track;
class TMS_Kalman;

class TMS_Reco_Truth_Override {
public:
  static bool FindTracksUsingTruth(std::vector<TMS_Hit> &hits, TMS_Event &event,
                                   std::vector<TMS_Track> &tracks);
  // static bool std::vector<TMS_Hit>& hits, std::vector<std::vector<TMS_Hit>>&
  // tracks
  static std::vector<std::pair<int, int>>
  GetAllUniqueParticles(std::vector<TMS_Hit> &hits);
  static std::vector<std::vector<TMS_Hit>>
  GetSeparateHitsByParticle(std::vector<TMS_Hit> &hits,
                            std::vector<std::pair<int, int>> particles);
  static TMS_Kalman RunKalmanFilter(std::vector<TMS_Hit> hits, double charge);
  static TMS_Track RunTrackReco(std::vector<TMS_Hit> hits);
  static bool CanKalmanProcess(std::vector<TMS_Hit> &hits);
  static double GetEnergyForOccupancy(std::vector<TMS_Hit> &hits);
  static void AdjustHitInformation(std::vector<TMS_Hit> &hits);
};

#endif
