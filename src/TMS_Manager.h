#ifndef __TMS_MANAGER_H__
#define __TMS_MANAGER_H__

#include <iostream>
#include "toml.hpp"

// Just a global parameter manager
class TMS_Manager {

  public:
    static TMS_Manager& GetInstance() {
      static TMS_Manager Instance;
      return Instance;
    }

    void SetFileName(std::string file) { 
      std::cout << "Setting global edep-sim filename to " << file << std::endl;
      Filename = file; 
    };
    std::string GetFileName() { return Filename; };

    int Get_Reco_MinHits() { return _RECO_MinHits; };

    int Get_Reco_DBSCAN_MinPoints() { return _RECO_DBSCAN_MinPoints; };
    double Get_Reco_DBSCAN_Epsilon() { return _RECO_DBSCAN_Epsilon; };

    int Get_Reco_HOUGH_MaxHough() { return _RECO_HOUGH_MaxHough; };
    double Get_Reco_HOUGH_MinInterp() { return _RECO_HOUGH_MinInterp; };
    double Get_Reco_HOUGH_MaxInterp() { return _RECO_HOUGH_MaxInterp; };
    double Get_Reco_HOUGH_MinSlope() { return _RECO_HOUGH_MinSlope; };
    double Get_Reco_HOUGH_MaxSlope() { return _RECO_HOUGH_MaxSlope; };
    double Get_Reco_HOUGH_HitMult() { return _RECO_HOUGH_HitMult; };

  private:
    TMS_Manager();
    TMS_Manager(TMS_Manager const &) = delete;
    void operator=(TMS_Manager const &) = delete;
    ~TMS_Manager() {};

    std::string Filename;

    int _RECO_MinHits;

    int _RECO_DBSCAN_MinPoints;
    double _RECO_DBSCAN_Epsilon;

    int _RECO_HOUGH_MaxHough;
    double _RECO_HOUGH_MinInterp;
    double _RECO_HOUGH_MaxInterp;
    double _RECO_HOUGH_MinSlope;
    double _RECO_HOUGH_MaxSlope;
    double _RECO_HOUGH_HitMult;

};

#endif
