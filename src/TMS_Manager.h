#ifndef __TMS_MANAGER_H__
#define __TMS_MANAGER_H__

#include <iostream>
#include <string>
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
    int Get_Reco_DBSCAN_PreDBNeighbours() { return _RECO_DBSCAN_PreDBNeighbours; };
    double Get_Reco_DBSCAN_PreDBDistance() { return _RECO_DBSCAN_PreDBDistance; };

    int Get_Reco_HOUGH_MaxHough() { return _RECO_HOUGH_MaxHough; };
    double Get_Reco_HOUGH_MinInterp() { return _RECO_HOUGH_MinInterp; };
    double Get_Reco_HOUGH_MaxInterp() { return _RECO_HOUGH_MaxInterp; };
    double Get_Reco_HOUGH_MinSlope() { return _RECO_HOUGH_MinSlope; };
    double Get_Reco_HOUGH_MaxSlope() { return _RECO_HOUGH_MaxSlope; };
    int Get_Reco_HOUGH_NSlope() { return _RECO_HOUGH_NSlope; };
    int Get_Reco_HOUGH_NInter() { return _RECO_HOUGH_NInter; };
    double Get_Reco_HOUGH_HitMult() { return _RECO_HOUGH_HitMult; };
    double Get_Reco_HOUGH_MergeTracks() { return _RECO_HOUGH_MergeTracks; };
    bool Get_Reco_HOUGH_RunAStar() { return _RECO_HOUGH_RunAStar; };
    bool Get_Reco_HOUGH_FirstCluster() { return _RECO_HOUGH_FirstCluster; };
    double Get_Reco_HOUGH_MinDist() { return _RECO_HOUGH_MinDist; };

    bool Get_Reco_EXTRAPOLATION_Extrapolation() { return _RECO_EXTRAPOLATION_Extrapolation; };
    int Get_Reco_EXTRAPOLATION_ExtrapolateDist() { return _RECO_EXTRAPOLATION_ExtrapolateDist; };
    int Get_Reco_EXTRAPOLATION_ExtrapolateLimit() { return _RECO_EXTRAPOLATION_ExtrapolateLimit; };
    int Get_Reco_EXTRAPOLATION_NumBarsEnd() { return _RECO_EXTRAPOLATION_NumBarsEnd; };
    int Get_Reco_EXTRAPOLATION_NumBarsStart() { return _RECO_EXTRAPOLATION_NumBarsStart; };

    int Get_Reco_TRACKMATCH_PlaneLimit() { return _RECO_TRACKMATCH_PlaneLimit; };
    int Get_Reco_TRACKMATCH_BarLimit() { return _RECO_TRACKMATCH_BarLimit; };
    int Get_Reco_TRACKMATCH_TimeLimit() { return _RECO_TRACKMATCH_TimeLimit; };
    int Get_Reco_TRACKMATCH_XTimeLimit() { return _RECO_TRACKMATCH_XTimeLimit; };
    float Get_Reco_TRACKMATCH_YAnchor() { return _RECO_TRACKMATCH_YAnchor; };
    double Get_Reco_TRACKMATCH_TiltAngle() { return _RECO_TRACKMATCH_TiltAngle; };
    float Get_Reco_TRACKMATCH_YDifference() { return _RECO_TRACKMATCH_YDifference; };
    int Get_Reco_TRACKMATCH_DirectionDistance() { return _RECO_TRACKMATCH_DirectionDistance; };

    bool Get_Reco_ASTAR_IsGreedy() { return _RECO_ASTAR_IsGreedy; };
    std::string Get_Reco_ASTAR_CostMetric() { return _RECO_ASTAR_CostMetric; };

    int Get_Reco_STOPPING_nLastHits() { return _RECO_STOPPING_nLastHits; };
    double Get_Reco_STOPPING_EnergyCut() { return _RECO_STOPPING_EnergyCut; };

    std::string Get_Reco_TrackMethod() { return _RECO_TRACK_METHOD; };
    bool Get_Reco_Clustering() { return _RECO_CLUSTERING; };

    bool Get_Reco_Kalman_Run() { return _RECO_KALMAN_RUN; };
    double Get_Reco_Kalman_Assumed_Charge() { return _RECO_KALMAN_ASSUMED_CHARGE; };

    bool Get_LightWeight_Truth() { return _TRUTH_LIGHTWEIGHT; };

    bool Get_Reco_TRACKSMOOTHING_UseTrackSmoothing() { return _RECO_TRACKSMOOTHING_UseTrackSmoothing; };
    std::string Get_Reco_TRACKSMOOTHING_TrackSmoothingStrategy() { return _RECO_TRACKSMOOTHING_TrackSmoothingStrategy; };
    double Get_Reco_TRACKSMOOTHING_MaxYDistanceBetweenUVTransitionPoints() { return _RECO_TRACKSMOOTHING_MaxYDistanceBetweenUVTransitionPoints; };
    double Get_Reco_TRACKSMOOTHING_UncertaintyGoodDirection() { return _RECO_TRACKSMOOTHING_UncertaintyGoodDirection; };
    double Get_Reco_TRACKSMOOTHING_UncertaintyBadDirection() { return _RECO_TRACKSMOOTHING_UncertaintyBadDirection; };
    double Get_Reco_TRACKSMOOTHING_UncertaintyForUVTransitionPoints() { return _RECO_TRACKSMOOTHING_UncertaintyForUVTransitionPoints; };

    bool Get_Reco_TIME_RunTimeSlicer() { return _RECO_TIME_RunTimeSlicer; };
    bool Get_Reco_TIME_RunSimpleTimeSlicer() { return _RECO_TIME_RunSimpleTimeSlicer; };
    double Get_RECO_TIME_TimeSlicerThresholdStart() { return _RECO_TIME_TimeSlicerThresholdStart; };
    double Get_RECO_TIME_TimeSlicerThresholdEnd() { return _RECO_TIME_TimeSlicerThresholdEnd; };
    double Get_RECO_TIME_TimeSlicerSliceUnit() { return _RECO_TIME_TimeSlicerSliceUnit; };
    int Get_RECO_TIME_TimeSlicerEnergyWindowInUnits() { return _RECO_TIME_TimeSlicerEnergyWindowInUnits; };
    int Get_RECO_TIME_TimeSlicerMinimumSliceWidthInUnits() { return _RECO_TIME_TimeSlicerMinimumSliceWidthInUnits; };
    double Get_RECO_TIME_TimeSlicerMaxTime() { return _RECO_TIME_TimeSlicerMaxTime; };
    
    double Get_RECO_CALIBRATION_EnergyCalibration() { return _RECO_CALIBRATION_EnergyCalibration; };
    
    double Get_FIDUCIAL_TMS_START_X() { return _FIDUCIAL_TMS_START_X; };
    double Get_FIDUCIAL_TMS_START_Y() { return _FIDUCIAL_TMS_START_Y; };
    double Get_FIDUCIAL_TMS_START_Z() { return _FIDUCIAL_TMS_START_Z; };
    double Get_FIDUCIAL_TMS_END_X() { return _FIDUCIAL_TMS_END_X; };
    double Get_FIDUCIAL_TMS_END_Y() { return _FIDUCIAL_TMS_END_Y; };
    double Get_FIDUCIAL_TMS_END_Z() { return _FIDUCIAL_TMS_END_Z; };
    double Get_ACTIVE_LAR_START_X() { return _ACTIVE_LAR_START_X; };
    double Get_ACTIVE_LAR_START_Y() { return _ACTIVE_LAR_START_Y; };
    double Get_ACTIVE_LAR_START_Z() { return _ACTIVE_LAR_START_Z; };
    double Get_ACTIVE_LAR_END_X() { return _ACTIVE_LAR_END_X; };
    double Get_ACTIVE_LAR_END_Y() { return _ACTIVE_LAR_END_Y; };
    double Get_ACTIVE_LAR_END_Z() { return _ACTIVE_LAR_END_Z; };

    double Get_LAR_FIDUCIAL_Z_CUT() { return _LAR_FIDUCIAL_Z_CUT; };
    double Get_LAR_FIDUCIAL_XY_CUT() { return _LAR_FIDUCIAL_XY_CUT; };
    double Get_LAR_OUTER_SHELL_THICKNESS() { return _LAR_OUTER_SHELL_THICKNESS; };
    double Get_ND_PHYSICS_MUON_LAR_SHELL_CUT_ENERGY() { return _ND_PHYSICS_MUON_LAR_SHELL_CUT_ENERGY; };

    bool Get_DrawPDF() { return _APPLICATIONS_DrawPDF; };
    int Get_MaximumNEvents() { return _APPLICATIONS_MaximumNEvents; };

    double Get_Geometry_YMIDDLE() { return _GEOMETRY_YMIDDLE; };
    
    double Get_Nersc_Spill_Period() { return _NERSC_SPILL_PERIOD; };
    void Set_Nersc_Spill_Period(double value) { _NERSC_SPILL_PERIOD = value; };

  private:
    TMS_Manager();
    TMS_Manager(TMS_Manager const &) = delete;
    void operator=(TMS_Manager const &) = delete;
    ~TMS_Manager() {};

    std::string Filename;

    int _RECO_MinHits;

    int _RECO_DBSCAN_MinPoints;
    double _RECO_DBSCAN_Epsilon;
    int _RECO_DBSCAN_PreDBNeighbours;
    double _RECO_DBSCAN_PreDBDistance;

    int _RECO_HOUGH_MaxHough;
    double _RECO_HOUGH_MinInterp;
    double _RECO_HOUGH_MaxInterp;
    double _RECO_HOUGH_MinSlope;
    double _RECO_HOUGH_MaxSlope;
    double _RECO_HOUGH_NSlope;
    double _RECO_HOUGH_NInter;
    double _RECO_HOUGH_HitMult;
    bool _RECO_HOUGH_MergeTracks;
    bool _RECO_HOUGH_RunAStar;
    bool _RECO_HOUGH_FirstCluster;
    double _RECO_HOUGH_MinDist;

    bool _RECO_EXTRAPOLATION_Extrapolation;
    int _RECO_EXTRAPOLATION_ExtrapolateDist;
    int _RECO_EXTRAPOLATION_ExtrapolateLimit;
    int _RECO_EXTRAPOLATION_NumBarsEnd;
    int _RECO_EXTRAPOLATION_NumBarsStart;

    int _RECO_TRACKMATCH_PlaneLimit;
    int _RECO_TRACKMATCH_BarLimit;
    int _RECO_TRACKMATCH_TimeLimit;
    int _RECO_TRACKMATCH_XTimeLimit;
    float _RECO_TRACKMATCH_YAnchor;
    double _RECO_TRACKMATCH_TiltAngle;
    float _RECO_TRACKMATCH_YDifference;
    int _RECO_TRACKMATCH_DirectionDistance;

    bool _RECO_TRACKSMOOTHING_UseTrackSmoothing;
    std::string _RECO_TRACKSMOOTHING_TrackSmoothingStrategy;
    double _RECO_TRACKSMOOTHING_MaxYDistanceBetweenUVTransitionPoints;
    double _RECO_TRACKSMOOTHING_UncertaintyGoodDirection;
    double _RECO_TRACKSMOOTHING_UncertaintyBadDirection;
    double _RECO_TRACKSMOOTHING_UncertaintyForUVTransitionPoints;

    bool _RECO_ASTAR_IsGreedy;
    std::string _RECO_ASTAR_CostMetric;

    bool _RECO_STOPPING_nLastHits;
    double _RECO_STOPPING_EnergyCut;

    std::string _RECO_TRACK_METHOD;
    bool _RECO_CLUSTERING;
    
    double _NERSC_SPILL_PERIOD;

    bool _RECO_KALMAN_RUN; // Whether we run Kalman filter or no
    double _RECO_KALMAN_ASSUMED_CHARGE; //set the assumed charge of the track in the kalman filter

    // Lightweight trajectory saving (ignore small trajectories and gammas)
    bool _TRUTH_LIGHTWEIGHT;
    
    bool _RECO_TIME_RunTimeSlicer;
    bool _RECO_TIME_RunSimpleTimeSlicer;
    double _RECO_TIME_TimeSlicerThresholdStart;
    double _RECO_TIME_TimeSlicerThresholdEnd;
    double _RECO_TIME_TimeSlicerSliceUnit;
    int _RECO_TIME_TimeSlicerEnergyWindowInUnits;
    int _RECO_TIME_TimeSlicerMinimumSliceWidthInUnits;
    double _RECO_TIME_TimeSlicerMaxTime;
    
    double _RECO_CALIBRATION_EnergyCalibration;
    
    double _FIDUCIAL_TMS_START_X;
    double _FIDUCIAL_TMS_START_Y;
    double _FIDUCIAL_TMS_START_Z;
    double _FIDUCIAL_TMS_END_X;
    double _FIDUCIAL_TMS_END_Y;
    double _FIDUCIAL_TMS_END_Z;
    double _ACTIVE_LAR_START_X;
    double _ACTIVE_LAR_START_Y;
    double _ACTIVE_LAR_START_Z;
    double _ACTIVE_LAR_END_X;
    double _ACTIVE_LAR_END_Y;
    double _ACTIVE_LAR_END_Z;
    
    double _LAR_FIDUCIAL_Z_CUT;
    double _LAR_FIDUCIAL_XY_CUT;
    double _LAR_OUTER_SHELL_THICKNESS;
    double _ND_PHYSICS_MUON_LAR_SHELL_CUT_ENERGY;

    bool _APPLICATIONS_DrawPDF;
    int _APPLICATIONS_MaximumNEvents;

    double _GEOMETRY_YMIDDLE;
};

#endif
