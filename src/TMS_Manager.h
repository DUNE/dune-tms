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
    std::string GetFileName() const { return Filename; };

    int Get_Reco_MinHits() const { return _RECO_MinHits; };

    int Get_Reco_DBSCAN_MinPoints() const { return _RECO_DBSCAN_MinPoints; };
    double Get_Reco_DBSCAN_Epsilon() const { return _RECO_DBSCAN_Epsilon; };
    int Get_Reco_DBSCAN_PreDBNeighbours() const { return _RECO_DBSCAN_PreDBNeighbours; };
    double Get_Reco_DBSCAN_PreDBDistance() const { return _RECO_DBSCAN_PreDBDistance; };

    int Get_Reco_HOUGH_MaxHough() const { return _RECO_HOUGH_MaxHough; };
    double Get_Reco_HOUGH_MinInterp() const { return _RECO_HOUGH_MinInterp; };
    double Get_Reco_HOUGH_MaxInterp() const { return _RECO_HOUGH_MaxInterp; };
    double Get_Reco_HOUGH_MinSlope() const { return _RECO_HOUGH_MinSlope; };
    double Get_Reco_HOUGH_MaxSlope() const { return _RECO_HOUGH_MaxSlope; };
    int Get_Reco_HOUGH_NSlope() const { return _RECO_HOUGH_NSlope; };
    int Get_Reco_HOUGH_NInter() const { return _RECO_HOUGH_NInter; };
    double Get_Reco_HOUGH_HitMult() const { return _RECO_HOUGH_HitMult; };
    double Get_Reco_HOUGH_MergeTracks() const { return _RECO_HOUGH_MergeTracks; };
    bool Get_Reco_HOUGH_RunAStar() const { return _RECO_HOUGH_RunAStar; };
    bool Get_Reco_HOUGH_FirstCluster() const { return _RECO_HOUGH_FirstCluster; };
    double Get_Reco_HOUGH_MinDist() const { return _RECO_HOUGH_MinDist; };

    bool Get_Reco_EXTRAPOLATION_Extrapolation() const { return _RECO_EXTRAPOLATION_Extrapolation; };
    int Get_Reco_EXTRAPOLATION_ExtrapolateDist() const { return _RECO_EXTRAPOLATION_ExtrapolateDist; };
    int Get_Reco_EXTRAPOLATION_ExtrapolateLimit() const { return _RECO_EXTRAPOLATION_ExtrapolateLimit; };
    int Get_Reco_EXTRAPOLATION_NumBarsEnd() const { return _RECO_EXTRAPOLATION_NumBarsEnd; };
    int Get_Reco_EXTRAPOLATION_NumBarsStart() const { return _RECO_EXTRAPOLATION_NumBarsStart; };

    int Get_Reco_TRACKMATCH_PlaneLimit() const { return _RECO_TRACKMATCH_PlaneLimit; };
    int Get_Reco_TRACKMATCH_BarLimit() const { return _RECO_TRACKMATCH_BarLimit; };
    int Get_Reco_TRACKMATCH_TimeLimit() const { return _RECO_TRACKMATCH_TimeLimit; };
    int Get_Reco_TRACKMATCH_XTimeLimit() const { return _RECO_TRACKMATCH_XTimeLimit; };
    float Get_Reco_TRACKMATCH_YAnchor() const { return _RECO_TRACKMATCH_YAnchor; };
    double Get_Reco_TRACKMATCH_TiltAngle() const { return _RECO_TRACKMATCH_TiltAngle; };
    float Get_Reco_TRACKMATCH_YDifference() const { return _RECO_TRACKMATCH_YDifference; };
    int Get_Reco_TRACKMATCH_DirectionDistance() const { return _RECO_TRACKMATCH_DirectionDistance; };

    bool Get_Reco_ASTAR_IsGreedy() const { return _RECO_ASTAR_IsGreedy; };
    std::string Get_Reco_ASTAR_CostMetric() const { return _RECO_ASTAR_CostMetric; };

    int Get_Reco_STOPPING_nLastHits() const { return _RECO_STOPPING_nLastHits; };
    double Get_Reco_STOPPING_EnergyCut() const { return _RECO_STOPPING_EnergyCut; };

    std::string Get_Reco_TrackMethod() const { return _RECO_TRACK_METHOD; };
    bool Get_Reco_Clustering() const { return _RECO_CLUSTERING; };

    bool Get_Reco_Kalman_Run() const { return _RECO_KALMAN_RUN; };
    double Get_Reco_Kalman_Assumed_Charge() const { return _RECO_KALMAN_ASSUMED_CHARGE; };
    double Get_Reco_Kalman_Use_Outlier_Rejection() const { return _RECO_KALMAN_USE_OUTLIER_REJECTION; };
    double Get_Reco_Kalman_Outlier_Rejection_Chi2_Threshold() const { return _RECO_KALMAN_OUTLIER_REJECTION_CHI2_THRESHOLD; };
    // New: Outlier handling and adaptive gating
    double Get_Reco_Kalman_Outlier_EarlyStepMultiplier() const { return _RECO_KALMAN_OUTLIER_EARLYSTEP_MULTIPLIER; };
    int    Get_Reco_Kalman_Outlier_EarlyStepsCount() const     { return _RECO_KALMAN_OUTLIER_EARLYSTEPS_COUNT; };
    double Get_Reco_Kalman_Outlier_MaxScale() const            { return _RECO_KALMAN_OUTLIER_MAXSCALE; };
    double Get_Reco_Kalman_Outlier_InnovationTraceCoeff() const{ return _RECO_KALMAN_OUTLIER_INNOVATIONTRACE_COEFF; };
    int    Get_Reco_Kalman_Outlier_ResetThreshold() const      { return _RECO_KALMAN_OUTLIER_RESET_THRESHOLD; };
    double Get_Reco_Kalman_Outlier_NudgeAlpha() const          { return _RECO_KALMAN_OUTLIER_NUDGE_ALPHA; };
    bool   Get_Reco_Kalman_Outlier_NudgeSlopes() const         { return _RECO_KALMAN_OUTLIER_NUDGE_SLOPES; };
    double Get_Reco_Kalman_Outlier_NudgeSlopesAlpha() const    { return _RECO_KALMAN_OUTLIER_NUDGE_SLOPES_ALPHA; };
    // New: Absolute residual gating and early ramp
    bool   Get_Reco_Kalman_Outlier_UseAbsoluteResidual() const { return _RECO_KALMAN_OUTLIER_USE_ABS_RESID; };
    double Get_Reco_Kalman_Outlier_AbsoluteResidualMax() const { return _RECO_KALMAN_OUTLIER_ABS_RESID_MAX; };
    double Get_Reco_Kalman_Outlier_AbsoluteResidualMaxX() const { return _RECO_KALMAN_OUTLIER_ABS_RESID_MAX_X; };
    double Get_Reco_Kalman_Outlier_AbsoluteResidualMaxY() const { return _RECO_KALMAN_OUTLIER_ABS_RESID_MAX_Y; };
    bool   Get_Reco_Kalman_Outlier_EarlyRampEnable() const     { return _RECO_KALMAN_OUTLIER_EARLY_RAMP_ENABLE; };
    int    Get_Reco_Kalman_Outlier_EarlyRampSteps() const      { return _RECO_KALMAN_OUTLIER_EARLY_RAMP_STEPS; };

    // Augmentation tuning
    int    Get_Reco_Kalman_Augment_MaximumPlaneGap() const      { return _RECO_KALMAN_AUGMENT_MAXIMUM_PLANE_GAP; };
    int    Get_Reco_Kalman_Augment_FreezeSteps() const          { return _RECO_KALMAN_AUGMENT_FREEZE_STEPS; };
    double Get_Reco_Kalman_Augment_MaxDeflectionPerStep() const { return _RECO_KALMAN_AUGMENT_MAX_DEFLECTION_PER_STEP; };
    double Get_Reco_Kalman_Augment_MomentumWindowZ() const      { return _RECO_KALMAN_AUGMENT_MOMENTUM_WINDOW_Z; };
    double Get_Reco_Kalman_Augment_FinalPlaneRelaxFactor() const { return _RECO_KALMAN_AUGMENT_FINAL_PLANE_RELAX_FACTOR; };
    double Get_Reco_Kalman_Augment_FinalPlaneRelaxAbsX() const  { return _RECO_KALMAN_AUGMENT_FINAL_PLANE_RELAX_ABS_X; };
    double Get_Reco_Kalman_Augment_FinalPlaneRelaxAbsY() const  { return _RECO_KALMAN_AUGMENT_FINAL_PLANE_RELAX_ABS_Y; };
    


    bool Get_LightWeight_Truth() const { return _TRUTH_LIGHTWEIGHT; };

    bool Get_Reco_TRACKSMOOTHING_UseTrackSmoothing() const { return _RECO_TRACKSMOOTHING_UseTrackSmoothing; };
    std::string Get_Reco_TRACKSMOOTHING_TrackSmoothingStrategy() const { return _RECO_TRACKSMOOTHING_TrackSmoothingStrategy; };
    double Get_Reco_TRACKSMOOTHING_MaxYDistanceBetweenUVTransitionPoints() const { return _RECO_TRACKSMOOTHING_MaxYDistanceBetweenUVTransitionPoints; };
    double Get_Reco_TRACKSMOOTHING_UncertaintyGoodDirection() const { return _RECO_TRACKSMOOTHING_UncertaintyGoodDirection; };
    double Get_Reco_TRACKSMOOTHING_UncertaintyBadDirection() const { return _RECO_TRACKSMOOTHING_UncertaintyBadDirection; };
    double Get_Reco_TRACKSMOOTHING_UncertaintyForUVTransitionPoints() const { return _RECO_TRACKSMOOTHING_UncertaintyForUVTransitionPoints; };

    bool Get_Reco_TIME_RunTimeSlicer() const { return _RECO_TIME_RunTimeSlicer; };
    bool Get_Reco_TIME_RunSimpleTimeSlicer() const { return _RECO_TIME_RunSimpleTimeSlicer; };
    double Get_RECO_TIME_TimeSlicerThresholdStart() const { return _RECO_TIME_TimeSlicerThresholdStart; };
    double Get_RECO_TIME_TimeSlicerThresholdEnd() const { return _RECO_TIME_TimeSlicerThresholdEnd; };
    double Get_RECO_TIME_TimeSlicerSliceUnit() const { return _RECO_TIME_TimeSlicerSliceUnit; };
    int Get_RECO_TIME_TimeSlicerEnergyWindowInUnits() const { return _RECO_TIME_TimeSlicerEnergyWindowInUnits; };
    int Get_RECO_TIME_TimeSlicerMinimumSliceWidthInUnits() const { return _RECO_TIME_TimeSlicerMinimumSliceWidthInUnits; };
    double Get_RECO_TIME_TimeSlicerMaxTime() const { return _RECO_TIME_TimeSlicerMaxTime; };
    
    double Get_RECO_CALIBRATION_EnergyCalibration() const { return _RECO_CALIBRATION_EnergyCalibration; };
    
    double Get_FIDUCIAL_TMS_START_X() const { return _FIDUCIAL_TMS_START_X; };
    double Get_FIDUCIAL_TMS_START_Y() const { return _FIDUCIAL_TMS_START_Y; };
    double Get_FIDUCIAL_TMS_START_Z() const { return _FIDUCIAL_TMS_START_Z; };
    double Get_FIDUCIAL_TMS_END_X() const { return _FIDUCIAL_TMS_END_X; };
    double Get_FIDUCIAL_TMS_END_Y() const { return _FIDUCIAL_TMS_END_Y; };
    double Get_FIDUCIAL_TMS_END_Z() const { return _FIDUCIAL_TMS_END_Z; };
    double Get_ACTIVE_LAR_START_X() const { return _ACTIVE_LAR_START_X; };
    double Get_ACTIVE_LAR_START_Y() const { return _ACTIVE_LAR_START_Y; };
    double Get_ACTIVE_LAR_START_Z() const { return _ACTIVE_LAR_START_Z; };
    double Get_ACTIVE_LAR_END_X() const { return _ACTIVE_LAR_END_X; };
    double Get_ACTIVE_LAR_END_Y() const { return _ACTIVE_LAR_END_Y; };
    double Get_ACTIVE_LAR_END_Z() const { return _ACTIVE_LAR_END_Z; };

    double Get_LAR_FIDUCIAL_DOWNSTREAM_Z_CUT() const { return _LAR_FIDUCIAL_DOWNSTREAM_Z_CUT; };
    double Get_LAR_FIDUCIAL_XY_CUT() const { return _LAR_FIDUCIAL_XY_CUT; };
    double Get_LAR_OUTER_SHELL_THICKNESS() const { return _LAR_OUTER_SHELL_THICKNESS; };
    double Get_ND_PHYSICS_MUON_LAR_SHELL_CUT_ENERGY() const { return _ND_PHYSICS_MUON_LAR_SHELL_CUT_ENERGY; };

    bool Get_DrawPDF() const { return _APPLICATIONS_DrawPDF; };
    int Get_MaximumNEvents() const { return _APPLICATIONS_MaximumNEvents; };

    double Get_Geometry_YMIDDLE() const { return _GEOMETRY_YMIDDLE; };
    
    double Get_Nersc_Spill_Period() const { return _NERSC_SPILL_PERIOD; };
    void Set_Nersc_Spill_Period(double value) { _NERSC_SPILL_PERIOD = value; };
    
    std::string Get_GEOMETRY_GitTag() const { return _GEOMETRY_git_tag; };
    std::string Get_GEOMETRY_GitBranch() const { return _GEOMETRY_git_branch; };
    std::string Get_GEOMETRY_GitCommit() const { return _GEOMETRY_git_commit; };
    
    int Get_GEOMETRY_NumberOfScintillatorPlanes() const { return _GEOMETRY_NumberOfScintillatorPlanes; };
    int Get_GEOMETRY_NumberOfSteelPlatesThin() const { return _GEOMETRY_NumberOfSteelPlatesThin; };
    int Get_GEOMETRY_NumberOfSteelPlatesThick() const { return _GEOMETRY_NumberOfSteelPlatesThick; };
    int Get_GEOMETRY_NumberOfSteelPlatesDouble() const { return _GEOMETRY_NumberOfSteelPlatesDouble; };

    int Get_Meta_Version_Major() const { return _META_version_major; };
    int Get_Meta_Version_Minor() const { return _META_version_minor; };
    int Get_Meta_Version_Patch() const { return _META_version_patch; };

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
    bool _RECO_KALMAN_USE_OUTLIER_REJECTION;
    double _RECO_KALMAN_OUTLIER_REJECTION_CHI2_THRESHOLD;
    double _RECO_KALMAN_OUTLIER_EARLYSTEP_MULTIPLIER;
    int    _RECO_KALMAN_OUTLIER_EARLYSTEPS_COUNT;
    double _RECO_KALMAN_OUTLIER_MAXSCALE;
    double _RECO_KALMAN_OUTLIER_INNOVATIONTRACE_COEFF;
    int    _RECO_KALMAN_OUTLIER_RESET_THRESHOLD;
    double _RECO_KALMAN_OUTLIER_NUDGE_ALPHA;
    bool   _RECO_KALMAN_OUTLIER_NUDGE_SLOPES;
    double _RECO_KALMAN_OUTLIER_NUDGE_SLOPES_ALPHA;
    bool   _RECO_KALMAN_OUTLIER_USE_ABS_RESID;
    double _RECO_KALMAN_OUTLIER_ABS_RESID_MAX;
    double _RECO_KALMAN_OUTLIER_ABS_RESID_MAX_X;
    double _RECO_KALMAN_OUTLIER_ABS_RESID_MAX_Y;
    bool   _RECO_KALMAN_OUTLIER_EARLY_RAMP_ENABLE;
    int    _RECO_KALMAN_OUTLIER_EARLY_RAMP_STEPS;

    // Augmentation tuning (must be provided in TOML)
    int    _RECO_KALMAN_AUGMENT_MAXIMUM_PLANE_GAP;
    int    _RECO_KALMAN_AUGMENT_FREEZE_STEPS;
    double _RECO_KALMAN_AUGMENT_MAX_DEFLECTION_PER_STEP;
    double _RECO_KALMAN_AUGMENT_MOMENTUM_WINDOW_Z;
    double _RECO_KALMAN_AUGMENT_FINAL_PLANE_RELAX_FACTOR;
    double _RECO_KALMAN_AUGMENT_FINAL_PLANE_RELAX_ABS_X;
    double _RECO_KALMAN_AUGMENT_FINAL_PLANE_RELAX_ABS_Y;


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
    
    double _LAR_FIDUCIAL_DOWNSTREAM_Z_CUT;
    double _LAR_FIDUCIAL_XY_CUT;
    double _LAR_OUTER_SHELL_THICKNESS;
    double _ND_PHYSICS_MUON_LAR_SHELL_CUT_ENERGY;

    bool _APPLICATIONS_DrawPDF;
    int _APPLICATIONS_MaximumNEvents;

    double _GEOMETRY_YMIDDLE;
    
    std::string _GEOMETRY_git_tag;
    std::string _GEOMETRY_git_branch;
    std::string _GEOMETRY_git_commit;
    int _GEOMETRY_NumberOfScintillatorPlanes;
    int _GEOMETRY_NumberOfSteelPlatesThin;
    int _GEOMETRY_NumberOfSteelPlatesThick;
    int _GEOMETRY_NumberOfSteelPlatesDouble;
    int _META_version_major;
    int _META_version_minor;
    int _META_version_patch;
};

#endif
