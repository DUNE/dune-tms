#ifndef __TMS_READOUT_MANAGER_H__
#define __TMS_READOUT_MANAGER_H__

#include <iostream>
#include <string>
#include "toml.hpp"

// Just a global parameter manager
class TMS_Readout_Manager {

  public:
    static TMS_Readout_Manager& GetInstance() {
      static TMS_Readout_Manager Instance;
      return Instance;
    }

    void SetFileName(std::string file) { 
      std::cout << "Setting global edep-sim filename to " << file << std::endl;
      Filename = file; 
    };
    std::string GetFileName() { return Filename; };
    
    double Get_Sim_Readout_ReadoutTime() { return _SIM_READOUT_ReadoutTime; };
    double Get_Sim_Readout_Deadtime() { return _SIM_READOUT_Deadtime; };
    double Get_Sim_Readout_ZombieTime() { return _SIM_READOUT_ZombieTime; };
    double Get_Sim_Readout_ReadoutNoise() { return _SIM_READOUT_ReadoutNoise; };
    double Get_Sim_Readout_PedestalSubtractionThreshold() { return _SIM_READOUT_PedestalSubtractionThreshold; };
    
    bool Get_Sim_Optical_ShouldSimulatePoisson() { return _SIM_OPTICAL_ShouldSimulatePoisson; };
    double Get_Sim_Optical_BirksConstant() { return _SIM_OPTICAL_BirksConstant; };
    bool Get_Sim_Optical_ShouldSimulateFiberLengths() { return _SIM_OPTICAL_ShouldSimulateFiberLengths; };
    double Get_Sim_Optical_WSFAttenuationLength() { return _SIM_OPTICAL_WSFAttenuationLength; };
    double Get_Sim_Optical_WSFLengthMultiplier() { return _SIM_OPTICAL_WSFLengthMultiplier; };
    double Get_Sim_Optical_WSFEndReflectionEff() { return _SIM_OPTICAL_WSFEndReflectionEff; };
    double Get_Sim_Optical_AdditionalFiberLength() { return _SIM_OPTICAL_AdditionalFiberLength; };
    double Get_Sim_Optical_AdditionalFiberAttenuationLength() { return _SIM_OPTICAL_AdditionalFiberAttenuationLength; };
    double Get_Sim_Optical_AdditionalFiberCouplingEff() { return _SIM_OPTICAL_AdditionalFiberCouplingEff; };
    double Get_Sim_Optical_ReadoutCouplingEff() { return _SIM_OPTICAL_ReadoutCouplingEff; };
    double Get_Sim_Optical_LightYield() { return _SIM_OPTICAL_LightYield; };
    
    double Get_Sim_Noise_DarkNoiseRate() { return _SIM_NOISE_DarkNoiseRate; };
    double Get_Sim_Noise_DarkNoiseMinPE() { return _SIM_NOISE_DarkNoiseMinPE; };


  private:
    TMS_Readout_Manager();
    TMS_Readout_Manager(TMS_Readout_Manager const &) = delete;
    void operator=(TMS_Readout_Manager const &) = delete;
    ~TMS_Readout_Manager() {};

    std::string Filename;
  
  // Parameters for simulating readout, used in TMS_Event::SimulateDeadtime().
  double _SIM_READOUT_ReadoutTime;
  double _SIM_READOUT_Deadtime;
  double _SIM_READOUT_ZombieTime;
  double _SIM_READOUT_ReadoutNoise;
  double _SIM_READOUT_PedestalSubtractionThreshold;
  
  // Parameters related to optical model, used in TMS_Event::SimulateOpticalModel()
  bool _SIM_OPTICAL_ShouldSimulatePoisson;
  double _SIM_OPTICAL_BirksConstant;
  bool _SIM_OPTICAL_ShouldSimulateFiberLengths;
  double _SIM_OPTICAL_WSFAttenuationLength;
  double _SIM_OPTICAL_WSFLengthMultiplier;
  double _SIM_OPTICAL_WSFEndReflectionEff;
  double _SIM_OPTICAL_AdditionalFiberLength;
  double _SIM_OPTICAL_AdditionalFiberAttenuationLength;
  double _SIM_OPTICAL_AdditionalFiberCouplingEff;
  double _SIM_OPTICAL_ReadoutCouplingEff;
  double _SIM_OPTICAL_LightYield;
  
  double _SIM_NOISE_DarkNoiseRate;
  double _SIM_NOISE_DarkNoiseMinPE;
};

#endif
