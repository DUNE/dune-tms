#include "TMS_Readout_Manager.h"

TMS_Readout_Manager::TMS_Readout_Manager() {

  if (!std::getenv("TMS_DIR")) {
    std::cerr << "Need ${TMS_DIR} environment set for reconstruction, please export TMS_DIR" << std::endl;
    throw; 
  }

  std::string filename;
  if (std::getenv("TMS_SIM_TOML") != NULL) filename = std::string(std::getenv("TMS_SIM_TOML"));
  else filename = std::string(std::getenv("TMS_DIR"))+"/config/TMS_Readout_Default_Config.toml";

  std::cout << "Creating TMS Readout Manager instance using TOML: " << filename << std::endl;

  // Read the TOML file
  const auto data = toml::parse(filename);
  
  // Parameters for simulating readout, used in TMS_Event::SimulateDeadtime().
  _SIM_READOUT_ReadoutTime = toml::find<double>(data, "Sim", "Readout", "ReadoutTime");
  _SIM_READOUT_Deadtime = toml::find<double>(data, "Sim", "Readout", "Deadtime");
  _SIM_READOUT_ZombieTime = toml::find<double>(data, "Sim", "Readout", "ZombieTime");
  _SIM_READOUT_ReadoutNoise = toml::find<double>(data, "Sim", "Readout", "ReadoutNoise");
  _SIM_READOUT_PedestalSubtractionThreshold = toml::find<double>(data, "Sim", "Readout", "PedestalSubtractionThreshold");
  
  // Parameters related to readout, used in TMS_Event::SimulateOpticalModel()
  _SIM_OPTICAL_ShouldSimulatePoisson = toml::find<bool>(data, "Sim", "Optical", "ShouldSimulatePoisson");
  _SIM_OPTICAL_BirksConstant = toml::find<double>(data, "Sim", "Optical", "BirksConstant");
  
  _SIM_OPTICAL_ShouldSimulateFiberLengths = toml::find<bool>(data, "Sim", "Optical", "ShouldSimulateFiberLengths");
  _SIM_OPTICAL_WSFAttenuationLength = toml::find<double>(data, "Sim", "Optical", "WSFAttenuationLength");
  _SIM_OPTICAL_WSFLengthMultiplier = toml::find<double>(data, "Sim", "Optical", "WSFLengthMultiplier");
  _SIM_OPTICAL_WSFEndReflectionEff = toml::find<double>(data, "Sim", "Optical", "WSFEndReflectionEff");
  _SIM_OPTICAL_AdditionalFiberLength = toml::find<double>(data, "Sim", "Optical", "AdditionalFiberLength");
  _SIM_OPTICAL_AdditionalFiberAttenuationLength = toml::find<double>(data, "Sim", "Optical", "AdditionalFiberAttenuationLength");
  _SIM_OPTICAL_AdditionalFiberCouplingEff = toml::find<double>(data, "Sim", "Optical", "AdditionalFiberCouplingEff");
  
  _SIM_OPTICAL_ReadoutCouplingEff = toml::find<double>(data, "Sim", "Optical", "ReadoutCouplingEff");
  _SIM_OPTICAL_LightYield = toml::find<double>(data, "Sim", "Optical", "LightYield");
  
  _SIM_NOISE_DarkNoiseRate = toml::find<double>(data, "Sim", "Noise", "DarkNoiseRate");
  _SIM_NOISE_DarkNoiseMinPE = toml::find<double>(data, "Sim", "Noise", "DarkNoiseMinPE");
  
  
  
  
}
