#include "TMS_Readout_Manager.h"

#include <cmath>
#include <stdexcept>

TMS_Readout_Manager::TMS_Readout_Manager() {

  if (!std::getenv("TMS_DIR")) {
    std::cerr << "Need ${TMS_DIR} environment set for reconstruction, please export TMS_DIR" << std::endl;
    throw std::runtime_error("TMS_DIR environment variable not set");
  }

  std::string filename;
  if (std::getenv("TMS_READOUT_TOML") != NULL) filename = std::string(std::getenv("TMS_READOUT_TOML"));
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
  // Keys added with the response-element pipeline are optional (defaults = shipped values) so older
  // custom readout configs keep working with the flag off
  _SIM_OPTICAL_PoissonAfterAttenuation = toml::find_or<bool>(data, "Sim", "Optical", "PoissonAfterAttenuation", false);
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

  _SIM_DETSIM_UseResponseElements = toml::find<bool>(data, "Sim", "DetSim", "UseResponseElements");
  _SIM_DETSIM_DepositBinLength = toml::find_or<double>(data, "Sim", "DetSim", "DepositBinLength", 1.0);
  _SIM_DETSIM_PassageMaxGap = toml::find_or<double>(data, "Sim", "DetSim", "PassageMaxGap", 0.5);
  // Optional with the shipped values, so older custom readout configs load with the flag off
  _SIM_DETSIM_FrontEndTimingMode = toml::find_or<bool>(data, "Sim", "DetSim", "FrontEndTimingMode", false);

  _SIM_FRONTEND_FastShaperPeakingTime = toml::find_or<double>(data, "Sim", "FrontEnd", "FastShaperPeakingTime", 15.0);
  _SIM_FRONTEND_FastShaperOrder = toml::find_or<int>(data, "Sim", "FrontEnd", "FastShaperOrder", 2);
  _SIM_FRONTEND_DiscriminatorThreshold = toml::find_or<double>(data, "Sim", "FrontEnd", "DiscriminatorThreshold", 1.7);
  _SIM_FRONTEND_DiscriminatorStopThreshold = toml::find_or<double>(data, "Sim", "FrontEnd", "DiscriminatorStopThreshold", -1.0);
  _SIM_FRONTEND_ToTMax = toml::find_or<double>(data, "Sim", "FrontEnd", "ToTMax", 255.5);
  _SIM_FRONTEND_TDCStep = toml::find_or<double>(data, "Sim", "FrontEnd", "TDCStep", 0.5);
  _SIM_FRONTEND_TimingJitter = toml::find_or<double>(data, "Sim", "FrontEnd", "TimingJitter", 0.1);
  _SIM_FRONTEND_SiPMCrosstalkProbability = toml::find_or<double>(data, "Sim", "FrontEnd", "SiPMCrosstalkProbability", 0.03);
  _SIM_FRONTEND_SiPMPixels = toml::find_or<int>(data, "Sim", "FrontEnd", "SiPMPixels", 1584);
  _SIM_FRONTEND_EnergyFromToT = toml::find_or<bool>(data, "Sim", "FrontEnd", "EnergyFromToT", false);

  // Settings that would make the timing-mode simulation loop forever or divide by zero fail here
  // with a clear message, instead of hanging on the first event
  if (_SIM_DETSIM_UseResponseElements && _SIM_DETSIM_FrontEndTimingMode) {
    if (_SIM_FRONTEND_FastShaperOrder <= 0) throw std::runtime_error("Sim.FrontEnd.FastShaperOrder must be positive (the pulse tail never decays otherwise)");
    if (!(_SIM_FRONTEND_FastShaperPeakingTime > 0)) throw std::runtime_error("Sim.FrontEnd.FastShaperPeakingTime must be positive");
    if (!(_SIM_FRONTEND_DiscriminatorThreshold > 0)) throw std::runtime_error("Sim.FrontEnd.DiscriminatorThreshold must be positive");
    // A positive stop threshold above the start threshold would end every pulse the moment it starts
    // (repeated crossings, zero ToT); <= 0 keeps the documented "same as the start threshold" fallback
    if (_SIM_FRONTEND_DiscriminatorStopThreshold > _SIM_FRONTEND_DiscriminatorThreshold) throw std::runtime_error("Sim.FrontEnd.DiscriminatorStopThreshold must not exceed DiscriminatorThreshold (use <= 0 for the same value)");
    // pixels * (1 - exp(-q / pixels)) blows up (or is undefined) for a non-positive pixel count
    if (_SIM_FRONTEND_SiPMPixels <= 0) throw std::runtime_error("Sim.FrontEnd.SiPMPixels must be positive (the saturation curve is undefined otherwise)");
    if (!(_SIM_FRONTEND_ToTMax > 0)) throw std::runtime_error("Sim.FrontEnd.ToTMax must be positive");
    if (!(_SIM_FRONTEND_TDCStep > 0)) throw std::runtime_error("Sim.FrontEnd.TDCStep must be positive");
    if (!(_SIM_FRONTEND_SiPMCrosstalkProbability >= 0 && _SIM_FRONTEND_SiPMCrosstalkProbability < 1)) {
      throw std::runtime_error("Sim.FrontEnd.SiPMCrosstalkProbability must be in [0, 1) (the crosstalk charge loop never ends otherwise)");
    }
  }

  _SIM_TIMING_ScintillatorDecayTime = toml::find_or<double>(data, "Sim", "Timing", "ScintillatorDecayTime", 3.0);
  _SIM_TIMING_WLSDecayTime = toml::find_or<double>(data, "Sim", "Timing", "WLSDecayTime", 20.0);
  _SIM_TIMING_FiberRefractiveIndex = toml::find_or<double>(data, "Sim", "Timing", "FiberRefractiveIndex", 1.5);
  _SIM_TIMING_ElectronicTimeNoise = toml::find_or<double>(data, "Sim", "Timing", "ElectronicTimeNoise", 1.0);
  _SIM_TIMING_MaxTimingPhotons = toml::find_or<int>(data, "Sim", "Timing", "MaxTimingPhotons", 300);

  // Settings that would divide by zero, draw from an undefined distribution or throw on the first
  // event fail here with a clear message instead. The Sim.Timing keys are used by the legacy timing
  // model too, so they are checked whichever pipeline runs; the deposit settings only matter for
  // the response-element pipeline.
  if (!(_SIM_TIMING_ScintillatorDecayTime > 0) || !std::isfinite(_SIM_TIMING_ScintillatorDecayTime)) throw std::runtime_error("Sim.Timing.ScintillatorDecayTime must be positive and finite");
  if (!(_SIM_TIMING_WLSDecayTime > 0) || !std::isfinite(_SIM_TIMING_WLSDecayTime)) throw std::runtime_error("Sim.Timing.WLSDecayTime must be positive and finite");
  if (!(_SIM_TIMING_FiberRefractiveIndex > 0) || !std::isfinite(_SIM_TIMING_FiberRefractiveIndex)) throw std::runtime_error("Sim.Timing.FiberRefractiveIndex must be positive and finite");
  if (!(_SIM_TIMING_ElectronicTimeNoise >= 0) || !std::isfinite(_SIM_TIMING_ElectronicTimeNoise)) throw std::runtime_error("Sim.Timing.ElectronicTimeNoise must be zero or positive and finite");
  if (_SIM_TIMING_MaxTimingPhotons < 1) throw std::runtime_error("Sim.Timing.MaxTimingPhotons must be at least 1");
  if (_SIM_DETSIM_UseResponseElements) {
    if (!(_SIM_DETSIM_DepositBinLength > 0) || !std::isfinite(_SIM_DETSIM_DepositBinLength)) throw std::runtime_error("Sim.DetSim.DepositBinLength must be positive and finite");
    if (!(_SIM_DETSIM_PassageMaxGap >= 0) || !std::isfinite(_SIM_DETSIM_PassageMaxGap)) throw std::runtime_error("Sim.DetSim.PassageMaxGap must be zero or positive");
  }
}
