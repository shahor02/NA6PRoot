// NA6PCCopyright

#ifndef NA6P_BEAM_PARAMS_H_
#define NA6P_BEAM_PARAMS_H_

#include "ConfigurableParam.h"
#include "ConfigurableParamHelper.h"

class NA6PBeam;
class TParticle;

struct NA6PBeamParam : public na6p::conf::ConfigurableParamHelper<NA6PBeamParam> {
  std::string particle = "Lead";
  float A = 208.f;
  float Z = 82.f;
  float PDG = 1000822080;
  float energyPerNucleon = 40.f; // GeV
  float meanX = 0.f;
  float meanY = 0.f;
  float meanSlopeX = 0.f;
  float meanSlopeY = 0.f;
  float sigX = 0.05f;
  float sigY = 0.05f;
  float sigSlopeX = 0.f;
  float sigSlopeY = 0.f;
  float corrX = 0.f; // correlation between sigX and sigSlopeX
  float corrY = 0.f; // correlation between sigY and sigSlopeY

  void generate(NA6PBeam& beam) const;
  void generate(TParticle& part, float z = 0.f) const;

  NA6PParamDef(NA6PBeamParam, "beam");
};

#endif
