// NA6PCCopyright


#ifndef NA6P_GENCUT_PARAMS_H_
#define NA6P_GENCUT_PARAMS_H_

#include "ConfigurableParam.h"
#include "ConfigurableParamHelper.h"

// generic cuts for generator

struct NA6PGenCutParam : public na6p::conf::ConfigurableParamHelper<NA6PGenCutParam>
{
  size_t maxTrailsPerParticle = 1000000;
  float etaMin = 0.f;
  float etaMax = 6.f;
  float ptMin = 0.f;
  float ptMax = 4.f;
  float phiMin = 0.f;
  float phiMax = 6.2831853f;

  NA6PParamDef(NA6PGenCutParam, "gencut");
};

#endif
