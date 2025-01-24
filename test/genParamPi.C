#if !defined(__CINT__) || defined(__MAKECINT__)
#include "NA6PGenParam.h"
#include "NA6PBeamParam.h"
#include <TMath.h>
#endif

// example of genparam generator for pi-
NA6PGenerator* genParamPi()
{
  auto gen = new NA6PGenParam("pi-");
  gen->setPDGCode(-211);
  // gen->setNTracks(50);
  gen->setdNdY(50); // this will override eventual setNTracks

  gen->setTransIsPt(false);                           // acknowledge that the function is for mT
  gen->setFunTrans("x*TMath::Exp(-x*x/[0])", 0., 2.); // mt scaling, but we provide pT range
  std::vector<float> parTrans{0.16};
  gen->setParametersTrans(parTrans);
  //
  const auto& beamParam = NA6PBeamParam::Instance();
  const float ycm = beamParam.getYCM();
  gen->setLongIsY(true); //  acknowledge that the function is for Y
  gen->setFunLong("TMath::Exp( -TMath::Power( (x-[0])/[1], 2.))", 1., 4.);
  std::vector<float> parY{ycm, 1.0};
  gen->setParametersLong(parY);

  return gen;
}
