// NA6PCCopyright

#include "NA6PBeamParam.h"
#include "MiscUtils.h"

O2ParamImpl(NA6PBeamParam);

#include "NA6PBeam.h"
#include "TParticle.h"

void NA6PBeamParam::generate(NA6PBeam& beam) const
{
  auto xsx = MiscUtils::genCorrelatedPair(sigX, sigSlopeX, corrX);
  auto ysy = MiscUtils::genCorrelatedPair(sigY, sigSlopeY, corrY);
  beam.setX(meanX + xsx.first);
  beam.setSlopeX(meanSlopeX + xsx.second);
  beam.setY(meanY + ysy.first);
  beam.setSlopeY(meanSlopeY + ysy.second);
}

void NA6PBeamParam::generate(TParticle& part, float z) const
{
  NA6PBeam beam;
  generate(beam);
  double etot = A*energyPerNucleon;
  double mass = Z*0.9383 + (A-Z)*0.9396;
  double ptot = std::sqrt(etot*etot - mass*mass);
  double normI = 1./std::sqrt(1. + beam.getSlopeX()*beam.getSlopeX() + beam.getSlopeY()*beam.getSlopeY());
  part.SetMomentum(ptot*beam.getSlopeX()*normI, ptot*beam.getSlopeY()*normI, ptot*normI, etot);
  part.SetProductionVertex(beam.getX(), beam.getY(), z, 0.);
}
