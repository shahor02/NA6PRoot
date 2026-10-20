// NA6PCCopyright

#include "NA6PGenBeam.h"
#include "NA6PBeamParam.h"
#include "NA6PBeam.h"
#include "NA6PLayoutParam.h"
#include "NA6PMCStack.h"
#include <TParticle.h>

void NA6PGenBeam::init()
{
  NA6PGenerator::init();
  const auto& param = NA6PLayoutParam::Instance();
  if (param.nTargets < 1) {
    LOGP(warn, "No targets defined, the beam particle will be put at Z = {}", mZPut);
  } else {
    mZPut = param.shiftTargets[2] + param.posTargetZ[0] - param.thicknessTarget[0] / 2; // in the beginning of the upstream target
    LOGP(info, "The beam particle will be put in front of the 1st target, Z = {}", mZPut);
  }
}

void NA6PGenBeam::generatePrimaryVertex(int /*maxTrials*/)
{
  // just create the beam particle upstream from the target, define dummy vertex at that point and push the beam particle.
  const auto& beamPar = NA6PBeamParam::Instance();
  TParticle beam;
  beamPar.generate(beam, mZPut);
  setOrigin(beam.Vx(), beam.Vy(), beam.Vz());
  auto mcH = mStack->getEventHeader();
  mcH->setVX(beam.Vx());
  mcH->setVY(beam.Vy());
  mcH->setVZ(beam.Vz());
  mStack->setPVGenerated(true);
  int dummy = 0;
  getStack()->PushTrack(true, -1, beamPar.PDG, beam.Px(), beam.Py(), beam.Pz(), beam.Energy(),
                        beam.Vx(), beam.Vy(), beam.Vz(), 0., 0., 0., 0., TMCProcess::kPPrimary, dummy, 1., 0);
}

void NA6PGenBeam::generate()
{
  if (getStack()->isPVGenerated()) {
    LOGP(fatal, "NA6PGenBeam should not see externally generated vertex");
  }
  generatePrimaryVertex();
  auto mcHead = getStack()->getEventHeader();
  static std::string info = fmt::format("{}_pdg{}x{}", getName(), NA6PBeamParam::Instance().PDG, 1);
  mcHead->getGenHeaders().emplace_back(1, 0, mcHead->getNPrimaries(), 0, info);
  mcHead->incNPrimaries(1);
}
