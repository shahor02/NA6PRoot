// NA6PCCopyright

#include "NA6PTarget.h"
#include "NA6PDetector.h"
#include "NA6PTGeoHelper.h"
#include "NA6PLayoutParam.h"
#include "NA6PBeamParam.h"
#include "NA6PBeam.h"

#include <TGeoVolume.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TRandom.h>
#include <fairlogger/Logger.h>

void NA6PTarget::createMaterials()
{
  auto& matPool = NA6PTGeoHelper::instance().getMatPool();
  std::string nameM;
  nameM = addName("Lead");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 207.19, 82, 11.35);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kGray + 1);
  }
}

void NA6PTarget::createGeometry(TGeoVolume* world)
{
  const auto& param = NA6PLayoutParam::Instance();
  const auto& beamParam = NA6PBeamParam::Instance();
  const float Avogadro = 6.022e23;
  const float BarnToCM = 1.e-24;
  createMaterials();

  // Dimensions
  float boxMargin = 1.f;
  float boxDZ = std::abs(param.posTargetZ[param.nVerTelPlanes - 1] - param.posTargetZ[0]) + 2 * boxMargin;

  // Container
  auto* tgtBox = new TGeoTube("TgtBoxShape", 0, 2 * boxMargin, boxDZ / 2);
  TGeoVolume* tgtContainer = new TGeoVolume("TgtContainer", tgtBox, NA6PTGeoHelper::instance().getMedium("Air"));

  mTgtLambda.resize(param.nTargets);
  mTgtProb.resize(param.nTargets);
  double ltoL = 0;
  float zoffs = param.posTargetZ[0] + boxDZ / 2 - boxMargin; // offset to be added due to the placement of targets to the tgt box
  for (int i = 0; i < param.nTargets; ++i) {
    auto* targetShape = new TGeoTube(Form("TgtShape_%zu", i), 0, param.radTarget[i], param.thicknessTarget[i] / 2);
    const auto* med = NA6PTGeoHelper::instance().getMedium(addName(param.medTarget[i]));
    const auto* mat = med->GetMaterial();
    TGeoVolume* target = new TGeoVolume(Form("TgtVol_%zu", i), targetShape, med);
    auto* targetTransform = new TGeoCombiTrans(param.posTargetX[i], param.posTargetY[i], param.posTargetZ[i] + zoffs, NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0));
    target->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName(param.medTarget[i])));
    tgtContainer->AddNode(target, composeNonSensorVolID(i + 1), targetTransform);
    if (param.xsectionTarget[i] <= 0.) {
      LOGP(fatal, "The cross-section in barn for {} target {} is not defined, the beam is {}", param.medTarget[i], i, beamParam.particle);
    }
    mTgtLambda[i] = 1. / (param.xsectionTarget[i] * BarnToCM * Avogadro / mat->GetA() * mat->GetDensity());
    float l2Lambda = param.thicknessTarget[i] / mTgtLambda[i];
    mTgtProb[i] = 1. - std::exp(-l2Lambda);
    LOGP(info, "Added {:.3f} cm {} target: interaction probability {:.3f} for {} beam",
         param.thicknessTarget[i], param.medTarget[i], mTgtProb[i], beamParam.particle);
    ltoL += l2Lambda;
  }
  LOGP(info, "Total interaction probability on the targets is {:.3f}", 1. - std::exp(-ltoL));

  auto* vtTransform = new TGeoCombiTrans(param.shiftTargets[0],
                                         param.shiftTargets[1],
                                         param.shiftTargets[2] + (param.posTargetZ[0] + boxDZ / 2 - boxMargin),
                                         NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0));
  world->AddNode(tgtContainer, composeNonSensorVolID(0), vtTransform);
}

bool NA6PTarget::generateVertex(float& x, float& y, float& z, int maxTrials) const
{
  const auto& beamPar = NA6PBeamParam::Instance();
  const auto& layoutPar = NA6PLayoutParam::Instance();

  NA6PBeam beam;
  int ntrials = 0, ntrialsTouch = 0;

  while (ntrials < maxTrials) {
    beamPar.generate(beam);
    for (int it = 0; it < layoutPar.nTargets; it++) {
      // Calculate particle position at the Z position of the target center and check if the particle intersects the cylinder
      float xTgt = layoutPar.posTargetX[it] + layoutPar.shiftTargets[0];
      float yTgt = layoutPar.posTargetY[it] + layoutPar.shiftTargets[1];
      float zTgt = layoutPar.posTargetZ[it] + layoutPar.shiftTargets[2];
      float dX = beam.getX(zTgt) - xTgt;
      float dY = beam.getX(zTgt) - yTgt;
      float d2 = dX * dX + dY * dY;
      if (d2 <= layoutPar.radTarget[it] * layoutPar.radTarget[it]) {
        ntrialsTouch++;
        auto r = gRandom->Rndm();
        if (r < mTgtProb[it]) {
          z = zTgt - 0.5f * layoutPar.thicknessTarget[it] - mTgtLambda[it] * std::log(1. - r);
          x = beam.getX(z);
          y = beam.getY(z);
          LOGP(info, "generate on tgt {} : {} {} {} Ztgt = {}", it, x, y, z, zTgt);
          return true; // Interaction occurred, exit the function
        }
      }
    }
    ntrials++;
  }
  LOGP(fatal, "Failed to generate primary vertex after shooting {} beam particles, targets were traversed {} times, check beam and target settings",
       ntrials, ntrialsTouch);
  return false;
}
