// NA6PCCopyright

#include "NA6PTarget.h"
#include "NA6PDetector.h"
#include "NA6PTGeoHelper.h"
#include "NA6PLayoutParam.h"

#include <TGeoVolume.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TColor.h>
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

  createMaterials();

  // Dimensions
  float boxMargin = 1.f;
  float boxDZ = std::abs(param.posTargetZ[param.nVerTelPlanes - 1] - param.posTargetZ[0]) + 2 * boxMargin;

  // Container
  auto* tgtBox = new TGeoTube("TgtBoxShape", 0, 2 * boxMargin, boxDZ / 2);
  TGeoVolume* tgtContainer = new TGeoVolume("TgtContainer", tgtBox, NA6PTGeoHelper::instance().getMedium("Air"));

  float zoffs = param.posTargetZ[0] + boxDZ / 2 - boxMargin; // offset to be added due to the placement of targets to the tgt box
  for (int i = 0; i < param.nTargets; ++i) {
    auto* targetShape = new TGeoTube(Form("TgtShape_%zu", i), 0, param.radTarget[i], param.thicknessTarget[i] / 2);
    TGeoVolume* target = new TGeoVolume(Form("TgtVol_%zu", i), targetShape, NA6PTGeoHelper::instance().getMedium(addName(param.medTarget[i])));
    auto* targetTransform = new TGeoCombiTrans(param.posTargetX[i], param.posTargetY[i], param.posTargetZ[i] + zoffs, NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0));
    target->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName(param.medTarget[i])));
    tgtContainer->AddNode(target, composeNonSensorVolID(i + 1), targetTransform);
  }

  auto* vtTransform = new TGeoCombiTrans(param.shiftTargets[0],
                                         param.shiftTargets[1],
                                         param.shiftTargets[2] + (param.posTargetZ[0] + boxDZ / 2 - boxMargin),
                                         NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0));
  world->AddNode(tgtContainer, composeNonSensorVolID(0), vtTransform);
}
