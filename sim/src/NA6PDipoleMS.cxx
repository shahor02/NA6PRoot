// NA6PCCopyright
#include "NA6PDipoleMS.h"
#include "NA6PDetector.h"
#include "NA6PTGeoHelper.h"
#include "NA6PLayoutParam.h"

#include <TGeoVolume.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoParaboloid.h>
#include <TGeoTrd1.h>
#include <TGeoTrd2.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TColor.h>
#include <fairlogger/Logger.h>

void NA6PDipoleMS::createMaterials()
{
  auto& matPool = NA6PTGeoHelper::instance().getMatPool();
  std::string nameM;
  nameM = addName("Iron");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 55.845, 26, 7.874);
    NA6PTGeoHelper::instance().addMedium(nameM,"", kGray);
  }
  nameM = addName("Copper");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 63.546, 29, 8.96); // Copper density
    NA6PTGeoHelper::instance().addMedium(nameM, "" , kRed - 7);
  }
}

void NA6PDipoleMS::createGeometry(TGeoVolume *world)
{
  const auto& param = NA6PLayoutParam::Instance();

  createMaterials();

  // Define dimensions and positions
  float extDX = 440.f; // cm, external size
  float extDY = 400.f; // cm, external size
  float extDZ = 130.f; // cm, external size
  float appertureDX = 320.f; // cm, hole X apperture
  float appertureDY = 240.f; // cm, hole Y apperture

  // Iron box
  TGeoShape *ironBoxS = new TGeoBBox("IronBox", extDX/2, extDY/2, extDZ/2);
  // Apperture
  TGeoShape *apperBoxS = new TGeoBBox("ApperBox", appertureDX/2, appertureDY/2, extDZ);

  auto* sub = new TGeoSubtraction(ironBoxS, apperBoxS);
  auto* dipShape = new TGeoCompositeShape("DipoleMS", sub);
  auto* dipVol = new TGeoVolume("DipoleMS", dipShape, NA6PTGeoHelper::instance().getMedium(addName("Iron")));
  dipVol->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName("Iron")));

  // Add the full assembly to the world
  world->AddNode(dipVol, composeNonSensorVolID(0), new TGeoTranslation(param.posDipMS[0], param.posDipMS[1], param.posDipMS[2]));

}
