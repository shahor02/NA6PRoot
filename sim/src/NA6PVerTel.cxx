// NA6PCCopyright

#include "NA6PVerTel.h"
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

void NA6PVerTel::createMaterials()
{
  auto& matPool = NA6PTGeoHelper::instance().getMatPool();
  if (matPool.find("Silicon") == matPool.end()) {
    matPool["Silicon"] = new TGeoMaterial("Silicon", 28.09, 14, 2.33);
    NA6PTGeoHelper::instance().addMedium("Silicon");
  }
  if (matPool.find("CarbonFoam") == matPool.end()) {
    auto mixt = new TGeoMixture("CarbonFoam", 1, 0.5);
    mixt->AddElement(12.01, 6, 1.0); // Carbon-only mixture
    matPool["CarbonFoam"] = mixt;
    NA6PTGeoHelper::instance().addMedium("CarbonFoam");
  }
}

void NA6PVerTel::createGeometry(TGeoVolume *world)
{
  const auto& param = NA6PLayoutParam::Instance();
  
  createMaterials();

  // Dimensions

  float frameDX = 30.0f;
  float frameDY = 30.0f;
  float frameDZ = 0.50f;
  float frameHoleR = 0.424f;

  float pixChipHoleDX = 12.5f;
  float pixChipHoleDY = 13.0f;
  float dxyCut = 1.0f;

  float pixChipContainerDX = frameDX;
  float pixChipContainerDY = frameDY;
  float pixChipContainerDz = 0.05f;

  float pixChipDX = 14.69f;
  float pixChipDY = 14.69f;
  float pixChipDz = 50e-4f;
  float pixChipOffsX = 0.29f;
  float pixChipOffsY = 0.31f;

  float boxDZMargin = pixChipContainerDz + 0.5f;
  float boxDZ = param.posVerTelPlaneZ[param.nVerTelPlanes-1] - param.posVerTelPlaneZ[0] + 2*boxDZMargin;
  
  std::vector<float> chipHoleX = {pixChipDX/2 + dxyCut, -pixChipDX/2, -pixChipDX/2 - dxyCut, pixChipDX/2};
  std::vector<float> chipHoleY = {pixChipDY/2, pixChipDY/2 + dxyCut, -pixChipDY/2, -pixChipDY/2 - dxyCut};

  // Container
  auto *vtShape = new TGeoBBox("VTContainer", (frameDX+2.f)/2, (frameDY+2.f)/2, boxDZ/2);
  TGeoVolume *vtContainer = new TGeoVolume("VTContainer", vtShape, NA6PTGeoHelper::instance().getMedium("Air"));
  auto *vtTransform = new TGeoCombiTrans(param.shiftVerTel[0],
					 param.shiftVerTel[1],
					 param.shiftVerTel[2] +  (param.posVerTelPlaneZ[0] + boxDZ/2 - boxDZMargin),
					 NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0));
  world->AddNode(vtContainer, 1, vtTransform);

  // pixel station Frame with holes (box with subtracted holes)
  auto *pixStFrameBox = new TGeoBBox("PixStFrameBox", frameDX/2, frameDY/2, frameDZ/2);
  auto *beamPipeHole = new TGeoTube("PixStFrameBoxBPHole", 0, frameHoleR, frameDZ);
  auto *frameSubtraction = new TGeoSubtraction(pixStFrameBox, beamPipeHole);
  auto *pixStFrameShape = new TGeoCompositeShape("PixStFrameBoxHole0", frameSubtraction);
  for (size_t ii = 0; ii < chipHoleX.size(); ++ii) {
    auto *pixChipHole = new TGeoBBox("PixChipHole", pixChipHoleDX/2, pixChipHoleDY/2, frameDZ);
    auto *holeTransform = new TGeoTranslation(chipHoleX[ii], chipHoleY[ii], 0);
    frameSubtraction = new TGeoSubtraction(pixStFrameShape, pixChipHole, nullptr, holeTransform);
    pixStFrameShape = new TGeoCompositeShape(Form("PixStFrameBoxHole0%zu", ii), frameSubtraction);
  }
  TGeoVolume *pixStFrame = new TGeoVolume("PixStFrame", pixStFrameShape, NA6PTGeoHelper::instance().getMedium("CarbonFoam"));
  pixStFrame->SetLineColor(kRed+4);
      
  // Silicon Tracker Station
  auto *pixelStationShape = new TGeoBBox("PixelStationShape", pixChipContainerDX/2, pixChipContainerDY/2, pixChipContainerDz/2);
  auto *sensorShape = new TGeoBBox("SensorShape", pixChipDX/2, pixChipDY/2, pixChipDz/2);
  auto *pixelStationVol = new TGeoVolume("PixelStationVol", pixelStationShape, NA6PTGeoHelper::instance().getMedium("Air"));    
  TGeoVolume *pixelSensor = new TGeoVolume("PixelSensor", sensorShape, NA6PTGeoHelper::instance().getMedium("Silicon"));
  pixelSensor->SetLineColor(kGray-1);
  // place sensors
  std::vector<float> alpdx{pixChipDX/2+pixChipOffsX,  -pixChipDX/2+pixChipOffsY, -pixChipDX/2-pixChipOffsX, pixChipDX/2-pixChipOffsY};
  std::vector<float> alpdy{pixChipDY/2-pixChipOffsY, pixChipDY/2+pixChipOffsX, -pixChipDY/2+pixChipOffsY,  -pixChipDY/2-pixChipOffsX};

  for (size_t ii = 0; ii < alpdx.size(); ++ii) {
    auto *sensorTransform = new TGeoTranslation(alpdx[ii], alpdy[ii], 0);
    pixelStationVol->AddNode(pixelSensor, ii, sensorTransform);
  }
  // place frames + stations
  float zoffs = param.posVerTelPlaneZ[0] + boxDZ/2 - boxDZMargin; // offset to be added due to the placement of stations to the VT box
  for (int ll = 0; ll < param.nVerTelPlanes; ++ll) {
    auto *stationTransform = new TGeoCombiTrans(param.posVerTelPlaneX[ll], param.posVerTelPlaneY[ll], param.posVerTelPlaneZ[ll] - zoffs,
						NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0));    
    vtContainer->AddNode(pixelStationVol, ll, stationTransform);
    auto *frameTransform = new TGeoCombiTrans(param.posVerTelPlaneX[ll], param.posVerTelPlaneY[ll], param.posVerTelPlaneZ[ll]+ 0.5 * (frameDZ + pixChipDz) - zoffs,
					      NA6PTGeoHelper::rotAroundVector(0, 0.0, 0.0, 0.0));
    vtContainer->AddNode(pixStFrame, ll, frameTransform);
  } 
}

