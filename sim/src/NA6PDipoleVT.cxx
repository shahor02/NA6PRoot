// NA6PCCopyright
#include "NA6PDipoleVT.h"
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

void NA6PDipoleVT::createMaterials()
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
  nameM = addName("Air");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 2, 0.001);
    mixt->AddElement(new TGeoElement("N", "Nitrogen", 7, 14.01), 0.78);
    mixt->AddElement(new TGeoElement("O", "Oxygen", 8, 16.00), 0.22);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM);
  }
}

void NA6PDipoleVT::createGeometry(TGeoVolume *world)
{
  const auto& param = NA6PLayoutParam::Instance();

  createMaterials();

  // Define dimensions and positions
  float ironPoleR = 50.0f; // cm
  float ironPoleR2 = 60.0f; // cm
  float ironPoleH = 30.0f; // cm
  float ironPoleYPos = 20.0f; // cm

  float copperCoilRmin = ironPoleR2; // cm
  float copperCoilRmax = 85.0f; // cm
  float copperCoilH = 20.0f; // cm

  float ironBoxX = 62.5f; // cm
  float ironBoxY = 50.0f; // cm
  float ironBoxZ = 120.0f; // cm

  float ironTrapL = 100.0f; // cm
  float ironTrapW = 80.0f; // cm

  float magnetRoofL = 165.0f; // cm
  float magnetRoofW1 = 80.0f; // cm
  float magnetRoofW2 = 100.0f; // cm
  float magnetRoofH = 45.0f; // cm

  float ironRoofY = 90.f;
  float ironRoofL = 330.f;
  float ironRoofCutH = -40.f;
  float ironRoofW = 2.0*ironRoofL*(copperCoilRmax - ironBoxZ/2)/(ironRoofL-copperCoilRmax) + ironBoxZ; 

  // Iron Pole
  TGeoShape *ironPoleShape = new TGeoParaboloid("IronPole", ironPoleR, ironPoleR2, ironPoleH / 2);
  TGeoVolume *ironPole = new TGeoVolume("IronPole", ironPoleShape, NA6PTGeoHelper::instance().getMedium(addName("Iron")));
  ironPole->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName("Iron")));
  auto ironPoleComb = new TGeoCombiTrans(0, ironPoleYPos + ironPoleH/2, 0, NA6PTGeoHelper::rotAroundVector(1.0f, 0.0f, 0.0f, 270.0f));

  // Copper Coil
  TGeoShape *copperCoilShape = new TGeoTube("CopperCoilShape", copperCoilRmin, copperCoilRmax, copperCoilH / 2);
  TGeoVolume *copperCoil = new TGeoVolume("CopperCoil", copperCoilShape, NA6PTGeoHelper::instance().getMedium(addName("Copper")));
  copperCoil->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName("Copper")));
  auto copperCoilComb = new TGeoCombiTrans(0, ironPoleYPos + ironPoleH - copperCoilH/2, 0, NA6PTGeoHelper::rotAroundVector(1.0f, 0.0f, 0.0f, 90.0f));

  // Iron Box
  TGeoShape *ironBoxShape = new TGeoBBox("IronBox", ironBoxX/2, ironBoxY/2, ironBoxZ / 2);
  TGeoVolume *ironBox = new TGeoVolume("IronBox", ironBoxShape, NA6PTGeoHelper::instance().getMedium(addName("Iron")));
  ironBox->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName("Iron")));
  auto ironBoxComb = new TGeoCombiTrans(copperCoilRmax + ironTrapL + ironBoxX/2, ironBoxY/2, 0., new TGeoRotation());

  // Iron Trapezoid
  TGeoShape *ironTrapShape = new TGeoTrd2("IronTrap", ironBoxZ/2, ironTrapW/2, ironBoxY/2, ironBoxY/2, ironTrapL / 2);
  TGeoVolume *ironTrap = new TGeoVolume("IronTrap", ironTrapShape, NA6PTGeoHelper::instance().getMedium(addName("Iron")));
  ironTrap->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName("Iron")));
  auto ironTrapComb = new TGeoCombiTrans(copperCoilRmax + ironTrapL/2, ironBoxY/2, 0, NA6PTGeoHelper::rotAroundVector(0.0f, -1.0f, 0.0f, 90.0f));

  // Magnet Roof (solidFeTrapRoof)    
  TGeoShape *solidFeTrapRoof0 = new TGeoTrd2("solidFeTrapRoof0", ironRoofW/2.f, ironBoxZ/2., ironRoofY/2., ironRoofY/2., ironRoofL/2.);
  TGeoShape *solidRoofCut0 = new TGeoTubeSeg("solidRoofCut0", copperCoilRmax, 2.*copperCoilRmax, 2.*ironRoofY, 0.0, 180.0);
  // First subtraction
  auto trans1 = new TGeoCombiTrans("", 0.0, 0.0, copperCoilRmax - ironRoofL/2., NA6PTGeoHelper::rotAroundVector(-1., 0., 0., 90.f));
  auto sub1 = new TGeoSubtraction(solidFeTrapRoof0, solidRoofCut0, nullptr, trans1);
  TGeoCompositeShape *solidFeTrapRoof1 = new TGeoCompositeShape("solidFeTrapRoof1", sub1);
  // Second subtraction
  auto trans2 = new TGeoCombiTrans("", 0.0, ironRoofCutH, copperCoilRmax - ironRoofL/2., NA6PTGeoHelper::rotAroundVector(-1., 0., 0., 30.f));
  auto sub2 = new TGeoSubtraction(solidFeTrapRoof1, solidRoofCut0, nullptr, trans2);
  TGeoCompositeShape *solidFeTrapRoof2 = new TGeoCompositeShape("solidFeTrapRoof2", sub2);
  // Third subtraction
  auto trans3 = new TGeoCombiTrans("", 0.0, ironRoofCutH, ironRoofL/2. - copperCoilRmax, NA6PTGeoHelper::rotAroundVector(1., 0., 0., 45.f));
  auto sub3 = new TGeoSubtraction(solidFeTrapRoof2, solidRoofCut0, nullptr, trans3);
  TGeoCompositeShape *solidFeTrapRoof3 = new TGeoCompositeShape("solidFeTrapRoof3", sub3);
  // Fourth subtraction
  TGeoShape *solidRoofCut1 = new TGeoBBox("solidRoofCut1", 2*ironRoofW, ironRoofW/2., ironRoofW/2.);
  auto rot4 = NA6PTGeoHelper::rotAroundVector(1., 0., 0., 45.f);
  rot4->RotateZ(36.f);
  auto trans4 = new TGeoCombiTrans("", 0.0, -ironRoofW-ironRoofY/2, 0.0, rot4);
  auto sub4 = new TGeoSubtraction(solidFeTrapRoof3, solidRoofCut1, nullptr, trans4);
  TGeoCompositeShape *solidFeTrapRoof4 = new TGeoCompositeShape("solidFeTrapRoof4", sub4);
  // Fifth subtraction
  auto rot5 = NA6PTGeoHelper::rotAroundVector(1., 0., 0., 45.f);
  rot5->RotateZ(-36.f);
  auto trans5 = new TGeoCombiTrans("", 0.0, -ironRoofW-ironRoofY/2, 0.0, rot5);
  auto sub5 = new TGeoSubtraction(solidFeTrapRoof4, solidRoofCut1, nullptr, trans5);
  TGeoCompositeShape *solidFeTrapRoof = new TGeoCompositeShape("solidFeTrapRoof", sub5);
  // Final magnet roof volume
  TGeoVolume *magnetRoof = new TGeoVolume("MagnetRoof", solidFeTrapRoof, NA6PTGeoHelper::instance().getMedium(addName("Iron")));
  magnetRoof->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName("Iron")));    
  // Correct placement for magnet roof
  auto magnetRoofComb = new TGeoCombiTrans(ironRoofL/2 - copperCoilRmax, ironBoxY + ironRoofY/2, 0.0f, NA6PTGeoHelper::rotAroundVector(0.0f, 1.0f, 0.0f, 90.0f));

  // Half assembly for magnets
  TGeoVolumeAssembly *magnetHalfAssembly = new TGeoVolumeAssembly("MagnetHalfAssembly");

  // Place components into the half assembly
  magnetHalfAssembly->AddNode(ironPole, composeNonSensorVolID(1), ironPoleComb);
  magnetHalfAssembly->AddNode(copperCoil, composeNonSensorVolID(2), copperCoilComb);
  magnetHalfAssembly->AddNode(ironBox, composeNonSensorVolID(3), ironBoxComb);
  magnetHalfAssembly->AddNode(ironTrap, composeNonSensorVolID(4), ironTrapComb);
  magnetHalfAssembly->AddNode(magnetRoof, composeNonSensorVolID(5), magnetRoofComb);

  // Full assembly with two half assemblies
  TGeoVolumeAssembly *fMagnetAssembly = new TGeoVolumeAssembly("fMagnetAssembly");
  auto halfRot = new TGeoRotation();
  halfRot->RotateX(180);

  fMagnetAssembly->AddNode(magnetHalfAssembly, composeNonSensorVolID(6), new TGeoTranslation(0.0, 0.0, 0.0));
  fMagnetAssembly->AddNode(magnetHalfAssembly, composeNonSensorVolID(7), new TGeoCombiTrans(0.0, 0.0, 0.0, halfRot));

  // Add the full assembly to the world
  world->AddNode(fMagnetAssembly, composeNonSensorVolID(0), new TGeoTranslation(param.posDipIP[0], param.posDipIP[1], param.posDipIP[2]));
}
