// ROOT geometry creation translated from Geant4 materials and volumes
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
#include "TGeoTube.h"
#include "TGeoMatrix.h"
#include "TGeoCone.h"
#include "TGeoParaboloid.h"
#include "TGeoTrd1.h"
#include "TGeoTrd2.h"
#include "TGeoCompositeShape.h"
#include "TGeoBoolNode.h"
#include "TColor.h"
#include "TMath.h"
#include <cmath>
#include <unordered_map>
#include <string>

std::unordered_map<std::string, TGeoMaterial*> matPool;
std::unordered_map<std::string, TGeoMedium*> medPool;

void NA6PTGeoHelper::instance().addMedium(const std::string& medName, const std::string& matName = "")
{
  const auto& matN = matName.empty() ? medName : matName;
  if (medPool.find(medName) != medPool.end()) {
    throw std::runtime_error(Form("Medium %s was already created", medName.c_str()));
  }
  if (matPool.find(matN) == matPool.end()) {
    throw std::runtime_error(Form("Material %s does not exist", matN.c_str()));
  }
  medPool[medName] = new TGeoMedium(matN.c_str(), medPool.size(), matPool[matN]);
}

TGeoMedium* getMedium(const std::string& medName)
{
  auto m = medPool.find(medName);
  if (m == medPool.end()) {
    throw std::runtime_error(Form("Medium %s was not created", medName.c_str()));
  }
  return m->second;
}

TGeoRotation* NA6PTGeoHelper::rotAroundVector(float uX, float uY, float uZ, float ddelta)
{
  ddelta *= TMath::DegToRad();
  double sinDelta = std::sin(ddelta), cosDelta = std::cos(ddelta);
  double oneMinusCosDelta = 1.0 - cosDelta;
  double rmat[9] = {
    oneMinusCosDelta * uX * uX + cosDelta,
    oneMinusCosDelta * uX * uY - sinDelta * uZ,
    oneMinusCosDelta * uX * uZ + sinDelta * uY,
    oneMinusCosDelta * uY * uX + sinDelta * uZ,
    oneMinusCosDelta * uY * uY + cosDelta,
    oneMinusCosDelta * uY * uZ - sinDelta * uX,
    oneMinusCosDelta * uZ * uX - sinDelta * uY,
    oneMinusCosDelta * uZ * uY + sinDelta * uX,
    oneMinusCosDelta * uZ * uZ + cosDelta};
  auto rot = new TGeoRotation();
  rot->SetMatrix(rmat);
  return rot;
}

void CreateCommonMaterials()
{
  auto& matPool = NA6PTGeoHelper::instance().getMatPool();
  if (matPool.find("Air") == matPool.end()) {
    auto mixt = new TGeoMixture("Air", 2, 0.001);
    mixt->AddElement(new TGeoElement("N", "Nitrogen", 7, 14.01), 0.78);
    mixt->AddElement(new TGeoElement("O", "Oxygen", 8, 16.00), 0.22);
    matPool["Air"] = mixt;
    NA6PTGeoHelper::instance().addMedium("Air");
  }
}

void CreateDipoleVTMaterials()
{
  auto& matPool = NA6PTGeoHelper::instance().getMatPool();
  if (matPool.find("Iron") == matPool.end()) {
    matPool["Iron"] = new TGeoMaterial("Iron", 55.845, 26, 7.874); // Iron density
    NA6PTGeoHelper::instance().addMedium("Iron");
  }
  if (matPool.find("Copper") == matPool.end()) {
    matPool["Copper"] = new TGeoMaterial("Copper", 63.546, 29, 8.96); // Copper density
    NA6PTGeoHelper::instance().addMedium("Copper");
  }
}

void CreateTargetMaterials()
{
  auto& matPool = NA6PTGeoHelper::instance().getMatPool();
  if (matPool.find("Lead") == matPool.end()) {
    matPool["Lead"] = new TGeoMaterial("Lead", 207.19, 82, 11.35);
    NA6PTGeoHelper::instance().addMedium("Lead");
  }
}

void CreateVerTelMaterials()
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

void CreateTargetGeom(TGeoVolume* parent)
{
  CreateTargetMaterials();
  // Targets (Lead cylinders)
  float pbTargetDr0 = 0.3f;
  float pbTargetDr1 = 0.1f;
  float pbTargetDz = 0.15f;

  std::vector<float> targetZPositions = {-5.3, -5.15, -5.0, -4.85, -4.7};
  for (size_t i = 0; i < targetZPositions.size(); ++i) {
    float outerRadius = (i == 0) ? pbTargetDr0 : pbTargetDr1;
    auto* targetShape = new TGeoTube(Form("PbTargetShape_%zu", i), 0, outerRadius, pbTargetDz / 2);
    TGeoVolume* target = new TGeoVolume(Form("PbTarget_%zu", i), targetShape, NA6PTGeoHelper::instance().getMedium("Lead"));
    auto* targetTransform = new TGeoCombiTrans(0, targetZPositions[i], 0, NA6PTGeoHelper::rotAroundVector(1.0, 0.0, 0.0, 90.0));
    parent->AddNode(target, i, targetTransform);
  }
}

void CreateDipoleVTGeom(TGeoVolume* world)
{

  // Define materials for N6Magnets
  CreateDipoleVTMaterials();

  float posDipIP[3] = {0.f, 0.f, 0.f};

  // Define dimensions and positions
  float ironPoleR = 50.0f;    // cm
  float ironPoleR2 = 60.0f;   // cm
  float ironPoleH = 30.0f;    // cm
  float ironPoleYPos = 20.0f; // cm

  float copperCoilRmin = ironPoleR2; // cm
  float copperCoilRmax = 85.0f;      // cm
  float copperCoilH = 20.0f;         // cm

  float ironBoxX = 62.5f;  // cm
  float ironBoxY = 50.0f;  // cm
  float ironBoxZ = 120.0f; // cm

  float ironTrapL = 100.0f; // cm
  float ironTrapW = 80.0f;  // cm

  float magnetRoofL = 165.0f;  // cm
  float magnetRoofW1 = 80.0f;  // cm
  float magnetRoofW2 = 100.0f; // cm
  float magnetRoofH = 45.0f;   // cm

  float ironRoofY = 90.f;
  float ironRoofL = 330.f;
  float ironRoofCutH = -40.f;
  float ironRoofW = 2.0 * ironRoofL * (copperCoilRmax - ironBoxZ / 2) / (ironRoofL - copperCoilRmax) + ironBoxZ;

  // Iron Pole
  TGeoShape* ironPoleShape = new TGeoParaboloid("IronPole", ironPoleR, ironPoleR2, ironPoleH / 2);
  TGeoVolume* ironPole = new TGeoVolume("IronPole", ironPoleShape, NA6PTGeoHelper::instance().getMedium("Iron"));
  ironPole->SetLineColor(kGray);
  auto ironPoleComb = new TGeoCombiTrans(0, ironPoleYPos + ironPoleH / 2, 0, NA6PTGeoHelper::rotAroundVector(1.0f, 0.0f, 0.0f, 270.0f));

  // Copper Coil
  TGeoShape* copperCoilShape = new TGeoTube("CopperCoilShape", copperCoilRmin, copperCoilRmax, copperCoilH / 2);
  TGeoVolume* copperCoil = new TGeoVolume("CopperCoil", copperCoilShape, NA6PTGeoHelper::instance().getMedium("Copper"));
  copperCoil->SetLineColor(kRed - 10);
  auto copperCoilComb = new TGeoCombiTrans(0, ironPoleYPos + ironPoleH - copperCoilH / 2, 0, NA6PTGeoHelper::rotAroundVector(1.0f, 0.0f, 0.0f, 90.0f));

  // Iron Box
  TGeoShape* ironBoxShape = new TGeoBBox("IronBox", ironBoxX / 2, ironBoxY / 2, ironBoxZ / 2);
  TGeoVolume* ironBox = new TGeoVolume("IronBox", ironBoxShape, NA6PTGeoHelper::instance().getMedium("Iron"));
  ironBox->SetLineColor(kGray);
  auto ironBoxComb = new TGeoCombiTrans(copperCoilRmax + ironTrapL + ironBoxX / 2, ironBoxY / 2, 0., new TGeoRotation());

  // Iron Trapezoid
  TGeoShape* ironTrapShape = new TGeoTrd2("IronTrap", ironBoxZ / 2, ironTrapW / 2, ironBoxY / 2, ironBoxY / 2, ironTrapL / 2);
  TGeoVolume* ironTrap = new TGeoVolume("IronTrap", ironTrapShape, NA6PTGeoHelper::instance().getMedium("Iron"));
  ironTrap->SetLineColor(kGray);
  auto ironTrapComb = new TGeoCombiTrans(copperCoilRmax + ironTrapL / 2, ironBoxY / 2, 0, NA6PTGeoHelper::rotAroundVector(0.0f, -1.0f, 0.0f, 90.0f));

  // Magnet Roof (solidFeTrapRoof)
  TGeoShape* solidFeTrapRoof0 = new TGeoTrd2("solidFeTrapRoof0", ironRoofW / 2.f, ironBoxZ / 2., ironRoofY / 2., ironRoofY / 2., ironRoofL / 2.);
  TGeoShape* solidRoofCut0 = new TGeoTubeSeg("solidRoofCut0", copperCoilRmax, 2. * copperCoilRmax, 2. * ironRoofY, 0.0, 180.0);
  // First subtraction
  auto trans1 = new TGeoCombiTrans("", 0.0, 0.0, copperCoilRmax - ironRoofL / 2., NA6PTGeoHelper::rotAroundVector(-1., 0., 0., 90.f));
  auto sub1 = new TGeoSubtraction(solidFeTrapRoof0, solidRoofCut0, nullptr, trans1);
  TGeoCompositeShape* solidFeTrapRoof1 = new TGeoCompositeShape("solidFeTrapRoof1", sub1);
  // Second subtraction
  auto trans2 = new TGeoCombiTrans("", 0.0, ironRoofCutH, copperCoilRmax - ironRoofL / 2., NA6PTGeoHelper::rotAroundVector(-1., 0., 0., 30.f));
  auto sub2 = new TGeoSubtraction(solidFeTrapRoof1, solidRoofCut0, nullptr, trans2);
  TGeoCompositeShape* solidFeTrapRoof2 = new TGeoCompositeShape("solidFeTrapRoof2", sub2);
  // Third subtraction
  auto trans3 = new TGeoCombiTrans("", 0.0, ironRoofCutH, ironRoofL / 2. - copperCoilRmax, NA6PTGeoHelper::rotAroundVector(1., 0., 0., 45.f));
  auto sub3 = new TGeoSubtraction(solidFeTrapRoof2, solidRoofCut0, nullptr, trans3);
  TGeoCompositeShape* solidFeTrapRoof3 = new TGeoCompositeShape("solidFeTrapRoof3", sub3);
  // Fourth subtraction
  TGeoShape* solidRoofCut1 = new TGeoBBox("solidRoofCut1", 2 * ironRoofW, ironRoofW / 2., ironRoofW / 2.);
  auto rot4 = NA6PTGeoHelper::rotAroundVector(1., 0., 0., 45.f);
  rot4->RotateZ(36.f);
  auto trans4 = new TGeoCombiTrans("", 0.0, -ironRoofW - ironRoofY / 2, 0.0, rot4);
  auto sub4 = new TGeoSubtraction(solidFeTrapRoof3, solidRoofCut1, nullptr, trans4);
  TGeoCompositeShape* solidFeTrapRoof4 = new TGeoCompositeShape("solidFeTrapRoof4", sub4);
  // Fifth subtraction
  auto rot5 = NA6PTGeoHelper::rotAroundVector(1., 0., 0., 45.f);
  rot5->RotateZ(-36.f);
  auto trans5 = new TGeoCombiTrans("", 0.0, -ironRoofW - ironRoofY / 2, 0.0, rot5);
  auto sub5 = new TGeoSubtraction(solidFeTrapRoof4, solidRoofCut1, nullptr, trans5);
  TGeoCompositeShape* solidFeTrapRoof = new TGeoCompositeShape("solidFeTrapRoof", sub5);
  // Final magnet roof volume
  TGeoVolume* magnetRoof = new TGeoVolume("MagnetRoof", solidFeTrapRoof, NA6PTGeoHelper::instance().getMedium("Iron"));
  magnetRoof->SetLineColor(kGray);
  // Correct placement for magnet roof
  auto magnetRoofComb = new TGeoCombiTrans(ironRoofL / 2 - copperCoilRmax, ironBoxY + ironRoofY / 2, 0.0f, NA6PTGeoHelper::rotAroundVector(0.0f, 1.0f, 0.0f, 90.0f));

  // Half assembly for magnets
  TGeoVolumeAssembly* magnetHalfAssembly = new TGeoVolumeAssembly("MagnetHalfAssembly");

  // Place components into the half assembly
  magnetHalfAssembly->AddNode(ironPole, 1, ironPoleComb);
  magnetHalfAssembly->AddNode(copperCoil, 2, copperCoilComb);
  magnetHalfAssembly->AddNode(ironBox, 3, ironBoxComb);
  magnetHalfAssembly->AddNode(ironTrap, 4, ironTrapComb);
  magnetHalfAssembly->AddNode(magnetRoof, 5, magnetRoofComb);

  // Full assembly with two half assemblies
  TGeoVolumeAssembly* fMagnetAssembly = new TGeoVolumeAssembly("fMagnetAssembly");
  auto halfRot = new TGeoRotation();
  halfRot->RotateX(180);

  fMagnetAssembly->AddNode(magnetHalfAssembly, 1, new TGeoTranslation(0.0, 0.0, 0.0));
  fMagnetAssembly->AddNode(magnetHalfAssembly, 2, new TGeoCombiTrans(0.0, 0.0, 0.0, halfRot));

  // Add the full assembly to the world
  world->AddNode(fMagnetAssembly, 1, new TGeoTranslation(posDipIP[0], posDipIP[1], posDipIP[2]));
}

void CreateVerTelGeom(TGeoVolume* world)
{
  // Materials
  CreateVerTelMaterials();

  float posVTBox[3] = {0.f, 0.f, 0.f};

  // Dimensions
  float dipoleFieldR = 50.0f;
  float dipoleFieldY = 40.0f;

  float frameDX = 30.0f;
  float frameDY = 30.0f;
  float frameDz = 0.50f;
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

  std::vector<float> chipHoleX = {pixChipDX / 2 + dxyCut, -pixChipDX / 2, -pixChipDX / 2 - dxyCut, pixChipDX / 2};
  std::vector<float> chipHoleY = {pixChipDY / 2, pixChipDY / 2 + dxyCut, -pixChipDY / 2, -pixChipDY / 2 - dxyCut};

  std::vector<float> sensorZPositions = {7.1175, 15.1175, 20.1175, 25.1175, 38.1175};

  // IP Container
  auto* ipShape = new TGeoTube("IPContainer", 0, dipoleFieldR, dipoleFieldY / 2);
  TGeoVolume* ipContainer = new TGeoVolume("IPContainer", ipShape, NA6PTGeoHelper::instance().getMedium("Air"));
  auto* ipTransform = new TGeoCombiTrans(posVTBox[0], posVTBox[1], posVTBox[2], NA6PTGeoHelper::rotAroundVector(-1.0, 0.0, 0.0, 90.0));
  world->AddNode(ipContainer, 1, ipTransform);

  CreateTargetGeom(ipContainer);

  // pixel station Frame with holes (box with subtracted holes)
  auto* pixStFrameBox = new TGeoBBox("PixStFrameBox", frameDX / 2, frameDY / 2, frameDz / 2);
  auto* beamPipeHole = new TGeoTube("PixStFrameBoxBPHole", 0, frameHoleR, frameDz);
  auto* frameSubtraction = new TGeoSubtraction(pixStFrameBox, beamPipeHole);
  auto* pixStFrameShape = new TGeoCompositeShape("PixStFrameBoxHole0", frameSubtraction);
  for (size_t ii = 0; ii < chipHoleX.size(); ++ii) {
    auto* pixChipHole = new TGeoBBox("PixChipHole", pixChipHoleDX / 2, pixChipHoleDY / 2, frameDz);
    auto* holeTransform = new TGeoTranslation(chipHoleX[ii], chipHoleY[ii], 0);
    frameSubtraction = new TGeoSubtraction(pixStFrameShape, pixChipHole, nullptr, holeTransform);
    pixStFrameShape = new TGeoCompositeShape(Form("PixStFrameBoxHole0%zu", ii), frameSubtraction);
  }
  TGeoVolume* pixStFrame = new TGeoVolume("PixStFrame", pixStFrameShape, NA6PTGeoHelper::instance().getMedium("CarbonFoam"));
  pixStFrame->SetLineColor(kRed + 4);

  // Silicon Tracker Station
  auto* pixelStationShape = new TGeoBBox("PixelStationShape", pixChipContainerDX / 2, pixChipContainerDY / 2, pixChipContainerDz / 2);
  auto* sensorShape = new TGeoBBox("SensorShape", pixChipDX / 2, pixChipDY / 2, pixChipDz / 2);
  auto* pixelStationVol = new TGeoVolume("PixelStationVol", pixelStationShape, NA6PTGeoHelper::instance().getMedium("Air"));
  TGeoVolume* pixelSensor = new TGeoVolume("PixelSensor", sensorShape, NA6PTGeoHelper::instance().getMedium("Silicon"));
  pixelSensor->SetLineColor(kGray - 1);
  // place sensors
  std::vector<float> alpdx{pixChipDX / 2 + pixChipOffsX, -pixChipDX / 2 + pixChipOffsY, -pixChipDX / 2 - pixChipOffsX, pixChipDX / 2 - pixChipOffsY};
  std::vector<float> alpdy{pixChipDY / 2 - pixChipOffsY, pixChipDY / 2 + pixChipOffsX, -pixChipDY / 2 + pixChipOffsY, -pixChipDY / 2 - pixChipOffsX};

  for (size_t ii = 0; ii < alpdx.size(); ++ii) {
    auto* sensorTransform = new TGeoTranslation(alpdx[ii], alpdy[ii], 0);
    pixelStationVol->AddNode(pixelSensor, ii, sensorTransform);
  }
  // place frames + stations
  for (size_t ll = 0; ll < sensorZPositions.size(); ++ll) {
    auto* stationTransform = new TGeoCombiTrans(0, sensorZPositions[ll], 0, NA6PTGeoHelper::rotAroundVector(-1.0, 0.0, 0.0, 90.0));
    ipContainer->AddNode(pixelStationVol, ll, stationTransform);
    auto* frameTransform = new TGeoCombiTrans(0, sensorZPositions[ll] + 0.5 * (frameDz + pixChipDz), 0, NA6PTGeoHelper::rotAroundVector(-1.0, 0.0, 0.0, 90.0));
    ipContainer->AddNode(pixStFrame, ll, frameTransform);
  }
}

void CreateAbsorberMaterials()
{
  if (matPool.find("Iron") == matPool.end()) {
    matPool["Iron"] = new TGeoMaterial("Iron", 55.845, 26, 7.874);
    NA6PTGeoHelper::instance().addMedium("Iron");
  }
  if (matPool.find("BeO") == matPool.end()) {
    matPool["BeO"] = new TGeoMaterial("BeO", 9.012, 4, 3.01);
    NA6PTGeoHelper::instance().addMedium("BeO");
  }
  if (matPool.find("Graphite") == matPool.end()) {
    matPool["Graphite"] = new TGeoMaterial("Graphite", 12.01, 6, 2.267);
    NA6PTGeoHelper::instance().addMedium("Graphite");
  }
  if (matPool.find("Tungsten") == matPool.end()) {
    matPool["Tungsten"] = new TGeoMaterial("Tungsten", 183.84, 74, 19.3);
    NA6PTGeoHelper::instance().addMedium("Tungsten");
  }
}

void CreateMSMaterials()
{
  if (matPool.find("Scintillator") == matPool.end()) {
    matPool["Scintillator"] = new TGeoMaterial("Scintillator", 12.01, 6, 1.032);
    NA6PTGeoHelper::instance().addMedium("Scintillator");
  }
  if (matPool.find("Iron") == matPool.end()) {
    matPool["Iron"] = new TGeoMaterial("Iron", 55.845, 26, 7.874);
    NA6PTGeoHelper::instance().addMedium("Iron");
  }
  if (matPool.find("Silicon") == matPool.end()) {
    matPool["Silicon"] = new TGeoMaterial("Silicon", 28.09, 14, 2.33);
    NA6PTGeoHelper::instance().addMedium("Silicon");
  }
  if (matPool.find("BeO") == matPool.end()) {
    matPool["BeO"] = new TGeoMaterial("BeO", 9.012, 4, 3.01);
    NA6PTGeoHelper::instance().addMedium("BeO");
  }
  if (matPool.find("Graphite") == matPool.end()) {
    matPool["Graphite"] = new TGeoMaterial("Graphite", 12.01, 6, 2.267);
    NA6PTGeoHelper::instance().addMedium("Graphite");
  }
  if (matPool.find("Tungsten") == matPool.end()) {
    matPool["Tungsten"] = new TGeoMaterial("Tungsten", 183.84, 74, 19.3);
    NA6PTGeoHelper::instance().addMedium("Tungsten");
  }
  if (matPool.find("Air") == matPool.end()) {
    CreateCommonMaterials();
  }
}

void CreateAbsorberGeom(TGeoVolume* world)
{
  CreateAbsorberMaterials();

  // Dimensions and placements
  double AbsorberBeO_1_dx = 52.f;
  double AbsorberBeO_1_dy = 52.f;
  double AbsorberBeO_1_dz = 33.0f;
  double AbsorberBeO_1_dr = 2.47f;

  double AbsorberBeO_2_dx = 120.f;
  double AbsorberBeO_2_dy = 120.f;
  double AbsorberBeO_2_dz = 72.0f;
  double AbsorberBeO_2_dr = 5.48f;

  float AbsorberBeOH_1_pos = 71.5f;
  float AbsorberBeOH_2_pos = 124.0f;
  float AbsorberC_Wall_pos = 900.0f;
  float AbsorberC_1_dr = 11.f;
  float AbsorberC_1_dz_1X0 = 18.6f;
  double AbsorberC_1_dz = 130.2f;
  double AbsorberC_1_dx = 260.f;
  double AbsorberC_1_dy = 260.f;
  double AbsorberC_Wall_dz = 180.f;
  double AbsorberC_Wall_dr = 300.f;

  // Absorber BeO1
  auto* solidAbsorberBeO1 = new TGeoBBox("AbsorberBeO1", AbsorberBeO_1_dx / 2.0, AbsorberBeO_1_dy / 2.0, AbsorberBeO_1_dz / 2.0);
  auto* solidAbsorberHole1 = new TGeoTube("AbsorberHole1", 0.0, AbsorberBeO_1_dr, AbsorberBeO_1_dz / 2);
  auto* compositeBeO1 = new TGeoCompositeShape("AbsorberBeO1S", new TGeoSubtraction(solidAbsorberBeO1, solidAbsorberHole1));
  auto* logicAbsorberBeOH1 = new TGeoVolume("AbsorberBeO1Vol", compositeBeO1, NA6PTGeoHelper::instance().getMedium("BeO"));
  world->AddNode(logicAbsorberBeOH1, 1, new TGeoTranslation(0.0, 0.0, AbsorberBeOH_1_pos));
  logicAbsorberBeOH1->SetLineColor(kYellow);

  // Absorber BeO2
  auto* solidAbsorberBeO2 = new TGeoBBox("AbsorberBeO2", AbsorberBeO_2_dx / 2.0, AbsorberBeO_2_dy / 2.0, AbsorberBeO_2_dz / 2.0);
  auto* solidAbsorberHole2 = new TGeoTube("AbsorberHole2", 0.0, AbsorberBeO_2_dr, AbsorberBeO_2_dz / 2);
  auto* compositeBeO2 = new TGeoCompositeShape("AbsorberBeO2S", new TGeoSubtraction(solidAbsorberBeO2, solidAbsorberHole2));
  auto* logicAbsorberBeOH2 = new TGeoVolume("AbsorberBeO2Vol", compositeBeO2, NA6PTGeoHelper::instance().getMedium("BeO"));
  logicAbsorberBeOH2->SetLineColor(kOrange);
  world->AddNode(logicAbsorberBeOH2, 2, new TGeoTranslation(0.0, 0.0, AbsorberBeOH_2_pos));

  // Tungsten Plugs for BeO1 and BeO2
  auto* solidAbsorberPlugW1 = new TGeoTube("AbsorberPlugW1", 0.0, AbsorberBeO_1_dr, AbsorberBeO_1_dz / 2.0);
  auto* logicAbsorberPlugW1 = new TGeoVolume("AbsorberPlugW1Vol", solidAbsorberPlugW1, NA6PTGeoHelper::instance().getMedium("Tungsten"));
  logicAbsorberPlugW1->SetLineColor(kRed);
  world->AddNode(logicAbsorberPlugW1, 1, new TGeoTranslation(0.0, 0.0, AbsorberBeOH_1_pos));

  auto* solidAbsorberPlugW2 = new TGeoTube("AbsorberPlugW2", 0.0, AbsorberBeO_2_dr, AbsorberBeO_2_dz / 2);
  auto* logicAbsorberPlugW2 = new TGeoVolume("AbsorberPlugW2Vol", solidAbsorberPlugW2, NA6PTGeoHelper::instance().getMedium("Tungsten"));
  logicAbsorberPlugW2->SetLineColor(kMagenta);
  world->AddNode(logicAbsorberPlugW2, 2, new TGeoTranslation(0.0, 0.0, AbsorberBeOH_2_pos));

  // Graphite Absorber 1 Container
  auto* solidAbsorberC1_C = new TGeoBBox("AbsorberC1_C", AbsorberC_1_dx / 2.0, AbsorberC_1_dy / 2.0, AbsorberC_1_dz / 2.0);
  auto* logicAbsorberC1_C = new TGeoVolume("AbsorberC1_CVol", solidAbsorberC1_C, NA6PTGeoHelper::instance().getMedium("Air"));

  auto* solidAbsorberC1_1X0 = new TGeoBBox("AbsorberC1_1X0", AbsorberC_1_dx / 2.0, AbsorberC_1_dy / 2.0, AbsorberC_1_dz_1X0 / 2.0);
  auto* solidAbsorberC1_1X0Hole = new TGeoTube("AbsorberC1_1X0Hole", 0.0, AbsorberC_1_dr, AbsorberC_1_dz_1X0);
  auto* compositeC1_1X0 = new TGeoCompositeShape("AbsorberC1_1X0_S", new TGeoSubtraction(solidAbsorberC1_1X0, solidAbsorberC1_1X0Hole));
  auto* logicAbsorberC1_1X0 = new TGeoVolume("AbsorberC1_1X0Vol", compositeC1_1X0, NA6PTGeoHelper::instance().getMedium("Graphite"));

  auto* solidAbsorberPlugW_3_1X0 = new TGeoTube("AbsorberPlugW_3_1X0", 0.0, AbsorberC_1_dr, AbsorberC_1_dz_1X0 / 2.0);
  auto* logicAbsorberPlugW_3_1X0 = new TGeoVolume("AbsorberPlugW_3_1X0Vol", solidAbsorberPlugW_3_1X0, NA6PTGeoHelper::instance().getMedium("Tungsten"));
  logicAbsorberPlugW_3_1X0->SetLineColor(kRed);
  logicAbsorberC1_C->SetLineColor(kGray + 1);
  logicAbsorberC1_1X0->SetLineColor(kGreen);

  // Placement of 1X0 slices and plugs
  int nSlices = 7;
  for (int i = 0; i < nSlices; ++i) {
    double z = (i - 3) * AbsorberC_1_dz_1X0;
    logicAbsorberC1_C->AddNode(logicAbsorberC1_1X0, i, new TGeoTranslation(0.0, 0.0, z));
    logicAbsorberC1_C->AddNode(logicAbsorberPlugW_3_1X0, i, new TGeoTranslation(0.0, 0.0, z));
  }

  // Placement of Graphite Absorber Container
  double AbsorberC_1_pos = AbsorberBeOH_2_pos + AbsorberBeO_2_dz / 2.0 + AbsorberC_1_dz / 2.0;
  world->AddNode(logicAbsorberC1_C, 1, new TGeoTranslation(0.0, 0.0, AbsorberC_1_pos));

  // Absorber Wall
  auto* solidAbsorberC_Wall = new TGeoTube("AbsorberC_Wall", 0.0, AbsorberC_Wall_dr, AbsorberC_Wall_dz / 2.0);
  auto* logicAbsorberC_Wall = new TGeoVolume("AbsorberC_WallVol", solidAbsorberC_Wall, NA6PTGeoHelper::instance().getMedium("Graphite"));
  logicAbsorberC_Wall->SetLineColor(kBlue);
  world->AddNode(logicAbsorberC_Wall, 1, new TGeoTranslation(0.0, 0.0, AbsorberC_Wall_pos));
}
/*
void CreateMSGeom(TGeoVolume *world)
{
  // Create materials for N6MuTracker
  CreateMSMaterials();

  // Dimensions and placement values
}
*/

void CreateGeometry()
{
  // Initialize the geometry manager
  TGeoManager* geom = new TGeoManager("world", "ROOT Geometry from Geant4");
  CreateCommonMaterials();
  // Define the top-level world volume
  TGeoVolume* world = geom->MakeBox("World", NA6PTGeoHelper::instance().getMedium("Air"), 1000, 1000, 1000);
  geom->SetTopVolume(world);

  // Build subsystems
  CreateDipoleVTGeom(world);
  CreateVerTelGeom(world);
  CreateAbsorberGeom(world);
  //    CreateN6SiTracker(world);
  //    CreateN6TorMagnet(world);

  // Close the geometry
  geom->CloseGeometry();

  // Export to file
  geom->Export("geometry.root");
}

void test()
{
  CreateGeometry();
  gGeoManager->GetTopVolume()->Draw("ogl");
}
