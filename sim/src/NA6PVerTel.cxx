// NA6PCCopyright

#include "NA6PVerTel.h"
#include "NA6PDetector.h"
#include "NA6PTGeoHelper.h"
#include "NA6PLayoutParam.h"
#include "NA6PMCStack.h"

#include <TVirtualMC.h>
#include <TGeoVolume.h>
#include <TGeoNode.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoTorus.h>
#include <TGeoXtru.h>
#include <TGeoManager.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TColor.h>
#include <fairlogger/Logger.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

struct CoolingPipeParams {
  Double_t rInt;
  Double_t rExt;
  TGeoMedium* medSteel;
  TGeoMedium* medWater;
};

// Per-geometry-build index — reset before each call to BuildFrameWithPipes.
static Int_t gPipeNodeIdx = 0;

// Groove cut terms — reset before each call to BuildFrameWithPipes.
static std::vector<TString> gGrooveTerms;

static void AddStraightTube(TGeoVolume* assembly,
                            const CoolingPipeParams& p,
                            Double_t halfLen,
                            Double_t tx, Double_t ty,
                            Double_t phi = 90,
                            Double_t theta = 90,
                            Double_t pipeZ = 0.0,
                            Double_t frameHZ = 0.4,
                            Double_t fitClear = 0.01)
{
  Int_t idx = ++gPipeNodeIdx;
  TString tag = TString::Format("%04d", idx);

  auto* shSteel = new TGeoTube("shTubeS_" + tag, p.rInt, p.rExt, halfLen);
  auto* volSteel = new TGeoVolume("VolTubeS_" + tag, shSteel, p.medSteel);
  volSteel->SetLineColor(kGray + 2);

  auto* shWater = new TGeoTube("shTubeW_" + tag, 0, p.rInt, halfLen);
  auto* volWater = new TGeoVolume("VolTubeW_" + tag, shWater, p.medWater);
  volWater->SetLineColor(kBlue);

  auto* rot = new TGeoRotation("rotT_" + tag, phi, theta, 0);
  assembly->AddNode(volSteel, idx, new TGeoCombiTrans(tx, ty, pipeZ, rot));
  assembly->AddNode(volWater, idx, new TGeoCombiTrans(tx, ty, pipeZ, rot));

  auto* shFit = new TGeoTube("GrFit_" + tag, 0, p.rExt + fitClear, halfLen);
  auto* ctCutTube = new TGeoCombiTrans("GrFitT_" + tag, tx, ty, pipeZ, rot);
  ctCutTube->RegisterYourself();
  gGrooveTerms.push_back(TString::Format("%s:%s", shFit->GetName(), ctCutTube->GetName()));

  auto* shChan = new TGeoBBox("GrChan_" + tag,
                              p.rExt + fitClear,
                              frameHZ / 2. + 0.01,
                              halfLen + 0.01);
  auto* ctCutBox = new TGeoCombiTrans("GrChanT_" + tag,
                                      tx, ty, pipeZ + frameHZ / 2., rot);
  ctCutBox->RegisterYourself();
  gGrooveTerms.push_back(TString::Format("%s:%s", shChan->GetName(), ctCutBox->GetName()));
}

static void AddTorus(TGeoVolume* assembly,
                     const CoolingPipeParams& p,
                     Double_t rBend,
                     Double_t phiStart, Double_t phiRange,
                     Double_t tx, Double_t ty,
                     Double_t pipeZ = 0.0,
                     Double_t frameHZ = 0.4,
                     Double_t chanWidth = 0.125,
                     Double_t fitClear = 0.01)
{
  Int_t idx = ++gPipeNodeIdx;
  TString tag = TString::Format("%04d", idx);

  auto* shSteel = new TGeoTorus("shTorS_" + tag, rBend, p.rInt, p.rExt, phiStart, phiRange);
  auto* volSteel = new TGeoVolume("VolTorS_" + tag, shSteel, p.medSteel);
  volSteel->SetLineColor(kGray + 2);

  auto* shWater = new TGeoTorus("shTorW_" + tag, rBend, 0, p.rInt, phiStart, phiRange);
  auto* volWater = new TGeoVolume("VolTorW_" + tag, shWater, p.medWater);
  volWater->SetLineColor(kBlue);

  assembly->AddNode(volSteel, idx, new TGeoTranslation(tx, ty, pipeZ));
  assembly->AddNode(volWater, idx, new TGeoTranslation(tx, ty, pipeZ));

  auto* shFit = new TGeoTorus("GrFit_" + tag,
                              rBend,
                              0,
                              p.rExt + fitClear,
                              phiStart,
                              phiRange);
  auto* trFit = new TGeoTranslation("GrFitT_" + tag, tx, ty, pipeZ);
  trFit->RegisterYourself();
  gGrooveTerms.push_back(TString::Format("%s:%s", shFit->GetName(), trFit->GetName()));

  auto* shChan = new TGeoTubeSeg("GrChan_" + tag,
                                 rBend - chanWidth - fitClear,
                                 rBend + chanWidth + fitClear,
                                 frameHZ / 2. + 0.01,
                                 phiStart - 0.01,
                                 phiStart + phiRange + 0.01);
  auto* trChan = new TGeoTranslation("GrChanT_" + tag,
                                     tx, ty, pipeZ + frameHZ / 2.);
  trChan->RegisterYourself();
  gGrooveTerms.push_back(TString::Format("%s:%s", shChan->GetName(), trChan->GetName()));
}

static TGeoXtru* CreateRoundedRect(const char* name, double dx, double dy, double r, double dz)
{
  const int nSegments = 20;
  int nPoints = 4 * nSegments;
  double* x = new double[nPoints];
  double* y = new double[nPoints];
  for (int i = 0; i < nSegments; i++) {
    double phi = i * (TMath::Pi() / 2) / (nSegments - 1);
    x[i] = (dx - r) + r * std::cos(phi);
    y[i] = (dy - r) + r * std::sin(phi);
    x[i + nSegments] = -(dx - r) + r * std::cos(phi + TMath::Pi() / 2);
    y[i + nSegments] = (dy - r) + r * std::sin(phi + TMath::Pi() / 2);
    x[i + 2 * nSegments] = -(dx - r) + r * std::cos(phi + TMath::Pi());
    y[i + 2 * nSegments] = -(dy - r) + r * std::sin(phi + TMath::Pi());
    x[i + 3 * nSegments] = (dx - r) + r * std::cos(phi + 3 * TMath::Pi() / 2);
    y[i + 3 * nSegments] = -(dy - r) + r * std::sin(phi + 3 * TMath::Pi() / 2);
  }
  TGeoXtru* xtru = new TGeoXtru(2);
  xtru->SetName(name);
  xtru->DefinePolygon(nPoints, x, y);
  xtru->DefineSection(0, -dz);
  xtru->DefineSection(1, dz);
  delete[] x;
  delete[] y;
  return xtru;
}

static TGeoXtru* CreateStadium(const char* name, double halfLen, double halfH, double dz)
{
  const int nSemi = 30;
  int nPoints = 2 * nSemi;
  double* x = new double[nPoints];
  double* y = new double[nPoints];
  for (int i = 0; i < nSemi; i++) {
    double phi = -TMath::Pi() / 2 + i * TMath::Pi() / (nSemi - 1);
    x[i] = halfLen + halfH * std::cos(phi);
    y[i] = halfH * std::sin(phi);
  }
  for (int i = 0; i < nSemi; i++) {
    double phi = TMath::Pi() / 2 + i * TMath::Pi() / (nSemi - 1);
    x[nSemi + i] = -halfLen + halfH * std::cos(phi);
    y[nSemi + i] = halfH * std::sin(phi);
  }
  TGeoXtru* xtru = new TGeoXtru(2);
  xtru->SetName(name);
  xtru->DefinePolygon(nPoints, x, y);
  xtru->DefineSection(0, -dz);
  xtru->DefineSection(1, dz);
  delete[] x;
  delete[] y;
  return xtru;
}

// ============================================================
//  BuildFrameWithPipes
//
//  Returns a TGeoVolume assembly (air box) containing:
//    • the grooved Al frame composite shape
//    • the full cooling-pipe circuit (steel + water volumes)
// ============================================================
static TGeoVolume* BuildFrameWithPipes(const char* tag,
                                       const CoolingPipeParams& pipeP,
                                       TGeoMedium* medPlate,
                                       TGeoMedium* medAir,
                                       double frameHX, double frameHY,
                                       double frameHZ, double alCutHZ,
                                       double alCornerR,
                                       double alSlotR, double alSlotL)
{
  const Double_t rBend = 1.985;
  const Double_t pipeZ = 0.0;
  const Double_t fitClear = 0.01;
  const Double_t chanWidth = pipeP.rExt;

  const Double_t dzHor1 = 35.1 / 2.0;
  const Double_t dzVer1 = 24.2 / 2.0;
  const Double_t dzHor2 = 24.2 / 2.0;
  const Double_t dzVer2 = 15.7 / 2.0;
  const Double_t dzHor3 = 3.0 / 2.0;
  const Double_t dzVer3 = 8.7 / 2.0;
  const Double_t dzHor4 = 10.2 / 2.0;
  const Double_t dzVer4 = 10.2 / 2.0;
  const Double_t dzHor5 = 17.2 / 2.0;
  const Double_t dzVer5 = 1.8 / 2.0;
  const Double_t dzHor6 = 6.9 / 2.0;

  const Double_t tubeOffX = 10.9 / 2.0;
  const Double_t tubeOffY = 14.085;

  auto TX = [&](Double_t x0) { return x0 + tubeOffX; };
  auto TY = [&](Double_t y0) { return y0 + tubeOffY; };

  TString asmName = TString::Format("FrameAsm_%s", tag);
  TGeoVolume* assembly = new TGeoVolumeAssembly(asmName);
  assembly->SetVisibility(kFALSE);

  // ── Pipe layout ───────────────────────────────────────────────────
  AddStraightTube(assembly, pipeP, dzHor1,
                  TX(0), TY(0), 90, 90, pipeZ, frameHZ, fitClear);

  AddTorus(assembly, pipeP, rBend, 90, 90,
           TX(-dzHor1), TY(-rBend), pipeZ, frameHZ, chanWidth, fitClear);

  AddStraightTube(assembly, pipeP, dzVer1,
                  TX(-dzHor1 - rBend), TY(-rBend - dzVer1), 0, 90,
                  pipeZ, frameHZ, fitClear);

  AddTorus(assembly, pipeP, rBend, 180, 90,
           TX(-dzHor1), TY(-rBend - 2 * dzVer1), pipeZ, frameHZ, chanWidth, fitClear);

  AddStraightTube(assembly, pipeP, dzHor2,
                  TX(-(dzHor1 - dzHor2)), TY(-rBend - 2 * dzVer1 - rBend), 90, 90,
                  pipeZ, frameHZ, fitClear);

  AddTorus(assembly, pipeP, rBend, 270, 90,
           TX(2 * dzHor2 - dzHor1), TY(-rBend - 2 * dzVer1), pipeZ, frameHZ, chanWidth, fitClear);

  AddStraightTube(assembly, pipeP, dzVer2,
                  TX(2 * dzHor2 - dzHor1 + rBend), TY(-rBend - 2 * dzVer1 + dzVer2), 0, 90,
                  pipeZ, frameHZ, fitClear);

  AddTorus(assembly, pipeP, rBend, 0, 90,
           TX(2 * dzHor2 - dzHor1), TY(-rBend - 2 * dzVer1 + 2 * dzVer2),
           pipeZ, frameHZ, chanWidth, fitClear);

  AddStraightTube(assembly, pipeP, dzHor3,
                  TX(2 * dzHor2 - dzHor1 - dzHor3), TY(-2 * dzVer1 + 2 * dzVer2), 90, 90,
                  pipeZ, frameHZ, fitClear);

  AddTorus(assembly, pipeP, rBend, 90, 90,
           TX(2 * dzHor2 - dzHor1 - 2 * dzHor3), TY(-rBend - 2 * dzVer1 + 2 * dzVer2),
           pipeZ, frameHZ, chanWidth, fitClear);

  AddStraightTube(assembly, pipeP, dzVer3,
                  TX(2 * dzHor2 - dzHor1 - 2 * dzHor3 - rBend),
                  TY(-rBend - 2 * dzVer1 + 2 * dzVer2 - dzVer3), 0, 90,
                  pipeZ, frameHZ, fitClear);

  AddTorus(assembly, pipeP, rBend, 270, 90,
           TX(2 * dzHor2 - dzHor1 - 2 * dzHor3 - 2 * rBend),
           TY(-rBend - 2 * dzVer1 + 2 * dzVer2 - 2 * dzVer3),
           pipeZ, frameHZ, chanWidth, fitClear);

  AddStraightTube(assembly, pipeP, dzHor4,
                  TX(2 * dzHor2 - dzHor1 - 2 * dzHor3 - 2 * rBend - dzHor4),
                  TY(-2 * rBend - 2 * dzVer1 + 2 * dzVer2 - 2 * dzVer3), 90, 90,
                  pipeZ, frameHZ, fitClear);

  AddTorus(assembly, pipeP, rBend, 180, 90,
           TX(2 * dzHor2 - dzHor1 - 2 * dzHor3 - 2 * rBend - 2 * dzHor4),
           TY(-rBend - 2 * dzVer1 + 2 * dzVer2 - 2 * dzVer3),
           pipeZ, frameHZ, chanWidth, fitClear);

  AddStraightTube(assembly, pipeP, dzVer4,
                  TX(2 * dzHor2 - dzHor1 - 2 * dzHor3 - 3 * rBend - 2 * dzHor4),
                  TY(-rBend - 2 * dzVer1 + 2 * dzVer2 - 2 * dzVer3 + dzVer4), 0, 90,
                  pipeZ, frameHZ, fitClear);

  AddTorus(assembly, pipeP, rBend, 90, 90,
           TX(2 * dzHor2 - dzHor1 - 2 * dzHor3 - 2 * rBend - 2 * dzHor4),
           TY(-rBend - 2 * dzVer1 + 2 * dzVer2 - 2 * dzVer3 + 2 * dzVer4),
           pipeZ, frameHZ, chanWidth, fitClear);

  AddStraightTube(assembly, pipeP, dzHor5,
                  TX(2 * dzHor2 - dzHor1 - 2 * dzHor3 - 2 * rBend - 2 * dzHor4 + dzHor5),
                  TY(-2 * dzVer1 + 2 * dzVer2 - 2 * dzVer3 + 2 * dzVer4), 90, 90,
                  pipeZ, frameHZ, fitClear);

  AddTorus(assembly, pipeP, rBend, 270, 90,
           TX(2 * dzHor2 - dzHor1 - 2 * dzHor3 - 2 * rBend - 2 * dzHor4 + 2 * dzHor5),
           TY(-2 * dzVer1 + 2 * dzVer2 - 2 * dzVer3 + 2 * dzVer4 + rBend),
           pipeZ, frameHZ, chanWidth, fitClear);

  AddStraightTube(assembly, pipeP, dzVer5,
                  TX(2 * dzHor2 - dzHor1 - 2 * dzHor3 - 2 * rBend - 2 * dzHor4 + 2 * dzHor5 + rBend),
                  TY(-2 * dzVer1 + 2 * dzVer2 - 2 * dzVer3 + 2 * dzVer4 + rBend + dzVer5), 0, 90,
                  pipeZ, frameHZ, fitClear);

  AddTorus(assembly, pipeP, rBend, 90, 90,
           TX(2 * dzHor2 - dzHor1 - 2 * dzHor3 - 2 * rBend - 2 * dzHor4 + 2 * dzHor5 + 2 * rBend),
           TY(-2 * dzVer1 + 2 * dzVer2 - 2 * dzVer3 + 2 * dzVer4 + rBend + 2 * dzVer5),
           pipeZ, frameHZ, chanWidth, fitClear);

  AddStraightTube(assembly, pipeP, dzHor6,
                  TX(2 * dzHor2 - dzHor1 - 2 * dzHor3 - 2 * rBend - 2 * dzHor4 + 2 * dzHor5 + 2 * rBend + dzHor6),
                  TY(-2 * dzVer1 + 2 * dzVer2 - 2 * dzVer3 + 2 * dzVer4 + 2 * rBend + 2 * dzVer5), 90, 90,
                  pipeZ, frameHZ, fitClear);

  // ── Grooved Al frame composite shape ─────────────────────────────
  TString base = TString::Format("FrameBase_%s", tag);
  TString hCn = TString::Format("HoleC_%s", tag);
  TString hHn = TString::Format("HoleH_%s", tag);
  TString hVn = TString::Format("HoleV_%s", tag);
  TString tTopN = TString::Format("TTop_%s", tag);
  TString tBotN = TString::Format("TBot_%s", tag);
  TString tLN = TString::Format("TL_%s", tag);
  TString tRN = TString::Format("TR_%s", tag);

  new TGeoBBox(base, frameHX, frameHY, frameHZ);
  CreateRoundedRect(hCn, 6.5, 6.5, alCornerR, alCutHZ);
  CreateStadium(hHn, alSlotL, alSlotR, alCutHZ);
  CreateRoundedRect(hVn, 2.5, 5.5, alCornerR, alCutHZ);

  auto* tTop = new TGeoTranslation(tTopN, 0, 10.5, 0);
  tTop->RegisterYourself();
  auto* tBot = new TGeoTranslation(tBotN, 0, -10.5, 0);
  tBot->RegisterYourself();
  auto* tL = new TGeoTranslation(tLN, -10.5, -0.5, 0);
  tL->RegisterYourself();
  auto* tR = new TGeoTranslation(tRN, 10.5, -0.5, 0);
  tR->RegisterYourself();

  TString grooveExpr = "";
  for (const auto& term : gGrooveTerms) {
    grooveExpr += " + " + term;
  }

  TString alExpr = base + " - (" + hCn + " + " + hHn + ":" + tTopN + " + " + hHn + ":" + tBotN + " + " + hVn + ":" + tLN + " + " + hVn + ":" + tRN + grooveExpr + ")";

  TString shapeName = TString::Format("FrameShape_%s", tag);
  TString volName = TString::Format("Frame_%s", tag);

  auto* frameShape = new TGeoCompositeShape(shapeName, alExpr.Data());
  TGeoVolume* frame = new TGeoVolume(volName, frameShape, medPlate);
  frame->SetLineColor(kAzure - 9);

  assembly->AddNode(frame, 1, new TGeoTranslation(0, 0, 0));

  return assembly;
}

void NA6PVerTel::createMaterials()
{
  auto& matPool = NA6PTGeoHelper::instance().getMatPool();
  std::string nameM;
  nameM = addName("Silicon");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 28.09, 14, 2.33);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kCyan + 1);
  }
  nameM = addName("CarbonFoam");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 1, 0.5);
    mixt->AddElement(12.01, 6, 1.0);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM, "", kBlue - 6);
  }
  nameM = addName("CarbonFoamLight");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 1, 0.09);
    mixt->AddElement(12.01, 6, 1.0);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM, "", kBlue - 6);
  }
  nameM = addName("CarbonFiber");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 1, 1.91);
    mixt->AddElement(12.01, 6, 1.0);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM, "", kGray + 1);
  }
  nameM = addName("Air");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 2, 0.001);
    mixt->AddElement(new TGeoElement("N", "Nitrogen", 7, 14.01), 0.78);
    mixt->AddElement(new TGeoElement("O", "Oxygen", 8, 16.00), 0.22);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM);
  }
  nameM = addName("Aluminium");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 26.98, 13, 2.7);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kAzure - 9);
  }
  nameM = addName("Steel");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 55.84, 26, 7.87);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kGray + 2);
  }
  nameM = addName("Water");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 2, 1.0);
    mixt->AddElement(new TGeoElement("H", "Hydrogen", 1, 1.008), 2);
    mixt->AddElement(new TGeoElement("O", "Oxygen", 8, 16.00), 1);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM, "", kBlue);
  }
}

void NA6PVerTel::createGeometry(TGeoVolume* world)
{
  const auto& param = NA6PLayoutParam::Instance();

  createMaterials();

  // ── Pixel station dimensions ──────────────────────────────────────
  float pixChipContainerDX = 30.0f;
  float pixChipContainerDY = 30.0f;
  float pixChipContainerDz = 0.1f;
  float pixChipDz = 50e-4f;
  float pixChipOffsX = 0.29f;
  float pixChipOffsY = 0.31f;
  float carbonPlateDz = 400e-4f;

  // ── Frame dimensions ───────────────────────────────────────────
  const double frameHX = 19.0;
  const double frameHY = 16.0;
  const double frameHZ = param.coolingPlateThickness / 2.;
  const double alCutHZ = 0.5;
  const double alCornerR = 1.0;
  const double alSlotR = 2.5;
  const double alSlotL = 10.5;
  float pixChipDX = 13.5996f;
  const float pixChipDYFirstLayers = 5.8692f;
  const float pixChipDYOtherLayers = 13.6948f;
  const float frameGlueGap = 0.01f;

  // ── Cooling pipe parameters ───────────────────────────────────────
  CoolingPipeParams pipeP;
  pipeP.rInt = 0.100;
  pipeP.rExt = 0.125;
  pipeP.medSteel = NA6PTGeoHelper::instance().getMedium(addName("Steel"));
  pipeP.medWater = NA6PTGeoHelper::instance().getMedium(addName("Water"));

  TGeoMedium* medPlate = NA6PTGeoHelper::instance().getMedium(addName(param.medPlateVerTel));
  TGeoMedium* medAir = NA6PTGeoHelper::instance().getMedium(addName("Air"));

  // ── VT container ─────────────────────────────────────────────────
  float boxDZMargin = pixChipContainerDz + 0.5f;
  float boxDZ = param.posVerTelPlaneZ[param.nVerTelPlanes - 1] - param.posVerTelPlaneZ[0] + 2 * boxDZMargin + 2 * static_cast<float>(frameHZ);

  const float maxCarbonPlateDX = 2.f * frameHX;
  const float maxCarbonPlateDY = 2.f * frameHY;
  const float vtHalfX = TMath::Max(static_cast<float>(frameHX + 4.5), maxCarbonPlateDX / 2.f + 1.f);
  const float vtHalfY = TMath::Max(static_cast<float>(frameHY + 3.f), maxCarbonPlateDY / 2.f + 1.f);

  auto* vtShape = new TGeoBBox("VTContainer",
                               vtHalfX,
                               vtHalfY,
                               boxDZ / 2);
  TGeoVolume* vtContainer = new TGeoVolume("VTContainer", vtShape, medAir);

  auto* vtTransform = new TGeoCombiTrans(
    param.shiftVerTel[0],
    param.shiftVerTel[1],
    param.shiftVerTel[2] + (param.posVerTelPlaneZ[0] + boxDZ / 2 - boxDZMargin),
    NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0));
  world->AddNode(vtContainer, composeNonSensorVolID(0), vtTransform);

  // ── Per-plane loop ────────────────────────────────────────────────
  float zoffs = param.posVerTelPlaneZ[0] + boxDZ / 2 - boxDZMargin;
  float frameRelZ = static_cast<float>(pixChipDz / 2 + carbonPlateDz + frameHZ + frameGlueGap);

  for (int ll = 0; ll < param.nVerTelPlanes; ++ll) {

    // Per-layer chip size: first two planes use tall chips, the rest use standard
    float pixChipDYlayer = (ll < 2) ? pixChipDYFirstLayers : pixChipDYOtherLayers;
    float carbonPlateDXlayer = 2 * frameHX;
    float carbonPlateDYlayer = (ll < 2) ? 15.f : 2 * frameHY;
    float pixelStationHX = TMath::Max(pixChipContainerDX / 2.f, carbonPlateDXlayer / 2.f + 0.01f);
    float pixelStationHY = TMath::Max(pixChipContainerDY / 2.f, carbonPlateDYlayer / 2.f + 0.01f);

    // ── Pixel sensor station (per-plane) ──────────────────────────
    TString stationShName = TString::Format("PixelStationShape_Pl%d", ll);
    TString sensorShName = TString::Format("SensorShape_Pl%d", ll);
    TString stationVName = TString::Format("PixelStationVol_Pl%d", ll);
    TString sensorVName = TString::Format("PixelSensor_Pl%d", ll);

    auto* pixelStationShape = new TGeoBBox(stationShName,
                                           pixelStationHX,
                                           pixelStationHY,
                                           pixChipContainerDz / 2);
    auto* sensorShape = new TGeoBBox(sensorShName,
                                     pixChipDX / 2, pixChipDYlayer / 2, pixChipDz / 2);

    auto* pixelStationVol = new TGeoVolume(stationVName, pixelStationShape, medAir);
    TGeoVolume* pixelSensor = new TGeoVolume(sensorVName, sensorShape,
                                             NA6PTGeoHelper::instance().getMedium(addName("Silicon")));
    pixelSensor->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName("Silicon")));

    // ── Carbon-fiber plate with central beam hole (per-plane) ─────
    TString cpFullShName = TString::Format("CarbonPlateFullShape_Pl%d", ll);
    TString cpHoleName = TString::Format("CarbonPlateBeamHole_Pl%d", ll);
    TString cpCompShName = TString::Format("CarbonPlateWithHoleShape_Pl%d", ll);
    TString cpVolName = TString::Format("CarbonPlateWithHole_Pl%d", ll);

    auto* carbonplateFullShape = new TGeoBBox(cpFullShName,
                                              carbonPlateDXlayer / 2, carbonPlateDYlayer / 2, carbonPlateDz / 2);
    auto* beamHole = new TGeoBBox(cpHoleName, pixChipOffsX, pixChipOffsY, carbonPlateDz);
    auto* holeRemoval = new TGeoSubtraction(carbonplateFullShape, beamHole);
    auto* cbPlateWithHoleShape = new TGeoCompositeShape(cpCompShName, holeRemoval);
    TGeoVolume* cbPlate = new TGeoVolume(cpVolName, cbPlateWithHoleShape,
                                         NA6PTGeoHelper::instance().getMedium(addName("CarbonFiber")));
    cbPlate->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName("CarbonFiber")));

    // ── Place sensors + carbon plate into the station volume ──────
    std::vector<float> alpdx = {pixChipDX / 2 + pixChipOffsX, -pixChipDX / 2 + pixChipOffsX,
                                -pixChipDX / 2 - pixChipOffsX, pixChipDX / 2 - pixChipOffsX};
    std::vector<float> alpdy = {pixChipDYlayer / 2 - pixChipOffsY, pixChipDYlayer / 2 + pixChipOffsY,
                                -pixChipDYlayer / 2 + pixChipOffsY, -pixChipDYlayer / 2 - pixChipOffsY};
    for (size_t ii = 0; ii < alpdx.size(); ++ii) {
      pixelStationVol->AddNode(pixelSensor, composeSensorVolID(ii),
                               new TGeoTranslation(alpdx[ii], alpdy[ii], 0));
    }
    pixelStationVol->AddNode(cbPlate, composeNonSensorVolID(20),
                             new TGeoCombiTrans(0., 0.,
                                                pixChipDz / 2 + carbonPlateDz / 2,
                                                NA6PTGeoHelper::rotAroundVector(0, 0.0, 0.0, 0.0)));

    vtContainer->AddNode(pixelStationVol, composeNonSensorVolID(ll),
                         new TGeoCombiTrans(param.posVerTelPlaneX[ll],
                                            param.posVerTelPlaneY[ll],
                                            param.posVerTelPlaneZ[ll] - zoffs,
                                            NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0)));

    if (param.useCoolingPlate) {
      gPipeNodeIdx = 0;
      gGrooveTerms.clear();

      TString planeTag = TString::Format("Pl%d", ll);
      TGeoVolume* alAsm = BuildFrameWithPipes(planeTag.Data(),
                                              pipeP, medPlate, medAir,
                                              frameHX, frameHY,
                                              frameHZ, alCutHZ,
                                              alCornerR, alSlotR, alSlotL);

      vtContainer->AddNode(alAsm, composeNonSensorVolID(ll + 40),
                           new TGeoCombiTrans(param.posVerTelPlaneX[ll],
                                              param.posVerTelPlaneY[ll],
                                              param.posVerTelPlaneZ[ll] - zoffs + frameRelZ,
                                              NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0)));
    }
  }
}

void NA6PVerTel::setAlignableEntries()
{
  const auto& param = NA6PLayoutParam::Instance();
  int svolCnt = 0;
  for (int ll = 0; ll < param.nVerTelPlanes; ++ll) {
    for (int ii = 0; ii < 4; ++ii) {
      int id = getActiveID() * 100 + svolCnt;
      std::string nm = fmt::format("VT_Lr{}_Sens{}", ll, ii);
      std::string path = fmt::format("/World_1/VTContainer_{}/PixelStationVol_{}/PixelSensor_{}", composeNonSensorVolID(0), composeNonSensorVolID(ll), composeSensorVolID(ii));
      gGeoManager->SetAlignableEntry(nm.c_str(), path.c_str(), id);
      LOGP(info, "Adding {} {} as alignable sensor {}", nm, path, id);
      svolCnt++;
    }
  }
}

bool NA6PVerTel::stepManager(int volID)
{
  int sensID = NA6PModule::volID2SensID(volID);
  if (sensID < 0) {
    LOGP(fatal, "Non-sensor volID={} was provided to stepManager of {}", volID, getName());
  }
  auto mc = TVirtualMC::GetMC();
  bool startHit = false, stopHit = false;
  unsigned char status = 0;
  if (mc->IsTrackEntering()) {
    status |= NA6PBaseHit::kTrackEntering;
  }
  if (mc->IsTrackInside()) {
    status |= NA6PBaseHit::kTrackInside;
  }
  if (mc->IsTrackExiting()) {
    status |= NA6PBaseHit::kTrackExiting;
  }
  if (mc->IsTrackOut()) {
    status |= NA6PBaseHit::kTrackOut;
  }
  if (mc->IsTrackStop()) {
    status |= NA6PBaseHit::kTrackStopped;
  }
  if (mc->IsTrackAlive()) {
    status |= NA6PBaseHit::kTrackAlive;
  }

  if ((status & NA6PBaseHit::kTrackEntering) || (status & NA6PBaseHit::kTrackInside && !mTrackData.mHitStarted)) {
    startHit = true;
  } else if ((status & (NA6PBaseHit::kTrackExiting | NA6PBaseHit::kTrackOut | NA6PBaseHit::kTrackStopped))) {
    stopHit = true;
  }
  if (!startHit) {
    mTrackData.mEnergyLoss += mc->Edep();
  }
  if (!(startHit | stopHit)) {
    return false; // do nothing
  }
  auto stack = (NA6PMCStack*)mc->GetStack();
  if (startHit) {
    mTrackData.mEnergyLoss = 0.;
    mc->TrackMomentum(mTrackData.mMomentumStart);
    mc->TrackPosition(mTrackData.mPositionStart);
    mTrackData.mTrkStatusStart = status;
    mTrackData.mHitStarted = true;
  }
  if (stopHit) {
    TLorentzVector positionStop, momentumStop;
    mc->TrackMomentum(momentumStop);
    mc->TrackPosition(positionStop);
    int stationID(-1);
    mc->CurrentVolOffID(1, stationID);
    stationID = volID2NonSensID(stationID);

    int chipindex = NChipsPerStation * stationID + sensID;
    auto* p = addHit(stack->GetCurrentTrackNumber(), chipindex, mTrackData.mPositionStart.Vect(), positionStop.Vect(),
                     mTrackData.mMomentumStart.Vect(), momentumStop.Vect(), positionStop.T(),
                     mTrackData.mEnergyLoss, mTrackData.mTrkStatusStart, status);
    if (mVerbosity > 0) {
      LOGP(info, "{} Tr{} {}", getName(), stack->GetCurrentTrackNumber(), p->asString());
    }
    // register the points in TParticle
    stack->addHit(getActiveIDBit());
    return true;
  }
  return false;
}

NA6PVerTelHit* NA6PVerTel::addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos, const TVector3& startMom, const TVector3& endMom,
                                  float endTime, float eLoss, unsigned char startStatus, unsigned char endStatus)
{
  mHits.emplace_back(trackID, detID, startPos, endPos, startMom, endMom, endTime, eLoss, startStatus, endStatus);
  return &(mHits.back());
}

void NA6PVerTel::createHitsOutput(const std::string& outDir)
{
  auto nm = fmt::format("{}Hits{}.root", outDir, getName());
  mHitsFile = TFile::Open(nm.c_str(), "recreate");
  mHitsTree = new TTree(fmt::format("hits{}", getName()).c_str(), fmt::format("{} Hits", getName()).c_str());
  mHitsTree->Branch(getName().c_str(), &hHitsPtr);
  LOGP(info, "Will store {} hits in {}", getName(), nm);
}

void NA6PVerTel::closeHitsOutput()
{
  if (mHitsTree && mHitsFile) {
    mHitsFile->cd();
    mHitsTree->Write();
    delete mHitsTree;
    mHitsTree = 0;
    mHitsFile->Close();
    delete mHitsFile;
    mHitsFile = 0;
  }
}

void NA6PVerTel::writeHits(const std::vector<int>& remapping)
{
  int nh = mHits.size();
  for (int i = 0; i < nh; i++) {
    auto& h = mHits[i];
    if (remapping[h.getTrackID()] < 0) {
      LOGP(error, "Track {} hit {} in {} was not remapped!", h.getTrackID(), i, getName());
    }
    h.setTrackID(remapping[h.getTrackID()]);
  }
  if (mHitsTree) {
    mHitsTree->Fill();
  }
}
