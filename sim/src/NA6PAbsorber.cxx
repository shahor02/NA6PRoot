// NA6PCCopyright

#include "NA6PAbsorber.h"
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

void NA6PAbsorber::createMaterials()
{
  auto& matPool = NA6PTGeoHelper::instance().getMatPool();
  std::string nameM;
  nameM = addName("Iron");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 55.845, 26, 7.874);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kGray);
  }
  nameM = addName("BeO");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 9.012, 4, 3.01);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kYellow + 2);
  }
  nameM = addName("Graphite");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 12.01, 6, 2.267);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kGray + 2);
  }
  nameM = addName("Tungsten");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 183.84, 74, 19.3);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kBlack);
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

void NA6PAbsorber::createGeometry(TGeoVolume* world)
{
  const auto& param = NA6PLayoutParam::Instance();

  createMaterials();

  float zplace = param.posZStartAbsorber;
  if (param.nAbsorberSlices < param.nAbsorberSlicesWall) {
    LOGP(fatal, "Number of MuonWall slices {} exceeds total number of absorber slices {}", param.nAbsorberSlicesWall, param.nAbsorberSlices);
  }
  for (int isl = 0; isl < param.nAbsorberSlices; isl++) {
    if (isl == param.nAbsorberSlices - param.nAbsorberSlicesWall) { // special setting for last slice (muon wall)
      zplace = param.posMuonWall;
    }
    TGeoShape* slice = nullptr;
    auto slnm = fmt::format("AbsSl{}{}", isl, param.medAbsorber[isl]);
    if (param.dimYAbsorber[isl] > 0) {
      slice = new TGeoBBox((slnm + "SH").c_str(), param.dimXAbsorber[isl] / 2, param.dimYAbsorber[isl] / 2, param.thicknessAbsorber[isl] / 2);
    } else {
      slice = new TGeoTube((slnm + "SH").c_str(), 0.0, param.dimXAbsorber[isl] / 2, param.thicknessAbsorber[isl] / 2);
    }
    if (param.radPlug[isl] > 0 && !param.medAbsorberPlug[isl].empty()) {
      auto hole = new TGeoTube((slnm + "_H").c_str(), 0.0, param.radPlug[isl], param.thicknessAbsorber[isl]);
      slice = new TGeoCompositeShape((slnm + "_HS").c_str(), new TGeoSubtraction(slice, hole)); // slice with plug hole
      auto slnmP = fmt::format("AbsSlPlug{}{}", isl, param.medAbsorber[isl]);
      auto* plug = new TGeoTube((slnmP + "_SH").c_str(), param.radPlugInt[isl], param.radPlug[isl], param.thicknessAbsorber[isl] / 2);
      auto* plugVol = new TGeoVolume(slnmP.c_str(), plug, NA6PTGeoHelper::instance().getMedium(addName(param.medAbsorberPlug[isl])));
      plugVol->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName(param.medAbsorberPlug[isl])));
      world->AddNode(plugVol, composeNonSensorVolID(isl + 50), new TGeoTranslation(0.0, 0.0, zplace + param.thicknessAbsorber[isl] / 2));
    }
    auto sliceVol = new TGeoVolume(slnm.c_str(), slice, NA6PTGeoHelper::instance().getMedium(addName(param.medAbsorber[isl])));
    sliceVol->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName(param.medAbsorber[isl])));
    world->AddNode(sliceVol, composeNonSensorVolID(isl), new TGeoTranslation(0.0, 0.0, zplace + param.thicknessAbsorber[isl] / 2));
    zplace += param.thicknessAbsorber[isl];
  }
}
