// NA6PCCopyright

#include <ranges>
#include <numeric>
#include <fairlogger/Logger.h>
#include <TGeoNode.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>

#include "NA6PLayoutParam.h"
#include "NA6PVerTelHit.h"
#include "NA6PVerTelDigitizer.h"

ClassImp(NA6PVerTelDigitizer)

  void NA6PVerTelDigitizer::init(const char* filename, const char* geoname)
{
  const auto& param = NA6PLayoutParam::Instance();
  mNumberOfModules = param.nVerTelPlanes * kNModulesPerLayer;
  mModules.resize(mNumberOfModules);
  mMatrices.resize(mNumberOfModules);

  if (!gGeoManager) {
    TFile* f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
      LOGP(error, "Cannot open geometry file {}", filename);
      return;
    }
    gGeoManager = (TGeoManager*)f->Get(geoname);
    if (!gGeoManager) {
      LOGP(error, "No geometry with name {} found in file {}", geoname, filename);
      f->Close();
      return;
    }
    f->Close();
  }

  std::vector<bool> isMatrixLoaded;
  isMatrixLoaded.assign(mNumberOfModules, false);
  mModuleHalfX.assign(mNumberOfModules, 0.);
  mModuleHalfY.assign(mNumberOfModules, 0.);
  TGeoIterator next(gGeoManager->GetTopVolume());
  TGeoNode* node = nullptr;
  while ((node = next())) {
    TString name = node->GetName();
    if (name.BeginsWith("PixelSensor_")) {
      Int_t currentLevel = next.GetLevel();
      if (currentLevel <= 0)
        continue;
      TGeoVolume* vol = node->GetVolume();
      if (!vol) {
        LOGP(error, "Volume null pointer");
        continue;
      }
      TGeoShape* shape = vol->GetShape();
      if (!shape) {
        LOGP(error, "Shape null pointer");
        continue;
      }
      if (shape->IsA() != TGeoBBox::Class()) {
        LOGP(error, "Module shape is not a TGeoBBox");
        continue;
      }
      TGeoNode* motherNode = next.GetNode(currentLevel - 1);
      int layer = motherNode->GetNumber() % 10;
      int modNum = node->GetNumber() % 10;
      int modIndex = kNModulesPerLayer * layer + modNum;
      if (modIndex < 0 || modIndex >= mNumberOfModules) {
        LOGP(error, "Wrong module index {} (layer={}, modNum={})", modIndex, layer, modNum);
        continue;
      }
      mMatrices[modIndex] = *(next.GetCurrentMatrix());
      TGeoBBox* box = static_cast<TGeoBBox*>(shape);
      mModuleHalfX[modIndex] = static_cast<float>(box->GetDX());
      mModuleHalfY[modIndex] = static_cast<float>(box->GetDY());
      isMatrixLoaded[modIndex] = true;
    }
  }
  for (int im = 0; im < mNumberOfModules; ++im) {
    if (!isMatrixLoaded[im]) {
      LOGP(error, "Matrix not loaded for module {}", im);
    }
  }
}

void NA6PVerTelDigitizer::process(const std::vector<NA6PVerTelHit>& hits, int layer)
{
  int nHits = hits.size();
  std::vector<int> hitIdx(nHits);
  std::iota(std::begin(hitIdx), std::end(hitIdx), 0);
  // sort hits to improve memory access
  std::sort(hitIdx.begin(), hitIdx.end(), [&hits](auto lhs, auto rhs) {
    return hits[lhs].getDetectorID() < hits[rhs].getDetectorID();
  });
  for (int i : hitIdx | std::views::filter([&](int idx) {
                 if (layer < 0)
                   return true;
                 return detID2Layer(hits[idx].getDetectorID()) == layer;
               })) {
    processHit(hits[i]);
  }
}

void NA6PVerTelDigitizer::processHit(NA6PVerTelHit hit)
{
  auto modID = hit.getDetectorID();
  printf("modID = %d\n", modID);
  auto& mod = mModules[modID];
  if (mod.isDisabled()) {
    LOGP(info, "Skipping disabled module {}", modID);
    return;
  }
  auto& matrix = mMatrices[modID];

  double xyzGloS[3] = {hit.getXIn(), hit.getYIn(), hit.getZIn()};
  double xyzGloE[3] = {hit.getXOut(), hit.getYOut(), hit.getZOut()};
  printf("Global coordinates = %f %f %f\n", xyzGloS[0], xyzGloS[1], xyzGloS[2]);

  double xyzLocS[3], xyzLocE[3];
  matrix.MasterToLocal(xyzGloS, xyzLocS);
  xyzLocS[0] += mModuleHalfX[modID];
  xyzLocS[1] += mModuleHalfY[modID];
  printf("Local coordinates: %f %f %f\n", xyzLocS[0], xyzLocS[1], xyzLocS[2]);
}
