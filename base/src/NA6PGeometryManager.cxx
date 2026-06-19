// NA6PCCopyright

#include <fairlogger/Logger.h>
#include <TGeoNode.h>
#include <TGeoPhysicalNode.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TFile.h>
#include "NA6PLayoutParam.h"
#include "NA6PGeometryManager.h"

bool NA6PGeometryManager::loadGeometry(const char* filename, const char* geoname)
{
  if (mGeoLoaded)
    return true;

  if (!gGeoManager) {
    TFile* f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
      LOGP(error, "Cannot open geometry file {}", filename);
      return false;
    }
    gGeoManager = (TGeoManager*)f->Get(geoname);
    if (!gGeoManager) {
      LOGP(error, "No geometry with name {} found in file {}", geoname, filename);
      f->Close();
      return false;
    }
    f->Close();
  }

  int nAlignMod = gGeoManager->GetNAlignable();
  std::vector<bool> isMatrixLoaded;
  bool doIter = false;
  if (nAlignMod == 0) {
    LOGP(info, "No alignable volumes in the geometry, resort to volume names for Vertex Telescope");
    const auto& param = NA6PLayoutParam::Instance();
    nAlignMod = param.nVerTelPlanes * kNVTModulesPerLayer;
    doIter = true;
  } else {
    LOGP(info, "Load geometry info for {} alignable volumes", nAlignMod);
  }
  isMatrixLoaded.assign(nAlignMod, false);
  mMatrices.resize(nAlignMod);
  mModuleHalfX.assign(nAlignMod, 0.);
  mModuleHalfY.assign(nAlignMod, 0.);

  if (doIter) {
    TGeoIterator next(gGeoManager->GetTopVolume());
    TGeoNode* node = nullptr;
    while ((node = next())) {
      TString name = node->GetName();
      if (name.BeginsWith("PixelSensor_")) {
        Int_t currentLevel = next.GetLevel();
        if (currentLevel <= 0)
          continue;
        TGeoNode* motherNode = next.GetNode(currentLevel - 1);
        int layer = motherNode->GetNumber() % 10;
        int modNum = node->GetNumber() % 10;
        int modIndex = kNVTModulesPerLayer * layer + modNum;
        if (modIndex < 0 || modIndex >= nAlignMod) {
          LOGP(error, "Wrong module index {} (layer={}, modNum={})", modIndex, layer, modNum);
          continue;
        }
        mMatrices[modIndex] = *(next.GetCurrentMatrix());
        bool sizeOk = fillModuleSize(modIndex, node->GetVolume());
        if (!sizeOk)
          continue;
        isMatrixLoaded[modIndex] = true;
      }
    }
  } else {
    for (int jMod = 0; jMod < nAlignMod; ++jMod) {
      TGeoPNEntry* entry = gGeoManager->GetAlignableEntry(jMod);
      if (!entry)
        continue;
      if (entry->GetGlobalOrig()) {
        mMatrices[jMod] = *(entry->GetGlobalOrig());
      } else {
        LOGP(error, "Module {} has a null global matrix", jMod);
        continue;
      }
      const char* path = entry->GetPath();
      if (!gGeoManager->cd(path)) {
        LOGP(error, "Failed to navigate to geometry path: {}", path);
        continue;
      }
      bool sizeOk = fillModuleSize(jMod, gGeoManager->GetCurrentVolume());
      if (!sizeOk)
        continue;
      isMatrixLoaded[jMod] = true;
    }
  }
  mGeoLoaded = true;
  for (int jMod = 0; jMod < nAlignMod; ++jMod) {
    if (!isMatrixLoaded[jMod]) {
      LOGP(error, "Matrix not loaded for module {}", jMod);
      mGeoLoaded = false;
    } else {
      LOGP(info, "Module {} Txyz {} {} {} size {} {}",
           jMod, mMatrices[jMod].GetTranslation()[0], mMatrices[jMod].GetTranslation()[1], mMatrices[jMod].GetTranslation()[2],
           mModuleHalfX[jMod] * 2, mModuleHalfY[jMod] * 2);
    }
  }
  return mGeoLoaded;
}

bool NA6PGeometryManager::fillModuleSize(int jMod, TGeoVolume* vol)
{
  if (!vol) {
    LOGP(error, "Null pointer to volume");
    return false;
  }
  TGeoShape* shape = vol->GetShape();
  if (!shape) {
    LOGP(error, "Null pointer to shape");
    return false;
  }
  if (shape->IsA() != TGeoBBox::Class()) {
    LOGP(error, "Module shape is not a TGeoBBox");
    return false;
  }
  TGeoBBox* box = static_cast<TGeoBBox*>(shape);
  mModuleHalfX[jMod] = static_cast<float>(box->GetDX());
  mModuleHalfY[jMod] = static_cast<float>(box->GetDY());
  return true;
}
