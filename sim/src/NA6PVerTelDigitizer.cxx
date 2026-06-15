// NA6PCCopyright

#include <ranges>
#include <numeric>
#include <fairlogger/Logger.h>
#include <TGeoNode.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TFile.h>
#include <TTree.h>

#include "NA6PLayoutParam.h"
#include "NA6PVerTelHit.h"
#include "NA6PVerTelDigitizer.h"

void NA6PVerTelDigitizer::init(const char* filename, const char* geoname)
{
  const auto& param = NA6PLayoutParam::Instance();
  mNumberOfModules = param.nVerTelPlanes * kNModulesPerLayer;
  mModules.resize(mNumberOfModules);
  mMatrices.resize(mNumberOfModules);
  mThresholds.assign(mNumberOfModules * NA6PVerTelSegmentation::NXTiles * NA6PVerTelSegmentation::NYSensors, kDefaultThresholdEl);

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
  createDigitsOutput();
}

void NA6PVerTelDigitizer::createDigitsOutput()
{
  auto nm = fmt::format("Digits{}.root", getName());
  mDigitsFile = TFile::Open(nm.c_str(), "recreate");
  mDigitsTree = new TTree(fmt::format("digits{}", getName()).c_str(), fmt::format("{} Digits", getName()).c_str());
  mDigitsTree->Branch(getName().c_str(), &hDigitsPtr);
  LOGP(info, "Will store {} hits in {}", getName(), nm);
}

void NA6PVerTelDigitizer::closeDigitsOutput()
{
  if (mDigitsTree && mDigitsFile) {
    mDigitsFile->cd();
    mDigitsTree->Write();
    delete mDigitsTree;
    mDigitsTree = nullptr;
    mDigitsFile->Close();
    delete mDigitsFile;
    mDigitsFile = nullptr;
  }
}

void NA6PVerTelDigitizer::writeDigits()
{
  if (mDigitsTree) {
    mDigitsTree->Fill();
  }
  LOGP(info, "Saved {} digits in tree with {} entries", mDigits.size(), mDigitsTree->GetEntries());
}

void NA6PVerTelDigitizer::process(const std::vector<NA6PVerTelHit>& hits, int layer)
{
  clearDigits();
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
  finalizeDigits();
  writeDigits();
}

void NA6PVerTelDigitizer::getHitLocalCoord(NA6PVerTelHit hit, double xyzLocS[3], double xyzLocE[3])
{
  auto modID = hit.getDetectorID();
  auto& matrix = mMatrices[modID];

  double xyzGloS[3] = {hit.getXIn(), hit.getYIn(), hit.getZIn()};
  double xyzGloE[3] = {hit.getXOut(), hit.getYOut(), hit.getZOut()};

  matrix.MasterToLocal(xyzGloS, xyzLocS);
  xyzLocS[0] += mModuleHalfX[modID];
  xyzLocS[1] += mModuleHalfY[modID];
  matrix.MasterToLocal(xyzGloE, xyzLocE);
  xyzLocE[0] += mModuleHalfX[modID];
  xyzLocE[1] += mModuleHalfY[modID];
  if (modID % kNModulesPerLayer == 1) {
    // swap x
    xyzLocS[0] = mModuleHalfX[modID] * 2. - xyzLocS[0];
    xyzLocE[0] = mModuleHalfX[modID] * 2. - xyzLocE[0];
  } else if (modID % kNModulesPerLayer == 2) {
    // swap x and y
    xyzLocS[0] = mModuleHalfX[modID] * 2. - xyzLocS[0];
    xyzLocS[1] = mModuleHalfY[modID] * 2. - xyzLocS[1];
    xyzLocE[0] = mModuleHalfX[modID] * 2. - xyzLocE[0];
    xyzLocE[1] = mModuleHalfY[modID] * 2. - xyzLocE[1];
  } else if (modID % kNModulesPerLayer == 3) {
    // swap y
    xyzLocS[1] = mModuleHalfY[modID] * 2. - xyzLocS[1];
    xyzLocE[1] = mModuleHalfY[modID] * 2. - xyzLocE[1];
  }
}

void NA6PVerTelDigitizer::processHit(NA6PVerTelHit hit)
{
  auto modID = hit.getDetectorID();
  auto& mod = mModules[modID];
  if (mod.isDisabled()) {
    LOGP(info, "Skipping disabled module {}", modID);
    return;
  }
  double xyzLocS[3], xyzLocE[3];
  getHitLocalCoord(hit, xyzLocS, xyzLocE);
  // we assume for the time being that the charge is deposited uniformly along all the 50 um Si thickness
  // -> to be improved once we have a better knowledge of charge collection in MOSAIX
  double deltaX = xyzLocE[0] - xyzLocS[0];
  double deltaY = xyzLocE[1] - xyzLocS[1];
  float pitchX = NA6PVerTelSegmentation::ActiveDX / NA6PVerTelSegmentation::NColsPerTile;
  float pitchY = NA6PVerTelSegmentation::ActiveDYTile / NA6PVerTelSegmentation::NRowsPerTile;
  int nSteps = std::max(std::abs(deltaX) / pitchX, std::abs(deltaY) / pitchY);
  if (nSteps < 1) nSteps = 1;
  if (nSteps > 100) nSteps = 100;
  float chargePerStep = hit.getHitValue() * kGeVToEl / nSteps;
  float stepX = deltaX / nSteps;
  float stepY = deltaY / nSteps;
  float x = static_cast<float>(xyzLocS[0]) + 0.5f * stepX;
  float y = static_cast<float>(xyzLocS[1]) + 0.5f * stepY;
  for (int iStep = 0; iStep < nSteps; ++iStep) {
    UShort_t rsu, tile, row, col;
    bool digOk = mSegmentation.localToIndices(x, y, rsu, tile, row, col);
    if (digOk) {
      auto key = mod.getOrderingKey(rsu, tile, row, col);
      PreDigit* pd = mod.findDigit(key);
      if (!pd) {
        mod.addDigit(key, rsu, tile, row, col, chargePerStep, hit.getTrackID());
      } else { // there is already a digit at this slot, account as PreDigitExtra contribution
        pd->charge += chargePerStep;
        // to be added: treat multiple particles (getTrackID) in the same pre-digit
      }
    }
    x += stepX;
    y += stepY;
  }
}

void NA6PVerTelDigitizer::finalizeDigits()
{
  for (int jMod = 0; jMod < mNumberOfModules; ++jMod) {
    auto& mod = mModules[jMod];
    if (mod.isDisabled()) {
      LOGP(info, "Skipping disabled module {}", jMod);
      continue;
    }
    auto& buffer = mod.getPreDigits();
    if (buffer.empty()) {
      continue;
    }
    auto itBeg = buffer.begin();
    auto iter = itBeg;
    for (; iter != buffer.end(); ++iter) {
      const auto& preDig = iter->second;
      int jTile = jMod * NA6PVerTelSegmentation::NXTiles * NA6PVerTelSegmentation::NYSensors + NA6PVerTelSegmentation::NTilesPerRSU * preDig.pixID.rsu +  preDig.pixID.tile;
      if (preDig.charge >= mThresholds[jTile]) {
        mDigits.emplace_back(static_cast<uint16_t>(jMod), preDig.pixID, preDig.particleID);
        auto dig = mDigits.back();
      }
    }
    buffer.erase(itBeg, iter); // erase processed entries; iter == end() for now
  }
}
