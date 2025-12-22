// NA6PCCopyright

#include "Propagator.h"
#include "fairlogger/Logger.h"
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>

Propagator::Propagator(bool uninitialized)
{
  if (uninitialized) {
    return;
  }
  ///< construct checking if needed components were initialized
  updateField();
}

//____________________________________________________________
void Propagator::updateField()
{
  if (!mField) {
    mField = static_cast<MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    if (!mField) {
      LOG(fatal) << "Magnetic field is not initialized!";
    }
  }
}

bool Propagator::hasGeometryLoaded() const
{
  return gGeoManager != nullptr;
}

void Propagator::getMeanMaterialBudgetFromGeom(const float* startF, const float* endF, float* mparam) const
{
  // "mparam" - parameters used for the energy and multiple scattering
  //  corrections:
  //
  // mparam[0] - mean density: sum(x_i*rho_i)/sum(x_i) [g/cm3]
  // mparam[1] - equivalent rad length fraction: sum(x_i/X0_i) [adimensional]
  // mparam[2] - mean A: sum(x_i*A_i)/sum(x_i) [adimensional]
  // mparam[3] - mean Z: sum(x_i*Z_i)/sum(x_i) [adimensional]
  // mparam[4] - length: sum(x_i) [cm]
  // mparam[5] - Z/A mean: sum(x_i*Z_i/A_i)/sum(x_i) [adimensional]
  // mparam[6] - number of boundary crosses
  mparam[0] = 0;
  mparam[1] = 1;
  mparam[2] = 0;
  mparam[3] = 0;
  mparam[4] = 0;
  mparam[5] = 0;
  mparam[6] = 0;

  if (gGeoManager == nullptr) {
    LOGP(fatal, "No geometry was loaded, cannot extract  material corrections");
  }
  float tolerance = 1e-9f;
  float length = std::sqrt((endF[0] - startF[0]) * (endF[0] - startF[0]) +
                           (endF[1] - startF[1]) * (endF[1] - startF[1]) +
                           (endF[2] - startF[2]) * (endF[2] - startF[2]));
  mparam[4] = length;
  if (length < tolerance)
    return;
  float invlen = 1.f / length;
  double dir[3] = {dir[0] = (endF[0] - startF[0]) * invlen, dir[1] = (endF[1] - startF[1]) * invlen, dir[2] = (endF[2] - startF[2]) * invlen};
  double start[3] = {startF[0], startF[1], startF[2]};
  TGeoNode* currNode = gGeoManager->InitTrack(start, dir);
  float sumSteps = 0.0f;
  float minStep = 1e-4f; // 1 micron

  float bparam[6]; // total parameters
  float lparam[6]; // local parameters
  for (int i = 0; i < 6; i++)
    bparam[i] = 0;

  while (sumSteps < length) {
    float stepMax = length - sumSteps;

    // Find next boundary (or max step to end point)
    TGeoNode* nextNode = gGeoManager->FindNextBoundaryAndStep(stepMax, false);
    if (!nextNode)
      break;
    float snext = gGeoManager->GetStep();
    // Handle numerical zero-step
    if (snext < tolerance) {
      gGeoManager->Step(minStep);
      snext = gGeoManager->GetStep();
      if (snext < tolerance)
        break; // protection in case still zero
    }

    TGeoMedium* med = currNode->GetVolume()->GetMedium();
    if (!med) {
      sumSteps += snext;
      continue;
    }
    TGeoMaterial* mat = med->GetMaterial();
    if (!mat) {
      sumSteps += snext;
      continue;
    }
    lparam[0] = mat->GetDensity();
    lparam[1] = mat->GetRadLen();
    lparam[2] = mat->GetA();
    lparam[3] = mat->GetZ();
    lparam[4] = length;
    lparam[5] = lparam[3] / lparam[2];
    if (mat->IsMixture()) {
      TGeoMixture* mixture = (TGeoMixture*)mat;
      lparam[5] = 0.f;
      float sum = 0.f;
      for (int iel = 0; iel < mixture->GetNelements(); iel++) {
        sum += mixture->GetWmixt()[iel];
        lparam[5] += mixture->GetZmixt()[iel] * mixture->GetWmixt()[iel] / mixture->GetAmixt()[iel];
      }
      lparam[5] /= sum;
    }
    bparam[0] += snext * lparam[0];
    bparam[1] += snext / lparam[1];
    bparam[2] += snext * lparam[2];
    bparam[3] += snext * lparam[3];
    bparam[5] += snext * lparam[5];
    sumSteps += snext;
    mparam[6] += 1.;
    currNode = nextNode;
  }
  mparam[0] = bparam[0] / sumSteps;
  mparam[1] = bparam[1];
  mparam[2] = bparam[2] / sumSteps;
  mparam[3] = bparam[3] / sumSteps;
  mparam[4] = sumSteps;
  mparam[5] = bparam[5] / sumSteps;
}

bool Propagator::propagateToZ(NA6PTrackParCov& track, float zToGo, MatCorrType matCorr, float maxStep) const
{
  //----------------------------------------------------------------
  //
  // Propagates the track to the plane at z (cm)
  // and correcting for the crossed material.
  // maxStep  - maximal step for propagation
  // matCorr  - material correction type, it is up to the user to make sure the pointer is attached (if LUT is requested)
  //----------------------------------------------------------------
  auto dz = zToGo - track.getZ();
  int dir = dz > 0.f ? 1 : -1;
  std::array<float, 3> b{};
  while (std::abs(dz) > Epsilon) {
    auto step = std::min(std::abs(dz), maxStep);
    if (dir < 0) {
      step = -step;
    }
    auto z = track.getZ() + step;
    auto xyz0 = track.getXYZ();
    getFieldXYZ(xyz0, b);
    auto correct = [&track, &xyz0, dir, matCorr, this]() {
      bool res = true;
      if (matCorr != MatCorrType::USEMatCorrNONE) {
        auto xyz1 = track.getXYZ();
        auto mb = this->getMeanMaterial(xyz0, xyz1);
        if (!track.correctForMeanMaterial(mb.meanX2X0, dir > 0 ? -mb.meanXRho : mb.meanXRho)) {
          res = false;
        }
      }
      return res;
    };

    if (!track.propagateToZ(z, b)) {
      return false;
    }
    if (!correct()) {
      return false;
    }
    dz = zToGo - track.getZ();
  }
  track.setZ(zToGo);
  return true;
}

bool Propagator::propagateToZ(NA6PTrackPar& track, float zToGo, MatCorrType matCorr, float maxStep) const
{
  //----------------------------------------------------------------
  //
  // Propagates the track to the plane at z (cm)
  // and correcting for the crossed material.
  // maxStep  - maximal step for propagation
  // matCorr  - material correction type, it is up to the user to make sure the pointer is attached (if LUT is requested)
  //----------------------------------------------------------------
  auto dz = zToGo - track.getZ();
  int dir = dz > 0.f ? 1 : -1;
  std::array<float, 3> b{};
  while (std::abs(dz) > Epsilon) {
    auto step = std::min(std::abs(dz), maxStep);
    if (dir < 0) {
      step = -step;
    }
    auto z = track.getZ() + step;
    auto xyz0 = track.getXYZ();
    getFieldXYZ(xyz0, b);
    auto correct = [&track, &xyz0, dir, matCorr, this]() {
      bool res = true;
      if (matCorr != MatCorrType::USEMatCorrNONE) {
        auto xyz1 = track.getXYZ();
        auto mb = this->getMeanMaterial(xyz0, xyz1);
        if (!track.correctForMeanMaterial(mb.meanX2X0, dir > 0 ? -mb.meanXRho : mb.meanXRho)) {
          res = false;
        }
      }
      return res;
    };

    if (!track.propagateToZ(z, b)) {
      return false;
    }
    if (!correct()) {
      return false;
    }
    dz = zToGo - track.getZ();
  }
  track.setZ(zToGo);
  return true;
}
