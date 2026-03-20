// NA6PCCopyright

#include "Propagator.h"
#include "fairlogger/Logger.h"
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TFile.h>

Propagator::Propagator(bool uninitialized)
{
  if (uninitialized) {
    return;
  }
  ///< construct checking if needed components were initialized
  updateField();
}

//____________________________________________________________
bool Propagator::loadField()
{
  if (dynamic_cast<MagneticField*>(TGeoGlobalMagField::Instance()->GetField())) {
    LOGP(warn, "Global magnetic field was already loaded");
    return true;
  }
  auto prop = Instance(true);
  prop->mOwnField = std::make_unique<MagneticField>();
  prop->mOwnField->loadField();
  prop->mOwnField->setAsGlobalField();
  prop->updateField();
  return true;
}

//____________________________________________________________
bool Propagator::loadGeometry(const std::string& path, const std::string geomName)
{
  if (gGeoManager) {
    LOGP(warn, "Geometry was already loaded");
    return true;
  }
  auto f = TFile::Open(path.c_str());
  if (!f || f->IsZombie()) {
    LOGP(error, "Could not open geometry file {}", path);
    return false;
  }
  if (!f->Get(geomName.c_str())) {
    LOGP(error, "Could not load geometry {} from file {}", geomName, path);
  }
  f->Close();
  delete f;
  return gGeoManager != nullptr;
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

template <typename T>
MatBudget Propagator::getMeanMaterialBudgetFromGeom(const T* startF, const T* endF) const
{
  // "mparam" - parameters used for the energy and multiple scattering
  //  corrections:
  MatBudget mbTot;
  if (gGeoManager == nullptr) {
    LOGP(fatal, "No geometry was loaded, cannot extract  material corrections");
  }
  float tolerance = 1e-9f;
  float dx = endF[0] - startF[0], dy = endF[1] - startF[1], dz = endF[2] - startF[2];
  float length = std::sqrt(dx * dx + dy * dy + dz * dz);
  if (length < tolerance) {
    mbTot.meanL = length;
    return mbTot;
  }
  float invlen = 1.f / length;
  double dir[3] = {dx * invlen, dy * invlen, dz * invlen}, start[3] = {startF[0], startF[1], startF[2]};
  TGeoNode* currNode = gGeoManager->InitTrack(start, dir);
  float sumSteps = 0.0f, minStep = 1e-4f; // 1 micron

  while (sumSteps < length) {
    MatBudget mbLoc;
    float stepMax = length - sumSteps;

    // Find next boundary (or max step to end point)
    TGeoNode* nextNode = gGeoManager->FindNextBoundaryAndStep(stepMax, false);
    if (!nextNode) {
      break;
    }
    float snext = gGeoManager->GetStep();
    // Handle numerical zero-step
    if (snext < tolerance) {
      gGeoManager->Step(minStep);
      snext = gGeoManager->GetStep();
      if (snext < tolerance) {
        break; // protection in case still zero
      }
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

    mbTot.meanRho += snext * mat->GetDensity();
    mbTot.meanX2X0 += snext / mat->GetRadLen();
    mbTot.meanA += snext * mat->GetA();
    mbTot.meanZ += snext * mat->GetZ();
    if (mat->IsMixture()) {
      TGeoMixture* mixture = (TGeoMixture*)mat;
      float meanZ2A = 0.f, sum = 0.f;
      for (int iel = 0; iel < mixture->GetNelements(); iel++) {
        sum += mixture->GetWmixt()[iel];
        meanZ2A += mixture->GetZmixt()[iel] * mixture->GetWmixt()[iel] / mixture->GetAmixt()[iel];
      }
      mbTot.meanZ2A += snext * meanZ2A / sum;
    } else {
      mbTot.meanZ2A += snext * mat->GetZ() / mat->GetA();
    }
    sumSteps += snext;
    mbTot.nBoundaries += 1;
    currNode = nextNode;
  }

  float norm = 1.f / sumSteps;
  mbTot.meanRho *= norm;
  mbTot.meanA *= norm;
  mbTot.meanZ *= norm;
  mbTot.meanZ2A *= norm;
  mbTot.meanL = sumSteps;
  return mbTot;
}

bool Propagator::propagateToZ(NA6PTrackParCov& track, float zToGo, const Propagator::PropOpt& opt) const
{
  //----------------------------------------------------------------
  //
  // Propagates the track to the plane at z (cm)
  // and correcting for the crossed material.
  // opt.maxStep  - maximal step for propagation
  // opt.byOnly   - pure dipole propagation
  // opt.matCorr  - material correction type, it is up to the user to make sure the pointer is attached (if LUT is requested)
  //----------------------------------------------------------------
  auto dz = zToGo - track.getZ();
  int dir = dz > 0.f ? 1 : -1;
  std::array<float, 3> b{};
  NA6PTrackPar& linRef = opt.linRef ? *opt.linRef : track;
  while (std::abs(dz) > Epsilon) {
    auto step = std::min(std::abs(dz), opt.maxStep);
    if (dir < 0) {
      step = -step;
    }
    auto z = track.getZ() + step;
    auto xyz0 = linRef.getXYZ();
    auto correct = [&track, &linRef, &xyz0, dir, matCorr = opt.matCorr, this]() {
      bool res = true;
      if (matCorr != MatCorrType::USEMatCorrNONE) {
        auto xyz1 = linRef.getXYZ();
        auto mb = this->getMeanMaterial(xyz0, xyz1);
        if (!track.correctForMaterial(mb.meanX2X0, (dir > 0 ? -mb.meanRho : mb.meanRho) * mb.meanL, linRef)) {
          res = false;
        }
      }
      return res;
    };
    if (!(opt.byOnly ? track.propagateToZ(z, getBy(xyz0), opt.linRef) : track.propagateToZ(z, getFieldXYZ(xyz0), opt.linRef))) {
      return false;
    }
    if (opt.fixCorrelations) {
      track.fixCorrelations();
    }
    if (!correct()) {
      return false;
    }
    dz = zToGo - track.getZ();
  }
  track.setZ(zToGo);
  if (opt.linRef) {
    opt.linRef->setZ(zToGo);
  }
  return true;
}

bool Propagator::propagateToZ(NA6PTrackPar& track, float zToGo, const Propagator::PropOpt& opt) const
{
  //----------------------------------------------------------------
  //
  // Propagates the track to the plane at z (cm)
  // and correcting for the crossed material.
  // opt.maxStep  - maximal step for propagation
  // opt.byOnly   - pure dipole propagation
  // opt.matCorr  - material correction type, it is up to the user to make sure the pointer is attached (if LUT is requested)
  //----------------------------------------------------------------
  auto dz = zToGo - track.getZ();
  int dir = dz > 0.f ? 1 : -1;
  std::array<float, 3> b{};
  while (std::abs(dz) > Epsilon) {
    auto step = std::min(std::abs(dz), opt.maxStep);
    if (dir < 0) {
      step = -step;
    }
    auto z = track.getZ() + step;
    auto xyz0 = track.getXYZ();
    getFieldXYZ(xyz0, b);
    auto correct = [&track, &xyz0, dir, matCorr = opt.matCorr, this]() {
      bool res = true;
      if (matCorr != MatCorrType::USEMatCorrNONE) {
        auto xyz1 = track.getXYZ();
        auto mb = this->getMeanMaterial(xyz0, xyz1);
        if (!track.correctForELoss((dir > 0 ? -mb.meanRho : mb.meanRho)) * mb.meanL) {
          res = false;
        }
      }
      return res;
    };
    if (!(opt.byOnly ? track.propagateParamToZ(z, getBy(xyz0)) : track.propagateParamToZ(z, getFieldXYZ(xyz0)))) {
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

template <class T>
bool Propagator::propagatePCAToLine(T& track, float x0, float y0, float tolerance, const Propagator::PropOpt& opt) const
{
  // propagate the track to the point of closest aproach to the straight line x(z) = x0, y(z) = y0
  float zca = 0;
  int iterLeft = 2;
  while (iterLeft--) {
    if (!track.getZPCAToLine(x0, y0, getBy(track.getXYZ()), zca) || !propagateToZ(track, zca, opt)) {
      return false;
    }
    float dx = x0 - track.getX(), dy = y0 - track.getY();
    if (dx * dx + dy * dy < tolerance * tolerance) {
      break;
    }
  }
  return true;
}

template MatBudget Propagator::getMeanMaterialBudgetFromGeom<float>(const float* start, const float* end) const;
template MatBudget Propagator::getMeanMaterialBudgetFromGeom<double>(const double* start, const double* end) const;
template bool Propagator::propagatePCAToLine(NA6PTrackPar& track, float x0, float y0, float tolerance, const Propagator::PropOpt& opt) const;
template bool Propagator::propagatePCAToLine(NA6PTrackParCov& track, float x0, float y0, float tolerance, const Propagator::PropOpt& opt) const;
