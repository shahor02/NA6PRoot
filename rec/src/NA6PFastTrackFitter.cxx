// NA6PCCopyright

#include <fmt/format.h>
#include <fairlogger/Logger.h>
#include <TGeoManager.h>
#include <TFile.h>
#include <TSystem.h>
#include "NA6PTrack.h"
#include "MagneticField.h"
#include "NA6PFastTrackFitter.h"
#include "NA6PVerTelCluster.h"
#include "NA6PMuonSpecCluster.h"

const float NA6PFastTrackFitter::kMassP = 0.938f;
const float NA6PFastTrackFitter::kMassK = 0.4937f;
const float NA6PFastTrackFitter::kMassPi = 0.1396f;
const float NA6PFastTrackFitter::kMassMu = 0.1057f;
const float NA6PFastTrackFitter::kMassE = 0.0005f;
const float NA6PFastTrackFitter::kAlmostZero = 1.e-12f;
// const float NA6PFastTrackFitter::fgkFldEps = 1.e-4f;

NA6PFastTrackFitter::NA6PFastTrackFitter() : mNLayers{5},
                                             mMaxChi2Cl{10.f},
                                             mIsSeedSet{false},
                                             mSeedOption{kThreePointSeed},
                                             mCharge{1},
                                             mMass{kMassPi},
                                             mPropagateToPrimVert{false},
                                             mPrimVertZ{0.0},
                                             mIsPrimVertSet{false},
                                             mCorrectForMaterial{true},
                                             mClusters(5)
{
  if (TGeoGlobalMagField::Instance()->GetField() == nullptr) {
    auto magField = new MagneticField();
    magField->loadField();
    magField->setAsGlobalField();
  } else {
    LOGP(info, "NA6PFastTrackFitter: TGeoGlobalMagField already initialized");
  }
  mSeedPos[0] = mSeedPos[1] = mSeedPos[2] = 0.;
  mSeedMom[0] = mSeedMom[1] = mSeedMom[2] = 1.; // 1 GeV for default momentum seed
}

void NA6PFastTrackFitter::addCluster(int jLay, const NA6PBaseCluster& cl)
{
  if (jLay < 0 || jLay >= mNLayers) {
    LOGP(error, "Invalid layer index {}", jLay);
    return;
  }
  mClusters[jLay] = &cl;
}

void NA6PFastTrackFitter::setSeed(const float* pos, const float* mom, int charge)
{
  if (!pos || !mom) {
    LOGP(error, "Null pointers passer seed");
    return;
  }
  for (int j = 0; j < 3; j++) {
    mSeedPos[j] = pos[j];
    mSeedMom[j] = mom[j];
  }
  setCharge(charge);
  mIsSeedSet = true;
}

void NA6PFastTrackFitter::setParticleHypothesis(int pdg)
{
  pdg = std::abs(pdg);
  if (pdg == 211)
    mMass = kMassPi;
  else if (pdg == 321)
    mMass = kMassK;
  else if (pdg == 2212)
    mMass = kMassP;
  else if (pdg == 11)
    mMass = kMassE;
  else if (pdg == 13)
    mMass = kMassMu;
  else
    LOGP(info, "pdg code {} not valid (only pi, K, p, e, and mu are accepted)", pdg);
}

bool NA6PFastTrackFitter::loadGeometry(const char* filename, const char* geoname)
{
  if (gGeoManager) {
    LOGP(info, "NA6PFastTrackFitter: Geometry was already loaded");
    return true;
  }
  if (gSystem->Exec(Form("ls -l %s > /dev/null", filename)) != 0) {
    LOGP(error, "filename {} does not exist", filename);
    return false;
  }
  TFile* f = TFile::Open(filename);
  gGeoManager = (TGeoManager*)f->Get(geoname);
  f->Close();
  if (gGeoManager)
    return true;
  else {
    LOGP(error, "No geometry with name {} found in file {}", geoname, filename);
    return false;
  }
}

void NA6PFastTrackFitter::getMeanMaterialBudgetFromGeom(float* startF, float* endF, float* mparam) const
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
    LOGP(info, "No geometry was loaded, won't apply materical corrections");
    return;
  }
  float tolerance = 1e-9f;
  float length = std::sqrt((endF[0] - startF[0]) * (endF[0] - startF[0]) +
			   (endF[1] - startF[1]) * (endF[1] - startF[1]) +
			   (endF[2] - startF[2]) * (endF[2] - startF[2]));
  mparam[4] = length;
  if (length < tolerance)
    return;
  float invlen = 1.f / length;
  double dir[3] = { dir[0] = (endF[0] - startF[0]) * invlen, dir[1] = (endF[1] - startF[1]) * invlen, dir[2] = (endF[2] - startF[2]) * invlen };
  double start[3] = {startF[0],startF[1],startF[2]};
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

int NA6PFastTrackFitter::propagateToZ(NA6PTrack* trc, float zTo) const
{
  float zCurr = trc->getZ();
  int dir = (zTo - zCurr) > 0 ? 1 : -1;
  return propagateToZ(trc, zCurr, zTo, dir);
}

int NA6PFastTrackFitter::propagateToZ(NA6PTrack* trc, float zFrom, float zTo, int dir) const
{

  if (trc->getTrackExtParam().getAlpha() < 0.f)
    return -1;
  if (std::abs(trc->getZ() - zFrom) > 1.e-4f) {
    LOGP(fatal, "Track is not in the expected z position: {} vs {}", zFrom, trc->getZ());
    return -1;
  }
  if (dir * (zTo - zFrom) < 0.f) {
    LOGP(fatal, "Wrong coordinates: Dir:{} Zstart:{} Zend:{}", dir, zFrom, zTo);
    return -1;
  }
  float start[3], end[3];
  trc->getXYZ(start);
  NA6PTrack tmpTrc = *trc;
  tmpTrc.propagateToZBxByBz(zTo);
  tmpTrc.getXYZ(end);
  float p1[3], p2[3];
  trc->getXYZ(p1);
  tmpTrc.getXYZ(p2);
  float mparam[7] = {};
  if (mCorrectForMaterial)
    getMeanMaterialBudgetFromGeom(start, end, mparam);
  float xrho = mparam[0] * mparam[4];
  float x2x0 = mparam[1];
  float corrELoss = 1;
  if (!trc->propagateToZBxByBz(zTo, 1.f, x2x0, -dir * xrho * corrELoss))
    return 0;
  return 1;
}

bool NA6PFastTrackFitter::updateTrack(NA6PTrack* trc, const NA6PBaseCluster* cl) const
{
  // update track with measured cluster
  // propagate to cluster
  float meas[2] = {cl->getX(), cl->getY()};
  float measErr2[3] = {cl->getSigXX(), cl->getSigYX(), cl->getSigYY()};

  float chi2 = trc->getTrackExtParam().getPredictedChi2(meas, measErr2);
  if (chi2 > mMaxChi2Cl)
    return false; // chi2 is too large

  if (!trc->update(meas, measErr2)) {
    return false;
  }
  trc->addCluster(cl, 1, chi2);
  //
  return true;
}

void NA6PFastTrackFitter::computeSeed()
{
  // compute track seed from the 3 (or 2) outermost clusters
  int nClus = getNumberOfClusters();
  if (nClus < 2) {
    LOGP(error, "Cannot compute seed with {} clusters", nClus);
    return;
  }
  if (nClus == 2 && mSeedOption == kThreePointSeed) {
    LOGP(error, "Cannot compute seed with the 3-cluster option and only {} clusters -> resort to 2-point seed", nClus);
  }

  for (int jLay = mNLayers - 1; jLay >= 0; --jLay) {
    if (mClusters[jLay]) {
      double pos[3] = {mClusters[jLay]->getX(), mClusters[jLay]->getY(), mClusters[jLay]->getZ()}, bxyz[3] = {};
      TGeoGlobalMagField::Instance()->Field(pos, bxyz);
      bool useTwoPoint = (mSeedOption == kTwoPointSeed) || (nClus == 2) || (std::abs(bxyz[1]) < kAlmostZero);
      for (int kLay = jLay - 1; kLay >= 0; --kLay) {
        if (mClusters[kLay]) {
          auto ux = mClusters[jLay]->getX() - mClusters[kLay]->getX();
          auto uy = mClusters[jLay]->getY() - mClusters[kLay]->getY();
          auto uz = mClusters[jLay]->getZ() - mClusters[kLay]->getZ();
          auto norm = std::sqrt(ux * ux + uy * uy + uz * uz);
          if (norm > kAlmostZero) {
            ux /= norm;
            uy /= norm;
            uz /= norm;
          }
          if (useTwoPoint) {
            mSeedMom[0] = ux;
            mSeedMom[1] = uy;
            mSeedMom[2] = uz;
            mIsSeedSet = true;
            return;
          }
          for (int lLay = kLay - 1; lLay >= 0; --lLay) {
            if (mClusters[lLay]) {
              auto x1 = mClusters[lLay]->getX();
              auto z1 = mClusters[lLay]->getZ();
              auto x2 = mClusters[kLay]->getX();
              auto z2 = mClusters[kLay]->getZ();
              auto x3 = mClusters[jLay]->getX();
              auto z3 = mClusters[jLay]->getZ();
              // circle fit: compute center (cx, cz) and radius
              auto determ = 2.0f * (x1 * prec_t(z2 - z3) - z1 * (x2 - x3) + (x2 * z3 - x3 * z2));
              float cx = 0.f, cz = 0.f, radius = 1e6f;
              if (std::abs(determ) > kAlmostZero) {
                auto a1 = x1 * x1 + z1 * z1;
                auto a2 = x2 * x2 + z2 * z2;
                auto a3 = x3 * x3 + z3 * z3;
                cx = (a1 * (z2 - z3) - z1 * (a2 - a3) + (a2 * z3 - a3 * z2)) / determ;
                cz = (a1 * (x2 - x3) - x1 * (a2 - a3) + (x2 * a3 - x3 * a2)) / (-determ);
                radius = std::sqrt((x1 - cx) * (x1 - cx) + (z1 - cz) * (z1 - cz));
              }
              if (radius > 1e5) {
                // resort to 2-point seed
                mSeedMom[0] = ux;
                mSeedMom[1] = uy;
                mSeedMom[2] = uz;
                mIsSeedSet = true;
                return;
              }
              auto pt = 3.e-4f * std::abs(mCharge * bxyz[1]) * radius; // radius is in cm, By in kG, pt GeV/c
              auto nt = std::sqrt(ux * ux + uz * uz);
              if (nt < kAlmostZero)
                nt = 1.0f;
              mSeedMom[0] = pt * ux / nt;
              mSeedMom[1] = pt * uy / nt;
              mSeedMom[2] = pt * uz / nt;
              auto crossy = (x2 - x1) * (z3 - z2) - (z2 - z1) * (x3 - x2);
              int qSign = (crossy * bxyz[1] > 0.f) ? +1 : -1;
              if (mCharge * qSign < 0.f)
                mCharge = -mCharge;
              if (mSeedMom[2] < 0) { // enforce positive pz
                for (int j = 0; j < 3; ++j)
                  mSeedMom[j] = -mSeedMom[j];
              }
              mIsSeedSet = true;
              return;
            }
          }
        }
      }
    }
  }
}

void NA6PFastTrackFitter::printSeed() const
{
  if (mIsSeedSet) {
    LOGP(info, "Seed for tracking");
    LOGP(info, "  xSeed = {},   ySeed = {},  zSeed = {}", mSeedPos[0], mSeedPos[1], mSeedPos[2]);
    LOGP(info, "  pxSeed = {}, pySeed = {}, pzSeed = {}", mSeedMom[0], mSeedMom[1], mSeedMom[2]);
    LOGP(info, " charge = {}", mCharge);
  } else
    LOGP(info, "Seed not set");
}

NA6PTrack* NA6PFastTrackFitter::fitTrackPoints()
{

  int nClus = getNumberOfClusters();
  if (nClus < 2) {
    LOGP(error, "Cannot fit track with only {} clusters\n", nClus);
    return nullptr;
  }
  // Set seed from cluser in the outer layer if not set from outside
  if (!mIsSeedSet)
    computeSeed();
  if (!mIsSeedSet)
    LOGP(warn, "Track seed not computed properly, will run the fit w/o seed");
  NA6PTrack* currTr = new NA6PTrack();
  if (mIsSeedSet)
    currTr->init(mSeedPos, mSeedMom, mCharge);
  currTr->setMass(mMass);
  currTr->resetCovariance(-1);
  float zCurr = -999.;
  float zNext = zCurr;
  bool isGoodFit = true;
  bool propToNext = true;
  // construct bit mask
  uint clusterMask = 0;
  for (int jLay = 0; jLay < mNLayers; ++jLay) {
    if (mClusters[jLay])
      clusterMask |= (1 << jLay);
  }

  // track fit starts from the outer layer
  for (int jLay = mNLayers - 1; jLay >= 0; jLay--) {
    // update track with point in current layer
    if (mClusters[jLay]) {
      if (!updateTrack(currTr, mClusters[jLay])) {
        isGoodFit = false;
        break;
      } else {
        zCurr = mClusters[jLay]->getZ();
        bool isInnermostHit = (clusterMask & ((1 << jLay) - 1)) == 0;
        if (jLay == 0 || isInnermostHit) {
          if (mPropagateToPrimVert && mIsPrimVertSet)
            zNext = mPrimVertZ;
          else
            propToNext = false;
        } else
          zNext = mClusters[jLay - 1]->getZ();
      }
    }
    // propagate to next layer
    if (zCurr > 0 && propToNext) {
      float pxyzBef[3], pxyzAft[3];
      currTr->getPXYZ(pxyzBef);
      float momBef = currTr->getP();
      if (propagateToZ(currTr, zCurr, zNext, -1) != 1) {
        isGoodFit = false;
        break;
      }
    }
    zCurr = zNext;
  }

  if (!isGoodFit)
    return nullptr;
  return currTr;
}
