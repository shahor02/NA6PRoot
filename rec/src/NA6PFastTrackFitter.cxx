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
#include "Propagator.h"

const float NA6PFastTrackFitter::kAlmostZero = 1.e-12f;
// const float NA6PFastTrackFitter::fgkFldEps = 1.e-4f;

NA6PFastTrackFitter::NA6PFastTrackFitter() : mNLayers{5},
                                             mMaxChi2Cl{10.f},
                                             mIsSeedSet{false},
                                             mSeedOption{kThreePointSeed},
                                             mCharge{1},
                                             mPID{PID::Pion},
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
  mSeedPos[0] = mSeedPos[1] = mSeedPos[2] = 0.f;
  mSeedMom[0] = mSeedMom[1] = mSeedMom[2] = 1.f; // 1 GeV for default momentum seed
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
    mPID = PID::Pion;
  else if (pdg == 321)
    mPID = PID::Kaon;
  else if (pdg == 2212)
    mPID = PID::Proton;
  else if (pdg == 11)
    mPID = PID::Electron;
  else if (pdg == 13)
    mPID = PID::Muon;
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

bool NA6PFastTrackFitter::updateTrack(NA6PTrack* trc, const NA6PBaseCluster* cl) const
{
  // update track with measured cluster
  float chi2 = trc->getPredictedChi2(*cl);
  if (chi2 > mMaxChi2Cl)
    return false; // chi2 is too large

  if (!trc->update(*cl)) {
    return false;
  }
  trc->addCluster(cl, cl->getHitID(), chi2);
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
      mSeedPos[0] = mClusters[jLay]->getX();
      mSeedPos[1] = mClusters[jLay]->getY();
      mSeedPos[2] = mClusters[jLay]->getZ();
      float bxyz[3] = {};
      Propagator::Instance()->getFieldXYZ(mSeedPos, bxyz);
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
              auto determ = 2.0f * (x1 * (z2 - z3) - z1 * (x2 - x3) + (x2 * z3 - x3 * z2));
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
  // Set seed from cluster in the outer layer if not set from outside
  if (!mIsSeedSet)
    computeSeed();
  if (!mIsSeedSet)
    LOGP(warn, "Track seed not computed properly, will run the fit w/o seed");
  NA6PTrack* currTr = new NA6PTrack();
  if (mIsSeedSet)
    currTr->init(mSeedPos, mSeedMom, mCharge);
  currTr->setPID(mPID);
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
  bool first = true;
  int nclUpd = 0;
  for (int jLay = mNLayers - 1; jLay >= 0; jLay--) {
    // update track with point in current layer
    if (mClusters[jLay]) {
      if (first) {
        zCurr = mClusters[jLay]->getZ();
        currTr->setZ(zCurr);
        first = false;
      }
      if (!updateTrack(currTr, mClusters[jLay])) {
        isGoodFit = false;
        break;
      } else {
        bool isInnermostHit = (clusterMask & ((1 << jLay) - 1)) == 0;
        if (jLay == 0 || isInnermostHit) {
          if (mPropagateToPrimVert && mIsPrimVertSet)
            zNext = mPrimVertZ;
          else
            propToNext = false;
        } else
          zNext = mClusters[jLay - 1]->getZ();
      }
      nclUpd++;
    }
    // propagate to next layer
    if (zCurr > 0 && propToNext) {
      if (!Propagator::Instance()->propagateToZ(*currTr, zNext, {.fixCorrelations = nclUpd < 2})) {
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
