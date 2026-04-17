// NA6PCCopyright

#include <fmt/format.h>
#include <fairlogger/Logger.h>
#include <TGeoManager.h>
#include <TFile.h>
#include <TSystem.h>
#include "NA6PVertex.h"
#include "MagneticField.h"
#include "NA6PPhysConst.h"
#include "NA6PFastTrackFitter.h"
#include "NA6PVerTelCluster.h"
#include "NA6PMuonSpecCluster.h"
#include "NA6PLine.h"

NA6PFastTrackFitter::NA6PFastTrackFitter()
{
  if (TGeoGlobalMagField::Instance()->GetField() == nullptr) {
    Propagator::loadField();
  }
}

void NA6PFastTrackFitter::printClusters() const
{
  for (int jLay = 0; jLay < getNLayers(); ++jLay) {
    if (mClusters[jLay])
      LOGP(info, "Layer {} Cluster position = {} {} {}", jLay, mClusters[jLay]->getX(), mClusters[jLay]->getY(), mClusters[jLay]->getZ());
    else
      LOGP(info, "Layer {} No cluster", jLay);
  }
}

void NA6PFastTrackFitter::setSeed(const float* pos, const float* mom, int charge) // RSTODO
{
  if (!pos || !mom) {
    LOGP(error, "Null pointers passer seed");
    return;
  }
  mSeed.initParam(pos, mom, charge);
}

bool NA6PFastTrackFitter::constrainTrackToVertex(NA6PTrack& trc, const NA6PVertex& pv) const
{
  NA6PTrackParCov t = trc;
  if (!Propagator::Instance()->propagateToZ(t, pv.getZ(), mPropOpt)) {
    return false;
  }
  // create a pseudocluster for the vertex, use -1 for layer to distinguish it from detector hits
  NA6PBaseCluster vClu(pv.getX(), pv.getY(), pv.getZ(), 1, -1);
  auto [tx, ty] = t.getSlopesXY(); // project covariance matrix elements using track direction
  auto sxx = pv.getSigmaX2() + (tx * tx * pv.getSigmaZ2()) - (2.0f * tx * pv.getSigmaXZ());
  auto syy = pv.getSigmaY2() + (ty * ty * pv.getSigmaZ2()) - (2.0f * ty * pv.getSigmaYZ());
  auto sxy = pv.getSigmaXY() + (tx * ty * pv.getSigmaZ2()) - (tx * pv.getSigmaYZ()) - (ty * pv.getSigmaXZ());
  vClu.setErr(sxx, sxy, syy);
  auto chi2 = t.getPredictedChi2(vClu);
  if (chi2 > mMaxChi2Cl || !t.update(vClu)) {
    trc.getVertexConstrainedParam().invalidate();
    return false;
  }
  trc.setVertexConstrainedParam(t);
  return true;
}

int NA6PFastTrackFitter::getLayersForSeed(std::array<int, 3>& layForSeed)
{
  // define layers to be used in seed calculation depending on the options
  int nLayWithClus = countLayerWithClusters();
  if (nLayWithClus < 2) {
    LOGP(error, "Only {} layers with clusters, at least 2 needed for seeding", nLayWithClus);
    return 0;
  }
  layForSeed[2] = -1; // in case only 2 layers available
  if (nLayWithClus <= 3) {
    for (int i = 0; i < mNClusters; i++) {
      layForSeed[i] = mLayersWithClusters[i];
    }
    return nLayWithClus;
  }
  // select layers for seed determination
  if (mSeedOption == kInnermostAsSeed) {
    for (int i = 0; i < 3 && i < nLayWithClus; ++i) {
      layForSeed[i] = mLayersWithClusters[i];
    }
  } else if (mSeedOption == kOutermostAsSeed) {
    for (int i = 0; i < 3 && i < nLayWithClus; ++i) {
      layForSeed[i] = mLayersWithClusters[nLayWithClus - 1 - i];
    }
  } else { // kInMidOutAsSeed
    layForSeed[0] = mLayersWithClusters[0];
    layForSeed[2] = mLayersWithClusters[nLayWithClus - 1];
    if (layForSeed[0] == layForSeed[2]) {
      LOGP(error, "Innermost and outermost clusters coincide");
      return 0;
    }
    int mid = (layForSeed[0] + layForSeed[2]) / 2;
    // Find closest layer to midpoint
    int best = -1;
    int bestDist = 999;
    for (int i = 0; i < nLayWithClus; ++i) {
      int lay = mLayersWithClusters[i];
      if (lay == layForSeed[0] || lay == layForSeed[2]) {
        continue;
      }
      int dist = std::abs(lay - mid);
      if (dist < bestDist) {
        bestDist = dist;
        best = lay;
      }
    }
    layForSeed[1] = best;
    if (best == -1) {
      LOGP(warning, "Cannot find middle layer between {} and {}", layForSeed[0], layForSeed[2]);
    }
  }
  // count layers
  int nSeedLayers = 0;
  for (int i = 0; i < 3; ++i) {
    if (layForSeed[i] >= 0) {
      ++nSeedLayers;
    }
  }
  return nSeedLayers;
}

bool NA6PFastTrackFitter::computeSeed(int dir, std::array<int, 3>& layForSeed, NA6PTrackPar* seed)
{
  if (!seed) {
    seed = &mSeed;
  }
  seed->invalidate();
  int nSeedLayers = layForSeed[2] >= 0 ? 3 : 2;
  if (nSeedLayers == 2 && mSeedPoints == kThreePointSeed) {
    LOGP(error, "Cannot compute seed with the 3-layer option and only {} layers -> resort to 2-point seed", nSeedLayers);
  }
  auto getLayerInOrder = [&layForSeed, nSeedLayers, dir](int i) {
    return dir > 0 ? layForSeed[i] : layForSeed[nSeedLayers - 1 - i];
  };

  auto prop = Propagator::Instance();
  int jLay = getLayerInOrder(0);
  seed->setXYZ(mClusters[jLay]->getXYZ());
  float by = prop->getBy(mClusters[jLay]->getXYZ());
  bool useTwoPoint = (mSeedPoints == kTwoPointSeed) || (nSeedLayers == 2) || (std::abs(by) < NA6PTrackPar::kSmallBend);
  int kLay = getLayerInOrder(1);
  auto uvec = NA6PLine::getDiff(mClusters[jLay]->getXYZ(), mClusters[kLay]->getXYZ());
  float norm = NA6PLine::getNorm(uvec);
  if (norm > phys_const::kAlmost0F) {
    auto normI = uvec[2] < 0 ? -1.f / norm : 1. / norm; // Enforce track direction along +z
    for (int i = 0; i < 3; i++) {
      uvec[i] *= normI;
    }
  }
  auto setTwoPointKin = [&uvec, &seed]() {
    auto pxz = std::hypot(uvec[0], uvec[2]);
    seed->setTx(uvec[0] / pxz);
    seed->setTy(uvec[1] / pxz);
    seed->setQ2Pxz(1.f / pxz); // RSCHECK
    return true;
  };
  if (useTwoPoint) {
    return setTwoPointKin();
  }
  int lLay = getLayerInOrder(2);
  if (lLay < 0 || !mClusters[lLay]) {
    LOGP(error, "Third seed layer invalid");
    return false;
  }
  if (dir < 0) {
    std::swap(jLay, lLay); // swap lLay and jLay to order the layers from innermost to outermost
  }
  auto *clJ = mClusters[jLay], *clK = mClusters[kLay], *clL = mClusters[lLay];
  float x1 = clJ->getX(), z1 = clJ->getZ(), x2 = clK->getX(), z2 = clK->getZ(), x3 = clL->getX(), z3 = clL->getZ();
  if (mOptionForSeedB == kBatMidPoint) {
    by = prop->getBy(clK->getXYZ());         // get magnetic field in the middle point
  } else if (mOptionForSeedB == kMaximumB) { // use maximum magnetic field among the 3 hits
    by = std::max({prop->getBy(clJ->getXYZ()), prop->getBy(clK->getXYZ()), prop->getBy(clL->getXYZ())});
  } else { // integral of field
    constexpr int nSteps = 10;
    auto point = clJ->getXYZ();
    const auto &midPoint = clK->getXYZ(), &lastPoint = clL->getXYZ();
    std::array<float, 3> step = NA6PLine::getScaled(NA6PLine::getDiff(midPoint, lastPoint), 1.f / nSteps);
    float stepL = NA6PLine::getNorm(step), aveField = 0.f, totLength = 0.f;
    for (int js = 0; js <= nSteps; js++) {
      aveField += prop->getBy(point) * stepL;
      totLength += stepL;
      for (int k = 0; k < 3; k++) {
        point[k] += step[k];
      }
    }
    by = aveField / totLength;
  }
  // circle fit: compute center (cx, cz) and radius
  float determ = 2.0 * (x1 * (z2 - z3) + x2 * (z3 - z1) + x3 * (z1 - z2));
  float cx = 0.0, cz = 0.0, radius = 1e6;
  if (std::abs(determ) > phys_const::kAlmost0F) {
    float a1 = x1 * x1 + z1 * z1;
    float a2 = x2 * x2 + z2 * z2;
    float a3 = x3 * x3 + z3 * z3;
    cx = (a1 * (z2 - z3) + a2 * (z3 - z1) + a3 * (z1 - z2)) / determ;
    cz = (a1 * (x3 - x2) + a2 * (x1 - x3) + a3 * (x2 - x1)) / determ;
    radius = std::sqrt((x1 - cx) * (x1 - cx) + (z1 - cz) * (z1 - cz));
  }

  if (radius > 1e5) {
    return setTwoPointKin();
  }
  float pxzI = 1. / std::abs(NA6PTrackPar::kB2C * by * radius); // radius is in cm, By in kG, pxz GeV/c
  float ntI = 1.f / std::hypot(uvec[0], uvec[2]);
  float crossy = (x2 - x1) * (z3 - z2) - (z2 - z1) * (x3 - x2);
  int qSign = (crossy * by > 0) ? +1 : -1;
  seed->setTx(uvec[0] * ntI);
  seed->setTy(uvec[1] * ntI);
  seed->setQ2Pxz((crossy * by > 0) ? pxzI : -pxzI);
  return true;
}

bool NA6PFastTrackFitter::computeSeed(int dir, NA6PTrackPar* seed)
{
  // compute track seed from 3 (or 2) clusters
  // dir = 1 -> forward dir = -1 -> backward
  if (!seed) {
    seed = &mSeed;
  }
  std::array<int, 3> layForSeed = {-1, -1, -1};
  int origSeedOption = mSeedOption;
  if (dir == 1 && mSeedOption == kOutermostAsSeed) {
    mSeedOption = kInnermostAsSeed;
  }
  if (dir == -1 && mSeedOption == kInnermostAsSeed) {
    mSeedOption = kOutermostAsSeed;
  }
  int nSeedLayers = getLayersForSeed(layForSeed);
  if (nSeedLayers < 2) {
    seed->invalidate();
    LOGP(error, "Cannot compute seed with {} seeding layers", nSeedLayers);
    return false;
  }
  mSeedOption = origSeedOption;
  return computeSeed(dir, layForSeed, seed);
}

void NA6PFastTrackFitter::printSeed() const
{
  if (mSeed.isValid()) {
    LOGP(info, "Seed for tracking: {}", mSeed.asString());
  } else
    LOGP(info, "Seed not set");
}

/*
bool NA6PFastTrackFitter::fitTrackPoints(NA6PTrack& trackToFit, int dir, const NA6PTrackParCov* seed)
{
  if (mNClusters < 2) {
    LOGP(error, "Cannot fit track with only {} clusters\n", mNClusters);
    return false;
  }
  bool isGoodFit = true;
  auto prop = Propagator::Instance();
  trackToFit.reset();
  if (seed) {
    // Seed provided as a track from previous pass
    ((NA6PTrackParCov&)trackToFit) = *seed;
  } else {
    if (!mSeed.isValid()) { // Set seed from clusters in the outer (inner) layers if not set from outside
      dir > 0 ? computeSeedInner() : computeSeedOuter();
    }
    if (mSeed.isValid()) {
      ((NA6PTrackPar&)trackToFit) = mSeed;
    } else {
      LOGP(warn, "Track seed not computed properly, will run the fit w/o seed");
    }
  }
  trackToFit.setPID(mPID);
  trackToFit.resetCovariance(-1);
  int startLay = (dir > 0) ? 0 : getNLayers() - 1, endLay = (dir > 0) ? getNLayers() : -1, stepLay = (dir > 0) ? 1 : -1;
  int nClFit = 0;
  for (int jLay = startLay; jLay != endLay; jLay += stepLay) {
    if (mClusters[jLay]) {
      if (nClFit == 0) { // 1st cluster, just assign track to its position
        trackToFit.setXYZ(mClusters[jLay]->getXYZ());
      } else {
        if (!prop->propagateToZ(trackToFit, mClusters[jLay]->getZ(), mPropOpt)) {
          isGoodFit = false;
          break;
        }
      }
      if (!trackToFit.updateTrack(*mClusters[jLay], mMaxChi2Cl)) {
        isGoodFit = false;
        break;
      }
      nClFit++;
    }
  }
  if (isGoodFit && dir < 0 && mPropagateToPrimVert && mIsPrimVertSet && !prop->propagateToZ(trackToFit, mPrimVertZ, mPropOpt)) {
    isGoodFit = false;
  }
  return isGoodFit;
}
*/

//=========================================================
float NA6PFastTrackFitter::fitSeed(NA6PTrackParCov& seed, bool resetCovMat, int dir, NA6PTrackPar* linRef)
{
  // fit track (already seeded but not necessarily at the position of the 1st cluster to fit.
  // If the covariance matrix reset it requested, then the track position is imposed at the 1st cluster to fit (w/o propagating other track params)
  // and the covariance is reset. Otherwise (assuming that the fit is a continuation of the existing fit), the track is propagated to the 1st point
  // (including the cov. matrix propagation).
  // If linRef is provided, it is used for the KF linearization, otherwise the linearization is done wrt the track being updated.
  // The propagation is done with the PID of the provided seed (or that of linRef, if provided).
  //
  // In case of the sucessful fit return total chi2, otherwise negative number
  if (mNClusters < 2) {
    LOGP(error, "Cannot fit track with only {} clusters\n", mNClusters);
    return -1.f;
  }
  float chi2Tot = 0.f;
  auto prop = Propagator::Instance();
  int step = 1, startL = mMinLayerWithCl, stopL = mMaxLayerWithCl + step;
  if (dir < 0) {
    startL = mMaxLayerWithCl;
    step = -1;
    stopL = mMinLayerWithCl - 1;
  }
  if (resetCovMat) {
    auto cl = getCluster(startL);
    seed.setXYZ(cl->getXYZ());
    seed.resetCovariance(-1);
    auto prelQPerpI = linRef ? linRef->getParam(4) : seed.getParam(4);
    seed.setCovMatElem(4, 4, prelQPerpI * prelQPerpI * seed.getCovMatElem(4, 4));
    if (linRef) {
      if (!prop->propagateToZ(*linRef, cl->getZ(), mPropOpt)) { // linRef must be at the same position as the seed
        return -1.f;
      }
      linRef->setXYZ(cl->getXYZ());
    }
  }
  mPropOpt.linRef = linRef;

  for (int il = startL; il != stopL; il += step) {
    if (!mClusters[il]) {
      continue;
    }
    const auto& cl = *mClusters[il];
    float chi2 = 0;
    if (!prop->propagateToZ(seed, cl.getZ(), mPropOpt) ||
        (chi2 = seed.getPredictedChi2(cl)) > mMaxChi2Cl ||
        !seed.update(cl)) {
      chi2Tot = -1;
      break;
    }
    chi2Tot += chi2;
  }

  if (dir < 0 && mPropagateToPrimVert && mIsPrimVertSet && !prop->propagateToZ(seed, mPrimVertZ, mPropOpt)) {
    return -1.f;
  }

  mPropOpt.linRef = nullptr;
  return chi2Tot;
}

void NA6PFastTrackFitter::addClustersToTrack(NA6PTrack& tr)
{
  for (int i = mMinLayerWithCl; i <= mMaxLayerWithCl; i++) {
    const auto* cl = getCluster(i);
    if (cl) {
      tr.addCluster(cl, cl->getClusterIndex());
    }
  }
}
