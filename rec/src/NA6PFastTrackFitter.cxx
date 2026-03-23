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

ClassImp(NA6PFastTrackFitter)

  NA6PFastTrackFitter::NA6PFastTrackFitter()
{
  if (TGeoGlobalMagField::Instance()->GetField() == nullptr) {
    Propagator::loadField();
  }
  mClusters.resize(mNLayers);
}

void NA6PFastTrackFitter::addCluster(int jLay, const NA6PBaseCluster& cl)
{
  if (jLay < 0 || jLay >= mNLayers) {
    LOGP(error, "Invalid layer index {}", jLay);
    return;
  }
  mClusters[jLay] = &cl;
}

void NA6PFastTrackFitter::printClusters() const
{
  LOGP(info, "N layers = {}", mNLayers);
  for (int jLay = 0; jLay < mNLayers; ++jLay) {
    if (mClusters[jLay])
      LOGP(info, "Layer {} Cluster position = {} {} {}", jLay, mClusters[jLay]->getX(), mClusters[jLay]->getY(), mClusters[jLay]->getZ());
    else
      LOGP(info, "Layer {} No cluster", jLay);
  }
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
    trc.setStatusConstrained(false);
    return false;
  }
  trc.setVertexConstrainedParam(t);
  trc.setStatusConstrained(true);
  return true;
}

int NA6PFastTrackFitter::getLayersForSeed(std::array<int, 3>& layForSeed) const
{
  // define layers to be used in seed calculation depending on the options

  int layWithClus[kMaxLayers];
  int nLayWithClus = 0;
  for (int i = 0; i < mNLayers; ++i) {
    if (mClusters[i]) {
      layWithClus[nLayWithClus++] = i;
    }
  }
  if (nLayWithClus < 2) {
    LOGP(error, "Array with clusters has less than 2 clusters {}", nLayWithClus);
    return 0;
  }
  // select layers for seed determination
  if (mSeedOption == kInnermostAsSeed) {
    for (int i = 0; i < 3 && i < nLayWithClus; ++i) {
      layForSeed[i] = layWithClus[i];
    }
  } else if (mSeedOption == kOutermostAsSeed) {
    for (int i = 0; i < 3 && i < nLayWithClus; ++i) {
      layForSeed[i] = layWithClus[nLayWithClus - 1 - i];
    }
  } else { // kInMidOutAsSeed
    layForSeed[0] = layWithClus[0];
    layForSeed[2] = layWithClus[nLayWithClus - 1];
    if (layForSeed[0] == layForSeed[2]) {
      LOGP(error, "Innermost and outermost clusters coincide");
      return 0;
    }
    int mid = (layForSeed[0] + layForSeed[2]) / 2;
    // Find closest layer to midpoint
    int best = -1;
    int bestDist = 999;
    for (int i = 0; i < nLayWithClus; ++i) {
      int lay = layWithClus[i];
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

int NA6PFastTrackFitter::sortLayersForSeed(std::array<int, 3>& layForSeed, int dir) const
{
  // Sort according to fitting direction
  // dir = 1 -> forward (innermost first)
  // dir = -1 -> backward (outermost first)
  // Returns the number of valid layers after sorting (2 or 3)
  int tmpArr[3];
  int valid = 0;
  for (int i = 0; i < 3; ++i) {
    if (layForSeed[i] >= 0 && layForSeed[i] < mNLayers) {
      tmpArr[valid++] = layForSeed[i];
    }
  }
  if (valid < 2) {
    LOGP(error, "Less than 2 clusters found in different layers. Seed cannot be computed");
    return 0;
  }
  dir > 0 ? std::sort(tmpArr, tmpArr + valid) : std::sort(tmpArr, tmpArr + valid, std::greater<int>());
  for (int i = 0; i < 3; ++i) {
    layForSeed[i] = (i < valid) ? tmpArr[i] : -1;
  }
  return valid;
}

void NA6PFastTrackFitter::computeSeed(int dir, std::array<int, 3>& layForSeed)
{
  int nClus = getNumberOfClusters();
  if (nClus < 2) {
    LOGP(error, "Cannot compute seed with {} clusters", nClus);
    return;
  }
  if (nClus == 2 && mSeedPoints == kThreePointSeed) {
    LOGP(error, "Cannot compute seed with the 3-cluster option and only {} clusters -> resort to 2-point seed", nClus);
  }
  auto prop = Propagator::Instance();
  int nSeedLayers = sortLayersForSeed(layForSeed, dir);
  int jLay = layForSeed[0];
  if (jLay < 0 || !mClusters[jLay]) {
    LOGP(error, "First seed layer invalid");
    return;
  }
  mSeedPos[0] = mClusters[jLay]->getX();
  mSeedPos[1] = mClusters[jLay]->getY();
  mSeedPos[2] = mClusters[jLay]->getZ();
  float by = prop->getBy(mSeedPos);
  bool useTwoPoint = (mSeedPoints == kTwoPointSeed) || (nClus == 2) || (nSeedLayers == 2) || (std::abs(by) < NA6PTrackPar::kSmallBend);
  int kLay = layForSeed[1];
  if (kLay < 0 || !mClusters[kLay]) {
    LOGP(error, "Second seed layer invalid");
    return;
  }
  float ux = mClusters[jLay]->getX() - mClusters[kLay]->getX();
  float uy = mClusters[jLay]->getY() - mClusters[kLay]->getY();
  float uz = mClusters[jLay]->getZ() - mClusters[kLay]->getZ();
  float norm = std::sqrt(ux * ux + uy * uy + uz * uz);
  if (norm > phys_const::kAlmost0F) {
    // Enforce track direction along +z
    if (uz < 0) {
      norm = -norm;
    }
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
  int lLay = layForSeed[2];
  if (lLay < 0 || !mClusters[lLay]) {
    LOGP(error, "Third seed layer invalid");
    return;
  }
  if (dir < 0) {
    // swap lLay and jLay to order the layers from innermost to outermost
    int tmp = jLay;
    jLay = lLay;
    lLay = tmp;
  }
  auto* clJ = mClusters[jLay];
  auto* clK = mClusters[kLay];
  auto* clL = mClusters[lLay];
  float x1 = clJ->getX();
  float z1 = clJ->getZ();
  float x2 = clK->getX();
  float z2 = clK->getZ();
  float x3 = clL->getX();
  float z3 = clL->getZ();
  if (mOptionForSeedB == kBatMidPoint) {
    // get magnetic field in the middle point
    float midPoint[3] = {clK->getX(), clK->getY(), clK->getZ()};
    by = prop->getBy(midPoint);
  } else if (mOptionForSeedB == kMaximumB) { // use maximum magnetic field among the 3 hits
    by = std::max({prop->getBy(clJ->getXYZ()), prop->getBy(clK->getXYZ()), prop->getBy(clL->getXYZ())});
  } else { // integral of field
    const auto& firstPoint = clJ->getXYZ();
    const auto& midPoint = clK->getXYZ();
    const auto& lastPoint = clL->getXYZ();
    std::array<float, 3> step, nextPoint;
    int nSteps = 10;
    for (int jc = 0; jc < 3; jc++) {
      step[jc] = (midPoint[jc] - firstPoint[jc]) / (float)nSteps;
    }
    float stepL = std::sqrt(step[0] * step[0] + step[1] * step[1] + step[2] * step[2]);
    float aveField = 0.f;
    float totLength = 0.f;
    for (int js = 0; js <= nSteps; js++) {
      for (int jc = 0; jc < 3; jc++) {
        nextPoint[jc] = firstPoint[jc] + js * step[jc];
      }
      aveField += prop->getBy(nextPoint) * stepL;
      totLength += stepL;
    }
    for (int jc = 0; jc < 3; jc++) {
      step[jc] = (lastPoint[jc] - midPoint[jc]) / nSteps;
    }
    stepL = std::sqrt(step[0] * step[0] + step[1] * step[1] + step[2] * step[2]);
    for (int js = 0; js <= nSteps; js++) {
      for (int jc = 0; jc < 3; jc++) {
        nextPoint[jc] = midPoint[jc] + js * step[jc];
      }
      aveField += prop->getBy(nextPoint) * stepL;
      totLength += stepL;
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
    // resort to 2-point seed
    mSeedMom[0] = ux;
    mSeedMom[1] = uy;
    mSeedMom[2] = uz;
    mIsSeedSet = true;
    return;
  }
  float pxz = 3.e-4 * std::abs(mCharge * by) * radius; // radius is in cm, By in kG, pxz GeV/c
  float nt = std::sqrt(ux * ux + uz * uz);
  if (nt < phys_const::kAlmost0F)
    nt = 1.0;
  mSeedMom[0] = pxz * ux / nt;
  mSeedMom[1] = pxz * uy / nt;
  mSeedMom[2] = pxz * uz / nt;
  float crossy = (x2 - x1) * (z3 - z2) - (z2 - z1) * (x3 - x2);
  int qSign = (crossy * by > 0) ? +1 : -1;
  // Ensure charge sign consistent with curvature direction in B-field
  if (mCharge * qSign < 0) {
    mCharge = -mCharge;
  }
  if (mSeedMom[2] < 0) { // enforce positive pz
    for (int j = 0; j < 3; ++j) {
      mSeedMom[j] = -mSeedMom[j];
    }
  }
  mIsSeedSet = true;
  return;
}

void NA6PFastTrackFitter::computeSeed(int dir)
{
  // compute track seed from 3 (or 2) clusters
  // dir = 1 -> forward dir = -1 -> backward

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
    LOGP(error, "Cannot compute seed with {} seeding layers", nSeedLayers);
    return;
  }
  mSeedOption = origSeedOption;
  computeSeed(dir, layForSeed);
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

bool NA6PFastTrackFitter::fitTrackPoints(NA6PTrack& trackToFit, int dir, const NA6PTrackParCov* seed)
{
  int nClus = getNumberOfClusters();
  if (nClus < 2) {
    LOGP(error, "Cannot fit track with only {} clusters\n", nClus);
    return false;
  }
  bool isGoodFit = true;
  auto prop = Propagator::Instance();
  trackToFit.reset();
  if (seed) {
    // Seed provided as a track from previous pass
    ((NA6PTrackParCov&)trackToFit) = *seed;
    mSeedPos = seed->getXYZ();  // RSTODO why do we need this?
    mSeedMom = seed->getPXYZ(); // RSTODO use seed as opt.linRef !!!!
    mIsSeedSet = true;
  } else {
    if (!mIsSeedSet) { // Set seed from clusters in the outer (inner) layers if not set from outside
      dir > 0 ? computeSeedInner() : computeSeedOuter();
    }
    if (mIsSeedSet) {
      trackToFit.init(mSeedPos, mSeedMom, mCharge); // RSTODO why? the seed was already assigned
    } else {
      LOGP(warn, "Track seed not computed properly, will run the fit w/o seed");
    }
  }
  trackToFit.setPID(mPID);
  trackToFit.resetCovariance(-1);
  int startLay = (dir > 0) ? 0 : mNLayers - 1, endLay = (dir > 0) ? mNLayers : -1, stepLay = (dir > 0) ? 1 : -1;
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

//=========================================================
float NA6PFastTrackFitter::fitSeed(NA6PTrackParCov& seed, bool resetCovMat, int dir, NA6PTrackPar* linRef)
{
  std::vector<const NA6PBaseCluster*> clusv;
  clusv.reserve(mNLayers); // RSTODO consider to have this as a class member filled once
  for (int i = 0; i < mNLayers; i++) {
    if (mClusters[i]) {
      clusv.push_back(mClusters[i]);
    }
  }
  return fitSeed(seed, clusv, resetCovMat, dir, linRef);
}

float NA6PFastTrackFitter::fitSeed(NA6PTrackParCov& seed, std::vector<const NA6PBaseCluster*> clusters, bool resetCovMat, int dir, NA6PTrackPar* linRef)
{
  // fit track (already seeded but not necessarily at the position of the 1st cluster to fit.
  // If the covariance matrix reset it requested, then the track position is imposed at the 1st cluster to fit (w/o propagating other track params)
  // and the covariance is reset. Otherwise (assuming that the fit is a continuation of the existing fit), the track is propagated to the 1st point
  // (including the cov. matrix propagation).
  // If linRef is provided, it is used for the KF linearization, otherwise the linearization is done wrt the track being updated.
  // The propagation is done with the PID of the provided seed (or that of linRef, if provided).
  //
  // The cluster must be provided in the increasing Z order
  //
  // In case of the sucessful fit return total chi2, otherwise negative number
  int nClus = clusters.size();
  if (nClus < 2) {
    LOGP(error, "Cannot fit track with only {} clusters\n", nClus);
    return -1.f;
  }
  float chi2Tot = 0.f;
  auto prop = Propagator::Instance();
  int startC = 0, stopC = nClus, step = 1;
  if (dir < 0) {
    startC = nClus - 1;
    stopC = -1;
    step = -1;
  }
  if (resetCovMat) {
    seed.setXYZ(clusters[startC]->getXYZ());
    seed.resetCovariance(-1);
    if (linRef) {
      if (!prop->propagateToZ(*linRef, clusters[startC]->getZ(), mPropOpt)) { // linRef must be at the same position as the seed
        return -1.f;
      }
      linRef->setXYZ(clusters[startC]->getXYZ());
    }
  }
  mPropOpt.linRef = linRef;

  int prevLr = clusters[startC]->getLayer() - step; // to check clusters ordering
  for (int ic = startC; ic != stopC; ic += step) {
    const auto& cl = *mClusters[ic];
    if (cl.getLayer() <= prevLr) {
      LOGP(fatal, "Clusters are not in increasing Z order");
    }
    prevLr = cl.getLayer();
    float chi2 = 0;
    if (!prop->propagateToZ(seed, cl.getZ(), mPropOpt) ||
        (chi2 = seed.getPredictedChi2(cl)) > mMaxChi2Cl ||
        !seed.update(cl)) {
      chi2Tot = -1;
      break;
    }
    chi2Tot += chi2;
  }
  mPropOpt.linRef = nullptr;
  return chi2Tot;
}
