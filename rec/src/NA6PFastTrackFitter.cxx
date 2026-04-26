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

int NA6PFastTrackFitter::getLayersForSeed(int dir, std::array<int, 3>& layForSeed)
{
  // define layers to be used in seed calculation depending on the options, in the order of track fit direction (dir>0 : forward)
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
    if (dir < 0) {
      std::swap(layForSeed[0], layForSeed[nLayWithClus - 1]);
    }
    return nLayWithClus;
  }
  // select layers for seed determination
  if (mSeedOption == kEdgeClusters) {
    if (dir > 0) {
      for (int i = 0; i < 3 && i < nLayWithClus; ++i) {
        layForSeed[i] = mLayersWithClusters[i];
      }
    } else {
      for (int i = 0; i < 3 && i < nLayWithClus; ++i) {
        layForSeed[i] = mLayersWithClusters[nLayWithClus - 1 - i];
      }
    }
  } else { // max lever arm
    layForSeed[0] = mLayersWithClusters[0];
    layForSeed[2] = mLayersWithClusters[nLayWithClus - 1];
    int mid = (layForSeed[0] + layForSeed[2]) / 2, best = -1, bestDist = 999; // Find closest layer to midpoint
    for (int i = 1; i < nLayWithClus - 1; ++i) {
      int dist = std::abs(mLayersWithClusters[i] - mid);
      if (dist < bestDist) {
        bestDist = dist;
        best = mLayersWithClusters[i];
      }
    }
    layForSeed[1] = best;
    if (dir < 0) {
      std::swap(layForSeed[0], layForSeed[2]);
    }
  }
  return 3;
}

std::pair<double, double> NA6PFastTrackFitter::getFieldMomenta(const std::array<float, 3>& pos, const std::array<float, 3>& df, const int nSteps) const
{
  // calcualte By 1st and 2nd momenta between 2 points defined by the start position and the difference end-start.
  const auto vStep = NA6PLine::getScaled(df, 1.f / nSteps);
  auto posFQuery = NA6PLine::getSum(pos, NA6PLine::getScaled(vStep, 0.5f));
  auto prop = Propagator::Instance();
  double mom0 = 0, mom1 = 0.;
  for (int i = 0; i < nSteps; i++) {
    float by = prop->getBy(posFQuery);
    mom0 += by;
    mom1 += by * posFQuery[2];
    NA6PLine::add(posFQuery, vStep);
  }
  return {mom0 * vStep[2], mom1 * vStep[2]};
}

bool NA6PFastTrackFitter::computeSeed(int dir, std::array<int, 3>& layForSeed, NA6PTrackPar* seed)
{
  // layers are provided in increasing order if dir>0 (forward fit) and decreasing if dir<0 (backward fit), seed is defined at 1st provided cluster position
  if (!seed) {
    seed = &mSeed;
  }
  int nSeedLayers = layForSeed[2] >= 0 ? 3 : 2;
  const auto &xyz0 = mClusters[layForSeed[0]]->getXYZ(), xyz1 = mClusters[layForSeed[1]]->getXYZ();
  seed->setXYZ(xyz0);
  const auto dpos01 = NA6PLine::getDiff(xyz0, xyz1); // c1 - c0
  if (nSeedLayers == 2) {                            // straight line seed only
    auto pxzI = (dir > 0 ? 1.f : -1.f) / std::hypot(dpos01[0], dpos01[2]);
    seed->setTx(dpos01[0] * pxzI);
    seed->setTy(dpos01[1] * pxzI);
    seed->setQ2Pxz(1.f / mMostProbableP);
    return true;
  }
  const auto& xyz2 = mClusters[layForSeed[2]]->getXYZ();
  /*
    Extraction of seed parameters in parabolic approximation for track propagation in Z-dependent By field (XZ bending):
    x(z) = x_ref + bx_ref * (z - z_ref) + beta * ( z * M0(z_ref,z) - M1(z_ref,z) )
    y(z) = y_ref + by_ref * (z - z_ref)
    with beta = (q/pXZ) * k
    M0_01, M1_01 : signed field moments from z0 to z1
    M0_02, M1_02 : signed field moments from z0 to z2

    Derivation: start from the stepwise transport equations in the XZ plane, with z_i ordered along the trajectory and constant step size
    dz = z_i - z_{i-1}, with the local state at step i is (x_i, bx_i), where the slope bx_i = dx/dz at z = z_i.

    The recursive propagation from z_{i-1} to state at z_{i} is:
    x_i  = x_{i-1} + bx_{i-1} * dz + (beta / 2) * B(z_{i-1}) * dz^2
    bx_i = bx_{i-1}+ beta * B(z_{i-1}) * dz
    with beta = (q / p_XZ) * k.

    Iterate from the reference step 0 up to step n:
    bx_1 = bx_0 + beta * B(z_0) * dz
    bx_2 = bx_1 + beta * B(z_1) * dz = bx_0 + beta * dz * [ B(z_0) + B(z_1) ]
    ...
    =>  bx_n = bx_0 + beta * dz * sum_{j=0}^{n-1} B(z_j)                                          | (1)

    Similarly, iterating position x_i from 0 to n we get:
    x_n = x_0 + dz * sum_{i=1}^{n} bx_{i-1} + (beta/2) * dz^2 * sum_{i=1}^{n} B(z_{i-1})
    Substituting bx_{i-1} by (1) and reordering:
    =>  x_n = x_0 + bx_0 * (z_n - z_0) + beta * dz^2 * sum_{j=0}^{n-1} [ n - j - 1/2 ] * B(z_j)   | (2)

    In the limit of small dz swith to integrals == moments of the field:
    x(z) = x_0 + bx_0 * (z - z_0) + beta * [ z * M0(z_0, z) - M1(z_0, z) ]                        | (3)
    with
    M0(z_0, z) = integral_{z_0}^{z} B(s) ds
    M1(z_0, z) = integral_{z_0}^{z} s * B(s) ds
  */
  const auto dpos12 = NA6PLine::getDiff(xyz1, xyz2), dpos02 = NA6PLine::getDiff(xyz0, xyz2);
  ;                                                                                    // c2 - c1 and c2 - c0
  const auto M01 = getFieldMomenta(xyz0, dpos01), M12 = getFieldMomenta(xyz1, dpos12); // field moments

  // Field-dependent regressors: F_i = z_i * M0(0,i) - M1(0,i)
  const auto F1 = xyz1[2] * M01.first - M01.second, F2 = xyz2[2] * (M01.first + M12.first) - (M01.second + M12.second);
  const auto dz1f2 = dpos01[2] * F2, dz2f1 = dpos02[2] * F1; // (z1-z0)*F2 - (z2-z0)*F1
  const auto det = dz1f2 - dz2f1, scale = std::fabs(dz1f2) + std::fabs(dz2f1) + 1.0, detI = 1. / det;
  auto bx_ref = (dpos01[0] * F2 - dpos02[0] * F1) * detI;
  float beta = (dpos01[2] * dpos02[0] - dpos02[2] * dpos01[0]) * detI;
  float qOverPXZ = beta / NA6PTrackPar::kB2C;

  // YZ fit: y_i = y0 + by0 * (z - z0)
  // by_0 = [ (z1-z0)*(y1-y0) + (z2-z0)*(y2-y0) ] / [ d1^2 + d2^2 ]
  float denY = dpos01[2] * dpos01[2] + dpos02[2] * dpos02[2];
  float by_ref = (dpos01[2] * dpos01[1] + dpos02[2] * dpos02[1]) / denY;
  seed->setQ2Pxz(qOverPXZ);
  // convert slopes bx = px/pz and by = py/pz to tx = px/pxz and ty = py/pxz
  seed->setTx(bx_ref / std::sqrt(1.f + bx_ref * bx_ref));
  seed->setTy(by_ref * seed->getCosPsi());
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
  int nSeedLayers = getLayersForSeed(dir, layForSeed);
  if (nSeedLayers < 2) {
    seed->invalidate();
    LOGP(error, "Cannot compute seed with {} seeding layers", nSeedLayers);
    return false;
  }
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
