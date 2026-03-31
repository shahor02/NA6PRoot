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
#include "Seeder.h"

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

bool NA6PFastTrackFitter::constrainTrackToVertex(NA6PTrack& trc, const NA6PVertex& pv)
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
  mChi2Buffer.push_back(chi2);
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
  if (nSeedLayers == 2) {                              // straight line seed only
    const auto dpos01 = NA6PLine::getDiff(xyz0, xyz1); // c1 - c0
    seed->setXYZ(xyz0);
    auto pxzI = (dir > 0 ? 1.f : -1.f) / std::hypot(dpos01[0], dpos01[2]);
    seed->setTx(dpos01[0] * pxzI);
    seed->setTy(dpos01[1] * pxzI);
    seed->setQ2Pxz(1.f / mMostProbableP);
    return true;
  }
  const auto& xyz2 = mClusters[layForSeed[2]]->getXYZ();

  auto mcTruth = getMCTruthStatus(); // RSREM
  LOGP(debug, "seed Dir{} L:{}/{}/{} MC status: TSame={} TPID={}", dir, layForSeed[0], layForSeed[1], layForSeed[2], mcTruth.second, mcTruth.first);

  return Seeder::create(*seed, dir > 0, xyz0, xyz1, xyz2, mSeedImprovePrec);
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

float NA6PFastTrackFitter::fitSeed(NA6PTrackParCov& seed, bool resetCovMat, int dir, bool useLinRef)
{
  // fit track (already seeded but not necessarily at the position of the 1st cluster to fit.
  // If the covariance matrix reset it requested, then the track position is imposed at the 1st cluster to fit (w/o propagating other track params)
  // and the covariance is reset. Otherwise (assuming that the fit is a continuation of the existing fit), the track is propagated to the 1st point
  // (including the cov. matrix propagation).
  // If useLinRef is true, the clone of the seed is used for the KF linearization, otherwise the linearization is done wrt the track being updated.
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
    seed.resetCovariance(NA6PTrackParCov::MaxErrSelRescale);
  }
  NA6PTrackPar linRefLoc(seed), *linRef = useLinRef ? &linRefLoc : nullptr;
  mPropOpt.linRef = useLinRef ? &linRefLoc : nullptr;

  auto mcTruth = getMCTruthStatus(); // RSREM

  int cntCl = 0;
  for (int il = startL; il != stopL; il += step) {
    if (!mClusters[il]) {
      continue;
    }
    const auto& cl = *mClusters[il];
    float chi2 = 0;
    if (!prop->propagateToZ(seed, cl.getZ(), mPropOpt)) {
      LOGP(debug, "DBG fitSeed fail: propagate il={} z={} truthSame={} truthPID={} state={}", il, cl.getZ(), mcTruth.second, mcTruth.first, seed.asString());
      chi2Tot = -1;
      break;
    }
    if (mcTruth.second) {
      LOGP(debug, "DBGR Before update {}/{}: {}", cntCl, mNClusters, seed.asString());
    }
    chi2 = seed.getPredictedChi2(cl);
    mChi2Buffer.push_back(chi2);
    if (mcTruth.second) {
      LOGP(debug, "DBGR chi2={} for cluster {}", chi2, cl.asString());
    }
    if (chi2 > mMaxChi2Cl) {
      const auto dx = seed.getX() - cl.getX();
      const auto dy = seed.getY() - cl.getY();
      const auto exx = seed.getSigmaX2() + cl.getSigXX();
      const auto eyx = seed.getSigmaYX() + cl.getSigYX();
      const auto eyy = seed.getSigmaY2() + cl.getSigYY();
      const auto det = exx * eyy - eyx * eyx;
      LOGP(debug, "DBG fitSeed fail chi2 at {}/{} il={} chi2={:.3f} max={} truthSame={} truthPID={} clPID={} clXYZ=({:.3e},{:.3e},{:.3e}) state={}", cntCl, mNClusters,
           il, chi2, mMaxChi2Cl, mcTruth.second, mcTruth.first, cl.getParticleID(), cl.getX(), cl.getY(), cl.getZ(), seed.asString());
      LOGP(debug, "DBG predChi2 new: frame=labXY il={} truthSame={} truthPID={} clPID={} dx={} dy={} trXY=({},{}) clXY=({},{}) trCov=({:.3e},{:.3e},{:.3e}) clCov=({:.3e},{:.3e},{:.3e}) sumcov=({:.3e},{:.3e},{:.3e}) det={} chi2={}",
           il, mcTruth.second, mcTruth.first, cl.getParticleID(), dx, dy, seed.getX(), seed.getY(), cl.getX(), cl.getY(), seed.getSigmaX2(), seed.getSigmaYX(), seed.getSigmaY2(),
           cl.getSigXX(), cl.getSigYX(), cl.getSigYY(), exx, eyx, eyy, det, chi2);
      chi2Tot = -1;
      break;
    } else if (!seed.update(cl)) {
      LOGP(debug, "DBG fitSeed fail: update il={} chi2={:.3f} truthSame={} truthPID={} clPID={} clXYZ=({:.3e},{:.3e},{:.3e}) state={}",
           il, chi2, mcTruth.second, mcTruth.first, cl.getParticleID(), cl.getX(), cl.getY(), cl.getZ(), seed.asString());
      chi2Tot = -1;
      break;
    } else if (mcTruth.second) {
      LOGP(debug, "DBG fitSeed OK chi2 at {}/{} il={} chi2={:.3f}/{:.3f} max={} truthSame={} truthPID={} clPID={} clXYZ=({:.3e},{:.3e},{:.3e}) state={}", cntCl, mNClusters,
           il, chi2, chi2 + chi2Tot, mMaxChi2Cl, mcTruth.second, mcTruth.first, cl.getParticleID(), cl.getX(), cl.getY(), cl.getZ(), seed.asString());
    }
    if (mcTruth.second) {
      LOGP(debug, "DBGR After update {}/{}: {}", cntCl, mNClusters, seed.asString());
    }
    chi2Tot += chi2;
    cntCl++;
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
      tr.addCluster(cl);
    }
  }
}

std::pair<int, bool> NA6PFastTrackFitter::getMCTruthStatus()
{
  std::vector<std::pair<int, int>> occurrences;
  for (int jl = mMinLayerWithCl; jl <= mMaxLayerWithCl; ++jl) {
    if (!mClusters[jl]) {
      continue;
    }
    const int pid = mClusters[jl]->getParticleID();
    bool found{false};
    for (int i = 0; i < (int)occurrences.size(); i++) {
      if (pid == occurrences[i].first) {
        occurrences[i].second++;
        found = true;
        break;
      }
    }
    if (!found) {
      occurrences.emplace_back(pid, 1);
    }
  }
  std::sort(std::begin(occurrences), std::end(occurrences), [](auto e1, auto e2) { return e1.second > e2.second; });
  return {occurrences.front().second == mNClusters, occurrences.front().first};
}
