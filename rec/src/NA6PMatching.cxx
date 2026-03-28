// NA6PCCopyright

#include "NA6PMatching.h"
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <unordered_map>
#include <fairlogger/Logger.h>

#include <TMatrixD.h>
#include <TVectorD.h>

ClassImp(NA6PMatching)

  NA6PMatching::NA6PMatching() : NA6PReconstruction("Matching")
{
}

NA6PMatching::NA6PMatching(const char* recparfile,
                           const char* geofile,
                           const char* geoname) : NA6PReconstruction("Matching")
{
  mGeoFilName = geofile;
  mGeoObjName = geoname;
  mRecoParFilName = recparfile;
  initMatching();
}

bool NA6PMatching::initMatching()
{
  NA6PReconstruction::init(mGeoFilName.c_str(), mGeoObjName.c_str());
  mTrackFitter = new NA6PFastTrackFitter();
  mTrackFitter->loadGeometry(mGeoFilName.c_str(), mGeoObjName.c_str());
  mTrackFitter->enableMaterialCorrections();
  mTrackFitter->setPropagateToPrimaryVertex(false);
  mTrackFitter->setNLayers(11);
  mTrackFitter->setParticleHypothesis(PID::Muon); // muon mass hypothesis
  mTrackFitter->setUseIntegralBForSeed();
  configureFromRecoParam(mRecoParFilName);
  createTracksOutput();
  return true;
}

void NA6PMatching::configureFromRecoParam(const std::string& filename)
{
  if (filename != "") {
    na6p::conf::ConfigurableParamHelper<NA6PRecoParam>::updateFromFile(filename);
  }
  const auto& param = NA6PRecoParam::Instance();

  if (param.mtIsZMatchingSet)
    setZMatching(param.mtZMatching);
  setMCMatching(param.mtMCMatching);
  setMinVTHits(param.mtMinVTHits);
  setMinMSHits(param.mtMinMSHits);
  setMinTrackP(param.mtMinTrackP);
  setMaxChi2Match(param.mtMaxChi2Match);
  setMaxChi2Refit(param.mtMaxChi2Refit);
  setPMatchWindow(param.mtPMatchWindow);
  setPropagateTracksToPrimaryVertex(param.mtPropagateMatchedTracksToPV);
  setDoOutwardInwardFit(param.mtDoOutwardInwardFit);
}

void NA6PMatching::createTracksOutput()
{
  auto nm = fmt::format("Tracks{}.root", getName());
  mMatchedTrackFile = TFile::Open(nm.c_str(), "recreate");
  mMatchedTrackTree = new TTree(fmt::format("tracks{}", getName()).c_str(), fmt::format("{} Tracks", getName()).c_str());
  mMatchedTrackTree->Branch(getName().c_str(), &hMatchedTrackPtr);
  LOGP(info, "Will store {} tracks in {}", getName(), nm);
}

void NA6PMatching::writeTracks()
{
  if (mMatchedTrackTree) {
    mMatchedTrackTree->Fill();
    LOGP(info, "Saved {} tracks in tree with {} entries", mMatchedTracks.size(), mMatchedTrackTree->GetEntries());
  }
}

void NA6PMatching::closeTracksOutput()
{
  if (mMatchedTrackTree && mMatchedTrackFile) {
    mMatchedTrackFile->cd();
    mMatchedTrackTree->Write();
    delete mMatchedTrackTree;
    mMatchedTrackTree = nullptr;
    mMatchedTrackFile->Close();
    delete mMatchedTrackFile;
    mMatchedTrackFile = nullptr;
  }
}

void NA6PMatching::setVerTelTracks(std::vector<NA6PTrack>& tracks)
{
  hVerTelTrackPtr = &tracks;
}

void NA6PMatching::setMuonSpecTracks(std::vector<NA6PTrack>& tracks)
{
  hMuonSpecTrackPtr = &tracks;
}

void NA6PMatching::setMuonSpecClusters(std::vector<NA6PMuonSpecCluster>& clusters)
{
  hMuonSpecClusPtr = &clusters;
  for (size_t i = 0; i < clusters.size(); ++i) {
    clusters[i].setClusterIndex(static_cast<int>(i));
  }
}

void NA6PMatching::setVerTelClusters(std::vector<NA6PVerTelCluster>& clusters)
{
  hVerTelClusPtr = &clusters;
  for (size_t i = 0; i < clusters.size(); ++i) {
    clusters[i].setClusterIndex(static_cast<int>(i));
  }
}

void NA6PMatching::sortVTTracksByP(std::vector<NA6PTrack>& tracks)
{
  std::sort(tracks.begin(), tracks.end(),
            [](const NA6PTrack& a, const NA6PTrack& b) {
              return a.getP() > b.getP();
            });
}

double NA6PMatching::computeChi2(const double* par1, const double* cov1,
                                 const double* par2, const double* cov2)
{
  TVectorD delta(5);
  for (int i = 0; i < 5; ++i) {
    delta[i] = par1[i] - par2[i];
  }

  TMatrixDSym V(5);
  int idx = 0;
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      V(i, j) = cov1[idx] + cov2[idx];
      ++idx;
    }
  }

  TMatrixDSym Vinv = V.Invert();
  if (!V.IsValid())
    return 1e9;

  return Vinv.Similarity(delta);
}

void NA6PMatching::propToZMatching(std::vector<NA6PTrack>& tracks, double z, bool outer)
{
  for (auto& track : tracks) {
    if (abs(track.getOuterParam().getZ() - z) < 1e-6) {
      continue;
    }
    if (outer) { // RSTODO results of the propagation is not checked
      Propagator::Instance()->propagateToZ(track.getOuterParam(), z, mTrackFitter->getPropOpt());
    } else {
      Propagator::Instance()->propagateToZ(track, z, mTrackFitter->getPropOpt());
    }
  }
};

void NA6PMatching::addClustersToFitter(const NA6PTrack& trk, const auto* clusPtr)
{
  uint32_t cmap = trk.getClusterMap();
  for (int lr = 0; lr < NA6PTrack::kMaxLr; ++lr) {
    if (cmap & (1u << lr)) {
      mTrackFitter->addCluster((*clusPtr)[trk.getClusterIndex(lr)]);
    }
  }
}

std::tuple<int, bool, double> NA6PMatching::findBestChi2Match(
  const NA6PTrack& msTrack,
  const std::vector<size_t>& validVerTelIndices)
{
  const int msCharge = msTrack.getCharge(), msPartId = msTrack.getParticleID();
  const float msP = msTrack.getP();
  // RSTODO why using P? The natural would be to use getQ2Pxz: measured momentum in the bending plane, or at least getQ2P()
  auto lower = std::partition_point(validVerTelIndices.begin(), validVerTelIndices.end(),
                                    [&](size_t vIdx) { return (*hVerTelTrackPtr)[vIdx].getOuterParam().getP() > msP + mPMatchWindow; });

  auto upper = std::partition_point(lower, validVerTelIndices.end(),
                                    [&](size_t vIdx) { return (*hVerTelTrackPtr)[vIdx].getOuterParam().getP() >= msP - mPMatchWindow; });

  float bestChi2 = mMaxChi2Match;
  int bestIdx = -1;
  bool isTrueMatch = false;

  for (auto it = lower; it != upper; ++it) {
    size_t vIdx = *it;
    const auto& vtTrack = (*hVerTelTrackPtr)[vIdx];
    if (vtTrack.getCharge() != msCharge) {
      continue;
    }
    // RSTODO this is quite expensive operation, consider doing preliminare rejections based on the abs. differences or pulls
    const auto chi2 = msTrack.getPredictedChi2(vtTrack.getOuterParam());
    if (chi2 < bestChi2) {
      bestChi2 = chi2;
      bestIdx = vIdx;
      isTrueMatch = (msPartId > 0 && vtTrack.getParticleID() == msPartId);
    }
  }

  return {bestIdx, isTrueMatch, bestChi2};
}

std::unordered_map<int, int> NA6PMatching::buildMCMatchingIndex()
{
  std::unordered_map<int, int> vtByPid;
  vtByPid.reserve(hVerTelTrackPtr->size());

  for (size_t i = 0; i < hVerTelTrackPtr->size(); ++i) {
    const auto& vt = (*hVerTelTrackPtr)[i];
    int pid = vt.getParticleID();
    if (pid <= 0 || vt.getNHits() < mMinVTHits) {
      continue;
    }
    // Prefer track with more hits if there are duplicates with same particle ID
    auto it = vtByPid.find(pid);
    if (it == vtByPid.end() || vt.getNHits() > (*hVerTelTrackPtr)[it->second].getNHits()) {
      vtByPid[pid] = static_cast<int>(i);
    }
  }
  return vtByPid;
}

std::vector<size_t> NA6PMatching::prefilterVerTelTracks()
{
  std::vector<size_t> validIndices;
  validIndices.reserve(hVerTelTrackPtr->size());
  sortVTTracksByP(*hVerTelTrackPtr);

  for (size_t vIdx = 0; vIdx < hVerTelTrackPtr->size(); ++vIdx) {
    const auto& track = (*hVerTelTrackPtr)[vIdx];
    if (track.getNHits() >= mMinVTHits && track.getP() >= mMinTrackP) {
      validIndices.push_back(vIdx);
    }
  }
  return validIndices;
}

bool NA6PMatching::fitAndStoreMatchedTrack(const NA6PTrack& vtTrk,
                                           const NA6PTrack& msTrk,
                                           int particleId,
                                           double matchChi2)
{
  mTrackFitter->cleanupAndStartFit();

  addClustersToFitter(vtTrk, hVerTelClusPtr);
  addClustersToFitter(msTrk, hMuonSpecClusPtr);

  NA6PTrack matchedTrack;
  NA6PTrack seed;
  NA6PTrack* seedPtr = nullptr;

  if (mDoOutwardInwardFit) {
    // RSTODO here again, we copy the whole track just for the refit, consider using only TrackParCov part
    NA6PTrack vtSeed = vtTrk;
    if (!mTrackFitter->fitTrackPointsOutward(seed, &vtSeed)) {
      LOGP(warn, "Fit of matched track failed, skipping");
      return false;
    }
    seedPtr = &seed;
  }

  if (!mTrackFitter->fitTrackPointsInward(matchedTrack, seedPtr)) {
    LOGP(warn, "Refit of matched track failed, skipping");
    return false;
  }

  matchedTrack.setParticleID(particleId);
  matchedTrack.setMatchChi2(matchChi2);

  // RSTODO is not this done internally in the mTrackFitter?
  if (mPropagateTracksToPrimaryVertex && mPrimaryVertex) {
    Propagator::Instance()->propagateToZ(matchedTrack, mPrimaryVertex->getZ(), mTrackFitter->getPropOpt());
  }

  mMatchedTracks.push_back(std::move(matchedTrack));
  return true;
}

void NA6PMatching::runMCMatching()
{
  auto vtByPid = buildMCMatchingIndex();
  for (const auto& msTrack : *hMuonSpecTrackPtr) {
    if (msTrack.getNHits() < mMinMSHits)
      continue;

    int msPartId = msTrack.getParticleID();
    if (msPartId <= 0)
      continue;

    auto it = vtByPid.find(msPartId);
    if (it == vtByPid.end())
      continue;

    int vtIdx = it->second;
    const auto& vtTrack = (*hVerTelTrackPtr)[vtIdx];
    fitAndStoreMatchedTrack(vtTrack, msTrack, msPartId, 0.0);
  }
}

void NA6PMatching::runDataMatching()
{
  if (mIsZMatchingSet) {
    propToZMatching(*hVerTelTrackPtr, mZMatching, true);
    propToZMatching(*hMuonSpecTrackPtr, mZMatching);
  }

  auto validVerTelIndices = prefilterVerTelTracks();

  for (const auto& msTrack : *hMuonSpecTrackPtr) {
    if (msTrack.getNHits() < mMinMSHits)
      continue;

    auto [vtIdx, isTrueMatch, bestChi2] = findBestChi2Match(msTrack, validVerTelIndices);
    if (vtIdx < 0)
      continue;

    const auto& vtTrack = (*hVerTelTrackPtr)[vtIdx];
    int matchedPartID = isTrueMatch ? msTrack.getParticleID() : -std::abs(msTrack.getParticleID());
    fitAndStoreMatchedTrack(vtTrack, msTrack, matchedPartID, bestChi2);
  }
}

void NA6PMatching::runMatching()
{
  if (!mIsInitialized) {
    LOGP(error, "Magnetic field and geometry not initialized");
    return;
  }
  clearTracks();
  mMatchedTracks.reserve(hMuonSpecTrackPtr ? hMuonSpecTrackPtr->size() : 0);
  LOGP(info, "Process event with nVTTracks {} nMSTracks {}, primary vertex in z = {} cm",
       hVerTelTrackPtr ? hVerTelTrackPtr->size() : 0,
       hMuonSpecTrackPtr ? hMuonSpecTrackPtr->size() : 0,
       mPrimaryVertex ? mPrimaryVertex->getZ() : -999.);

  if (mMCMatching) {
    runMCMatching();
  } else {
    runDataMatching();
  }

  writeTracks();
}
