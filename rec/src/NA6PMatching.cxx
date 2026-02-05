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

bool NA6PMatching::init(const char* filename, const char* geoname)
{
  NA6PReconstruction::init(filename, geoname);
  mTrackFitter = new NA6PFastTrackFitter();
  mTrackFitter->loadGeometry(filename, geoname);
  mTrackFitter->enableMaterialCorrections();
  mTrackFitter->setPropagateToPrimaryVertex(true);
  mTrackFitter->setNLayers(11);
  mTrackFitter->setParticleHypothesis(13); // muon mass hypothesis
  mTrackFitter->setMaxChi2Cl(mMaxChi2Refit);
  return true;
}

void NA6PMatching::configureFromRecoParam(const std::string& filename)
{
  if (filename != "") {
    na6p::conf::ConfigurableParamHelper<NA6PRecoParam>::updateFromFile(filename);
  }
  const auto& param = NA6PRecoParam::Instance();

  if (param.isZMatchingSet)
    setZMatching(param.zMatching);
  setMCMatching(param.mcMatching);
  setMinVTHits(param.minVTHits);
  setMinMSHits(param.minMSHits);
  setMinTrackP(param.minTrackP);
  setMaxChi2Match(param.maxChi2Match);
  setMaxChi2Refit(param.maxChi2Refit);
  setPMatchWindow(param.pMatchWindow);
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
    if (abs(track.getZLabOuter() - z) < 1e-6)
      continue;
    if (outer) {
      mTrackFitter->propagateToZOuter(&track, z);
    } else
      mTrackFitter->propagateToZ(&track, z);
  }
};

void NA6PMatching::addClustersToFitter(const NA6PTrack& trk, const auto* clusPtr)
{
  uint32_t cmap = trk.getClusterMap();
  for (int lr = 0; lr < NA6PTrack::kMaxLr; ++lr) {
    if (cmap & (1u << lr)) {
      mTrackFitter->addCluster(lr, (*clusPtr)[trk.getClusterIndex(lr)]);
    }
  }
}


std::tuple<int, bool, double> NA6PMatching::findBestChi2Match(
  const NA6PTrack& msTrack,
  const std::vector<size_t>& validVerTelIndices)
{
  const auto msPar = msTrack.getTrackExtParam().getParameter();
  const auto msCov = msTrack.getTrackExtParam().getCovariance();
  const int msCharge = msTrack.getCharge();
  const float msP = msTrack.getP();
  const int msPartId = msTrack.getParticleID();

  auto lower = std::partition_point(validVerTelIndices.begin(), validVerTelIndices.end(),
    [&](size_t vIdx) { return mVerTelTracks[vIdx].getPOuter() > msP + mPMatchWindow; });
  
  auto upper = std::partition_point(lower, validVerTelIndices.end(),
    [&](size_t vIdx) { return mVerTelTracks[vIdx].getPOuter() >= msP - mPMatchWindow; });
  

  float bestChi2 = mMaxChi2Match;
  int bestIdx = -1;
  bool isTrueMatch = false;

  for (auto it = lower; it != upper; ++it) {
    size_t vIdx = *it;
    const auto& vtTrack = mVerTelTracks[vIdx];

    if (vtTrack.getCharge() * msCharge < 0)
      continue;

    const auto vtPar = vtTrack.getOuterParam().getParameter();
    const auto vtCov = vtTrack.getOuterParam().getCovariance();
    const double chi2 = computeChi2(msPar, msCov, vtPar, vtCov);

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
  vtByPid.reserve(mVerTelTracks.size());

  for (size_t i = 0; i < mVerTelTracks.size(); ++i) {
    const auto& vt = mVerTelTracks[i];
    int pid = vt.getParticleID();
    if (pid <= 0)
      continue;

    // Prefer track with more hits if there are duplicates with same particle ID
    auto it = vtByPid.find(pid);
    if (it == vtByPid.end() || vt.getNHits() > mVerTelTracks[it->second].getNHits()) {
      vtByPid[pid] = static_cast<int>(i);
    }
  }
  return vtByPid;
}

std::vector<size_t> NA6PMatching::prefilterVerTelTracks()
{
  std::vector<size_t> validIndices;
  validIndices.reserve(mVerTelTracks.size());
  sortVTTracksByP(mVerTelTracks);

  for (size_t vIdx = 0; vIdx < mVerTelTracks.size(); ++vIdx) {
    const auto& track = mVerTelTracks[vIdx];
    if (track.getNHits() >= mMinVTHits && track.getP() >= mMinTrackP) {
      validIndices.push_back(vIdx);
    }
  }
  return validIndices;
}

bool NA6PMatching::fitAndStoreMatchedTrack(const NA6PTrack& vtTrk, const NA6PTrack& msTrk, int particleId, double matchChi2)
{
  mTrackFitter->cleanupAndStartFit();

  addClustersToFitter(vtTrk, &mVerTelClusters);
  addClustersToFitter(msTrk, hMuonSpecClusPtr);

  std::unique_ptr<NA6PTrack> matchedTrack(mTrackFitter->fitTrackPoints());
  if (!matchedTrack) {
    LOGP(warn, "Refit of matched track failed, skipping");
    return false;
  }

  matchedTrack->setParticleID(particleId);
  matchedTrack->setMatchChi2(matchChi2);
  mTrackFitter->propagateToZ(matchedTrack.get(), mPrimaryVertex.Z());
  mMatchedTracks.push_back(*matchedTrack);

  return true;
}


void NA6PMatching::runMCMatching()
{
  auto vtByPid = buildMCMatchingIndex();

  for (const auto& msTrack : mMuonSpecTracks) {
    if (msTrack.getNHits() < mMinMSHits)
      continue;

    int msPartId = msTrack.getParticleID();
    if (msPartId <= 0)
      continue;

    auto it = vtByPid.find(msPartId);
    if (it == vtByPid.end())
      continue;

    int vtIdx = it->second;
    const auto& vtTrack = mVerTelTracks[vtIdx];
    fitAndStoreMatchedTrack(vtTrack, msTrack, msPartId, 0.0);
  }
}

void NA6PMatching::runDataMatching()
{
  if (mIsZMatchingSet) {
    propToZMatching(mVerTelTracks, mZMatching, true);
    propToZMatching(mMuonSpecTracks, mZMatching);
  }

  auto validVerTelIndices = prefilterVerTelTracks();

  for (const auto& msTrack : mMuonSpecTracks) {
    if (msTrack.getNHits() < mMinMSHits)
      continue;

    auto [vtIdx, isTrueMatch, bestChi2] = findBestChi2Match(msTrack, validVerTelIndices);
    if (vtIdx < 0)
      continue;

    const auto& vtTrack = mVerTelTracks[vtIdx];
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
  mMatchedTracks.reserve(mMuonSpecTracks.size());

  if (mMCMatching) {
    runMCMatching();
  } else {
    runDataMatching();
  }

  writeTracks();
}