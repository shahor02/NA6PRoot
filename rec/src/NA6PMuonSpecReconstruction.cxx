// NA6PCCopyright

#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <fairlogger/Logger.h>
#include "NA6PMuonSpecModularHit.h"
#include "NA6PTrackerCA.h"
#include "NA6PVertexerTracklets.h"
#include "NA6PMuonSpecReconstruction.h"
#include "NA6PLayoutParam.h"

ClassImp(NA6PMuonSpecReconstruction)

  NA6PMuonSpecReconstruction::NA6PMuonSpecReconstruction() : NA6PReconstruction("MuonSpec")
{
  initTracker();
}

NA6PMuonSpecReconstruction::~NA6PMuonSpecReconstruction()
{
  // move here since the NA6PTrackerCA was fwd-declared
}

bool NA6PMuonSpecReconstruction::initTracker()
{
  mMSTracker = std::make_unique<NA6PTrackerCA>();
  mMSTracker->configureFromRecoParamMS();
  const auto& param = NA6PRecoParam::Instance();
  if (param.msDoTrackMSTrackletMID) {
    mMSTracker->setNLayers(4);
    mMSTracker->setStartLayer(param.vtNLayers);
    for (int j = 0; j < mMSTracker->getNIterations(); ++j) {
      if (mMSTracker->getMinimumNumberOfClusters(j) > 4) {
        mMSTracker->setMinimumNumberOfClusters(j, 4);
      }
    }
    mMSTracker->setDoOutwardPropagation(true);
    mMSTracker->setZForOutwardPropagation(param.msZForMSMIDmatch);
  }
  createTracksOutput();
  return true;
}

void NA6PMuonSpecReconstruction::createClustersOutput()
{
  useOwnedClusterStorage(); // to write into mClusters
  auto nm = fmt::format("Clusters{}.root", getName());
  mClusFile = TFile::Open(nm.c_str(), "recreate");
  mClusTree = new TTree(fmt::format("clusters{}", getName()).c_str(), fmt::format("{} Clusters", getName()).c_str());
  mClusTree->Branch(getName().c_str(), &hClusPtr);
  LOGP(info, "Will store {} clusters in {}", getName(), nm);
}

void NA6PMuonSpecReconstruction::writeClusters()
{
  if (mClusTree) {
    mClusTree->Fill();
  }
  LOGP(info, "Saved {} clusters in tree with {} entries", mClusters.size(), mClusTree->GetEntries());
}

void NA6PMuonSpecReconstruction::closeClustersOutput()
{
  if (mClusTree && mClusFile) {
    mClusFile->cd();
    mClusTree->Write();
    delete mClusTree;
    mClusTree = nullptr;
    mClusFile->Close();
    delete mClusFile;
    mClusFile = nullptr;
  }
}

void NA6PMuonSpecReconstruction::setClusters(std::vector<NA6PMuonSpecCluster>& clusters)
{
  hClusPtr = &clusters;
  // Assign a transient index based on position in the vector
  // to be used for storing the cluster indices in the track
  for (size_t i = 0; i < clusters.size(); ++i) {
    clusters[i].setClusterIndex(static_cast<int>(i));
  }
}

void NA6PMuonSpecReconstruction::hitsToRecPoints(const std::vector<NA6PMuonSpecModularHit>& hits)
{
  int nHits = hits.size();
  const auto& layout = NA6PLayoutParam::Instance();

  for (int jHit = 0; jHit < nHits; ++jHit) {
    const auto& hit = hits[jHit];
    double x = hit.getX();
    double y = hit.getY();
    double z = hit.getZ();
    double ex2clu = 5.e-4;
    double ey2clu = 5.e-4;
    if (mCluResX > 0) {
      x = gRandom->Gaus(hit.getX(), mCluResX);
      ex2clu = mCluResX * mCluResX;
    }
    if (mCluResY > 0) {
      y = gRandom->Gaus(hit.getY(), mCluResY);
      ey2clu = mCluResY * mCluResY;
    }
    // very rough cluster size settings
    int clusiz = 1;
    int nDet = hit.getDetectorID();
    int idPart = hit.getTrackID();
    // Set layer based on Z position - find closest plane

    float minDist = std::abs(z - layout.posMSPlaneZ[0]);
    int layer = layout.nVerTelPlanes;

    for (int i = 1; i < layout.nMSPlanes; ++i) {
      float dist = std::abs(z - layout.posMSPlaneZ[i]);
      if (dist < minDist) {
        minDist = dist;
        layer = i + layout.nVerTelPlanes;
      }
    }
    mClusters.emplace_back(x, y, z, clusiz, layer);
    auto& clu = mClusters.back();
    clu.setErr(ex2clu, 0., ey2clu);
    clu.setDetectorID(nDet);
    clu.setParticleID(idPart);
    clu.setHitID(jHit);
  }
}

//____________________________________________________________________________________

void NA6PMuonSpecReconstruction::createTracksOutput()
{
  auto nm = fmt::format("Tracks{}.root", getName());
  mTrackFile = TFile::Open(nm.c_str(), "recreate");
  mTrackTree = new TTree(fmt::format("tracks{}", getName()).c_str(), fmt::format("{} Tracks", getName()).c_str());
  mTrackTree->Branch(getName().c_str(), &hTrackPtr);
  LOGP(info, "Will store {} tracks in {}", getName(), nm);
}

void NA6PMuonSpecReconstruction::writeTracks()
{
  if (mTrackTree) {
    mTrackTree->Fill();
  }
  LOGP(info, "Saved {} tracks in tree with {} entries", mTracks.size(), mTrackTree->GetEntries());
}

void NA6PMuonSpecReconstruction::closeTracksOutput()
{
  if (mTrackTree && mTrackFile) {
    mTrackFile->cd();
    mTrackTree->Write();
    delete mTrackTree;
    mTrackTree = nullptr;
    mTrackFile->Close();
    delete mTrackFile;
    mTrackFile = nullptr;
  }
}

void NA6PMuonSpecReconstruction::runTracking()
{
  clearTracks();
  const auto& param = NA6PRecoParam::Instance();
  mMSTracker->findTracks(getClusters(), mPrimaryVertex);
  if (param.msDoTrackMSTrackletMID == false) {
    mTracks = mMSTracker->getTracks();
    for (auto& t : mTracks) {
      t.setStatusMS(NA6PTrack::kFullMSMID);
    }
  } else {
    runMSTrackMIDTrackletMatching();
  }
  writeTracks();
}

void NA6PMuonSpecReconstruction::runMSTrackMIDTrackletMatching()
{
  const auto& param = NA6PRecoParam::Instance();
  const auto& layout = NA6PLayoutParam::Instance();
  auto& clusters = getClusters();

  std::vector<NA6PTrack> trks = mMSTracker->getTracks();
  int nTrks = trks.size();
  int firstMID = layout.nVerTelPlanes + layout.nMSPlanes - 2;
  std::vector<std::pair<NA6PMuonSpecCluster, NA6PMuonSpecCluster>> trkltsMID = mMSTracker->findTracklets(firstMID, firstMID + 1, clusters, mPrimaryVertex);
  int nTrklets = trkltsMID.size();
  NA6PFastTrackFitter* fitter = mMSTracker->getTrackFitter();

  // create lookup table of cluster indices
  std::vector<int> clusterLookup(clusters.size(), -1);
  for (size_t i = 0; i < clusters.size(); ++i) {
    int originalID = clusters[i].getClusterIndex();
    if (originalID >= 0 && (size_t)originalID < clusters.size())
      clusterLookup[originalID] = static_cast<int>(i);
  }

  // matching loop
  std::vector<MatchCandidate> candidates;
  for (int jT = 0; jT < nTrks; jT++) {
    NA6PTrack tr = trks[jT];
    for (int jS = 0; jS < nTrklets; jS++) {
      const auto& tracklet = trkltsMID[jS];
      NA6PMuonSpecCluster clu1 = tracklet.first;
      NA6PMuonSpecCluster clu2 = tracklet.second;
      float xyzClu1[3] = {clu1.getX(), clu1.getY(), clu1.getZ()};
      float xyzClu2[3] = {clu2.getX(), clu2.getY(), clu2.getZ()};

      if (xyzClu1[2] > xyzClu2[2]) {
        float tmp[3] = {xyzClu1[0], xyzClu1[1], xyzClu1[2]};
        for (int j = 0; j < 3; ++j) {
          xyzClu1[j] = xyzClu2[j];
          xyzClu2[j] = tmp[j];
        }
      }
      float zToProp = xyzClu1[2];
      float dirSegm[3];
      float norm = 0.f;
      for (int j = 0; j < 3; ++j) {
        dirSegm[j] = xyzClu2[j] - xyzClu1[j];
        norm += dirSegm[j] * dirSegm[j];
      }
      norm = std::sqrt(norm);
      for (int j = 0; j < 3; ++j)
        dirSegm[j] /= norm;
      auto trOuter = tr.getOuterParam();
      if (!Propagator::Instance()->propagateToZ(trOuter, zToProp, fitter->getPropOpt())) {
        continue;
      }
      auto xyzTr = trOuter.getXYZ<float>();
      auto pxyzTr = trOuter.getPXYZ<float>();
      norm = 0.f;
      for (int j = 0; j < 3; ++j)
        norm += pxyzTr[j] * pxyzTr[j];
      norm = std::sqrt(norm);
      for (int j = 0; j < 3; ++j)
        pxyzTr[j] /= norm;
      float distXY = std::sqrt((xyzClu1[0] - xyzTr[0]) * (xyzClu1[0] - xyzTr[0]) + (xyzClu1[1] - xyzTr[1]) * (xyzClu1[1] - xyzTr[1]));
      float costh = dirSegm[0] * pxyzTr[0] + dirSegm[1] * pxyzTr[1] + dirSegm[2] * pxyzTr[2];
      if (distXY < param.msMaxDistTrackMSTrackletMID && costh > param.msMinCosThetaTrackMSTrackletMID) {
        costh = std::clamp(costh, -1.f, 1.f);
        float score = distXY / 2. + (1 - costh) / 1.e-4;
        candidates.push_back({jT, jS, score});
      }
    }
  }
  std::sort(candidates.begin(), candidates.end(),
            [](auto& a, auto& b) { return a.score < b.score; });

  std::vector<bool> trackUsed(nTrks, false);
  std::vector<bool> segmentUsed(nTrklets, false);
  std::vector<int> trackMatch(nTrks, -1);

  for (auto& c : candidates) {
    if (trackUsed[c.track])
      continue;
    if (segmentUsed[c.segment])
      continue;
    trackMatch[c.track] = c.segment;
    trackUsed[c.track] = true;
    segmentUsed[c.segment] = true;
  }

  for (int jT = 0; jT < nTrks; jT++) {
    NA6PTrack tr = trks[jT];
    int jS = trackMatch[jT];
    if (jS < 0) {
      // track not matched to MID: save with specific status flag
      tr.setStatusMS(NA6PTrack::kMSNotMatchedToMID);
      mTracks.push_back(tr);
      continue;
    }
    const auto& tracklet = trkltsMID[jS];
    NA6PMuonSpecCluster clu1 = tracklet.first;
    NA6PMuonSpecCluster clu2 = tracklet.second;
    fitter->cleanupAndStartFit();
    NA6PTrack outTr = tr;
    static_cast<NA6PTrackParCov&>(outTr) = tr.getOuterParam();
    outTr.setChi2(tr.getChi2Outer());
    float zToProp = param.msZForMSMIDmatch;
    if (!Propagator::Instance()->propagateToZ(outTr, zToProp, fitter->getPropOpt())) {
      // track not refitted: save with specific status flag
      tr.setStatusMS(NA6PTrack::kMSMatchedToMIDnotRefitted);
      mTracks.push_back(tr);
      continue;
    }
    fitter->addCluster(clu1);
    fitter->addCluster(clu2);
    outTr.resetClusters();
    if (!fitter->fitSeedOutward(outTr, false)) {
      // track not refitted: save with specific status flag
      tr.setStatusMS(NA6PTrack::kMSMatchedToMIDnotRefitted);
      mTracks.push_back(tr);
      continue;
    }

    NA6PTrack refitInw = outTr;
    // add MS clusters for refit inward
    for (int jl = 0; jl < 4; ++jl) {
      int originalID = tr.getClusterIndex(jl + param.vtNLayers);
      if (originalID >= 0 && (size_t)originalID < clusterLookup.size()) {
        int jNewPos = clusterLookup[originalID];
        if (jNewPos < 0) {
          LOGP(error, "Cluster originalID={} not found for track {} layer {}", originalID, jT, jl);
          continue;
        }
        const auto& cl = clusters[jNewPos];
        fitter->addCluster(cl);
      }
    }
    float chi2Refit = fitter->fitSeedInward(refitInw, true);
    if (chi2Refit < 0.f) {
      // track not refitted: save with specific status flag
      tr.setStatusMS(NA6PTrack::kMSMatchedToMIDnotRefitted);
      mTracks.push_back(tr);
      continue;
    }
    refitInw.setChi2(chi2Refit);
    refitInw.setOuterParam(outTr);
    refitInw.resetClusters();
    fitter->addClustersToTrack(refitInw);
    if (param.vtDoConstrainedTrack) {
      fitter->constrainTrackToVertex(refitInw, *mPrimaryVertex);
    } else {
      refitInw.getVertexConstrainedParam().invalidate();
    }
    refitInw.setStatusMS(NA6PTrack::kMSMatchedToMIDRefitted);
    mTracks.push_back(refitInw);
  }
}
