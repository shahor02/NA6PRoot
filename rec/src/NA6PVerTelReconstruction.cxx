// NA6PCCopyright

#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <fairlogger/Logger.h>
#include "NA6PVerTelHit.h"
#include "NA6PTrackerCA.h"
#include "NA6PVerTelReconstruction.h"

ClassImp(NA6PVerTelReconstruction)

bool NA6PVerTelReconstruction::init(const char* filename, const char* geoname)
{
  NA6PReconstruction::init(filename, geoname);
  mVTTracker = new NA6PTrackerCA();
  mVTTracker->configureFromRecoParam();
  createTracksOutput();
  return true;
}

void NA6PVerTelReconstruction::createClustersOutput()
{
  auto nm = fmt::format("Clusters{}.root", getName());
  mClusFile = TFile::Open(nm.c_str(), "recreate");
  mClusTree = new TTree(fmt::format("clusters{}", getName()).c_str(), fmt::format("{} Clusters", getName()).c_str());
  mClusTree->Branch(getName().c_str(), &hClusPtr);
  LOGP(info, "Will store {} clusters in {}", getName(), nm);
}

void NA6PVerTelReconstruction::writeClusters()
{
  if (mClusTree) {
    mClusTree->Fill();
  }
  LOGP(info, "Saved {} clusters in tree with {} entries", mClusters.size(),mClusTree->GetEntries());
}

void NA6PVerTelReconstruction::closeClustersOutput()
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

void NA6PVerTelReconstruction::hitsToRecPoints(const std::vector<NA6PVerTelHit>& vtHits)
{
  int nHits = vtHits.size();
  for(int jHit = 0; jHit < nHits; ++jHit) {
    const auto& hit = vtHits[jHit];
    double x = hit.getX();
    double y = hit.getY();
    double z = hit.getZ();
    double ex2clu = 5.e-4;
    double ey2clu = 5.e-4;
    if (mCluRes > 0) {
      x = gRandom->Gaus(hit.getX(),mCluRes);
      y = gRandom->Gaus(hit.getY(),mCluRes);
      ex2clu = mCluRes*mCluRes;
      ey2clu = mCluRes*mCluRes;
    }
    double eloss = hit.getHitValue();
    // very rough cluster size settings
    int clusiz = 2;
    if (eloss > 2.e-5 && eloss < 5.e-5) clusiz = 3;
    else if (eloss > 5.e-5) clusiz = 4;
    int nDet=hit.getDetectorID();
    int idPart=hit.getTrackID();    
    mClusters.emplace_back(x, y, z, clusiz);
    auto& clu = mClusters.back();
    clu.setErr(ex2clu,0.,ey2clu);
    clu.setDetectorID(nDet);
    clu.setParticleID(idPart);
    clu.setHitID(jHit);    
  }
}

void NA6PVerTelReconstruction::createTracksOutput()
{
  auto nm = fmt::format("Tracks{}.root", getName());
  mTrackFile = TFile::Open(nm.c_str(), "recreate");
  mTrackTree = new TTree(fmt::format("tracks{}", getName()).c_str(), fmt::format("{} Tracks", getName()).c_str());
  mTrackTree->Branch(getName().c_str(), &hTrackPtr);
  LOGP(info, "Will store {} tracks in {}", getName(), nm);
}

void NA6PVerTelReconstruction::writeTracks()
{
  if (mTrackTree) {
    mTrackTree->Fill();
  }
  LOGP(info, "Saved {} tracks in tree with {} entries", mTracks.size(),mTrackTree->GetEntries());
}

void NA6PVerTelReconstruction::closeTracksOutput()
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

void NA6PVerTelReconstruction::runTracking()
{
  if (!mIsInitialized){
    LOGP(error,"Magnetic field and geometry not initialized");
    return;
  }
  clearTracks();
  mVTTracker->findTracks(mClusters,mPrimaryVertex);
  mTracks = mVTTracker->getTracks();
  writeTracks();
}
