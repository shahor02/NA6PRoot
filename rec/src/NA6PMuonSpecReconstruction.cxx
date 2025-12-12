// NA6PCCopyright

#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <fairlogger/Logger.h>
#include "NA6PMuonSpecModularHit.h"
#include "NA6PTrackerCA.h"
#include "NA6PVertexerTracklets.h"
#include "NA6PMuonSpecReconstruction.h"

ClassImp(NA6PMuonSpecReconstruction)

  NA6PMuonSpecReconstruction::NA6PMuonSpecReconstruction() : NA6PReconstruction("MuonSpec")
{
}

bool NA6PMuonSpecReconstruction::init(const char* filename, const char* geoname)
{
  NA6PReconstruction::init(filename, geoname);
  mMSTracker = new NA6PTrackerCA();
  mMSTracker->configureFromRecoParam();
  mMSTracker->setNLayers(6); 
  createTracksOutput();
  return true;
}

void NA6PMuonSpecReconstruction::createClustersOutput()
{
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

void NA6PMuonSpecReconstruction::hitsToRecPoints(const std::vector<NA6PMuonSpecModularHit>& hits)
{
  int nHits = hits.size();
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
    double eloss = hit.getHitValue();
    // very rough cluster size settings
    int clusiz = 1;
    int nDet = hit.getDetectorID();
    int idPart = hit.getTrackID();
    mClusters.emplace_back(x, y, z, clusiz);
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
  if (!mIsInitialized) {
    LOGP(error, "Magnetic field and geometry not initialized");
    return;
  }
  clearTracks();
  mMSTracker->findTracks(mClusters, mPrimaryVertex);
  mTracks = mMSTracker->getTracks();
  writeTracks();
}
