// NA6PCCopyright

#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <fairlogger/Logger.h>
#include "NA6PVerTelHit.h"
#include "NA6PTrackerCA.h"
#include "NA6PVertexerTracklets.h"
#include "NA6PVerTelReconstruction.h"

ClassImp(NA6PVerTelReconstruction)

  NA6PVerTelReconstruction::NA6PVerTelReconstruction() : NA6PReconstruction("VerTel"),
                                                         mGeoFilName{"geometry.root"},
                                                         mGeoObjName{"NA6P"},
                                                         mRecoParFilName{""}
{
}

NA6PVerTelReconstruction::NA6PVerTelReconstruction(const char* recparfile,
                                                   const char* geofile,
                                                   const char* geoname) : NA6PReconstruction("VerTel"),
                                                                          mGeoFilName{geofile},
                                                                          mGeoObjName{geoname},
                                                                          mRecoParFilName{recparfile}
{
  initAll();
}

bool NA6PVerTelReconstruction::initAll()
{
  bool retV = initVertexer();
  bool retC = initTracker();
  return (retV & retC);
}

bool NA6PVerTelReconstruction::initVertexer()
{
  if (!mVTTrackletVertexer)
    mVTTrackletVertexer = new NA6PVertexerTracklets();
  if (mRecoParFilName == "")
    LOGP(info, "Initializing vertexer with default parameters");
  else
    LOGP(info, "Initializing vertexer from file {}", mRecoParFilName.c_str());
  mVTTrackletVertexer->configureFromRecoParam(mRecoParFilName.c_str());
  createVerticesOutput();
  return true;
}

bool NA6PVerTelReconstruction::initTracker()
{
  NA6PReconstruction::init(mGeoFilName.c_str(), mGeoObjName.c_str());
  if (!mVTTracker)
    mVTTracker = new NA6PTrackerCA();
  mVTTracker->setNLayers(5);
  mVTTracker->setStartLayer(0);
  if (mRecoParFilName == "")
    LOGP(info, "Initializing tracker with default parameters");
  else
    LOGP(info, "Initializing tracker from file {}", mRecoParFilName.c_str());
  mVTTracker->configureFromRecoParamVT(mRecoParFilName.c_str());
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
  LOGP(info, "Saved {} clusters in tree with {} entries", mClusters.size(), mClusTree->GetEntries());
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

void NA6PVerTelReconstruction::setClusters(const std::vector<NA6PVerTelCluster>& clusters)
{
  mClusters = clusters;
  // Assign a transient index based on position in the vector
  // to be used for storing the cluster indices in the track
  for (size_t i = 0; i < mClusters.size(); ++i) {
    mClusters[i].setClusterIndex(static_cast<int>(i));
  }
}

void NA6PVerTelReconstruction::hitsToRecPoints(const std::vector<NA6PVerTelHit>& hits)
{
  int nHits = hits.size();
  for (int jHit = 0; jHit < nHits; ++jHit) {
    const auto& hit = hits[jHit];
    double x = hit.getX();
    double y = hit.getY();
    double z = hit.getZ();
    double ex2clu = 5.e-4;
    double ey2clu = 5.e-4;
    if (mCluRes > 0) {
      x = gRandom->Gaus(hit.getX(), mCluRes);
      y = gRandom->Gaus(hit.getY(), mCluRes);
      ex2clu = mCluRes * mCluRes;
      ey2clu = mCluRes * mCluRes;
    }
    double eloss = hit.getHitValue();
    // very rough cluster size settings
    int clusiz = 2;
    if (eloss > 2.e-5 && eloss < 5.e-5)
      clusiz = 3;
    else if (eloss > 5.e-5)
      clusiz = 4;
    int nDet = hit.getDetectorID();
    int idPart = hit.getTrackID();
    int layer = nDet / 4;
    mClusters.emplace_back(x, y, z, clusiz, layer);
    auto& clu = mClusters.back();
    clu.setErr(ex2clu, 0., ey2clu);
    clu.setDetectorID(nDet);
    clu.setParticleID(idPart);
    clu.setHitID(jHit);
  }
}

//____________________________________________________________________________________

void NA6PVerTelReconstruction::createVerticesOutput()
{
  auto nm = fmt::format("Vertices{}.root", getName());
  mVertexFile = TFile::Open(nm.c_str(), "recreate");
  mVertexTree = new TTree(fmt::format("vertices{}", getName()).c_str(), fmt::format("{} Vertices", getName()).c_str());
  mVertexTree->Branch(getName().c_str(), &hVerticesPtr);
  LOGP(info, "Will store {} vertices in {}", getName(), nm);
}

void NA6PVerTelReconstruction::writeVertices()
{
  if (mVertexTree) {
    mVertexTree->Fill();
  }
  LOGP(info, "Saved {} vertices in tree with {} entries", mVertices.size(), mVertexTree->GetEntries());
}

void NA6PVerTelReconstruction::closeVerticesOutput()
{
  if (mVertexTree && mVertexFile) {
    mVertexFile->cd();
    mVertexTree->Write();
    delete mVertexTree;
    mVertexTree = nullptr;
    mVertexFile->Close();
    delete mVertexFile;
    mVertexFile = nullptr;
  }
}

void NA6PVerTelReconstruction::runVertexerTracklets()
{
  if (!mVTTrackletVertexer) {
    LOGP(info, "Tracklet vertexer not initialized, will call initialization");
    initVertexer();
  }
  clearVertices();
  mVTTrackletVertexer->findVertices(mClusters, mVertices);
  writeVertices();
}
//____________________________________________________________________________________

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
  LOGP(info, "Saved {} tracks in tree with {} entries", mTracks.size(), mTrackTree->GetEntries());
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
  if (!mIsInitialized) {
    LOGP(error, "Magnetic field and geometry not initialized");
    return;
  }
  if (!mVTTracker) {
    LOGP(info, "Tracker not initialized, will call default initialization");
    initTracker();
  }
  clearTracks();
  mVTTracker->setDoOutwardPropagation(true);
  mVTTracker->setPropagateTracksToPrimaryVertex(true);
  mVTTracker->findTracks(mClusters, mPrimaryVertex);
  mTracks = mVTTracker->getTracks();
  writeTracks();
}
