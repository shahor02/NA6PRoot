// NA6PCCopyright

// Minimal matching reconstruction driver

#ifndef NA6P_MATCHING_H
#define NA6P_MATCHING_H

#include <string>
#include <Rtypes.h>
#include <TVector3.h>

#include "NA6PReconstruction.h"
#include "NA6PVerTelCluster.h"
#include "NA6PMuonSpecCluster.h"
#include "NA6PFastTrackFitter.h"

#include "NA6PTrack.h"
#include "NA6PVertex.h"
class TH1D;
class TH2D;

class TFile;
class TTree;
class NA6PFastTrackFitter;

class NA6PMatching : public NA6PReconstruction
{
 public:
  NA6PMatching();
  ~NA6PMatching() override = default;

  void useMCMatching() { mMCMatching = true; }
  void useChi2Matching() { mMCMatching = false; }

  void setMaxChi2Match(double v) { mMaxChi2Match = v; }
  void setMaxChi2Refit(double v) { mMaxChi2Refit = v; }
  void setMinTrackP(double p) { mMinTrackP = p; }
  void setMinVTHits(int n) { mMinVTHits = n; }
  void setMinMSHits(int n) { mMinMSHits = n; }
  void setPMatchWindow(double w) { mPMatchWindow = w; }
  void setZMatching(double z)
  {
    mZMatching = z;
    mIsZMatchingSet = true;
  }

  void setPrimaryVertexPosition(double x, double y, double z)
  {
    mPrimaryVertex.SetXYZ(x, y, z);
  }

  void setVerTelTracks(const std::vector<NA6PTrack>& tracks)
  {
    mVerTelTracks = tracks;
  }

  void setMuonSpecTracks(const std::vector<NA6PTrack>& tracks)
  {
    mMuonSpecTracks = tracks;
  }

  void setVerTelClusters(const std::vector<NA6PVerTelCluster>& clusters)
  {
    mVerTelClusters = clusters;
  }

  void setMuonSpecClusters(const std::vector<NA6PMuonSpecCluster>& clusters)
  {
    mMuonSpecClusters = clusters;
  }

  bool init(const char* filename, const char* geoname = "NA6P") override;

  void createTracksOutput() override;
  void clearTracks() override
  {
    mMatchedTracks.clear();
  }
  void writeTracks() override;
  void closeTracksOutput() override;
  void sortVTTracksByP(std::vector<NA6PTrack>& tracks);
  double computeChi2(const double* par1, const double* cov1,
                     const double* par2, const double* cov2);
  void propToZMatching(std::vector<NA6PTrack>& tracks, double z, bool outer = false);
  void runMatching();
  
  bool fitAndStoreMatchedTrack(const NA6PTrack& vtTrk, const NA6PTrack& msTrk, int particleId, double matchChi2);

 private:
  void addClustersToFitter(const NA6PTrack& trk, const auto* clusPtr);
  void runMCMatching();
  void runDataMatching();
  std::unordered_map<int, int> buildMCMatchingIndex();

  std::tuple<int, bool, double> findBestChi2Match(
      const NA6PTrack& msTrack,
      const std::vector<size_t>& validVerTelIndices);

  std::vector<size_t> prefilterVerTelTracks();

  std::vector<NA6PVerTelCluster> mVerTelClusters, *hVerTelClusPtr = &mVerTelClusters;         // vector of clusters
  std::vector<NA6PMuonSpecCluster> mMuonSpecClusters, *hMuonSpecClusPtr = &mMuonSpecClusters; // vector of clusters                                                     // cluster resolution, cm (for fast simu)
  TVector3 mPrimaryVertex{0.0, 0.0, 0.0};                                                     // primary vertex position
  double mZMatching = 38.1175;
  bool mIsZMatchingSet = false;
  bool mMCMatching = true;
  int mMinVTHits = 5;
  int mMinMSHits = 6;
  double mMinTrackP = 2.0;
  double mMaxChi2Match = std::numeric_limits<double>::max();
  double mMaxChi2Refit = std::numeric_limits<double>::max();
  double mPMatchWindow = 3; // GeV/c

  NA6PFastTrackFitter* mTrackFitter = nullptr;
  std::vector<NA6PTrack> mVerTelTracks, *hVerTelTrackPtr = &mVerTelTracks;       // Vertex telescope tracks
  std::vector<NA6PTrack> mMuonSpecTracks, *hMuonSpecTrackPtr = &mMuonSpecTracks; // Muon spectrometer tracks

  std::vector<NA6PTrack> mMatchedTracks, *hMatchedTrackPtr = &mMatchedTracks; // Matched tracks

  TFile* mMatchedTrackFile = nullptr; // file with Matched tracks
  TTree* mMatchedTrackTree = nullptr; // tree of Matched tracks

  ClassDefNV(NA6PMatching, 1);
};

#endif