// NA6PCCopyright

// Based on:
// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef NA6P_VERTEXER_TRACKLETS_H
#define NA6P_VERTEXER_TRACKLETS_H

#include <Rtypes.h>

class NA6PVertex;

// Vertex reconstruction from tracklets in the telescope

// structures for temporary objects used in track finding

struct TrackletForVertex {
  int startingLayer;
  int firstClusterIndex;
  int secondClusterIndex;
  float tanL;
  float phi;
  float pxpz;
  float pypz;
  bool isSignal;
};

struct TracklIntersection {
  float zeta;
  float sigmazeta;
  float tanl;
  int firstClusterIndex;
  int secondClusterIndex;
};

class NA6PVertexerTracklets
{
 public:
  enum { kYZ,
         kRZ,
         k3D,
         kXZ };
  enum { kNoWeight,
         kTanL,
         kSigma };
  enum { kMaxPileupVertices = 5,
         kMaxBinsForPeakFind = 1000 };
  enum { kKDE,
         kHistoPeak };
  enum { kStandardKDE,
         kAdaptiveKDE };

  NA6PVertexerTracklets();
  ~NA6PVertexerTracklets() = default;

  void configurePeakFinding(float zmin = -20.0, float zmax = 5., int nbins = 250)
  {
    mZMin = zmin;
    mZMax = zmax;
    mNBinsForPeakFind = nbins;
    mZBinWidth = (mZMax - mZMin) / mNBinsForPeakFind;
    mHistIntersec.assign(nbins, 0);
  }
  void configureKDE(float width, int nbins = 500)
  {
    mNGridKDE = nbins;
    mKDEBandwidth = width;
  }

  // Settters
  void setNLayersVT(int n) { mNLayersVT = n; }
  void setLayerToStart(int lay) { mLayerToStart = lay; }
  void setBeamX(float x) { mBeamX = x; }
  void setBeamY(float y) { mBeamY = y; }
  void setRecoInYZPlane() { mRecoType = kYZ; }
  void setRecoInXZPlane() { mRecoType = kXZ; }
  void setRecoInRZPlane() { mRecoType = kRZ; }
  void setRecoIn3D() { mRecoType = k3D; }
  void setMaxDCAxy(float v) { mMaxDCAxy = v; }
  void setMaxDeltaThetaTracklet(float v) { mMaxDeltaThetaTracklet = v; }
  void setMaxDeltaPhiTracklet(float v) { mMaxDeltaPhiTracklet = v; }
  void setMaxDeltaTanLamInOut(float v) { mMaxDeltaTanLamInOut = v; }
  void setMaxDeltaPhiInOut(float v) { mMaxDeltaPhiInOut = v; }
  void setMaxDeltaPxPzInOut(float v) { mMaxDeltaPxPzInOut = v; }
  void setMaxDeltaPyPzInOut(float v) { mMaxDeltaPyPzInOut = v; }
  void setUseHistoForPeakFinding() { mMethod = kHistoPeak; }
  void setUseKDEForPeakFinding() { mMethod = kKDE; }
  void setUseNoWeight() { mWeightedMeanOption = kNoWeight; }
  void setUseTanLWeight() { mWeightedMeanOption = kTanL; }
  void setUseSigmaWeight() { mWeightedMeanOption = kSigma; }
  void setZRange(float zmin, float zmax)
  {
    mZMin = zmin;
    mZMax = zmax;
    if (mNBinsForPeakFind > 0)
      configurePeakFinding(mZMin, mZMax, mNBinsForPeakFind);
  }
  void setNBinsForPeakFind(int nbins)
  {
    mNBinsForPeakFind = nbins;
    configurePeakFinding(mZMin, mZMax, mNBinsForPeakFind);
  }
  void setPeakWidthBins(int bins) { mPeakWidthBins = bins; }
  void setMinCountsInPeak(int counts) { mMinCountsInPeak = counts; }
  void setUseStandardKDE() { mKDEOption = kStandardKDE; }
  void setUseAdaptiveKDE() { mKDEOption = kAdaptiveKDE; }
  void setVerbosity(bool opt = true) { mVerbose = opt; }
  void configureFromRecoParam(const std::string& filename = "");
  void printConfiguration() const;

  // methods for vertex calculation
  void findVertices(std::vector<NA6PBaseCluster>& cluArr,
                    std::vector<NA6PVertex>& vertices);

  void sortClustersByLayerAndEta(std::vector<NA6PBaseCluster>& cluArr,
                                 std::vector<int>& firstIndex,
                                 std::vector<int>& lastIndex);
  void sortTrackletsByLayerAndIndex(std::vector<TrackletForVertex>& tracklets,
                                    std::vector<int>& firstIndex,
                                    std::vector<int>& lastIndex);
  void computeLayerTracklets(const std::vector<NA6PBaseCluster>& cluArr,
                             const std::vector<int>& firstIndex,
                             const std::vector<int>& lastIndex,
                             std::vector<TrackletForVertex>& tracklets);
  void selectTracklets(const std::vector<TrackletForVertex>& tracklets,
                       const std::vector<int>& firstIndex,
                       const std::vector<int>& lastIndex,
                       const std::vector<NA6PBaseCluster>& cluArr,
                       std::vector<TrackletForVertex>& selTracklets);
  void filterOutUsedTracklets(const std::vector<TrackletForVertex>& tracklets,
                              std::vector<TrackletForVertex>& usableTracklets);
  void computeIntersections(const std::vector<TrackletForVertex>& selTracklets,
                            const std::vector<NA6PBaseCluster>& cluArr,
                            std::vector<TracklIntersection>& zIntersec);
  bool findVertexHistoPeak(std::vector<TracklIntersection>& zIntersec,
                           std::vector<NA6PVertex>& vertices);

  float gaussKernel(float u)
  {
    return std::exp(-0.5 * u * u) / std::sqrt(2.0 * M_PI);
  }
  void resetClusters(int nClus)
  {
    mIsClusterUsed.resize(nClus);
    for (int jClu = 0; jClu < nClus; jClu++)
      mIsClusterUsed[jClu] = false;
  }
  bool findVertexKDE(const std::vector<TracklIntersection>& zIntersec,
                     std::vector<NA6PVertex>& vertices);
  void printStats(const std::vector<TrackletForVertex>& candidates,
                  const std::vector<NA6PBaseCluster>& cluArr,
                  const std::string& label);
  short getMethodForPeakFinding() const { return mMethod; }

 private:
  int mNLayersVT = 5;                    // number of layers in the VT
  int mLayerToStart = 0;                 // innermost layer used in tracklets
  std::vector<bool> mIsClusterUsed = {}; // flag for used clusters
  float mBeamX = 0.;                     // beam transverse coordindates
  float mBeamY = 0.;                     // beam transverse coordindates
  float mMaxDeltaThetaTracklet = 0.6;    // selections for tracklet building
  float mMaxDeltaPhiTracklet = 0.05;     // selections for tracklet building
  float mMaxDeltaTanLamInOut = 1.;       // selections for tracklet validation
  float mMaxDeltaPhiInOut = 0.2;         // selections for tracklet validation
  float mMaxDeltaPxPzInOut = 99.;        // selections for tracklet validation
  float mMaxDeltaPyPzInOut = 0.005;      // selections for tracklet validation
  short mRecoType = kYZ;                 // method to compute tracklet intersections with beam axis
  float mMaxDCAxy = 0.25;                // selection for tracklet intersection, cm
  short mMethod = kKDE;                  // method for peak finding (KDE vs histo)
  int mWeightedMeanOption = kNoWeight;   // option for weigthed mean
  float mZMin = -20.0;                   // z range, min, cm
  float mZMax = 5.;                      // z range, max, cm
  float mZWindowWidth = 1.25;            // window around peak, cm
  int mNBinsForPeakFind = 250;           // 0.1 cm per bin
  float mZBinWidth = 0.1;                // bin width (recalculated from n bins)
  std::vector<int> mHistIntersec;        // histogram for the peak finding method
  int mPeakWidthBins = 3;                // number of bins for integrating histo peak
  int mMinCountsInPeak = 3;              // minimum entries in histo peak
  int mKDEOption = kStandardKDE;         // option for using uncertainties in KDE sigma
  int mNGridKDE = 500;                   // number of points to sample the KDE
  float mKDEBandwidth = 0.5;             // smoothing parameter, cm
  bool mVerbose = false;

  ClassDefNV(NA6PVertexerTracklets, 1);
};

#endif
