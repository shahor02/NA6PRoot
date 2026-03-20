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

#ifndef NA6P_FAST_TRACK_FITTER_H
#define NA6P_FAST_TRACK_FITTER_H

#include <string>
#include <Rtypes.h>
#include "NA6PBaseCluster.h"
#include "NA6PTrack.h"
#include "Propagator.h"

// Fast track fit based on Kalman filter

class TGeoManager;
class NA6PTrack;
class NA6PVertex;

class NA6PFastTrackFitter
{

 public:
  enum { kTwoPointSeed = 0,
         kThreePointSeed = 1 };
  enum { kOutermostAsSeed = 0,
         kInnermostAsSeed = 1,
         kInMidOutAsSeed = 2 };
  enum { kBatMidPoint = 0,
         kMaximumB = 1,
         kIntegralB = 2 };

  static constexpr int kMaxLayers = 20;

  NA6PFastTrackFitter();
  ~NA6PFastTrackFitter(){};

  // setters for configurable parameters
  void setMaxChi2Cl(float v = 10) { mMaxChi2Cl = v; }
  void setNLayers(int n)
  {
    mNLayers = n;
    mClusters.clear();
    mClusters.resize(n);
  }
  void setParticleHypothesis(int pdg) { mPID.setFromPDG(pdg); }

  void enableMaterialCorrections() { mPropOpt.matCorr = Propagator::MatCorrType::USEMatCorrTGeo; }
  void disableMaterialCorrections() { mPropOpt.matCorr = Propagator::MatCorrType::USEMatCorrNONE; }
  void setMaxPropagationStep(float v) { mPropOpt.maxStep = std::max(0.1f, std::abs(v)); }
  void setUseByPropagation(bool v = true) { mPropOpt.byOnly = v; }

  void setPropagateToPrimaryVertex(bool opt = true) { mPropagateToPrimVert = opt; }
  void setPrimaryVertexZ(float zvert)
  {
    mPrimVertZ = zvert;
    mIsPrimVertSet = true;
  }
  void setSeed(const float* pos, const float* mom, int charge = 1);
  void setCharge(int ch)
  {
    if (ch != 0)
      mCharge = ch;
  };
  void unsetSeed() { mIsSeedSet = false; }
  void setSeedFromTwoHits() { mSeedPoints = kTwoPointSeed; }
  void setSeedFromThreeHits() { mSeedPoints = kThreePointSeed; }
  void setSeedFromOutermostHits() { mSeedOption = kOutermostAsSeed; }
  void setSeedFromInnermostHits() { mSeedOption = kInnermostAsSeed; }
  void setSeedFromInMidOutHits() { mSeedOption = kInMidOutAsSeed; }
  void setUseBatMidPointForSeed() { mOptionForSeedB = kBatMidPoint; }
  void setUseMaximumBForSeed() { mOptionForSeedB = kMaximumB; }
  void setUseIntegralBForSeed() { mOptionForSeedB = kIntegralB; }
  int getLayersForSeed(std::array<int, 3>& layForSeed) const;
  int sortLayersForSeed(std::array<int, 3>& layForSeed, int dir) const;
  void computeSeed(int dir, std::array<int, 3>& layForSeed);
  void computeSeed(int dir = -1);
  void computeSeedOuter() { computeSeed(-1); }
  void computeSeedInner() { computeSeed(1); }
  void printClusters() const;
  void printSeed() const;
  const auto& getSeedMomentum() const { return mSeedMom; }
  const auto& getSeedPosition() const { return mSeedPos; }
  int getCharge() const { return mCharge; }
  void addCluster(int jLay, const NA6PBaseCluster& cl);
  void resetClusters()
  {
    std::fill(mClusters.begin(), mClusters.end(), nullptr);
  }
  void cleanupAndStartFit()
  {
    resetClusters();
    unsetSeed();
  }

  bool loadGeometry(const std::string& filename = "geometry.root", const std::string geoname = "NA6P") { return Propagator::loadGeometry(filename, geoname); }

  bool fitTrackPoints(NA6PTrack& trackToFit, int dir = -1, const NA6PTrackParCov* seed = nullptr);
  bool fitTrackPointsInward(NA6PTrack& trackToFit, const NA6PTrackParCov* seed = nullptr) { return fitTrackPoints(trackToFit, -1, seed); }
  bool fitTrackPointsOutward(NA6PTrack& trackToFit, const NA6PTrackParCov* seed = nullptr) { return fitTrackPoints(trackToFit, 1, seed); }

  bool constrainTrackToVertex(NA6PTrack& trc, const NA6PVertex& pv) const;
  const auto& getPropOpt() const { return mPropOpt; }

 protected:
  int mNLayers = 5;                   // number of active
  float mMaxChi2Cl = 10.;             // max cluster-track chi2
  bool mIsSeedSet = false;            // flag for set seed
  int mSeedOption = kOutermostAsSeed; // seed option (see enum)
  int mSeedPoints = kThreePointSeed;  // number of hits used for seed
  int mOptionForSeedB = kBatMidPoint; // option for B field usage in seed

  std::array<float, 3> mSeedPos{};                 // seed for track position
  std::array<float, 3> mSeedMom = {0.f, 0.f, 1.f}; // seed for track momentum
  int mCharge = 1;                    // track charge for seed

  bool mPropagateToPrimVert = false;  // flag for propagation to primary vertex
  bool mIsPrimVertSet = false;        // flag for presence of prim vert z
  float mPrimVertZ = 0.0;             // primary vertex z

  Propagator::PropOpt mPropOpt{}; // propagator options

  PID mPID;                                      // PID to impose
  std::vector<const NA6PBaseCluster*> mClusters; // array with clusters

 protected:
  int getNumberOfClusters() const
  {
    int nClus = 0;
    for (int jLay = 0; jLay < mNLayers; ++jLay)
      if (mClusters[jLay])
        nClus++;
    return nClus;
  }

  ClassDefNV(NA6PFastTrackFitter, 1);
};

#endif
