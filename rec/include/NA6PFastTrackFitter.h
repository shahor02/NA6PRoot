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

// Fast track fit based on Kalman filter

class TGeoManager;
class NA6PTrack;

class NA6PFastTrackFitter
{

 public:
  enum { kTwoPointSeed = 0,
         kThreePointSeed = 1 };

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
  void setParticleHypothesis(int pdg);
  void enableMaterialCorrections() { mCorrectForMaterial = true; }
  void disableMaterialCorrections() { mCorrectForMaterial = false; }
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
  void setSeedFromTwoOutermostHits() { mSeedOption = kTwoPointSeed; }
  void setSeedFromThreeOutermostHits() { mSeedOption = kThreePointSeed; }
  void computeSeed();
  void printSeed() const;

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

  bool loadGeometry(const char* filename = "geometry.root", const char* geoname = "NA6P");

  NA6PTrack* fitTrackPoints();
  bool updateTrack(NA6PTrack* trc, const NA6PBaseCluster* cl) const;
  int propagateToZ(NA6PTrack* trc, float zFrom, float zTo, int dir) const;
  int propagateToZ(NA6PTrack* trc, float zTo) const;
  void getMeanMaterialBudgetFromGeom(float* start, float* end, float* mparam) const;

  static const float kMassP;
  static const float kMassK;
  static const float kMassPi;
  static const float kMassMu;
  static const float kMassE;
  static const float kAlmostZero;

 protected:
  int mNLayers = 5;                  // number of active
  float mMaxChi2Cl = 10.f;           // max cluster-track chi2
  bool mIsSeedSet = false;           // flag for set seed
  int mSeedOption = kThreePointSeed; // seed option (see enum)
  float mSeedPos[3];                // seed for track position
  float mSeedMom[3];                // seed for track momentum
  int mCharge = 1;                   // track charge for seed
  float mMass = kMassPi;            // mass hypothesis for particle
  bool mPropagateToPrimVert = false; // flag for propagation to primary vertex
  float mPrimVertZ = 0.0f;           // primary vertex z
  bool mIsPrimVertSet = false;       // flag for presence of prim vert z
  bool mCorrectForMaterial = true;   // flag for material corrections

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
