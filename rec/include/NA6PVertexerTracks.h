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

#ifndef NA6P_VERTEXER_TRACKS_H
#define NA6P_VERTEXER_TRACKS_H

#include <vector>
#include <NA6PLine.h>
#include <NA6PTrack.h>
#include <NA6PVertex.h>

struct TrackVF {
  
  NA6PLine mLine;           // straight line representation of the track
  float mSig2XI = 0.f;      // XX component of inverse cov.matrix
  float mSig2YI = 0.f;      // YY component of inverse cov.matrix
  float mSigXYI = 0.f;      // XY component of inverse cov.matrix
  float mWeightHisto = 0.f; // weight for histogram seeder 

  TrackVF() = default;
  TrackVF(const NA6PTrack& src) {
    // NB: src must already be propagated to its DCA to the beam axis
    // before constructing TrackVF — call propagateToDCABeamAxis() first
    double xyz[3], pxyz[3];
    src.getXYZ(xyz);
    src.getPXYZ(pxyz);
    mLine = NA6PLine::fromPointAndDirection(xyz, pxyz);
    float sxx = src.getSigmaX2(), syy = src.getSigmaY2(), sxy = src.getSigmaXY();
    float cx = mLine.mCosinesDirector[0];
    float cy = mLine.mCosinesDirector[1];
    float cz = mLine.mCosinesDirector[2];
    float iczcz = 1.f / (cz * cz);
    float szz = (cx * cx * sxx + cy * cy * syy + 2.f * cx * cy * sxy) * iczcz;
    float det = sxx * syy - sxy * sxy;
    if (det <= 1e-20f) {
      mWeightHisto = -1;
      return;
    }
    auto detI = 1. / det;
    mSig2XI = sxx * detI;
    mSig2YI = syy * detI;
    mSigXYI = -sxy * detI;
    mWeightHisto = (szz > 0.f) ? 1.f / szz : 0.f;
  }
  bool isValid() const { return mWeightHisto >= 0.f; }
  std::array<float,3> residual(const float vtxPos[3]) const
  {
    // vector from closest point on line to vtxPos
    auto comps = mLine.getDCAComponents(vtxPos);
    return {comps[0], comps[3], comps[5]}; // dx, dy, dz components
  }
  ClassDefNV(TrackVF, 1);
};

struct VertexSeed {
  float x = 0.f, y = 0.f, z = 0.f;
  double wghSum = 0.;    // sum of tracks weights
  double wghChi2 = 0.;   // sum of tracks weighted chi2's
  double cxx = 0., cyy = 0., czz = 0., cxy = 0., cxz = 0., cyz = 0., cx0 = 0., cy0 = 0., cz0 = 0.; // elements of lin.equation matrix
  float scaleSigma2 = 1.;  // scaling parameter on top of Tukey param
  float scaleSigma2Prev = 1.;
  float maxScaleSigma2Tested = 0.;
  float scaleSig2ITuk2I = 0; // inverse squared Tukey parameter scaled by scaleSigma2
  int nScaleSlowConvergence = 0;
  int nScaleIncrease = 0;
  int nIterations = 0;

  void setScale(float scale2, float tukey2I)
  {
    scaleSigma2Prev = scaleSigma2;
    scaleSigma2 = scale2;
    scaleSig2ITuk2I = tukey2I / scale2;
  }

  void resetForNewIteration()
  {
    wghSum = 0.;
    wghChi2 = 0.;
    cxx = cyy = czz = cxy = cxz = cyz = cx0 = cy0 = cz0 = 0.;
  }

  VertexSeed() = default;
  VertexSeed(const NA6PVertex& src)
  {
    x = src.getX();
    y = src.getY();
    z = src.getZ();
  }

};

class NA6PVertexerTracks
{

 public:
  
  enum class FitStatus : int { Failure,
                               PoolEmpty,
                               NotEnoughTracks,
                               IterateFurther,
                               OK };
  
  NA6PVertexerTracks();
  ~NA6PVertexerTracks() = default;

  void setBeamX(float x) { mBeamX = x; }
  void setBeamY(float y) { mBeamY = y; }
  void configurePeakFinding(float zmin = -20.0, float zmax = 5., int nbins = 250)
  {
    mZMin = zmin;
    mZMax = zmax;
    mNBinsForPeakFind = nbins;
    mZBinWidth = (mZMax - mZMin) / mNBinsForPeakFind;
    mHistZ.assign(nbins, 0.f);
  }
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
  void createTracksPool(const std::vector<NA6PTrack>& tracks);
  void buildAndFillHistoZ();
  int findPeakBin();
  int findVertices();

  
 private:
  std::vector<TrackVF> mTracksPool;  ///< tracks in internal representation used for vertexing
  float mBeamX = 0.;                 // beam transverse coordindates
  float mBeamY = 0.;                 // beam transverse coordindates
  float mMaxDCA = 0.1;               // cut on DCA of track to beam line (cm)
  float mZMin = -20.0;               // z range, min, cm
  float mZMax = 5.;                  // z range, max, cm
  int mNBinsForPeakFind = 250;       // 0.1 cm per bin
  float mZBinWidth = 0.1;            // bin width 
  std::vector<float> mHistZ;         // histogram for the peak finding method
  std::vector<int>   mFilledBinsZ;   // indices of non-empty bins
  int mMaxVerticesPerCluster = 5;    // one vertex per target
  int mMaxTrialsPerCluster = 2;      //
  
  ClassDefNV(NA6PVertexerTracks, 1);
};

#endif
