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
  enum { kUsed,
         kNoVtx = -1,
         kDiscarded = kNoVtx - 1 };

  NA6PLine mLine;      // straight line representation of the track
  float mSig2XI = 0.f; // XX component of inverse cov.matrix
  float mSig2YI = 0.f; // YY component of inverse cov.matrix
  float mSig2ZI = 0.f; // ZZ component of inverse cov matrix
  float mSigXYI = 0.f; // XY component of inverse cov matrix
  float mSigXZI = 0.f; // XZ component of inverse cov matrix
  float mSigYZI = 0.f; // YZ component of inverse cov matrix
  int trackIndex = -1; // track index
  float wgh = 0.f;     ///< track weight wrt current vertex seed
  int vtxID = kNoVtx;  // assigned vertex

  TrackVF() = default;
  TrackVF(const NA6PTrack& src, int id)
  {
    // NB: src must already be propagated to its DCA to the beam axis
    // before constructing TrackVF — call propagateToDCABeamAxis() first
    trackIndex = id;
    double xyz[3], pxyz[3];
    src.getXYZ(xyz);
    src.getPXYZ(pxyz);
    mLine = NA6PLine::fromPointAndDirection(xyz, pxyz);
    float sxx = src.getSigmaX2(), syy = src.getSigmaY2(), sxy = src.getSigmaXY();
    float cx = mLine.mCosinesDirector[0];
    float cy = mLine.mCosinesDirector[1];
    float cz = mLine.mCosinesDirector[2];
    float pt2 = cx * cx + cy * cy;
    float szz = (pt2 > 1e-6f) ? 0.5f * (sxx + syy) * (cz * cz) / pt2 : 1e10f;
    float det = sxx * syy - sxy * sxy;
    if (det <= 1e-20f) {
      mSig2ZI = -1.f;
      return;
    }
    float detI = 1.f / det;
    mSig2XI = syy * detI;
    mSig2YI = sxx * detI;
    mSigXYI = -sxy * detI;
    // z weight from error propagation of transverse cov along track direction
    // xz and yz correlations neglected (mSigXZI = mSigYZI = 0)
    mSig2ZI = (szz > 0.f) ? 1.f / szz : 0.f;
  }
  bool isValid() const { return mSig2ZI > 0.f; }
  bool canUse() const { return vtxID == kNoVtx; }
  bool canAssign() const { return wgh > 0. && vtxID == kNoVtx; }

  std::array<float, 3> getResiduals(const std::array<float, 3>& vtxPos) const
  {
    // vector from closest point on line to vtxPos
    auto comps = mLine.getDCAComponents(vtxPos);
    return {comps[0], comps[3], comps[5]}; // dx, dy, dz components
  }

  float evalChi2ToVertex(const std::array<float, 3>& vtxPos) const
  {
    auto res = getResiduals(vtxPos); // track-vertex residuals and chi2
    float dx = res[0], dy = res[1];
    return evalChi2ToVertex(dx, dy);
  }
  float evalChi2ToVertex(float dx, float dy) const
  {
    constexpr float NDOF2I = 0.5f;
    float chi2T = dx * dx * mSig2XI + 2.f * dx * dy * mSigXYI + dy * dy * mSig2YI;
    chi2T *= NDOF2I;
    return chi2T;
  }

  ClassDefNV(TrackVF, 1);
};

struct VertexSeed {
  float x = 0.f, y = 0.f, z = 0.f;
  double wghSum = 0.;                                                                              // sum of tracks weights
  double wghChi2 = 0.;                                                                             // sum of tracks weighted chi2's
  double cxx = 0., cyy = 0., czz = 0., cxy = 0., cxz = 0., cyz = 0., cx0 = 0., cy0 = 0., cz0 = 0.; // elements of lin.equation matrix
  float scaleSigma2 = 1.;                                                                          // scaling parameter on top of Tukey param
  float scaleSigma2Prev = 1.;
  float maxScaleSigma2Tested = 0.;
  float scaleSig2ITuk2I = 0; // inverse squared Tukey parameter scaled by scaleSigma2
  int nScaleSlowConvergence = 0;
  int nScaleIncrease = 0;
  int nIterations = 0;
  float chi2 = 0.f;      ///< final vertex chi2
  int nContributors = 0; ///< number of tracks contributing to current iteration

  void setChi2(float c) { chi2 = c; }
  float getChi2() const { return chi2; }
  int getNContributors() const { return nContributors; }
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
    nContributors = 0;
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
  void setTukey(float t) { mTukey2I = t > 0.f ? 1.f / (t * t) : 1.f / (kDefTukey * kDefTukey); }
  void setMaxVerticesPerCluster(int n) { mMaxVerticesPerCluster = n; }
  void setMaxTrialsPerCluster(int n) { mMaxTrialsPerCluster = n; }
  void setInitScaleSigma2(float s) { mInitScaleSigma2 = s; }
  void setMaxIterations(int n) { mMaxIterations = n; }
  void setMinTracksPerVtx(int n) { mMinTracksPerVtx = n; }
  void setMinScale2(float s) { mMinScale2 = s; }
  void setMaxScale2(float s) { mMaxScale2 = s; }
  void setMaxChi2Mean(float c) { mMaxChi2Mean = c; }
  void setMaxNScaleIncreased(int n) { mMaxNScaleIncreased = n; }
  void setSlowConvergenceFactor(float f) { mSlowConvergenceFactor = f; }
  void setMaxNScaleSlowConvergence(int n) { mMaxNScaleSlowConvergence = n; }
  void setAcceptableScale2(float s) { mAcceptableScale2 = s; }
  void setUpscaleFactor(float f) { mUpscaleFactor = f; }

  void setVerbosity(bool opt = true) { mVerbose = opt; }

  void createTracksPool(const std::vector<NA6PTrack>& tracks);
  void buildAndFillHistoZ();
  int findPeakBin();
  int findVertices(std::vector<NA6PVertex>& vertices);
  bool fitVertex(float zSeed, NA6PVertex& vtx);
  FitStatus fitIteration(VertexSeed& vtxSeed);
  FitStatus evalIterations(VertexSeed& vtxSeed) const;
  bool upscaleSigma(VertexSeed& vtxSeed) const;
  bool solveVertex(VertexSeed& vtxSeed) const;
  void accountTrack(TrackVF& trc, VertexSeed& vtxSeed) const;
  void finalizeVertex(NA6PVertex& vtx, std::vector<NA6PVertex>& vertices);

 private:
  std::vector<TrackVF> mTracksPool;               ///< tracks in internal representation used for vertexing
  float mBeamX = 0.;                              ///< beam transverse coordindates
  float mBeamY = 0.;                              ///< beam transverse coordindates
  float mMaxDCA = 0.05;                           ///< cut on DCA of track to beam line (cm)
  float mZMin = -20.0;                            ///< z range, min, cm
  float mZMax = 5.;                               ///< z range, max, cm
  int mNBinsForPeakFind = 250;                    ///< 0.1 cm per bin
  float mZBinWidth = 0.1;                         ///< bin width
  std::vector<float> mHistZ;                      ///< histogram for the peak finding method
  std::vector<int> mFilledBinsZ;                  ///< indices of non-empty bins
  int mMaxVerticesPerCluster = 5;                 ///< one vertex per target
  int mMaxTrialsPerCluster = 100;                 //
  static constexpr float kAlmost0F = 1e-7f;       ///< tiny float
  //  static constexpr float kScaleStability = 0.1f;  ///< tolerance for scale change
  static constexpr float kDefTukey = 5.0f;        ///< def.value for tukey constant
  float mTukey2I = 1.f / (kDefTukey * kDefTukey); ///< 1./[Tukey parameter]^2
  float mInitScaleSigma2 = 10.f;                  ///< scaling parameter on top of Tukey param
  int mMaxIterations = 20;                        ///< max iterations per vertex fit
  int mMinTracksPerVtx = 2;                       ///< minimum number of tracks per vertex
  float mMinScale2 = 1.;                          ///< min scaling factor^2
  float mMaxScale2 = 50.;                         ///< max slaling factor^2
  float mMaxChi2Mean = 30.;                       ///< max mean chi2 of vertex to accept
  int mMaxNScaleIncreased = 5;                    ///< max number of scaling-non-decreasing iterations
  float mSlowConvergenceFactor = 0.5;             ///< consider convergence as slow if ratio new/old scale2 exceeds it
  int mMaxNScaleSlowConvergence = 5;              ///< max number of weak scaling decrease iterations
  float mAcceptableScale2 = 4.;                   ///< if below this factor, try to refit with minScale2
  float mUpscaleFactor = 9.;                      ///< factor for upscaling if not candidate is found
  bool mVerbose = false;                          ///< verbosity flag

  ClassDefNV(NA6PVertexerTracks, 1);
};

#endif
