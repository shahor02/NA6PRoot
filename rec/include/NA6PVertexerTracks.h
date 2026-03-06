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

class NA6PVertexerTracks
{

 public:
  void setBeamX(float x) { mBeamX = x; }
  void setBeamY(float y) { mBeamY = y; }
  void createTracksPool(const std::vector<NA6PTrack>& tracks);

 private:
  std::vector<TrackVF> mTracksPool;  ///< tracks in internal representation used for vertexing
  float mBeamX = 0.;                 // beam transverse coordindates
  float mBeamY = 0.;                 // beam transverse coordindates
  float mMaxDCA = 0.1;               // cut on DCA of track to beam line (cm)
  
  ClassDefNV(NA6PVertexerTracks, 1);
};

#endif
