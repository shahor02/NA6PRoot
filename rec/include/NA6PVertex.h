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

#ifndef NA6P_VERTEX_H
#define NA6P_VERTEX_H

#include <string>
#include <Rtypes.h>
#include <Math/Point3D.h>

// Basic vertex class

class NA6PVertex
{
 public:
  enum CovElems : int { kCovXX,
                        kCovXY,
                        kCovYY,
                        kCovXZ,
                        kCovYZ,
                        kCovZZ };
  static constexpr int kNCov = 6;

  enum vertTypes { kTrackletPrimaryVertex,
                   kTrackPrimaryVertex,
                   kBaseVertex };

  NA6PVertex() = default;
  NA6PVertex(const ROOT::Math::XYZPointF& pos, const std::array<float, kNCov>& cov, int nCont, float chi2);
  NA6PVertex(const float* xyz, int nCont);
  NA6PVertex(const NA6PVertex&) = default;
  NA6PVertex& operator=(const NA6PVertex&) = default;
  virtual ~NA6PVertex() = default;

  void init(const float* xyz, const float* cov, int nCont, float chi2);
  void setX(float x) { mPos.SetX(x); }
  void setY(float y) { mPos.SetY(y); }
  void setZ(float z) { mPos.SetZ(z); }
  void setXYZ(float x, float y, float z)
  {
    setX(x);
    setY(y);
    setZ(z);
  }
  void setPos(const ROOT::Math::XYZPointF& p) { mPos = p; }

  void setSigmaX2(float v) { mCov[kCovXX] = v; }
  void setSigmaY2(float v) { mCov[kCovYY] = v; }
  void setSigmaZ2(float v) { mCov[kCovZZ] = v; }
  void setSigmaXY(float v) { mCov[kCovXY] = v; }
  void setSigmaXZ(float v) { mCov[kCovXZ] = v; }
  void setSigmaYZ(float v) { mCov[kCovYZ] = v; }
  void setSigmaX(float val) { setSigmaX2(val * val); }
  void setSigmaY(float val) { setSigmaY2(val * val); }
  void setSigmaZ(float val) { setSigmaZ2(val * val); }

  void setCov(float sxx, float sxy, float syy, float sxz, float syz, float szz)
  {
    setSigmaX2(sxx);
    setSigmaY2(syy);
    setSigmaZ2(szz);
    setSigmaXY(sxy);
    setSigmaXZ(sxz);
    setSigmaYZ(syz);
  }
  void setCov(const std::array<float, kNCov>& cov) { mCov = cov; }
  void setNContributors(ushort v) { mNContributors = v; }
  void addContributor() { mNContributors++; }
  void setChi2(float v) { mChi2 = v; }
  void setVertexType(short t) { mVertexType = t; }

  float getX() const { return mPos.X(); }
  float getY() const { return mPos.Y(); }
  float getZ() const { return mPos.Z(); }
  float getSigmaX2() const { return mCov[kCovXX]; }
  float getSigmaY2() const { return mCov[kCovYY]; }
  float getSigmaZ2() const { return mCov[kCovZZ]; }
  float getSigmaXY() const { return mCov[kCovXY]; }
  float getSigmaXZ() const { return mCov[kCovXZ]; }
  float getSigmaYZ() const { return mCov[kCovYZ]; }
  float getSigmaX() const { return std::sqrt(getSigmaX2()); }
  float getSigmaY() const { return std::sqrt(getSigmaY2()); }
  float getSigmaZ() const { return std::sqrt(getSigmaZ2()); }
  std::array<float, kNCov> getCov() const { return mCov; }
  ROOT::Math::XYZPointF getXYZ() const { return mPos; }
  int getNContributors() const { return mNContributors; }
  float getChi2() const { return mChi2; }
  short getVertexType() const { return mVertexType; }

  virtual void print() const;
  std::string asString() const;

 protected:
  ROOT::Math::XYZPointF mPos{0., 0., 0.};
  std::array<float, kNCov> mCov{};
  float mChi2 = 0.0;
  int mNContributors = 0;
  short mVertexType = kBaseVertex;

  ClassDefNV(NA6PVertex, 1);
};

#endif
