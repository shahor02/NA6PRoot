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

#ifndef NA6P_TRACK_H
#define NA6P_TRACK_H

#include <string>
#include <Rtypes.h>
#include "NA6PTrackParCov.h"

class NA6PBaseCluster;

class NA6PTrack : public NA6PTrackParCov
{
 public:
  enum { kNDOF = 5,
         kMaxLr = 16 };

  NA6PTrack();
  NA6PTrack(const float* xyz, const float* pxyz, int sign, float errLoose = -2);
  NA6PTrack& operator=(const NA6PTrack&) = default;
  NA6PTrack(const NA6PTrack&) = default;
  ~NA6PTrack() = default;

  void reset();
  void setOuterParam(const NA6PTrackParCov& p) { mOuter = p; }
  void setVertexConstrainedParam(const NA6PTrackParCov& p) { mConstrained = p; }
  void setStatusRefitInward(bool status) { mStatusRefitInward = status; }
  void setStatusConstrained(bool status) { mStatusConstrained = status; }
  bool getStatusRefitInward() const { return mStatusRefitInward; }
  bool getStatusConstrained() const { return mStatusConstrained; }

  bool updateTrack(const NA6PBaseCluster& cl, float maxChi2);

  float getChi2() const { return mChi2VT + mChi2MS; }
  float getChi2VT() const { return mChi2VT; }
  float getChi2MS() const { return mChi2MS; }
  float getMatchChi2() const { return mMatchChi2; }

  const NA6PTrackParCov& getOuterParam() const { return mOuter; }
  NA6PTrackParCov& getOuterParam() { return mOuter; }

  const NA6PTrackParCov& getVertexConstrainedParam() const { return mConstrained; }
  NA6PTrackParCov& getVertexConstrainedParam() { return mConstrained; }

  int getNVTHits() const { return mNClustersVT; }
  int getNMSHits() const { return mNClustersMS; }
  int getNTRHits() const { return mNClustersTR; }
  int getNHits() const { return mNClusters; }

  uint32_t getClusterMap() const { return mClusterMap; }
  int getClusterIndex(int lr) const { return lr < kMaxLr ? mClusterIndices[lr] : -1; }
  int getParticleLabel(int lr) const { return lr < kMaxLr ? mClusterPartID[lr] : -2; }
  int getParticleID() const { return mParticleID; }
  int getCAIteration() const { return mCAIteration; }

  float getChi2VTOuter() const { return mChi2VTOuter; }
  float getChi2MSOuter() const { return mChi2MSOuter; }
  float getChi2VTRefit() const { return mChi2VTRefit; }
  float getChi2MSRefit() const { return mChi2MSRefit; }
  float getNormChi2() const { return mNClusters < 3 ? 0 : (mChi2VT + mChi2MS) / ((mNClusters << 1) - kNDOF); }

  void setMatchChi2(float chi2) { mMatchChi2 = chi2; }
  void setChi2VT(float chi2) { mChi2VT = chi2; }
  void setChi2MS(float chi2) { mChi2MS = chi2; }
  void setChi2VTOuter(float chi2) { mChi2VTOuter = chi2; }
  void setChi2MSOuter(float chi2) { mChi2MSOuter = chi2; }
  void setChi2VTRefit(float chi2) { mChi2VTRefit = chi2; }
  void setChi2MSRefit(float chi2) { mChi2MSRefit = chi2; }
  void setParticleLabel(int idx, int lr)
  {
    if (lr < kMaxLr)
      mClusterPartID[lr] = idx;
  }
  void setClusterIndex(int idx, int lr)
  {
    if (lr < kMaxLr)
      mClusterIndices[lr] = idx;
  }
  void setParticleID(int idx) { mParticleID = idx; }
  void setCAIteration(int iter) { mCAIteration = iter; }

  template <typename ClusterType>
  void addCluster(const ClusterType* clu, int cluIndex, float chi2);

  void print() const;
  std::string asString() const;

 protected:
  NA6PTrackParCov mOuter{};                  // parametrization for outward fit
  NA6PTrackParCov mConstrained{};            // parametrization with vertex constrain
  std::array<int, kMaxLr> mClusterIndices{}; // cluster indices
  std::array<int, kMaxLr> mClusterPartID{};  // particle ID (per cluster) // RSTOD why this is needed? This info must be available from the cluster indices

  float mMatchChi2 = 0.f;                    // total chi2
  float mChi2VT = 0.f;                       // total chi2 VT
  float mChi2MS = 0.f;                       // total chi2 MS
  float mChi2VTOuter = 0.f;                  // total chi2 VT outward fit
  float mChi2MSOuter = 0.f;                  // total chi2 MS outward fit
  float mChi2VTRefit = 0.f;                  // total chi2 VT inward refit
  float mChi2MSRefit = 0.f;                  // total chi2 MS inward refit
  uint32_t mClusterMap = 0;                  // pattern of clusters per layer
  int mNClusters = 0;                        // total hits
  int mNClustersVT = 0;                      // total VT hits
  int mNClustersMS = 0;                      // total MS hits
  int mNClustersTR = 0;                      // total TR hits
  int mParticleID = -1;                      // particle ID (MC truth)
  int mCAIteration = -1;                     //! CA iteration (for debug)
  bool mStatusRefitInward = false;           // boolean for status of inward refit
  bool mStatusConstrained = false;           // boolean for status of contrained track

 private:
  ClassDefNV(NA6PTrack, 1)
};

#endif
