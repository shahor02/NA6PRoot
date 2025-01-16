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

#ifndef NA6P_BASE_HIT_H
#define NA6P_BASE_HIT_H

#include <string>
#include <TVector3.h>
#include <Rtypes.h>

// Mother class of all hit classes

class NA6PBaseHit
{
 public:
  enum HitStatus_t {
    kTrackEntering = 0x1,
    kTrackInside = 0x1 << 1,
    kTrackExiting = 0x1 << 2,
    kTrackOut = 0x1 << 3,
    kTrackStopped = 0x1 << 4,
    kTrackAlive = 0x1 << 5
  };

  NA6PBaseHit() = default;
  NA6PBaseHit(int trackid, short did, const TVector3& xyzIn, const TVector3& xyzOut, const TVector3& momIn, float time, float val, uint8_t statusStart, uint8_t statusEnd);

  auto getTrackID() const { return mTrackID; }
  void setTrackID(int id) { mTrackID = id; }

  // getting the cartesian coordinates
  auto getXIn() const { return mPosIn[0]; }
  auto getYIn() const { return mPosIn[1]; }
  auto getZIn() const { return mPosIn[2]; }
  auto* getPosIn() const { return mPosIn; }

  auto getXOut() const { return mPosOut[0]; }
  auto getYOut() const { return mPosOut[1]; }
  auto getZOut() const { return mPosOut[2]; }
  auto* getPosOut() const { return mPosOut; }

  auto getX() const { return 0.5 * (mPosIn[0] + mPosOut[0]); }
  auto getY() const { return 0.5 * (mPosIn[1] + mPosOut[1]); }
  auto getZ() const { return 0.5 * (mPosIn[2] + mPosOut[2]); }

  auto* getMomIn() const { return mMomIn; }
  auto getPX() const { return mMomIn[0]; }
  auto getPY() const { return mMomIn[1]; }
  auto getPZ() const { return mMomIn[2]; }

  auto getHitValue() const { return mHitValue; }
  auto getTime() const { return mTime; }
  auto getDetectorID() const { return mDetectorID; }

  auto getStatusEnd() const { return mTrackStatusEnd; }
  auto getStatusStart() const { return mTrackStatusStart; }

  Bool_t isEntering() const { return mTrackStatusEnd & kTrackEntering; }
  Bool_t isInside() const { return mTrackStatusEnd & kTrackInside; }
  Bool_t isExiting() const { return mTrackStatusEnd & kTrackExiting; }
  Bool_t isOut() const { return mTrackStatusEnd & kTrackOut; }
  Bool_t isStopped() const { return mTrackStatusEnd & kTrackStopped; }
  Bool_t isAlive() const { return mTrackStatusEnd & kTrackAlive; }

  Bool_t isEnteringStart() const { return mTrackStatusStart & kTrackEntering; }
  Bool_t isInsideStart() const { return mTrackStatusStart & kTrackInside; }
  Bool_t isExitingStart() const { return mTrackStatusStart & kTrackExiting; }
  Bool_t isOutStart() const { return mTrackStatusStart & kTrackOut; }
  Bool_t isStoppedStart() const { return mTrackStatusStart & kTrackStopped; }
  Bool_t isAliveStart() const { return mTrackStatusStart & kTrackAlive; }

  // modifiers
  void setTime(float time) { mTime = time; }
  void setHitValue(float val) { mHitValue = val; }
  void setDetectorID(short detID) { mDetectorID = detID; }

  void setPosIn(const TVector3& v) { setPosIn(v.X(), v.Y(), v.Z()); }
  void setPosIn(float x, float y, float z)
  {
    setXIn(x);
    setYIn(y);
    setZIn(z);
  }
  void setXIn(float x) { mPosIn[0] = x; }
  void setYIn(float y) { mPosIn[1] = y; }
  void setZIn(float z) { mPosIn[2] = z; }

  void setPosOut(const TVector3& v) { setPosOut(v.X(), v.Y(), v.Z()); }
  void setPosOut(float x, float y, float z)
  {
    setXOut(x);
    setYOut(y);
    setZOut(z);
  }
  void setXOut(float x) { mPosOut[0] = x; }
  void setYOut(float y) { mPosOut[1] = y; }
  void setZOut(float z) { mPosOut[2] = z; }

  void setMomIn(const TVector3& v)
  {
    mMomIn[0] = v.X();
    mMomIn[1] = v.Y();
    mMomIn[2] = v.Z();
  }

  void print() const;
  std::string asString() const;

 protected:
  float mPosIn[3] = {};          // cartesian position of Hit at the entrance
  float mPosOut[3] = {};         // cartesian position of Hit at the exit
  float mMomIn[3] = {};          // momentum at entrance
  float mTime = 0.f;             // time of flight
  float mHitValue = 0.;          // hit value
  int mTrackID = 0;              // track_id
  short mDetectorID = 0;         // the detector/sensor id
  uint8_t mTrackStatusEnd = 0;   // MC status flag at exit
  uint8_t mTrackStatusStart = 0; // MC status at starting point

  ClassDefNV(NA6PBaseHit, 1);
};

#endif
