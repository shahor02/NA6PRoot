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

#ifndef NA6P_BASE_CLUSTER_H
#define NA6P_BASE_CLUSTER_H

#include <string>
#include <Rtypes.h>

// Mother class of all cluster classes

class NA6PBaseCluster
{
 public:

  NA6PBaseCluster() = default;
  NA6PBaseCluster(float x, float y, float z, int clusiz, int layer);
  NA6PBaseCluster(const NA6PBaseCluster&) = default;
  NA6PBaseCluster& operator=(const NA6PBaseCluster&) = default;
  virtual ~NA6PBaseCluster() {}

  int getLayer() const { return mLayer; }
  int   getClusterSize() const { return mCluSiz; }
  // Lab frame coordinates
  auto getXLab()       const {return mPos[0];}
  auto getYLab()       const {return mPos[1];}
  auto getZLab()       const {return mPos[2];}
  // Tracking frame coordinates (right-handed):
  // X_tracking = Z_lab, Y_tracking = X_lab, Z_tracking = Y_lab
  auto getXTF()        const {return mPos[2];}
  auto getYTF()        const {return mPos[0];}
  auto getZTF()        const {return mPos[1];}
  auto getX()          const {return getXLab();}
  auto getY()          const {return getYLab();}
  auto getZ()          const {return getZLab();}
  auto getSigYY()      const { return mSigYY; }
  auto getSigYZ()      const { return mSigYZ; }
  auto getSigZZ()      const { return mSigZZ; }
  auto getDetectorID() const { return mDetectorID; }
  auto getParticleID() const { return mParticleID; }
  auto getHitID()      const { return mHitID; }

  void setDetectorID(int id) { mDetectorID = id; }
  void setParticleID(int id) { mParticleID = id; }
  void setHitID(int id) { mHitID = id; }
  void setPos(float x, float y, float z){
    mPos[0] = x; mPos[1] = y; mPos[2] = z;
  }
  void setX(float v)    {mPos[0] = v;}
  void setY(float v)    {mPos[1] = v;}
  void setZ(float v)    {mPos[2] = v;}
  void setErr(float syy, float syz, float szz){
    mSigYY = syy; mSigYZ = syz; mSigZZ = szz;
  }

  virtual void print() const;
  std::string asString() const;
  
 protected:
  float mPos[3] = {};          // cartesian position of cluster in lab frame
  float mSigYY = 0.f;          // covariance matrix elements
  float mSigYZ = 0.f;          // covariance matrix elements
  float mSigZZ = 0.f;          // covariance matrix elements
  int   mCluSiz = 0;           // cluster size
  int   mParticleID = 0;       // particle ID in Kine tree (MC truth)
  int   mHitID = -1;           // hit ID (for test of hitsToRecPoints)
  short mDetectorID = 0;       // the detector/sensor id
  int8_t mLayer = -1;

  ClassDefNV(NA6PBaseCluster, 1);
};

#endif
