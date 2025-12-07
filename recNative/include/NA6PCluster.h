// NA6PCCopyright

#ifndef NA6P_CLUSTER_H
#define NA6P_CLUSTER_H

#include <Rtypes.h>

// base cluster

class NA6PCluster
{
 public:
  NA6PCluster() = default;
  NA6PCluster(float x, float y, float z, float sxx, float sxy, float syy, short sensID, short ui) : mX(x), mY(y), mZ(z), mSigXX(sxx), mSigYY(syy), mSigXY(sxy), mSensorID(sensID), mUserInfo(ui) {}
  auto getX() const { return mX; }
  auto getY() const { return mY; }
  auto getZ() const { return mZ; }
  auto getSigXX() const { return mSigXX; }
  auto getSigYY() const { return mSigYY; }
  auto getSigXY() const { return mSigXY; }
  auto getSensorID() const { return mSensorID; }
  auto getParticleID() const { return mUserInfo; }

  void setX(float v) { mX = v; }
  void setY(float v) { mY = v; }
  void setZ(float v) { mZ = v; }
  void setSigXX(float v) { mSigXX = v; }
  void setSigXY(float v) { mSigXY = v; }
  void setSigYY(float v) { mSigYY = v; }
  void setSensorID(short v) { mSensorID = v; }
  void setUserInfo(short v) { mUserInfo = v; }
  void print() const;
  std::string asString() const;

 protected:
  float mX = 0.f;
  float mY = 0.f;
  float mZ = 0.f;
  float mSigXX = 0.f;
  float mSigYY = 0.f;
  float mSigXY = 0.f;
  short mSensorID = -1;
  short mUserInfo = 0;

  ClassDefNV(NA6PCluster, 1);
};

#endif
