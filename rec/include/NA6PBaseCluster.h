// NA6PCCopyright

#ifndef NA6P_BASE_CLUSTER_H
#define NA6P_BASE_CLUSTER_H

#include <string>
#include <Rtypes.h>

// Mother class of all cluster classes

class NA6PBaseCluster
{
 public:

  NA6PBaseCluster(float x, float y, float z, int clusiz = 0, int layer = 0);
  NA6PBaseCluster(float x, float y, float z, float sxx, float syx, float syy, int clusiz = 0, int layer = 0);
  NA6PBaseCluster() = default;
  NA6PBaseCluster(const NA6PBaseCluster&) = default;
  NA6PBaseCluster& operator=(const NA6PBaseCluster&) = default;
  ~NA6PBaseCluster() {}

  int getLayer() const { return mLayer; }
  int getClusterSize() const { return mCluSiz; }
  auto getX()          const {return mX;}
  auto getY()          const {return mY;}
  auto getZ()          const {return mZ;}
  auto getSigXX()      const { return mSigXX; }
  auto getSigYX()      const { return mSigYX; }
  auto getSigYY()      const { return mSigYY; }
  
  auto getDetectorID() const { return mDetectorID; }
  auto getParticleID() const { return mParticleID; }
  auto getHitID()      const { return mHitID; }

  void setDetectorID(int id) { mDetectorID = id; }
  void setParticleID(int id) { mParticleID = id; }
  void setHitID(int id) { mHitID = id; }
  void setPos(float x, float y, float z){
    mX = x;
    mY = y;
    mZ = z;
  }
  void setX(float v)    {mX = v;}
  void setY(float v)    {mY = v;}
  void setZ(float v)    {mZ = v;}
  void setErr(float sxx, float syx, float syy){
    mSigXX = sxx;
    mSigYX = syx;
    mSigXX = syy;
  }
  
  void print() const;
  std::string asString() const;
  
 protected:
  float mX = 0.f;
  float mY = 0.f;
  float mZ = 0.f;
  float mSigXX = 0.f;
  float mSigYY = 0.f;
  float mSigYX = 0.f;
  int   mCluSiz = 0;           // cluster size
  int   mParticleID = 0;       // particle ID in Kine tree (MC truth)
  int   mHitID = -1;           // hit ID (for test of hitsToRecPoints)
  short mDetectorID = 0;       // the detector/sensor id
  int8_t mLayer = -1;

  ClassDefNV(NA6PBaseCluster, 1);
};

#endif
