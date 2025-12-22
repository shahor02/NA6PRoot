// NA6PCCopyright

#include "NA6PBaseCluster.h"
#include <fmt/format.h>
#include <fairlogger/Logger.h>

NA6PBaseCluster::NA6PBaseCluster(float x, float y, float z, int clusiz, int layer)
  : mX{x}, mY{y}, mZ{z}, mCluSiz{clusiz}, mLayer{int8_t(layer)}
{
}

NA6PBaseCluster::NA6PBaseCluster(float x, float y, float z, float sxx, float syx, float syy, int clusiz, int layer)
  : mX{x}, mY{y}, mZ{z}, mSigXX{sxx}, mSigYX{syx}, mSigYY{syy}, mCluSiz{clusiz}, mLayer{int8_t(layer)}
{
}

std::string NA6PBaseCluster::asString() const
{
  return fmt::format("Cluster: Det:{} XYZ:{:+.3e},{:+.3e},{:+.3e} Sig2:{:.3e},{:.3e},{:.3e} cluster size:{} ParticleID:{}",
                     mDetectorID, mX, mY, mZ, mSigXX, mSigYX, mSigYY, mCluSiz, mParticleID);
}

void NA6PBaseCluster::print() const
{
  LOGP(info, "{}", asString());
}
