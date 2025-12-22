// NA6PCCopyright

#include "NA6PBaseCluster.h"
#include <fmt/format.h>
#include <fairlogger/Logger.h>

NA6PBaseCluster::NA6PBaseCluster(float x, float y, float z, int clusiz, int layer)
  : mPos{x, y, z}, mCluSiz(clusiz), mLayer(layer)
{
}

std::string NA6PBaseCluster::asString() const
{
  return fmt::format("Cluster: Det:{} XYZlab:{:.3f},{:.3f},{:.3f} XYZtrackingframe:{:.3f},{:.3f},{:.3f} cluster size:{} ParticleID:{}",
                     mDetectorID, getXLab(), getYLab(), getZLab(), getXTF(), getYTF(), getZTF(), mCluSiz, mParticleID);
}

void NA6PBaseCluster::print() const
{
  LOGP(info, "{}", asString());
}
