// NA6PCCopyright

#include <fmt/format.h>
#include <fairlogger/Logger.h>
#include "NA6PVerTelCluster.h"
#include "NA6PMuonSpecCluster.h"
#include "NA6PTrack.h"
#include "NA6PLayoutParam.h"

//_______________________________________________________________________
NA6PTrack::NA6PTrack()
{
  mClusterIndices.fill(-1);
  mClusterPartID.fill(-2);
}

//_______________________________________________________________________
NA6PTrack::NA6PTrack(const float* xyz, const float* pxyz, int sign, float errLoose)
{
  mClusterIndices.fill(-1);
  mClusterPartID.fill(-2);
  init(xyz, pxyz, sign, errLoose);
}

//_______________________________________________________________________
void NA6PTrack::reset()
{
  mChi2 = 0;
  mChi2Outer = 0;
  mClusterMap = 0;
  mStatusRefitInward = false;
  mStatusConstrained = false;
  resetCovariance();
  mParticleID = -2;
  mClusterIndices.fill(-1);
  mClusterPartID.fill(-2);
  mNClusters = 0;
}

//_______________________________________________________________________
std::string NA6PTrack::asString() const
{
  auto pxyz = getPXYZ();
  return fmt::format("Track: Nclusters:{} chi2:{} chi2:{} pos:{:.4f},{:.4f},{:.4f} mom:{:.3f},{:.3f},{:.3f}",
                     mNClusters, mChi2, mChi2Outer, getX(), getY(), getZ(), pxyz[0], pxyz[1], pxyz[2]);
}

//_______________________________________________________________________
void NA6PTrack::print() const
{
  LOGP(info, "{}", asString());
}

//_______________________________________
template <typename ClusterType>
void NA6PTrack::addCluster(const ClusterType* clu)
{

  mNClusters++;
  int trackID = clu->getParticleID();
  int nLay = clu->getLayer();
  mClusterPartID[nLay] = trackID;
  mClusterIndices[nLay] = clu->getClusterIndex();
  mClusterMap |= (1 << nLay);
}

template void NA6PTrack::addCluster<NA6PBaseCluster>(const NA6PBaseCluster*);
template void NA6PTrack::addCluster<NA6PVerTelCluster>(const NA6PVerTelCluster*);
template void NA6PTrack::addCluster<NA6PMuonSpecCluster>(const NA6PMuonSpecCluster*);
