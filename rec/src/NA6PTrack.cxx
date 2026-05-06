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
  // RSTODO consider suppressing this method
  mMatchChi2 = 0;
  mChi2VT = 0;
  mChi2MS = 0;
  mChi2VTOuter = 0;
  mChi2MSOuter = 0;
  mChi2VTRefit = 0;
  mChi2MSRefit = 0;
  mClusterMap = 0;
  mStatusRefitInward = false;
  mStatusConstrained = false;
  resetCovariance();
  mParticleID = -2;
  mClusterIndices.fill(-1);
  mClusterPartID.fill(-2);
  mNClusters = mNClusters = mNClustersMS = mNClustersTR = 0;
}

//_______________________________________________________________________
std::string NA6PTrack::asString() const
{
  auto pxyz = getPXYZ();
  return fmt::format("Track: Nclusters:{} NVTclusters:{} NMSclusters:{} NTRclusters:{} match chi2:{} chi2VT:{} chi2MS:{} pos:{:.4f},{:.4f},{:.4f} mom:{:.3f},{:.3f},{:.3f}",
                     mNClusters, mNClustersVT, mNClustersMS, mNClustersTR, mMatchChi2, mChi2VT, mChi2MS, getX(), getY(), getZ(), pxyz[0], pxyz[1], pxyz[2]);
}

//_______________________________________________________________________
void NA6PTrack::print() const
{
  LOGP(info, "{}", asString());
}

//_______________________________________
template <typename ClusterType>
void NA6PTrack::addCluster(const ClusterType* clu, int cluIndex, float chi2)
{

  mNClusters++;
  int trackID = clu->getParticleID();

  int nLay = clu->getLayer();
  mClusterPartID[nLay] = trackID;
  mClusterIndices[nLay] = cluIndex;
  mClusterMap |= (1 << nLay);

  if (nLay < NA6PLayoutParam::Instance().nVerTelPlanes) {
    mNClustersVT++;
    mChi2VT += chi2;
  } else if (nLay >= NA6PLayoutParam::Instance().nVerTelPlanes + NA6PLayoutParam::Instance().nMSPlanes - 2) {
    mNClustersTR++;
    mChi2MS += chi2;
  } else {
    mNClustersMS++;
    mChi2MS += chi2;
  }
}

bool NA6PTrack::updateTrack(const NA6PBaseCluster& cl, float maxChi2)
{
  auto chi2 = getPredictedChi2(cl);
  if (chi2 > maxChi2 || !update(cl)) {
    return false; // chi2 is too large
  }
  addCluster(&cl, cl.getClusterIndex(), chi2);
  return true;
}

template void NA6PTrack::addCluster<NA6PBaseCluster>(const NA6PBaseCluster*, int, float);
template void NA6PTrack::addCluster<NA6PVerTelCluster>(const NA6PVerTelCluster*, int, float);
template void NA6PTrack::addCluster<NA6PMuonSpecCluster>(const NA6PMuonSpecCluster*, int, float);
