// NA6PCCopyright

#include <fmt/format.h>
#include <fairlogger/Logger.h>
#include <TGeoGlobalMagField.h>
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
NA6PTrack::NA6PTrack(const float* xyz, const float* pxyz, int sign, float errLoose) : NA6PTrackParCov(xyz, pxyz, sign, errLoose)
{
  // initialize arrays
  mClusterIndices.fill(-1);
  mClusterPartID.fill(-2);
}

//_______________________________________________________________________
void NA6PTrack::reset()
{
  mChi2 = 0;
  mClusterMap = 0;
  resetCovariance(-1);
  mParticleID = -2;
  mClusterIndices.fill(-1);
  mClusterPartID.fill(-2);
  mNClusters = mNClusters = mNClustersMS = mNClustersTR = 0;
}

//_______________________________________________________________________
std::string NA6PTrack::asString() const
{
  float pxyz[3];
  getPXYZ(pxyz);
  return fmt::format("Track: Nclusters:{} NVTclusters:{} NMSclusters:{} NTRclusters:{} chi2:{} chi2VT:{} chi2MS:{} pos:{:.4f},{:.4f},{:.4f} mom:{:.3f},{:.3f},{:.3f}",
                     mNClusters, mNClustersVT, mNClustersMS, mNClustersTR, mChi2, mChi2VT, mChi2MS, getX(), getY(), getZ(), pxyz[0], pxyz[1], pxyz[2]);
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
  mChi2 += chi2;
  int trackID = clu->getParticleID();

  int nLay = clu->getLayer();
  mClusterPartID[nLay] = trackID;
  mClusterIndices[nLay] = cluIndex;
  mClusterMap |= (1<<nLay);

  if (nLay < NA6PLayoutParam::Instance().nVerTelPlanes) {
    mNClustersVT++;
    mChi2VT += chi2;
  }
  else if (nLay >= NA6PLayoutParam::Instance().nVerTelPlanes + NA6PLayoutParam::Instance().nMSPlanes - 2) {
    mNClustersMS++;
    mChi2MS += chi2;
  }
  else {
    mNClustersTR++;
    mChi2MS += chi2;
  }
}

template void NA6PTrack::addCluster<NA6PBaseCluster>(const NA6PBaseCluster*, int, float);
template void NA6PTrack::addCluster<NA6PVerTelCluster>(const NA6PVerTelCluster*, int, float);
template void NA6PTrack::addCluster<NA6PMuonSpecCluster>(const NA6PMuonSpecCluster*, int, float);
