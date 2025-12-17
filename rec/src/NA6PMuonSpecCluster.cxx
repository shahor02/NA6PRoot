// NA6PCCopyright

#include "NA6PMuonSpecCluster.h"
#include "NA6PLayoutParam.h"
#include <cmath>

NA6PMuonSpecCluster::NA6PMuonSpecCluster(float x, float y, float z, int clusiz)
  : NA6PBaseCluster(x, y, z, clusiz)
{
  // Set layer based on Z position and layout parameters
  float zCluster = getZLab();
  const auto& layout = NA6PLayoutParam::Instance();
  const float dzWindow = 20.f; // half-width window around plane Z
  for (int i = 0; i < layout.nMSPlanes; ++i) {
    float zPlane = layout.posMSPlaneZ[i];
    if (zCluster > (zPlane - dzWindow) && zCluster < (zPlane + dzWindow)) {
      mLayer = i + layout.nVerTelPlanes;
      return;
    }
  }
}