// NA6PCCopyright

#include "NA6PVerTelCluster.h"
#include "NA6PLayoutParam.h"


NA6PVerTelCluster::NA6PVerTelCluster(float x, float y, float z, int clusiz, int nDet)
  : NA6PBaseCluster(x, y, z, clusiz)
{
	mDetectorID = nDet;
	mLayer = nDet / 4;
}