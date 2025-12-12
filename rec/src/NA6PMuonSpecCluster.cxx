// NA6PCCopyright

#include "NA6PMuonSpecCluster.h"
#include "NA6PLayoutParam.h"
#include <cmath>

int NA6PMuonSpecCluster::getLayer() const {
	float z = getZLab();
	const auto& layout = NA6PLayoutParam::Instance();
	const float dzWindow = 20.f; // half-width window around plane Z
	for (int i = 0; i < layout.nMSPlanes; ++i) {
		float zPlane = layout.posMSPlaneZ[i];
		if (z > (zPlane - dzWindow) && z < (zPlane + dzWindow)) {
			return i;
		}
	}
	return -1;
}

