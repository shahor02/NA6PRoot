// NA6PCCopyright

#include "NA6PVerTelCluster.h"
#include "NA6PLayoutParam.h"

int NA6PVerTelCluster::getLayer() const {
	float z = getZLab();
	const auto& layout = NA6PLayoutParam::Instance();
	const float dzWindow = 2.5f; // tighter window for VT planes
	for (int i = 0; i < layout.nVerTelPlanes; ++i) {
		float zPlane = layout.posVerTelPlaneZ[i];
		if (z > (zPlane - dzWindow) && z < (zPlane + dzWindow)) {
			return i;
		}
	}
	return -1;
}
// NA6PCCopyright

#include "NA6PVerTelCluster.h"
