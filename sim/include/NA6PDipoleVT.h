// NA6PCCopyright
#ifndef NA6P_DIPOLEVT_H_
#define NA6P_DIPOLEVT_H_

#include "NA6PModule.h"

class TGeoVolume;

class NA6PDipoleVT : public NA6PModule
{
 public:
  NA6PDipoleVT() : NA6PModule("DipoleVT") {}
  ~NA6PDipoleVT() override = default;
  void createMaterials() override;
  void createGeometry(TGeoVolume* world) override;
};

#endif
