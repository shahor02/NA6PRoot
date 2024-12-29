// NA6PCCopyright
#ifndef NA6P_DIPOLEIP_H_
#define NA6P_DIPOLEIP_H_

#include "NA6PModule.h"

class TGeoVolume;

class NA6PDipoleIP : public NA6PModule
{
 public:
  NA6PDipoleIP() : NA6PModule("DipoleIP") {}
  ~NA6PDipoleIP() override = default;
  void createMaterials() override;
  void createGeometry(TGeoVolume *world) override;
};

#endif
