// NA6PCCopyright
#ifndef NA6P_DIPOLEMS_H_
#define NA6P_DIPOLEMS_H_

#include "NA6PModule.h"

class TGeoVolume;

class NA6PDipoleMS : public NA6PModule
{
 public:
  NA6PDipoleMS() : NA6PModule("DipoleMS") {}
  ~NA6PDipoleMS() override = default;
  void createMaterials() override;
  void createGeometry(TGeoVolume* world) override;
};

#endif
