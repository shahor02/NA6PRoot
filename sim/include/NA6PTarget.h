// NA6PCCopyright

#ifndef NA6P_TARGET_H_
#define NA6P_TARGET_H_

#include "NA6PModule.h"

class TGeoVolume;

class NA6PTarget : public NA6PModule
{
 public:
  NA6PTarget() : NA6PModule("Target") {}
  ~NA6PTarget() override = default;
  void createMaterials() override;
  void createGeometry(TGeoVolume *world) override;
};


#endif
