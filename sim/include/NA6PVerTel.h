// NA6PCCopyright

#ifndef NA6P_VERTEL_H_
#define NA6P_VERTEL_H_

#include "NA6PModule.h"

class TGeoVolume;

class NA6PVerTel : public NA6PModule
{
 public:
  NA6PVerTel() : NA6PModule("VerTel") {}
  ~NA6PVerTel() override = default;
  void createMaterials() override;
  void createGeometry(TGeoVolume *world) override;
};


#endif
