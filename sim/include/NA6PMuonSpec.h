// NA6PCCopyright

#ifndef NA6P_MUONSPEC_H_
#define NA6P_MUONSPEC_H_

#include "NA6PModule.h"

class TGeoVolume;

class NA6PMuonSpec : public NA6PModule
{
 public:
  NA6PMuonSpec() : NA6PModule("MuonSpec") {}
  ~NA6PMuonSpec() override = default;
  void createMaterials() override;
  void createGeometry(TGeoVolume *world) override;
};

#endif
