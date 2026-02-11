// NA6PCCopyright
#ifndef NA6P_MAGNETS_GDML_H_
#define NA6P_MAGNETS_GDML_H_

#include "NA6PModule.h"
#include <string>

class TGeoVolume;

/// Module for importing magnet geometry from an external GDML file
/// (e.g. BothMagnets.gdml) and aligning it to layout positions.
class NA6PMagnetsGDML : public NA6PModule
{
 public:
  NA6PMagnetsGDML() : NA6PModule("MagnetsGDML") {}
  ~NA6PMagnetsGDML() override = default;

  void createMaterials() override;
  void createGeometry(TGeoVolume* world) override;

  // Public helper methods for GDML geometry processing
  void alignMagnetsToLayout(TGeoVolume* world);
  void assignMediaToGDMLVolumes();
  void harmoniseVolumeNames();
};

#endif
