// NA6PCCopyright
#ifndef NA6P_DETECTOR_H_
#define NA6P_DETECTOR_H_

#include <string>
#include <unordered_map>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoMatrix.h>

class NA6PModule;

// Helper for geometry builder
class NA6PDetector
{
 public:
  void createGeometry(const std::string& name = "NA6P");
  void createCommonMaterials();

  NA6PModule* getModule(const std::string& name) {
    auto m = mModules.find(name);
    return m != mModules.end() ? m->second : nullptr;
  }

  void addModule(NA6PModule* m);

 private:
  std::unordered_map<std::string, NA6PModule*> mModules;

};

#endif
