// NA6PCCopyright
#ifndef NA6P_DETECTOR_H_
#define NA6P_DETECTOR_H_

#include <string>
#include <vector>
#include <unordered_map>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoMatrix.h>

class NA6PModule;

// Helper for geometry builder
class NA6PDetector
{
 public:
  static NA6PDetector& instance()
  {
    static NA6PDetector inst;
    return inst;
  }

  // Switch between analytic (C++) and GDML-based magnet geometry.
  // Must be called BEFORE the first call to instance(), i.e. before
  // NA6PMC is constructed.
  static void setUseGDMLMagnets(bool v);

  void createGeometry(const std::string& name = "NA6P");
  void createCommonMaterials();

  auto getNModules() const { return (int)mModulesVec.size(); }
  auto getNActiveModules() const { return (int)mActiveModulesVec.size(); }

  NA6PModule* getModule(const std::string& name)
  {
    auto m = mModulesMap.find(name);
    return m != mModulesMap.end() ? m->second : nullptr;
  }

  NA6PModule* getModule(int id)
  {
    return mModulesVec[id];
  }

  NA6PModule* getActiveModule(int id)
  {
    return mActiveModulesVec[id];
  }

  void addModule(NA6PModule* m);

  void setVerbosity(int v);
  auto getVerbosity() const { return mVerbosity; }

  void createHitsOutput(const std::string& outDir = "");
  void closeHitsOutput();
  void writeHits(const std::vector<int>& remapping);

 private:
  NA6PDetector();

  std::unordered_map<std::string, NA6PModule*> mModulesMap;
  std::vector<NA6PModule*> mModulesVec;
  std::vector<NA6PModule*> mActiveModulesVec;
  int mVerbosity = 0;
};

#endif
