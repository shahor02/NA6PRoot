// NA6PCCopyright
#ifndef NA6P_MODULE_H_
#define NA6P_MODULE_H_

#include <string>

class TGeoVolume;

// Base (pasive or active) module of NA6P

class NA6PModule
{
 public:
  NA6PModule(std::string name) : mName(name) {}
  virtual ~NA6PModule() = default;
  virtual void createMaterials() = 0;
  virtual void createGeometry(TGeoVolume *base) = 0;

  bool isActive() const { return mActive; }
  void setActive(bool v) { mActive = v; }
  const std::string& getName() const { return mName; }
  
 protected:
  std::string mName{};
  bool mActive{false};
};

#endif
