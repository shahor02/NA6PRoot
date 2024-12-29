// NA6PCCopyright
#ifndef NA6P_TGEOHELPER_H_
#define NA6P_TGEOHELPER_H_

#include <TGeoMedium.h>
#include <TGeoMatrix.h>
#include <TGeoMatrix.h>
#include <unordered_map>
#include <string>


class NA6PTGeoHelper
{
 public:
  static NA6PTGeoHelper& instance() {
    static NA6PTGeoHelper inst;
    return inst;
  }

  auto& getMatPool() { return mMatPool; }
  auto& getMedPool() { return mMedPool; }

  void addMedium(const std::string& medName, const std::string& matName="");
  TGeoMedium* getMedium(const std::string& medName) const;
  static TGeoRotation* rotAroundVector(float uX, float uY, float uZ, float ddelta);
  
 private:
  NA6PTGeoHelper() = default;
  ~NA6PTGeoHelper() = default;
  
  std::unordered_map<std::string, TGeoMaterial*> mMatPool;
  std::unordered_map<std::string, TGeoMedium*> mMedPool;

};

#endif
