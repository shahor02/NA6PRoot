// NA6PCCopyright
#ifndef NA6P_TGEOHELPER_H_
#define NA6P_TGEOHELPER_H_

#include <TGeoMedium.h>
#include <TGeoMatrix.h>
#include <TGeoMatrix.h>
#include <TColor.h>
#include <unordered_map>
#include <string>

class NA6PTGeoHelper
{
 public:
  enum class EProc {
    kPAIR = 0,
    kCOMP,
    kPHOT,
    kPFIS,
    kDRAY,
    kANNI,
    kBREM,
    kHADR,
    kMUNU,
    kDCAY,
    kLOSS,
    kMULS,
    kCKOV,
    kLABS,
    kRAYL
  };

  /// cuts available
  enum class ECut {
    kCUTGAM = 0,
    kCUTELE,
    kCUTNEU,
    kCUTHAD,
    kCUTMUO,
    kBCUTE,
    kBCUTM,
    kDCUTE,
    kDCUTM,
    kPPCUTM,
    kTOFMAX
  };

  static NA6PTGeoHelper& instance()
  {
    static NA6PTGeoHelper inst;
    return inst;
  }

  auto& getMatPool() { return mMatPool; }
  auto& getMedPool() { return mMedPool; }
  const char* getMediumCutName(ECut cut) const;
  const char* getPhysicsProcessName(EProc process) const;

  void loadCutsAndProcessesFromFile(const std::string& fname);
  void addMedium(const std::string& medName, const std::string& matName = "", Color_t col = kGray);
  TGeoMedium* getMedium(const std::string& medName, bool fatalIfMissing = true) const;
  Color_t getMediumColor(const std::string& medName, bool fatalIfMissing = true) const;
  static TGeoRotation* rotAroundVector(float uX, float uY, float uZ, float ddelta);

 private:
  NA6PTGeoHelper() = default;
  ~NA6PTGeoHelper() = default;

  std::unordered_map<std::string, TGeoMaterial*> mMatPool;
  std::unordered_map<std::string, TGeoMedium*> mMedPool;
  std::unordered_map<std::string, Color_t> mColorPool;

  /// fixed names of cuts
  const static std::unordered_map<ECut, const char*> CutIDToName;
  /// fixed names of processes
  const static std::unordered_map<EProc, const char*> ProcessIDToName;
};

#endif
