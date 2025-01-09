// NA6PCCopyright

#ifndef NA6P_MCGEN_HEADER_H
#define NA6P_MCGEN_HEADER_H

#include <Rtypes.h>
#include <string>
#include <unordered_map>

class NA6PMCGenHeader
{
 public:
  NA6PMCGenHeader() = default;
  NA6PMCGenHeader(uint32_t np, uint32_t ns, uint32_t proffs, uint32_t secoffs, const std::string& info) : mInfo(info), mNPrimaries(np), mNSecondaries(ns), mPrimariesOffset(proffs), mSecondariesOffset(secoffs) {}

  uint32_t getNPrimaries() const { return mNPrimaries; }
  uint32_t getNSecondaries() const { return mNSecondaries; }
  uint32_t getPrimariesOffset() const { return mPrimariesOffset; }
  uint32_t getSecondariesOffset() const { return mSecondariesOffset; }
  const std::string& getInfo() const { return mInfo; }
  std::string& getInfo() { return mInfo; }
  const std::unordered_map<std::string, float>& getKeyVals() const { return mKeyVals; }
  std::unordered_map<std::string, float>& getKeyVals() { return mKeyVals; }
  bool hasKey(const std::string& s) const { return mKeyVals.find(s) != mKeyVals.end(); }
  float getValue(const std::string& k) const;
  int getNKeyVals() const { return static_cast<int>(mKeyVals.size()); }

  void setNPrimaries(uint32_t nPrimaries) { mNPrimaries = nPrimaries; }
  void setNSecondaries(uint32_t nSecondaries) { mNSecondaries = nSecondaries; }
  void setPrimariesOffset(uint32_t primariesOffset) { mPrimariesOffset = primariesOffset; }
  void setSecondariesOffset(uint32_t secondariesOffset) { mSecondariesOffset = secondariesOffset; }
  void setInfo(const std::string& info) { mInfo = info; }
  void setKeyVals(const std::unordered_map<std::string, float>& keyVals) { mKeyVals = keyVals; }
  void addKeyVal(const std::string& k, float v) { mKeyVals[k] = v; }

  void clear();
  void print(bool full = true) const;
  std::string asString(bool full = true) const;

 protected:
  uint32_t mNPrimaries = 0;                        // number of primaries
  uint32_t mNSecondaries = 0;                      // number of secondaries
  uint32_t mPrimariesOffset = 0;                   // 1st primary entry in the kine tree
  uint32_t mSecondariesOffset = 0;                 // 1st primary entry in the kine tree
  std::string mInfo{};                             // text info
  std::unordered_map<std::string, float> mKeyVals; // optional metadata

  ClassDefNV(NA6PMCGenHeader, 1);
};

#endif
