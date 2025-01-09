// NA6PCCopyright
#ifndef NA6P_MISCUTILS_H
#define NA6P_MISCUTILS_H

#include <string>
#include <Rtypes.h>
#include <utility>

// set of various static utils

class MiscUtils
{
 public:
  static void silenceStdOut(const std::string& warn = "");
  static void reviveStdOut(const std::string& warn = "");

  static std::pair<float, float> genCorrelatedPair(float sig0, float sig1, float corr);
  
  ClassDefNV(MiscUtils,1);
};

#endif
