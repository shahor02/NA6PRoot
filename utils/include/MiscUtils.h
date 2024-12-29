// NA6PCCopyright
#ifndef NA6P_MISCUTILS_H
#define NA6P_MISCUTILS_H

#include <string>
#include <Rtypes.h>

// set of various static utils

class MiscUtils
{
 public:
  static void silenceStdOut(const std::string& warn = "");
  static void reviveStdOut(const std::string& warn = "");

  ClassDefNV(MiscUtils,1);
};

#endif
