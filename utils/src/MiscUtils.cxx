// NA6PCCopyright
#include "MiscUtils.h"
#include <fairlogger/Logger.h>
#include <iostream>


void MiscUtils::silenceStdOut(const std::string& warn)
{
  if (!warn.empty()) {
    LOGP(warn, "Silencing stdout {}", warn);
  }
  std::cout.setstate(std::ios::failbit);
}

void MiscUtils::reviveStdOut(const std::string& warn)
{
  std::cout.clear();
}

