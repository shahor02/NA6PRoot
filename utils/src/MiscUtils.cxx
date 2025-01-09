// NA6PCCopyright
#include "MiscUtils.h"
#include <TRandom.h>
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

std::pair<float, float> MiscUtils::genCorrelatedPair(float sig0, float sig1, float corr)
{
  float a, b;
  gRandom->Rannor(a, b);
  try {
    a *= sig0;
    b *= sig1;
    b = corr * a + std::sqrt(1.f - corr * corr) * b;
  } catch (std::exception const& e) {
    LOGP(fatal, "Failed to generate correlated pair with sig0={} sig1={} corr={}, reason: {}", sig0, sig1, corr, e.what());
  }
  return {a, b};
}
