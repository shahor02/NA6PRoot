// NA6PCCopyright

#include "NA6PCluster.h"
#include <fairlogger/Logger.h>

std::string NA6PCluster::asString() const
{
  return fmt::format("sens:{} XYZ:[{:+.3f},{:+.3f},{:+.3f}] cov:[{:.3f},{:+.3f},{:.3f}] UserInfo: {}", mSensorID, mX, mY, mZ, mSigXX, mSigXY, mSigYY, mUserInfo);
}

void NA6PCluster::print() const
{
  LOGP(info, "{}", asString());
}
