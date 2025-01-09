// NA6PCCopyright

#include "NA6PBaseHit.h"
#include <fmt/format.h>
#include <fairlogger/Logger.h>

NA6PBaseHit::NA6PBaseHit(int trackid, short did, const TVector3& xyzIn, const TVector3& xyzOut, const TVector3& momIn, float time, float val, uint8_t statusStart, uint8_t statusEnd) : mTime(time), mHitValue(val), mTrackID(trackid), mDetectorID(did), mTrackStatusStart(statusStart), mTrackStatusEnd(statusEnd)
{
  setPosIn(xyzIn);
  setPosOut(xyzOut);
  setMomIn(momIn);
}

std::string NA6PBaseHit::asString() const
{
  return fmt::format("Hit: Det:{} Track:{} T:{} V:{} XYZIn:{:.3f},{:.3f},{:.3f} XYZOut:{:.3f},{:.3f},{:.3f} PIn:{:.3f},{:.3f},{:.3f} Status:{:B}:{:B}",
                     mDetectorID, mTrackID, mTime, mHitValue, getXIn(), getYIn(), getZIn(), getXOut(), getYOut(), getZOut(), getPX(), getPY(), getPZ(), mTrackStatusStart, mTrackStatusEnd);
}

void NA6PBaseHit::print() const
{
  LOGP(info, "{}", asString());
}
