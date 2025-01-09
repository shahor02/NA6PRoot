// NA6PCCopyright

#include "NA6PMCEventHeader.h"
#include <fairlogger/Logger.h>
#include <fmt/format.h>

std::string NA6PMCEventHeader::asString(bool genh) const
{
  int nh = getNGenHeaders();
  auto rep = fmt::format("Run:{} Ev:{} Vtx:[{:6.3}, {:6.3}, {:6.3}] NTracks:{} NPrim:{} | {} GenHeaders",
                         mRunNumber, mEventID, mVX, mVY, mVZ, mNTracks, mNPrimaries, nh);
  if (genh) {
    rep += ": ";
    for (int ih = 0; ih < nh; ih++) {
      rep += fmt::format("(H{}: {})", ih, getGenHeader(ih).asString());
    }
  }
  return rep;
}

void NA6PMCEventHeader::print(bool genh) const
{
  LOGP(info, "{}", asString(genh));
}

void NA6PMCEventHeader::clear()
{
  mGenHeaders.clear();
  mEventID = 0;
  mVX = mVY = mVZ = 0.f;
  mNTracks = mNPrimaries = 0;
}
