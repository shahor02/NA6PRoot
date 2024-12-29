// NA6PCCopyright

#include "NA6PMCGenHeader.h"
#include <fairlogger/Logger.h>
#include <fmt/format.h>

std::string NA6PMCGenHeader::asString(bool full) const
{
  std::string rep = fmt::format("Gen:{} NPrim:{}({}) NSec:{}({})", mInfo, mNPrimaries, mPrimariesOffset, mNSecondaries, mSecondariesOffset);
  if (full && mKeyVals.size()) {
    rep += " [";
    for (auto& v : mKeyVals) {
      rep += fmt::format("({} : {})", v.first, v.second);
    }
    rep += "]";
  }
  return rep;
}

void NA6PMCGenHeader::print(bool full) const
{
  LOGP(info, "{}", asString(full));
}

float NA6PMCGenHeader::getValue(const std::string& k) const
{
  auto ent = mKeyVals.find(k);
  return ent != mKeyVals.end() ? ent->second : 0.f;
}

void NA6PMCGenHeader::clear()
{
  mInfo = "";
  mKeyVals.clear();
  mNPrimaries = mNSecondaries = mPrimariesOffset = mSecondariesOffset = 0;
}
