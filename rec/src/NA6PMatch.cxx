// NA6PCCopyright

#include <fmt/format.h>
#include <fairlogger/Logger.h>
#include "NA6PMatch.h"

//_______________________________________________________________________
std::string NA6PMatch::asString() const
{
  return fmt::format("Match: VTid:{} MSid:{} chi2Match:{:.2f} chi2NormRefit:{:.2f} {}",
                     mIndexVT, mIndexMS, getChi2Match(), getChi2Norm(), ((NA6PTrackParCov*)this)->asString());
}

//_______________________________________________________________________
void NA6PMatch::print() const
{
  LOGP(info, "{}", asString());
}
