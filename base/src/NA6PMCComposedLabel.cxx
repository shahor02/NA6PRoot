// NA6PCCopyright

#include <fmt/format.h>
#include <fairlogger/Logger.h>
#include "NA6PMCComposedLabel.h"

//_____________________________________________
std::string NA6PMCComposedLabel::asString() const
{
  if (isValid()) {
    return fmt::format("[sourceID:{}/eventID:{}/fake:{}/trackID:{:6d}]", getSourceID(), getEventID(), isFake() ? 'Y' : 'N', getTrackID());
  }
  return isNoise() ? "[noise]" : "[unset]";
}

//_____________________________________________
void NA6PMCComposedLabel::print() const
{
  LOGP(info, "{}", asString());
}
