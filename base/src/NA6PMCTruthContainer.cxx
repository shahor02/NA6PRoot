#include <fmt/format.h>
#include <fairlogger/Logger.h>
#include <TBuffer.h>
#include "NA6PMCTruthContainer.h"

//_____________________________________________
std::string NA6PMCTruthContainer::asString() const
{
  return fmt::format("NA6PMCTruthContainer index = {} for {} elements(s), flat buffer size {}", getIndexedSize(), getNElements(), mStreamerData.size());
}

//_____________________________________________
void NA6PMCTruthContainer::print() const
{
  LOGP(info, "{}", asString());
}

//_____________________________________________
void NA6PMCTruthContainer::Streamer(TBuffer& R__b)
{
  if (R__b.IsReading()) {
    R__b.ReadClassBuffer(NA6PMCTruthContainer::Class(), this);
    inflate();
  } else {
    deflate();
    R__b.WriteClassBuffer(NA6PMCTruthContainer::Class(), this);
    inflate();
  }
}
