// NA6PCCopyright

#include "NA6PVerTelDigit.h"
#include <fmt/format.h>
#include <fairlogger/Logger.h>

NA6PVerTelDigit::NA6PVerTelDigit(UShort_t detID, UShort_t rsu, UShort_t tile, UShort_t row, UShort_t col, int pid) : mDetectorID(detID), mParticleID(pid)
{
  mPixID.rsu = rsu;
  mPixID.tile = tile;
  mPixID.row = row;
  mPixID.col = col;
}

std::string NA6PVerTelDigit::asString() const
{
  return fmt::format("Digit: Det:{} RSU:{} Tile:{} Row:{} Col:{} PartId:{}",
                     mDetectorID, getRSU(), getTile(), getRow(), getCol(), mParticleID);
}

void NA6PVerTelDigit::print() const
{
  LOGP(info, "{}", asString());
}
