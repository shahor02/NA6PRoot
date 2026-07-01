// NA6PCCopyright

#include "NA6PVerTelDigit.h"
#include <fmt/format.h>
#include <fairlogger/Logger.h>

NA6PVerTelDigit::NA6PVerTelDigit(uint16_t detID, const VTPixID& id) : mDetectorID(detID), mPixID(id)
{
}

NA6PVerTelDigit::NA6PVerTelDigit(uint16_t detID, uint32_t rsu, uint32_t tile, uint32_t row, uint32_t col) : mDetectorID(detID)
{
  mPixID.rsu = rsu;
  mPixID.tile = tile;
  mPixID.row = row;
  mPixID.col = col;
}

std::string NA6PVerTelDigit::asString() const
{
  return fmt::format("Digit: Det:{} RSU:{} Tile:{} Row:{} Col:{}",
                     mDetectorID, getRSU(), getTile(), getRow(), getCol());
}

void NA6PVerTelDigit::print() const
{
  LOGP(info, "{}", asString());
}
