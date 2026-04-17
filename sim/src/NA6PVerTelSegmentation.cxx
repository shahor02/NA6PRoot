// NA6PCCopyright

#include "NA6PVerTelSegmentation.h"

const float NA6PVerTelSegmentation::DeadXShort = 0.15;    // short end cap
const float NA6PVerTelSegmentation::DeadXLong = 0.45;     // long end cap   (readout)
const float NA6PVerTelSegmentation::DeadYBottom = 0.0525; // bottom dead zone
const float NA6PVerTelSegmentation::DeadYTop = 0.0525;    // top dead zone

const float NA6PVerTelSegmentation::DeadTopBotHalves = 0.012;                    // dead space between top and bottom halves
const float NA6PVerTelSegmentation::DeadXTile = 0.002;                           // dead space between tiles in X (apart from DeadXDataBackbone after each 3 tiles)
const float NA6PVerTelSegmentation::DeadXDataBackbone = 0.006;                   // dead space between segments
const float NA6PVerTelSegmentation::XSizeTot = 12.9996 + DeadXShort + DeadXLong; // 13.5996;
const float NA6PVerTelSegmentation::YSizeTot = 13.5898 + DeadYBottom + DeadYTop; // 13.6948; // readout side
const int NA6PVerTelSegmentation::NXTiles = 36;
const int NA6PVerTelSegmentation::NXSegments = 12; // group of 3 tiles
const int NA6PVerTelSegmentation::NYSensors = 7;
const float NA6PVerTelSegmentation::DYSens = YSizeTot / NYSensors;
const float NA6PVerTelSegmentation::DXSegment = (XSizeTot - DeadXShort - DeadXLong) / NXSegments;
const float NA6PVerTelSegmentation::DXTile = (DXSegment - DeadXDataBackbone) / 3;

int NA6PVerTelSegmentation::isInAcc(float x, float y) const
{

  float absOffX = std::abs(mOffsX) - mInterChipGap / 2;
  float absOffY = std::abs(mOffsY) - mInterChipGap / 2;
  if ((x < -absOffX && y < absOffY) ||
      (x > -absOffX && y < -absOffY)) {
    x = -x;
    y = -y;
  }

  x -= mOffsX;
  if (x < 0.) {
    y += mOffsY - mInterChipGap / 2;
    x += XSizeTot + mInterChipGap / 2; // relate to chip local left/bottom corner
  } else {
    y -= mOffsY + mInterChipGap / 2;
    y = YSizeTot - y;
    x = (XSizeTot + mInterChipGap / 2) - x; // relate to chip local left/bottom corner
  }

  if (x < 0 || x > XSizeTot || y < 0 || y > YSizeTot)
    return 0;

  if (x < DeadXLong || x > XSizeTot - DeadXShort || y < DeadYBottom || y > YSizeTot - DeadYTop)
    return -1;
  x -= DeadXLong;
  int sens = int(y / DYSens);
  y -= sens * DYSens;
  if (y < DeadYBottom || y > DYSens - DeadYTop)
    return -1;
  int topbot = int(y / (DYSens / 2));
  y -= topbot * (DYSens / 2);
  if (y < DeadTopBotHalves / 2)
    return -1;
  int segm = int(x / DXSegment);
  x -= segm * DXSegment;
  if (x > DXSegment - DeadXDataBackbone)
    return -1;
  int tile = int(x / DXTile);
  x -= tile * DXTile;
  if (x < DeadXTile)
    return -1;

  return 1;
}
