// NA6PCCopyright

#include "NA6PVerTelSegmentation.h"

void NA6PVerTelSegmentation::setStaggered(bool val)
{
  mStaggered = val;
  if (val) {
    mDeadXLongEff = DeadXLong + DeadXShort;
    mDeadXShortEff = 0.f;
    mDeadYTopEff = DeadYTop + DeadYBottom;
    mDeadYBottomEff = 0.f;
  } else {
    mDeadXLongEff = DeadXLong;
    mDeadXShortEff = DeadXShort;
    mDeadYTopEff = DeadYTop;
    mDeadYBottomEff = DeadYBottom;
  }
}

int NA6PVerTelSegmentation::isInAcc(float x, float y) const
{

  float absOffX = std::abs(mOffsX) - mInterChipGap / 2;
  float absOffY = std::abs(mOffsY) - mInterChipGap / 2;
  if ((x < -absOffX && y < absOffY) ||
      (x > -absOffX && y < -absOffY)) {
    // flip signs for Q3 and Q4
    x = -x;
    y = -y;
  }

  x -= mOffsX;
  if (x < 0.) { // Q2 + Q4
    y += mOffsY - mInterChipGap / 2;
    x += XSizeTot + mInterChipGap / 2; // relate to chip local left/bottom corner
  } else {                             // Q1 + Q3
    y -= mOffsY + mInterChipGap / 2;
    x = (XSizeTot + mInterChipGap / 2) - x; // relate to chip local left/bottom corner
  }

  if (x < 0 || x > XSizeTot || y < 0 || y > YSizeTot)
    return 0;

  if (x < mDeadXLongEff || x > XSizeTot - mDeadXShortEff || y < mDeadYBottomEff || y > YSizeTot - mDeadYTopEff)
    return -1;
  x -= mDeadXLongEff;
  int sens = int(y / DYSens);
  y -= sens * DYSens;
  if (y < mDeadYBottomEff || y > DYSens - mDeadYTopEff)
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
