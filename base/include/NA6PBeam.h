// NA6PCCopyright

#ifndef NA6P_BEAM_H_
#define NA6P_BEAM_H_

#include <Rtypes.h>

class NA6PBeam
{
 public:
  float getX() const { return mX; }
  float getY() const { return mY; }
  float getZ() const { return mZ; }
  float getSlopeX() const { return mSlopeX; }
  float getSlopeY() const { return mSlopeY; }
  float getX(float z) const { return mX + mSlopeX*(z - mZ); }
  float getY(float z) const { return mY + mSlopeY*(z - mZ); }

  void setX(float x) { mX = x; }
  void setY(float y) { mY = y; }
  void setZ(float z) { mZ = z; }  
  void setSlopeX(float slopeX) { mSlopeX = slopeX; }
  void setSlopeY(float slopeY) { mSlopeY = slopeY; }


  
 protected:
  float mX = 0.f;
  float mY = 0.f;
  float mZ = 0.f;
  float mSlopeX = 0.f;
  float mSlopeY = 0.f;

  ClassDefNV(NA6PBeam, 1);
};

#endif
