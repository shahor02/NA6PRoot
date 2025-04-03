// NA6PCCopyright
#ifndef _MAGFIELD_REGION_H
#define _MAGFIELD_REGION_H

#include <Rtypes.h>

class MagneticFieldRegion
{
 public:
  const std::string& getName() const { return mName; }
  float getRefX() const { return mRefPos[0]; }
  float getRefY() const { return mRefPos[1]; }
  float getRefZ() const { return mRefPos[2]; }
  float getXMin() const { return mXMin; }
  float getYMin() const { return mYMin; }
  float getZMin() const { return mZMin; }
  float getXMax() const { return mXMax; }
  float getYMax() const { return mYMax; }
  float getZMax() const { return mZMax; }

  void loadFlukaField(const std::string& filename);

  template <typename T = float>
  bool isInside(const T* xyz) const
  {
    return xyz[2] >= mBoxPos[2][0] && xyz[2] <= mBoxPos[2][1] &&
           xyz[1] >= mBoxPos[1][0] && xyz[1] <= mBoxPos[1][1] &&
           xyz[0] >= mBoxPos[0][0] && xyz[0] <= mBoxPos[0][1];
  }

  template <typename T = float>
  bool getField(const T* xyz, T* bxbybz) const
  {
    float z = xyz[2] - mRefPos[2];
    if (z < mZMin || z > mZMax) {
      return false;
    }
    float y = xyz[1] - mRefPos[1];
    if (y < mYMin || y > mYMax) {
      return false;
    }
    float x = xyz[0] - mRefPos[0];
    if (x < mXMin || x > mXMax) {
      return false;
    }
    // Find the indices of the surrounding grid points
    int i = findGridIndex(x, mXMin, mDXI, mNX1);
    int j = findGridIndex(y, mYMin, mDYI, mNY1);
    int k = findGridIndex(z, mZMin, mDZI, mNZ1);

    // Perform second-degree interpolation
    interpolateField(i, j, k, x, y, z, bxbybz);
    return true;
  }

  template <typename T = float>
  bool addField(const T* xyz, T* bxbybz) const
  {
    float z = xyz[2] - mRefPos[2];
    if (z < mZMin || z > mZMax) {
      return false;
    }
    float y = xyz[1] - mRefPos[1];
    if (y < mYMin || y > mYMax) {
      return false;
    }
    float x = xyz[0] - mRefPos[0];
    if (x < mXMin || x > mXMax) {
      return false;
    }
    // Find the indices of the surrounding grid points
    int i = findGridIndex(x, mXMin, mDXI, mNX1);
    int j = findGridIndex(y, mYMin, mDYI, mNY1);
    int k = findGridIndex(z, mZMin, mDZI, mNZ1);

    // Perform second-degree interpolation
    interpolateField(i, j, k, x, y, z, bxbybz);
    return true;
  }

  void setRefPosition(float x, float y, float z);
  void setName(const std::string& s) { mName = s; }

 private:
  void cacheValues(const std::string& line, std::vector<float>& cachev);
  void interpolateField(int i, int j, int k, float x, float y, float z, double* bxbybz) const;
  void interpolateFieldAdd(int i, int j, int k, float x, float y, float z, double* bxbybz) const;
  float getFieldComponent(int i, int j, int k, int dim) const
  {
    //    int index = i * ny * nz + j * nz + k;
    int index = i + mNX * (j + mNY * k);
    return mFieldData[3 * index + dim];
  }
  int findGridIndex(float value, float minValue, float stepInv, int numBins) const
  {
    auto ind = int((value - minValue) * stepInv);
    return ind < numBins ? (ind >= 0 ? ind : 0) : numBins - 1;
  }

  float mRefPos[3] = {}; // position of the field reference point
  float mBoxPos[3][2];   // box position in lab

  float mXMin = 0.f, mXMax = 0.f, mYMin = 0.f, mYMax = 0.f, mZMin = 0.f, mZMax = 0.f; // field extent wrt reference
  int mNX = 0, mNY = 0, mNZ = 0, mNX1 = 0, mNY1 = 0, mNZ1 = 0;
  float mDX = -1.f, mDY = -1.f, mDZ = -1.f;
  float mDXI = -1.f, mDYI = -1.f, mDZI = -1.f;
  std::vector<float> mFieldData;
  std::string mName{};

  ClassDefNV(MagneticFieldRegion, 1);
};

#endif
