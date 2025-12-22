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
  void loadOpera3DField(const std::string& filename);

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
    interpolateFieldAdd(i, j, k, x, y, z, bxbybz);
    return true;
  }

  void setRefPosition(float x, float y, float z);
  void setName(const std::string& s) { mName = s; }

 private:
  void cacheValues(const std::string& line, std::vector<float>& cachev);
  template <typename T = float>
  void interpolateField(int i, int j, int k, float x, float y, float z, T* bxbybz) const;
  template <typename T = float>
  void interpolateFieldAdd(int i, int j, int k, float x, float y, float z, T* bxbybz) const;
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

template <typename T>
void MagneticFieldRegion::interpolateField(int i, int j, int k, float x, float y, float z, T* bxbybz) const
{
  float x0 = mXMin + i * mDX;
  float y0 = mYMin + j * mDY;
  float z0 = mZMin + k * mDZ;
  float tx = (x - x0) * mDXI, tx1 = 1.f - tx;
  float ty = (y - y0) * mDYI, ty1 = 1.f - ty;
  float tz = (z - z0) * mDZI, tz1 = 1.f - tz;

  /*
  // parabolic interpolation
  // Perform trilinear interpolation
  for (int dim = 0; dim < 3; ++dim) {
    bxbybz[dim] = tx1 * ty1 * tz1 * getFieldComponent(i, j, k, dim) +
      tx * ty1 * tz1 * getFieldComponent(i + 1, j, k, dim) +
      tx1 * ty * tz1 * getFieldComponent(i, j + 1, k, dim) +
      tx * ty * tz1 * getFieldComponent(i + 1, j + 1, k, dim) +
      tx1 * ty1 * tz * getFieldComponent(i, j, k + 1, dim) +
      tx * ty1 * tz * getFieldComponent(i + 1, j, k + 1, dim) +
      tx1 * ty * tz * getFieldComponent(i, j + 1, k + 1, dim) +
      tx * ty * tz * getFieldComponent(i + 1, j + 1, k + 1, dim);
    }
  */
  // linear interpolation
  for (int dim = 0; dim < 3; ++dim) {
    float c00 = tx1 * getFieldComponent(i, j, k, dim) + tx * getFieldComponent(i + 1, j, k, dim);
    float c01 = tx1 * getFieldComponent(i, j, k + 1, dim) + tx * getFieldComponent(i + 1, j, k + 1, dim);
    float c10 = tx1 * getFieldComponent(i, j + 1, k, dim) + tx * getFieldComponent(i + 1, j + 1, k, dim);
    float c11 = tx1 * getFieldComponent(i, j + 1, k + 1, dim) + tx * getFieldComponent(i + 1, j + 1, k + 1, dim);

    float c0 = ty1 * c00 + ty * c10;
    float c1 = ty1 * c01 + ty * c11;

    bxbybz[dim] = tz1 * c0 + tz * c1;
  }
}

template <typename T>
void MagneticFieldRegion::interpolateFieldAdd(int i, int j, int k, float x, float y, float z, T* bxbybz) const
{
  float x0 = mXMin + i * mDX;
  float y0 = mYMin + j * mDY;
  float z0 = mZMin + k * mDZ;
  float tx = (x - x0) * mDXI, tx1 = 1.f - tx;
  float ty = (y - y0) * mDYI, ty1 = 1.f - ty;
  float tz = (z - z0) * mDZI, tz1 = 1.f - tz;

  /*
  // parabolic interpolation
  // Perform trilinear interpolation
  for (int dim = 0; dim < 3; ++dim) {
    bxbybz[dim] += tx1 * ty1 * tz1 * getFieldComponent(i, j, k, dim) +
      tx * ty1 * tz1 * getFieldComponent(i + 1, j, k, dim) +
      tx1 * ty * tz1 * getFieldComponent(i, j + 1, k, dim) +
      tx * ty * tz1 * getFieldComponent(i + 1, j + 1, k, dim) +
      tx1 * ty1 * tz * getFieldComponent(i, j, k + 1, dim) +
      tx * ty1 * tz * getFieldComponent(i + 1, j, k + 1, dim) +
      tx1 * ty * tz * getFieldComponent(i, j + 1, k + 1, dim) +
      tx * ty * tz * getFieldComponent(i + 1, j + 1, k + 1, dim);
    }
  */
  // linear interpolation
  for (int dim = 0; dim < 3; ++dim) {
    float c00 = tx1 * getFieldComponent(i, j, k, dim) + tx * getFieldComponent(i + 1, j, k, dim);
    float c01 = tx1 * getFieldComponent(i, j, k + 1, dim) + tx * getFieldComponent(i + 1, j, k + 1, dim);
    float c10 = tx1 * getFieldComponent(i, j + 1, k, dim) + tx * getFieldComponent(i + 1, j + 1, k, dim);
    float c11 = tx1 * getFieldComponent(i, j + 1, k + 1, dim) + tx * getFieldComponent(i + 1, j + 1, k + 1, dim);

    float c0 = ty1 * c00 + ty * c10;
    float c1 = ty1 * c01 + ty * c11;

    bxbybz[dim] += tz1 * c0 + tz * c1;
  }
}

#endif
