// NA6PCCopyright

#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include "StringUtils.h"
#include "fairlogger/Logger.h"
#include "MagneticFieldRegion.h"

void MagneticFieldRegion::loadFlukaField(const std::string& filename)
{
  // read field provided by MGNCREAT / MGNDATA fluka cards
  if (mName.empty()) {
    mName = filename;
  }
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file: " + filename);
  }
  std::string line;
  int seenMGNCREAT = 0;
  std::vector<float> MGNCREATcache;
  while (std::getline(file, line)) {
    if (line.find("MGNCREAT") != std::string::npos) {
      seenMGNCREAT = true;
      cacheValues(line, MGNCREATcache);
    } else if (line.find("MGNDATA") != std::string::npos) {
      if (seenMGNCREAT == 1) {
        if (MGNCREATcache.size() >= 13) {
          mXMin = MGNCREATcache[7];
          mYMin = MGNCREATcache[8];
          mZMin = MGNCREATcache[9];
          mXMax = MGNCREATcache[10];
          mYMax = MGNCREATcache[11];
          mZMax = MGNCREATcache[12];
          mNX = static_cast<int>(MGNCREATcache[4]);
          mNY = static_cast<int>(MGNCREATcache[5]);
          mNZ = static_cast<int>(MGNCREATcache[6]);
          mNX1 = mNX - 1;
          mNY1 = mNY - 1;
          mNZ1 = mNZ - 1;
          mDX = (mXMax - mXMin) / mNX1;
          mDY = (mYMax - mYMin) / mNY1;
          mDZ = (mZMax - mZMin) / mNZ1;
          mDXI = 1. / mDX;
          mDYI = 1. / mDY;
          mDZI = 1. / mDZ;
          seenMGNCREAT = 2;
          mFieldData.reserve(3 * mNX * mNY * mNZ);
        } else {
          LOGP(fatal, "MGNCREAT collected only {} values", MGNCREATcache.size());
        }
      }
      cacheValues(line, mFieldData);
    }
  }
  file.close();
  if (mFieldData.size() != size_t(3 * mNX * mNY * mNZ)) {
    LOGP(fatal, "Expected {} field points for {}x{}x{} grid, loaded {}", 3 * mNX * mNY * mNZ, mNX, mNY, mNZ, mFieldData.size());
  }
}

void MagneticFieldRegion::cacheValues(const std::string& line, std::vector<float>& cachev)
{
  std::istringstream iss(line);
  auto vstr = na6p::utils::Str::tokenize(line, ',', true, false);
  for (auto& s : vstr) {
    try {
      cachev.push_back(std::stof(s));
    } catch (...) {
      // ignore
    }
  }
}

void MagneticFieldRegion::interpolateField(int i, int j, int k, float x, float y, float z, double* bxbybz) const
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

void MagneticFieldRegion::interpolateFieldAdd(int i, int j, int k, float x, float y, float z, double* bxbybz) const
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

void MagneticFieldRegion::setRefPosition(float x, float y, float z)
{
  std::string rep;
  mRefPos[0] = x;
  mRefPos[1] = y;
  mRefPos[2] = z;
  mBoxPos[0][0] = mXMin + mRefPos[0];
  mBoxPos[0][1] = mXMax + mRefPos[0];
  rep += fmt::format("{:.1f}<X<{:.1f} ", mBoxPos[0][0], mBoxPos[0][1]);
  mBoxPos[1][0] = mYMin + mRefPos[1];
  mBoxPos[1][1] = mYMax + mRefPos[1];
  rep += fmt::format("{:.1f}<Y<{:.1f} ", mBoxPos[1][0], mBoxPos[1][1]);
  mBoxPos[2][0] = mZMin + mRefPos[2];
  mBoxPos[2][1] = mZMax + mRefPos[2];
  rep += fmt::format("{:.1f}<Z<{:.1f} ", mBoxPos[2][0], mBoxPos[2][1]);
  LOGP(info, "Adding Field {} in {}", mName, rep);
}
