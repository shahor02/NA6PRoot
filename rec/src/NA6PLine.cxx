// NA6PCCopyright

#include "NA6PLine.h"

NA6PLine::NA6PLine(const NA6PLine& source)
{
  for (int i{0}; i < 3; ++i) {
    mOriginPoint[i] = source.mOriginPoint[i];
    mCosinesDirector[i] = source.mCosinesDirector[i];
  }
}

NA6PLine::NA6PLine(const float firstPoint[3], const float secondPoint[3])
{
  for (int index{0}; index < 3; ++index) {
    mOriginPoint[index] = firstPoint[index];
    mCosinesDirector[index] = secondPoint[index] - firstPoint[index];
  }

  float inverseNorm{1.f / std::sqrt(mCosinesDirector[0] * mCosinesDirector[0] + mCosinesDirector[1] * mCosinesDirector[1] +
                                    mCosinesDirector[2] * mCosinesDirector[2])};

  for (int index{0}; index < 3; ++index) {
    mCosinesDirector[index] *= inverseNorm;
  }
}

std::array<float, 6> NA6PLine::getDCAComponents(const float point[3]) const
{
  std::array<float, 6> components{0., 0., 0., 0., 0., 0.};
  float cdelta{0.};
  for (int i{0}; i < 3; ++i) {
    cdelta -= mCosinesDirector[i] * (mOriginPoint[i] - point[i]);
  }

  components[0] = mOriginPoint[0] - point[0] + mCosinesDirector[0] * cdelta;
  components[3] = mOriginPoint[1] - point[1] + mCosinesDirector[1] * cdelta;
  components[5] = mOriginPoint[2] - point[2] + mCosinesDirector[2] * cdelta;
  components[1] = std::sqrt(components[0] * components[0] + components[3] * components[3]);
  components[2] = std::sqrt(components[0] * components[0] + components[5] * components[5]);
  components[4] = std::sqrt(components[3] * components[3] + components[5] * components[5]);

  return components;
}

bool NA6PLine::getClosestPoints(const NA6PLine& line1, const NA6PLine& line2,
                                float p1[3], float p2[3], float precision)
{

  if (areParallel(line1, line2, precision)) {
    p1[0] = line1.mOriginPoint[0];
    p1[1] = line1.mOriginPoint[1];
    p1[2] = line1.mOriginPoint[2];
    p2[0] = line2.mOriginPoint[0];
    p2[1] = line2.mOriginPoint[1];
    p2[2] = line2.mOriginPoint[2];
    return false;
  }

  // scalar products
  float w0[3];
  float uu{0.f}, vv{0.f}, uv{0.f}, w0u{0.f}, w0v{0.f};
  for (int jj = 0; jj < 3; ++jj) {
    w0[jj] = line1.mOriginPoint[jj] - line2.mOriginPoint[jj];
    uu += line1.mCosinesDirector[jj] * line1.mCosinesDirector[jj];
    vv += line2.mCosinesDirector[jj] * line2.mCosinesDirector[jj];
    uv += line1.mCosinesDirector[jj] * line2.mCosinesDirector[jj];
    w0u += w0[jj] * line1.mCosinesDirector[jj];
    w0v += w0[jj] * line2.mCosinesDirector[jj];
  }
  float denom = uu * vv - uv * uv;

  float s = (uv * w0v - vv * w0u) / denom;
  float t = (uu * w0v - uv * w0u) / denom;

  p1[0] = line1.mOriginPoint[0] + s * line1.mCosinesDirector[0];
  p1[1] = line1.mOriginPoint[1] + s * line1.mCosinesDirector[1];
  p1[2] = line1.mOriginPoint[2] + s * line1.mCosinesDirector[2];

  p2[0] = line2.mOriginPoint[0] + t * line2.mCosinesDirector[0];
  p2[1] = line2.mOriginPoint[1] + t * line2.mCosinesDirector[1];
  p2[2] = line2.mOriginPoint[2] + t * line2.mCosinesDirector[2];

  return true;
}

bool NA6PLine::getCrossingPoint(const NA6PLine& line1, const NA6PLine& line2, float p[3], const float precision)
{
  float p1[3], p2[3];
  bool retCode = getClosestPoints(line1, line2, p1, p2, precision);
  if (!retCode)
    return false;
  for (int jj = 0; jj < 3; ++jj)
    p[jj] = 0.5f * (p1[jj] + p2[jj]);
  return true;
}

bool NA6PLine::areParallel(const NA6PLine& firstLine, const NA6PLine& secondLine, const float precision)
{
  float crossProdX{firstLine.mCosinesDirector[1] * secondLine.mCosinesDirector[2] -
                   firstLine.mCosinesDirector[2] * secondLine.mCosinesDirector[1]};
  float module{std::abs(firstLine.mCosinesDirector[1] * secondLine.mCosinesDirector[2]) +
               std::abs(firstLine.mCosinesDirector[2] * secondLine.mCosinesDirector[1])};
  if (std::abs(crossProdX) > precision * module) {
    return false;
  }

  float crossProdY{-firstLine.mCosinesDirector[0] * secondLine.mCosinesDirector[2] +
                   firstLine.mCosinesDirector[2] * secondLine.mCosinesDirector[0]};
  module = std::abs(firstLine.mCosinesDirector[0] * secondLine.mCosinesDirector[2]) +
           std::abs(firstLine.mCosinesDirector[2] * secondLine.mCosinesDirector[0]);
  if (std::abs(crossProdY) > precision * module) {
    return false;
  }

  float crossProdZ = firstLine.mCosinesDirector[0] * secondLine.mCosinesDirector[1] -
                     firstLine.mCosinesDirector[1] * secondLine.mCosinesDirector[0];
  module = std::abs(firstLine.mCosinesDirector[0] * secondLine.mCosinesDirector[1]) +
           std::abs(firstLine.mCosinesDirector[1] * secondLine.mCosinesDirector[0]);
  if (std::abs(crossProdZ) > precision * module) {
    return false;
  }

  return true;
}
