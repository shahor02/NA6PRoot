// NA6PCCopyright

#include "NA6PLine.h"
#include "NA6PTrackPar.h"

NA6PLine::NA6PLine(const NA6PTrackPar& t)
{
  mOriginPoint = t.getXYZ();
  auto pxz2p = 1.f / t.getP2Pxz();
  mCosinesDirector = {t.getTx() * pxz2p, t.getTy() * pxz2p, t.getCosPsi() * pxz2p};
}

bool NA6PLine::getClosestPoints(const NA6PLine& line2, std::array<float, 3>& p1, std::array<float, 3>& p2, float precision) const
{
  if (areParallel(line2, precision)) {
    p1 = mOriginPoint;
    p2 = line2.mOriginPoint;
    return false;
  }

  // scalar products
  float uu{0.f}, vv{0.f}, uv{0.f}, w0u{0.f}, w0v{0.f};
  for (int jj = 0; jj < 3; ++jj) {
    auto w0 = mOriginPoint[jj] - line2.mOriginPoint[jj];
    uu += mCosinesDirector[jj] * mCosinesDirector[jj];
    vv += line2.mCosinesDirector[jj] * line2.mCosinesDirector[jj];
    uv += mCosinesDirector[jj] * line2.mCosinesDirector[jj];
    w0u += w0 * mCosinesDirector[jj];
    w0v += w0 * line2.mCosinesDirector[jj];
  }
  float denom = uu * vv - uv * uv;

  float s = (uv * w0v - vv * w0u) / denom;
  float t = (uu * w0v - uv * w0u) / denom;

  p1[0] = mOriginPoint[0] + s * mCosinesDirector[0];
  p1[1] = mOriginPoint[1] + s * mCosinesDirector[1];
  p1[2] = mOriginPoint[2] + s * mCosinesDirector[2];

  p2[0] = line2.mOriginPoint[0] + t * line2.mCosinesDirector[0];
  p2[1] = line2.mOriginPoint[1] + t * line2.mCosinesDirector[1];
  p2[2] = line2.mOriginPoint[2] + t * line2.mCosinesDirector[2];

  return true;
}

bool NA6PLine::getCrossingPoint(const NA6PLine& line2, std::array<float, 3>& p, const float precision) const
{
  std::array<float, 3> p1, p2;
  if (!getClosestPoints(line2, p1, p2, precision)) {
    return false;
  }
  for (int jj = 0; jj < 3; ++jj) {
    p[jj] = 0.5f * (p1[jj] + p2[jj]);
  }
  return true;
}

bool NA6PLine::areParallel(const NA6PLine& secondLine, float precision) const
{
  float crossProdX = mCosinesDirector[1] * secondLine.mCosinesDirector[2] - mCosinesDirector[2] * secondLine.mCosinesDirector[1];
  float modul = std::abs(mCosinesDirector[1] * secondLine.mCosinesDirector[2]) + std::abs(mCosinesDirector[2] * secondLine.mCosinesDirector[1]);
  if (std::abs(crossProdX) > precision * modul) {
    return false;
  }
  float crossProdY = -mCosinesDirector[0] * secondLine.mCosinesDirector[2] + mCosinesDirector[2] * secondLine.mCosinesDirector[0];
  modul = std::abs(mCosinesDirector[0] * secondLine.mCosinesDirector[2]) + std::abs(mCosinesDirector[2] * secondLine.mCosinesDirector[0]);
  if (std::abs(crossProdY) > precision * modul) {
    return false;
  }
  float crossProdZ = mCosinesDirector[0] * secondLine.mCosinesDirector[1] - mCosinesDirector[1] * secondLine.mCosinesDirector[0];
  modul = std::abs(mCosinesDirector[0] * secondLine.mCosinesDirector[1]) + std::abs(mCosinesDirector[1] * secondLine.mCosinesDirector[0]);
  if (std::abs(crossProdZ) > precision * modul) {
    return false;
  }
  return true;
}
