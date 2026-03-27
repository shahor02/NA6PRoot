// NA6PCCopyright

// Based on:
// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef NA6P_LINE_H
#define NA6P_LINE_H

#include <cmath>
#include <array>
#include <numeric>
#include <Rtypes.h>

class NA6PTrackPar; // fwd declaration only

// class for operations with straight lines (tracklets)

struct NA6PLine {
  NA6PLine() = default;
  NA6PLine(const NA6PLine& source) = default;
  NA6PLine& operator=(const NA6PLine& source) = default;
  NA6PLine(const NA6PTrackPar& t);

  template <typename T = float>
  NA6PLine(const std::array<T, 3>& start, const std::array<T, 3>& end);

  template <typename T = float>
  static NA6PLine fromTwoPoints(const std::array<T, 3>& p0, const std::array<T, 3>& p1)
  {
    return NA6PLine(p0, p1);
  }

  template <typename T = float>
  static NA6PLine fromPointAndDirection(const std::array<T, 3>& xyz, const std::array<T, 3>& dir);

  template <typename T1, typename T2>
  static std::array<T1, 3> getDiff(const std::array<T1, 3>& start, const std::array<T2, 3>& end)
  {
    return {end[0] - start[0], end[1] - start[1], end[2] - start[2]};
  }

  template <typename T>
  float getDistanceFromPoint(const std::array<T, 3>& point) const;

  template <typename T>
  std::array<float, 6> getDCAComponents(const std::array<T, 3>& point) const;

  float getDCA(const NA6PLine&, float precision = 1e-7f) const;

  bool getClosestPoints(const NA6PLine& line2, std::array<float, 3>& p1, std::array<float, 3>& p2, float precision = 1e-7f) const;

  bool getCrossingPoint(const NA6PLine& line2, std::array<float, 3>& p, float precision = 1e-7f) const;

  bool areParallel(const NA6PLine&, float precision = 1e-7f) const;

  bool isSameLine(const NA6PLine& other, float precision) const { return !areParallel(other, precision) ? false : getDistanceFromPoint(other.mOriginPoint) <= precision; }

  bool isEmpty() const { return (mOriginPoint[0] == 0.f && mOriginPoint[1] == 0.f && mOriginPoint[2] == 0.f) && (mCosinesDirector[0] == 0.f && mCosinesDirector[1] == 0.f && mCosinesDirector[2] == 0.f); }

  std::array<float, 3> mOriginPoint{};
  std::array<float, 3> mCosinesDirector{};

  ClassDefNV(NA6PLine, 1);
};

template <typename T>
inline NA6PLine::NA6PLine(const std::array<T, 3>& start, const std::array<T, 3>& end)
{
  *this = fromPointAndDirection(start, getDiff(start, end));
}

template <typename T>
inline NA6PLine NA6PLine::fromPointAndDirection(const std::array<T, 3>& xyz, const std::array<T, 3>& dir)
{
  T norm = std::hypot(dir[0], dir[1], dir[2]), inv = (norm > T(1e-12)) ? T(1) / norm : T(0);
  return {{xyz[0], xyz[1], xyz[2]}, {
                                      dir[0] * inv,
                                      dir[1] * inv,
                                      dir[2] * inv,
                                    }};
}

template <typename T>
inline float NA6PLine::getDistanceFromPoint(const std::array<T, 3>& point) const
{
  const auto dif = getDiff(mOriginPoint, point);
  const float t = std::inner_product(dif.begin(), dif.end(), mCosinesDirector.begin(), 0.f);
  const auto dxyz = getDiff(dif, std::array<float, 3>{t * mCosinesDirector[0], t * mCosinesDirector[1], t * mCosinesDirector[2]});
  return std::hypot(dxyz[0], dxyz[1], dxyz[2]);
}

template <typename T>
std::array<float, 6> NA6PLine::getDCAComponents(const std::array<T, 3>& point) const
{
  const auto dif = getDiff(mOriginPoint, point);
  const auto t = std::inner_product(dif.begin(), dif.end(), mCosinesDirector.begin(), 0.f);

  std::array<float, 6> components{dif[0] - t * mCosinesDirector[0], 0., 0.,
                                  dif[1] - t * mCosinesDirector[1], 0.,
                                  dif[2] - t * mCosinesDirector[2]};
  const float cc0 = components[0] * components[0], cc3 = components[3] * components[3], cc5 = components[5] * components[5];
  components[1] = std::sqrt(cc0 + cc3);
  components[2] = std::sqrt(cc0 + cc5);
  components[4] = std::sqrt(cc3 + cc5);
  return components;
}

inline float NA6PLine::getDCA(const NA6PLine& secondLine, const float precision) const
{
  std::array<float, 3> nxyz{mCosinesDirector[1] * secondLine.mCosinesDirector[2] - mCosinesDirector[2] * secondLine.mCosinesDirector[1],
                            -mCosinesDirector[0] * secondLine.mCosinesDirector[2] + mCosinesDirector[2] * secondLine.mCosinesDirector[0],
                            mCosinesDirector[0] * secondLine.mCosinesDirector[1] - mCosinesDirector[1] * secondLine.mCosinesDirector[0]};
  const float norm2 = std::inner_product(nxyz.begin(), nxyz.end(), nxyz.begin(), 0.f);
  if (norm2 <= precision * precision) {
    return getDistanceFromPoint(secondLine.mOriginPoint);
  }
  const auto dif = getDiff(secondLine.mOriginPoint, mOriginPoint);
  return std::abs(std::inner_product(dif.begin(), dif.end(), nxyz.begin(), 0.f)) / std::sqrt(norm2);
}

#endif // NA6P_LINE_H
