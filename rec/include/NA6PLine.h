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
#include <Rtypes.h>

// class for operations with straight lines (tracklets)

class NA6PLine
{
 public:
  NA6PLine() = default;
  NA6PLine(const NA6PLine& source);
  NA6PLine(const float firstPoint[3], const float secondPoint[3]);

  static float getDistanceFromPoint(const NA6PLine& line, const float point[3]);
  static float getDistanceFromPoint(const NA6PLine& line, const std::array<float, 3> point)
  {
    return getDistanceFromPoint(line, point.data());
  }
  float getDistanceFromPoint(const float point[3]) const
  {
    return getDistanceFromPoint(*this, point);
  }
  std::array<float, 6> getDCAComponents(const float point[3]) const;
  std::array<float, 6> getDCAComponents(const std::array<float, 3> point) const
  {
    return getDCAComponents(point.data());
  }
  static float getDCA(const NA6PLine&, const NA6PLine&, const float precision = 1e-7f);
  float getDCA(const NA6PLine& other, float precision = 1e-7f) const
  {
    return getDCA(*this, other, precision);
  }
  static bool getClosestPoints(const NA6PLine& line1, const NA6PLine& line2,
                               float p1[3], float p2[3], float precision = 1e-7f);
  static bool getCrossingPoint(const NA6PLine& line1, const NA6PLine& line2, float p[3], const float precision = 1e-7f);
  static bool areParallel(const NA6PLine&, const NA6PLine&, const float precision = 1e-7f);
  bool isSameLine(const NA6PLine& other, float precision) const
  {
    if (!areParallel(*this, other, precision))
      return false;
    return getDistanceFromPoint(*this, other.mOriginPoint) <= precision;
  }
  bool isEmpty() const { return (mOriginPoint[0] == 0.f && mOriginPoint[1] == 0.f && mOriginPoint[2] == 0.f) &&
                                (mCosinesDirector[0] == 0.f && mCosinesDirector[1] == 0.f && mCosinesDirector[2] == 0.f); }

  float mOriginPoint[3] = {0.0};
  float mCosinesDirector[3] = {0.0};

  ClassDefNV(NA6PLine, 1);
};

// static functions:

inline float NA6PLine::getDistanceFromPoint(const NA6PLine& line, const float point[3])
{
  const float dx = point[0] - line.mOriginPoint[0];
  const float dy = point[1] - line.mOriginPoint[1];
  const float dz = point[2] - line.mOriginPoint[2];
  const float d = (dx * line.mCosinesDirector[0]) + (dy * line.mCosinesDirector[1]) + (dz * line.mCosinesDirector[2]);

  const float vx = dx - (d * line.mCosinesDirector[0]);
  const float vy = dy - (d * line.mCosinesDirector[1]);
  const float vz = dz - (d * line.mCosinesDirector[2]);

  return std::sqrt(vx * vx + vy * vy + vz * vz);
}

inline float NA6PLine::getDCA(const NA6PLine& firstLine, const NA6PLine& secondLine, const float precision)
{
  const float nx = (firstLine.mCosinesDirector[1] * secondLine.mCosinesDirector[2]) -
                   (firstLine.mCosinesDirector[2] * secondLine.mCosinesDirector[1]);
  const float ny = -(firstLine.mCosinesDirector[0] * secondLine.mCosinesDirector[2]) +
                   (firstLine.mCosinesDirector[2] * secondLine.mCosinesDirector[0]);
  const float nz = (firstLine.mCosinesDirector[0] * secondLine.mCosinesDirector[1]) -
                   (firstLine.mCosinesDirector[1] * secondLine.mCosinesDirector[0]);
  const float norm2 = (nx * nx) + (ny * ny) + (nz * nz);

  if (norm2 <= precision * precision) {
    return getDistanceFromPoint(firstLine, secondLine.mOriginPoint);
  }

  const float dx = secondLine.mOriginPoint[0] - firstLine.mOriginPoint[0];
  const float dy = secondLine.mOriginPoint[1] - firstLine.mOriginPoint[1];
  const float dz = secondLine.mOriginPoint[2] - firstLine.mOriginPoint[2];
  const float triple = (dx * nx) + (dy * ny) + (dz * nz);

  return std::abs(triple) / std::sqrt(norm2);
}

#endif // NA6P_LINE_H
