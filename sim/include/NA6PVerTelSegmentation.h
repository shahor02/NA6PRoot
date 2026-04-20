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

#ifndef NA6P_VERTEL_SEGMENTATION_H
#define NA6P_VERTEL_SEGMENTATION_H

#include <Rtypes.h>

// define dead areas and segmentation of the VT modules

class NA6PVerTelSegmentation
{
 public:
  static constexpr float DeadXLong = 0.45;                            // right end cap
  static constexpr float DeadXShort = 0.15;                           // left end cap (readout)
  static constexpr float DeadYBottom = 0.0525;                        // bottom dead zone
  static constexpr float DeadYTop = 0.0525;                           // top dead zone
  static constexpr float DeadTopBotHalves = 0.012;                    // dead space between top and bottom halves
  static constexpr float DeadXTile = 0.002;                           // dead space between tiles in X
  static constexpr float DeadXDataBackbone = 0.006;                   // dead space between every segment (triplet of tiles)
  static constexpr float XSizeTot = 12.9996 + DeadXShort + DeadXLong; // 13.5996;
  static constexpr float YSizeTot = 13.5898 + DeadYBottom + DeadYTop; // 13.6948; // readout side
  static constexpr int NXTiles = 36;                                  // tiles (3 tiles per segment)
  static constexpr int NXSegments = 12;                               // group of 3 tiles
  static constexpr int NYSensors = 7;                                 // number of rows along y
  static constexpr float DXSegment = (XSizeTot - DeadXShort - DeadXLong) / NXSegments;
  static constexpr float DXTile = (DXSegment - DeadXDataBackbone) / 3;
  static constexpr float DYSens = YSizeTot / NYSensors;

  NA6PVerTelSegmentation() = default;

  void setOffsetX(float val) { mOffsX = val; }
  void setOffsetY(float val) { mOffsY = val; }
  void setInterChipGap(float val) { mInterChipGap = val; }
  void setStaggered(bool val);

  int isInAcc(float x, float y) const;

 protected:
  float mOffsX = 0.3;
  float mOffsY = -0.3;
  float mInterChipGap = 0.02;
  bool mStaggered = false;
  float mDeadXLongEff = DeadXLong;
  float mDeadXShortEff = DeadXShort;
  float mDeadYBottomEff = DeadYBottom;
  float mDeadYTopEff = DeadYTop;

  ClassDefNV(NA6PVerTelSegmentation, 1);
};

#endif
