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
  
  NA6PVerTelSegmentation() = default;

  void setOffsetX(float val) { mOffsX = val;}
  void setOffsetY(float val) { mOffsY = val;}
  void setInterChipGap(float val) { mInterChipGap = val;}
  
  int isInAcc(float x, float y) const;
  
 protected:
  static const float XSizeTot;
  static const float YSizeTot; // readout side
  static const float DeadXLong;  // right end cap
  static const float DeadXShort; // left end cap (readout)
  static const float DeadYBottom; // bottom dead zone
  static const float DeadYTop; // top dead zone
  static const float DeadTopBotHalves; // dead space between top and bottom halves
  static const float DeadXTile;  // dead space between tiles in X
  static const float DeadXDataBackbone; // dead space between every segment (triplet of tiles)
  static const int NXTiles;
  static const int NXSegments;
  static const int NYSensors;
  static const float DXTile;
  static const float DXSegment;
  static const float DYSens;
  
  float mOffsX = 0.3;
  float mOffsY = -0.3;
  float mInterChipGap = 0.02;
  
  ClassDefNV(NA6PVerTelSegmentation, 1);
};

#endif
