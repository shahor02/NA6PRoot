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

#ifndef NA6P_VERTEL_DIGIT_H
#define NA6P_VERTEL_DIGIT_H

#include <Rtypes.h>

// object for VT digits

struct VTPixID {

  static constexpr int kColBits = 10; // 0-161,  max 1023
  static constexpr int kRowBits = 10; // 0-459,  max 1023
  static constexpr int kTileBits = 4; // 0-11,   max 15
  static constexpr int kRsuBits = 6;  // 0-41,   max 63

  static_assert(kColBits + kRowBits + kTileBits + kRsuBits == 30,
                "VTPixID fields must fit in 32 bits");

  uint32_t col : kColBits = 0;
  uint32_t row : kRowBits = 0;
  uint32_t tile : kTileBits = 0;
  uint32_t rsu : kRsuBits = 0;

  uint32_t pack() const
  {
    return uint32_t(col) |
           (uint32_t(row) << kColBits) |
           (uint32_t(tile) << (kColBits + kRowBits)) |
           (uint32_t(rsu) << (kColBits + kRowBits + kTileBits));
  }

  static VTPixID unpack(uint32_t val)
  {
    VTPixID id;
    id.col = val & ((1 << kColBits) - 1);
    id.row = (val >> kColBits) & ((1 << kRowBits) - 1);
    id.tile = (val >> (kColBits + kRowBits)) & ((1 << kTileBits) - 1);
    id.rsu = (val >> (kColBits + kRowBits + kTileBits)) & ((1 << kRsuBits) - 1);
    return id;
  }

  ClassDefNV(VTPixID, 1);
};

class NA6PVerTelDigit
{
 public:
  NA6PVerTelDigit() = default;
  NA6PVerTelDigit(UShort_t detID, UShort_t rsu, UShort_t tile,
                  UShort_t row, UShort_t col, int pid = -1);

  UShort_t getDetectorID() const { return mDetectorID; }
  const VTPixID& getPixID() const { return mPixID; }
  UShort_t getRSU() const { return mPixID.rsu; }
  UShort_t getTile() const { return mPixID.tile; }
  UShort_t getRow() const { return mPixID.row; }
  UShort_t getCol() const { return mPixID.col; }
  int getParticleID() const { return mParticleID; }

  void setDetectorID(int id) { mDetectorID = id; }
  void setPixID(const VTPixID& id) { mPixID = id; }
  void setRSU(int id) { mPixID.rsu = id; }
  void setTile(int id) { mPixID.tile = id; }
  void setRow(int id) { mPixID.row = id; }
  void setCol(int id) { mPixID.col = id; }
  void setParticleID(int id) { mParticleID = id; }

  void print() const;
  std::string asString() const;

 protected:
  UShort_t mDetectorID = 0; // the detector/sensor id
  VTPixID mPixID;           // pixel identifier
  int mParticleID = -1;     // Particle ID

  ClassDefNV(NA6PVerTelDigit, 1);
};

#endif
