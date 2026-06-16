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
// NOTE: VTPixID uses C++ bit fields which ROOT cannot stream correctly.
// Do not access mPixID branches directly via TTree::Scan() or TTree::Draw().
// Always use the getters (getRSU(), getTile(), getRow(), getCol())
// to access pixel indices.

struct VTPixID {

  static constexpr int kColBits = 10; // 0-155,  max 1023
  static constexpr int kRowBits = 10; // 0-443,  max 1023
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
  NA6PVerTelDigit(uint16_t detID, const VTPixID& id, int pid = -1);
  NA6PVerTelDigit(uint16_t detID, uint32_t rsu, uint32_t tile,
                  uint32_t row, uint32_t col, int pid = -1);

  uint16_t getDetectorID() const { return mDetectorID; }
  const VTPixID& getPixID() const { return mPixID; }
  uint32_t getRSU() const { return mPixID.rsu; }
  uint32_t getTile() const { return mPixID.tile; }
  uint32_t getRow() const { return mPixID.row; }
  uint32_t getCol() const { return mPixID.col; }
  int getParticleID() const { return mParticleID; }

  void setDetectorID(uint16_t id) { mDetectorID = id; }
  void setPixID(const VTPixID& id) { mPixID = id; }
  void setRSU(uint32_t id) { mPixID.rsu = id; }
  void setTile(uint32_t id) { mPixID.tile = id; }
  void setRow(uint32_t id) { mPixID.row = id; }
  void setCol(uint32_t id) { mPixID.col = id; }
  void setParticleID(int id) { mParticleID = id; }

  void print() const;
  std::string asString() const;

 protected:
  uint16_t mDetectorID = 0; // the detector/sensor id
  VTPixID mPixID;           // pixel identifier
  int mParticleID = -1;     // Particle ID

  ClassDefNV(NA6PVerTelDigit, 1);
};

#endif
