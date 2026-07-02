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
// VTPixID stores rsu/tile/row/col packed into a single uint32_t,
// to ensure correct ROOT streaming.
// Always use the getters/setters to access pixel indices.

struct VTPixID {
  static constexpr int kColBits = 10; // 0-155,  max 1023
  static constexpr int kRowBits = 10; // 0-443,  max 1023
  static constexpr int kTileBits = 4; // 0-11,   max 15
  static constexpr int kRsuBits = 6;  // 0-41,   max 63
  static_assert(kColBits + kRowBits + kTileBits + kRsuBits == 30,
                "VTPixID fields must fit in 32 bits");

  uint32_t mPacked = 0; // the ONLY streamed member

  uint32_t getCol() const { return mPacked & ((1u << kColBits) - 1); }
  uint32_t getRow() const { return (mPacked >> kColBits) & ((1u << kRowBits) - 1); }
  uint32_t getTile() const { return (mPacked >> (kColBits + kRowBits)) & ((1u << kTileBits) - 1); }
  uint32_t getRsu() const { return (mPacked >> (kColBits + kRowBits + kTileBits)) & ((1u << kRsuBits) - 1); }

  void setCol(uint32_t v) { mPacked = (mPacked & ~(((1u << kColBits) - 1))) | (v & ((1u << kColBits) - 1)); }
  void setRow(uint32_t v) { mPacked = (mPacked & ~(((1u << kRowBits) - 1) << kColBits)) | ((v & ((1u << kRowBits) - 1)) << kColBits); }
  void setTile(uint32_t v) { mPacked = (mPacked & ~(((1u << kTileBits) - 1) << (kColBits + kRowBits))) | ((v & ((1u << kTileBits) - 1)) << (kColBits + kRowBits)); }
  void setRsu(uint32_t v) { mPacked = (mPacked & ~(((1u << kRsuBits) - 1) << (kColBits + kRowBits + kTileBits))) | ((v & ((1u << kRsuBits) - 1)) << (kColBits + kRowBits + kTileBits)); }

  ClassDefNV(VTPixID, 2);
};

class NA6PVerTelDigit
{
 public:
  NA6PVerTelDigit() = default;
  NA6PVerTelDigit(uint16_t detID, const VTPixID& id);
  NA6PVerTelDigit(uint16_t detID, uint32_t rsu, uint32_t tile,
                  uint32_t row, uint32_t col);

  uint16_t getDetectorID() const { return mDetectorID; }
  const VTPixID& getPixID() const { return mPixID; }
  uint32_t getRSU() const { return mPixID.getRsu(); }
  uint32_t getTile() const { return mPixID.getTile(); }
  uint32_t getRow() const { return mPixID.getRow(); }
  uint32_t getCol() const { return mPixID.getCol(); }

  void setDetectorID(uint16_t id) { mDetectorID = id; }
  void setPixID(const VTPixID& id) { mPixID = id; }
  void setRSU(uint32_t id) { mPixID.setRsu(id); }
  void setTile(uint32_t id) { mPixID.setTile(id); }
  void setRow(uint32_t id) { mPixID.setRow(id); }
  void setCol(uint32_t id) { mPixID.setCol(id); }

  void print() const;
  std::string asString() const;

 protected:
  uint16_t mDetectorID = 0; // the detector/sensor id
  VTPixID mPixID;           // pixel identifier

  ClassDefNV(NA6PVerTelDigit, 1);
};

#endif
