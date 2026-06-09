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

#ifndef NA6P_VERTEL_PREDIGITCONTAINER_H
#define NA6P_VERTEL_PREDIGITCONTAINER_H

#include <Rtypes.h>

// container of pre-digits per module

struct PreDigit {
  UShort_t rsu = 0;    ///< RSU index in the detector [0, 41]
  UShort_t tile = 0;   ///< Tile index in the RSU [0, 11]
  UShort_t row = 0;    ///< Pixel row in the tile [0, 459]
  UShort_t col = 0;    ///< Pixel column in the tile [0, 161]
  float charge = 0.f;  ///< collected charg in pixel
  int particleID = -1; ///< label of MC particle

  PreDigit(UShort_t rs = 0, UShort_t tl = 0, UShort_t rw = 0, UShort_t cl = 0, float ch = 0.f, int lbl = 0)
    : rsu(rs), tile(tl), row(rw), col(cl), charge(ch), particleID(lbl) {}

  ClassDefNV(PreDigit, 1);
};

class NA6PVerTelPreDigitContainer
{
 public:
  static constexpr int kColBits = 10; // 0-161,  max 1023
  static constexpr int kRowBits = 10; // 0-459,  max 1023
  static constexpr int kTileBits = 4; // 0-11,   max 15
  static constexpr int kRsuBits = 6;  // 0-41,   max 63
  static constexpr int kRofBits = 34; // reserved for future use

  static constexpr ULong64_t kColMask = (1ULL << kColBits) - 1;
  static constexpr ULong64_t kRowMask = (1ULL << kRowBits) - 1;
  static constexpr ULong64_t kTileMask = (1ULL << kTileBits) - 1;
  static constexpr ULong64_t kRsuMask = (1ULL << kRsuBits) - 1;
  static constexpr ULong64_t kRofMask = (1ULL << kRofBits) - 1;

  static constexpr int kRowShift = kColBits;
  static constexpr int kTileShift = kRowShift + kRowBits;
  static constexpr int kRsuShift = kTileShift + kTileBits;
  static constexpr int kRofShift = kRsuShift + kRsuBits;

  static_assert(kColBits + kRowBits + kTileBits + kRsuBits + kRofBits == 64,
                "Key bit fields must sum to exactly 64 bits");

  NA6PVerTelPreDigitContainer(UShort_t idx = 0) : mDetectorIndex(idx){};
  ~NA6PVerTelPreDigitContainer() = default;

  std::map<ULong64_t, PreDigit>& getPreDigits() { return mPreDigits; }
  UShort_t getDetectorIndex() const { return mDetectorIndex; }

  void setDetectorIndex(UShort_t ind) { mDetectorIndex = ind; }

  bool isEmpty() const { return mPreDigits.empty(); }
  static ULong64_t getOrderingKey(UShort_t rsu, UShort_t tile, UShort_t row, UShort_t col)
  {
    return (static_cast<ULong64_t>(rsu) << kRsuShift) |
           (static_cast<ULong64_t>(tile) << kTileShift) |
           (static_cast<ULong64_t>(row) << kRowShift) |
           static_cast<ULong64_t>(col);
  }
  static UShort_t key2col(ULong64_t key) { return static_cast<UShort_t>(key & kColMask); }
  static UShort_t key2row(ULong64_t key) { return static_cast<UShort_t>((key >> kRowShift) & kRowMask); }
  static UShort_t key2tile(ULong64_t key) { return static_cast<UShort_t>((key >> kTileShift) & kTileMask); }
  static UShort_t key2rsu(ULong64_t key) { return static_cast<UShort_t>((key >> kRsuShift) & kRsuMask); }
  static void unpackKey(ULong64_t key,
                        UShort_t& rsu, UShort_t& tile,
                        UShort_t& row, UShort_t& col)
  {
    col = key2col(key);
    row = key2row(key);
    tile = key2tile(key);
    rsu = key2rsu(key);
  }

  bool isDisabled() const { return mDisabled; }
  void disable(bool v) { mDisabled = v; }

 protected:
  UShort_t mDetectorIndex = 0;              ///< Detector index
  bool mDisabled = false;                   ///< flag to disable module
  std::map<ULong64_t, PreDigit> mPreDigits; ///< Map of fired pixels

  ClassDefNV(NA6PVerTelPreDigitContainer, 1);
};

#endif
