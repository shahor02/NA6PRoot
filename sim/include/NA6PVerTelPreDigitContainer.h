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
#include "NA6PVerTelDigit.h"
#include "NA6PMCComposedLabel.h"

// container of pre-digits per module

struct PreDigit {
  VTPixID pixID;                                             ///< pixel identifier from NA6PVerTelDigit
  float charge = 0.f;                                        ///< collected charg in pixel
  std::vector<std::pair<float, NA6PMCComposedLabel>> labels; ///< <energy, label>, sorted by energy before digitization

  PreDigit(UShort_t rs = 0, UShort_t tl = 0, UShort_t rw = 0, UShort_t cl = 0, float ch = 0.f, NA6PMCComposedLabel lbl = {}) : charge(ch)
  {
    pixID.setRsu(rs);
    pixID.setTile(tl);
    pixID.setRow(rw);
    pixID.setCol(cl);
    labels.emplace_back(ch, lbl);
  }

  void addContribution(float ch, const NA6PMCComposedLabel& lbl)
  {
    charge += ch;
    for (auto& p : labels) {
      if (p.second == lbl) {
        p.first += ch;
        return;
      }
    }
    labels.emplace_back(ch, lbl);
  }

  void sortLabelsByEnergy()
  {
    std::sort(labels.begin(), labels.end(),
              [](const auto& a, const auto& b) {
                return a.first > b.first;
              });
  }

  ClassDefNV(PreDigit, 1);
};

class NA6PVerTelPreDigitContainer
{
 public:
  static constexpr int kRofShift = 8 * sizeof(uint32_t); // shift by VTPixID size

  static_assert(VTPixID::kColBits + VTPixID::kRowBits +
                    VTPixID::kTileBits + VTPixID::kRsuBits <=
                  kRofShift,
                "VTPixID bit fields must fit within kRofShift in the 64-bit key");

  NA6PVerTelPreDigitContainer(UShort_t idx = 0) : mDetectorIndex(idx){};
  ~NA6PVerTelPreDigitContainer() = default;

  std::map<ULong64_t, PreDigit>& getPreDigits() { return mPreDigits; }
  PreDigit* findDigit(ULong64_t key);
  void addDigit(ULong64_t key, UShort_t rs, UShort_t tl, UShort_t rw, UShort_t cl, float ch, const NA6PMCComposedLabel& lbl);

  UShort_t getDetectorIndex() const { return mDetectorIndex; }
  void setDetectorIndex(UShort_t ind) { mDetectorIndex = ind; }

  bool isEmpty() const { return mPreDigits.empty(); }
  static ULong64_t getOrderingKey(UShort_t rsu, UShort_t tile, UShort_t row, UShort_t col)
  {
    VTPixID id;
    id.setCol(col);
    id.setRow(row);
    id.setTile(tile);
    id.setRsu(rsu);
    return static_cast<ULong64_t>(id.mPacked);
  }

  static UShort_t key2col(ULong64_t key) { return key2pixID(key).getCol(); }
  static UShort_t key2row(ULong64_t key) { return key2pixID(key).getRow(); }
  static UShort_t key2tile(ULong64_t key) { return key2pixID(key).getTile(); }
  static UShort_t key2rsu(ULong64_t key) { return key2pixID(key).getRsu(); }

  static void unpackKey(ULong64_t key,
                        UShort_t& rsu, UShort_t& tile,
                        UShort_t& row, UShort_t& col)
  {
    VTPixID id = key2pixID(key);
    rsu = id.getRsu();
    tile = id.getTile();
    row = id.getRow();
    col = id.getCol();
  }

  bool isDisabled() const { return mDisabled; }
  void disable(bool v) { mDisabled = v; }

 protected:
  UShort_t mDetectorIndex = 0;              ///< Detector index
  bool mDisabled = false;                   ///< flag to disable module
  std::map<ULong64_t, PreDigit> mPreDigits; ///< Map of fired pixels

 private:
  static VTPixID key2pixID(ULong64_t key)
  {
    VTPixID id;
    id.mPacked = static_cast<uint32_t>(key & 0xFFFFFFFF);
    return id;
  }

  ClassDefNV(NA6PVerTelPreDigitContainer, 1);
};

//_______________________________________________________________________
inline PreDigit* NA6PVerTelPreDigitContainer::findDigit(ULong64_t key)
{
  // finds the digit corresponding to global key
  auto digitentry = mPreDigits.find(key);
  return digitentry != mPreDigits.end() ? &(digitentry->second) : nullptr;
}
//_______________________________________________________________________
inline void NA6PVerTelPreDigitContainer::addDigit(ULong64_t key, UShort_t rs, UShort_t tl, UShort_t rw, UShort_t cl, float ch, const NA6PMCComposedLabel& lbl)
{
  mPreDigits.emplace(std::make_pair(key, PreDigit(rs, tl, rw, cl, ch, lbl)));
}

#endif
