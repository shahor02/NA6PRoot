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

/// \file PID.h
/// \brief particle ids, masses, names class definition
/// \author ruben.shahoyan@cern.ch

#ifndef NA6P_TRACK_PID_H_
#define NA6P_TRACK_PID_H_

#include <cstdint>

#include "Rtypes.h"
#include "NA6PPhysConst.h"

class PID
{
 public:
  // particle identifiers, continuos starting from 0
  typedef uint8_t ID;
  static constexpr ID Electron = 0;
  static constexpr ID Muon = 1;
  static constexpr ID Pion = 2;
  static constexpr ID Kaon = 3;
  static constexpr ID Proton = 4;

  static constexpr ID First = Electron;
  static constexpr ID Last = Proton;   ///< if extra IDs added, update this !!!
  static constexpr ID NIDs = Last + 1; ///< number of defined IDs
  static constexpr const char* sNames[NIDs + 1] = {"Electron", "Muon", "Pion", "Kaon", "Proton", nullptr};
  static constexpr float sMasses[NIDs] = {phys_const::MassElectron, phys_const::MassMuon, phys_const::MassPionCharged, phys_const::MassKaonCharged, phys_const::MassProton};
  static constexpr float sMasses2[NIDs] = {phys_const::MassElectron * phys_const::MassElectron, phys_const::MassMuon* phys_const::MassMuon, phys_const::MassPionCharged* phys_const::MassPionCharged, phys_const::MassKaonCharged* phys_const::MassKaonCharged, phys_const::MassProton* phys_const::MassProton};
  static constexpr float sMassesI[NIDs] = {1.f / phys_const::MassElectron, 1.f / phys_const::MassMuon, 1.f / phys_const::MassPionCharged, 1.f / phys_const::MassKaonCharged, 1.f / phys_const::MassProton};
  static constexpr float sMasses2Z[NIDs] = {phys_const::MassElectron, phys_const::MassMuon, phys_const::MassPionCharged, phys_const::MassKaonCharged, phys_const::MassProton};
  static constexpr float sCharges[NIDs] = {1, 1, 1, 1, 1};

  PID() = default;
  PID(ID id) : mID(id) {}
  PID(const char* name);
  PID(const PID& src) = default;
  PID& operator=(const PID& src) = default;

  ID getID() const { return mID; }
  operator ID() const { return getID(); }

  float getMass() const { return getMass(mID); }
  float getMass2() const { return getMass2(mID); }
  float getMassInv() const { return getMassInv(mID); }
  float getMass2Z() const { return getMass2Z(mID); }
  const char* getName() const { return getName(mID); }
  int getCharge() const { return getCharge(mID); }

  static float getMass(ID id) { return sMasses[id]; }
  static float getMass2(ID id) { return sMasses2[id]; }
  static float getMassInv(ID id) { return sMassesI[id]; }
  static float getMass2Z(ID id) { return sMasses2Z[id]; }
  static int getCharge(ID id) { return sCharges[id]; }
  static const char* getName(ID id) { return sNames[id]; }

 private:
  ID mID = Pion;

  // are 2 strings equal ? (trick from Giulio)
  static constexpr bool sameStr(char const* x, char const* y) { return !*x && !*y ? true : /* default */ (*x == *y && sameStr(x + 1, y + 1)); }
  static constexpr ID nameToID(char const* name, ID id) { return id > Last ? id : sameStr(name, sNames[id]) ? id
                                                                                                            : nameToID(name, id + 1); }

  ClassDefNV(PID, 1);
};

#endif
