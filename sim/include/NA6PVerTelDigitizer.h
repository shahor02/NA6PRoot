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

#ifndef NA6P_VERTEL_DIGITIZER_H
#define NA6P_VERTEL_DIGITIZER_H

#include "NA6PVerTelPreDigitContainer.h"
#include <Rtypes.h>
#include <TGeoMatrix.h>

// steers hits -> digits step

class NA6PVerTelDigitizer
{
 public:
  static constexpr int kNModulesPerLayer = 4;

  NA6PVerTelDigitizer() = default;
  ~NA6PVerTelDigitizer() = default;

  static int detID2Layer(int detID) { return detID / kNModulesPerLayer; }

  void init(const char* filename = "geometry.root", const char* geoname = "NA6P");
  void process(const std::vector<NA6PVerTelHit>& hits, int layer = -1);
  void processHit(NA6PVerTelHit hit);

 protected:
  int mNumberOfModules = 0;                          ///< number of modules
  std::vector<NA6PVerTelPreDigitContainer> mModules; ///< Array of module pre-digits containers
  std::vector<TGeoHMatrix> mMatrices{};              ///< local-to-global transforms
  std::vector<float> mModuleHalfX;                   ///< Module half length along x
  std::vector<float> mModuleHalfY;                   ///< Module half length along y
  ClassDefNV(NA6PVerTelDigitizer, 1);
};

#endif
