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

#ifndef NA6P_GEOMETRY_MANAGER_H
#define NA6P_GEOMETRY_MANAGER_H

#include <Rtypes.h>
#include <TGeoMatrix.h>

class TFile;
class TGeoVolume;

// helper class to interface to the geometry (martices, sizes) of sensors

class NA6PGeometryManager
{
 public:
  static constexpr int kNVTModulesPerLayer = 4;

  NA6PGeometryManager() = default;
  ~NA6PGeometryManager() = default;

  bool loadGeometry(const char* filename = "geometry.root", const char* geoname = "NA6P");
  bool isGeometryLoaded() const { return mGeoLoaded; }

  const TGeoHMatrix& getMatrix(int jMod) const { return mMatrices[jMod]; }
  float getModuleHalfX(int jMod) const { return mModuleHalfX[jMod]; }
  float getModuleHalfY(int jMod) const { return mModuleHalfY[jMod]; }
  float getModuleFullX(int jMod) const { return mModuleHalfX[jMod] * 2.0f; }
  float getModuleFullY(int jMod) const { return mModuleHalfY[jMod] * 2.0f; }

 private:
  bool fillModuleSize(int jMod, TGeoVolume* vol);

 protected:
  std::vector<TGeoHMatrix> mMatrices{}; ///< local-to-global transforms
  std::vector<float> mModuleHalfX;      ///< Module half length along x
  std::vector<float> mModuleHalfY;      ///< Module half length along y
  bool mGeoLoaded{false};               ///< flag for successful load of geo
  ClassDefNV(NA6PGeometryManager, 1);
};

#endif
