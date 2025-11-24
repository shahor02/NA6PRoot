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

#ifndef NA6P_VERTEL_RECONSTRUCTION_H
#define NA6P_VERTEL_RECONSTRUCTION_H

#include <string>
#include <Rtypes.h>

// Class to steer the VT reconstruction

#include "NA6PReconstruction.h"

class TFile;
class TTree;
class NA6PVerTelHit;
class NA6PBaseCluster;

class NA6PVerTelReconstruction : public NA6PReconstruction
{
 public:

  NA6PVerTelReconstruction() : NA6PReconstruction("VerTel") {}
  ~NA6PVerTelReconstruction() override = default;

  // methods to steer cluster reconstruction
  void createClustersOutput() override;
  void clearClusters() override { mClusters.clear(); }
  void writeClusters() override;
  void closeClustersOutput() override;
  // fast method to smear the hits bypassing digitization and cluster finder
  void setClusterSpaceResolution( double clures ) {mCluRes = clures;}
  void hitsToRecPoints(const std::vector<NA6PVerTelHit>& vtHits);

 private:
  std::vector<NA6PBaseCluster> mClusters, *hClusPtr = &mClusters;  // vector of clusters
  TFile* mClusFile = nullptr;                                      // file with clusters
  TTree* mClusTree = nullptr;                                      // tree of clusters
  double mCluRes = 5.e-4;                                          // cluster resolution, cm (for fast simu)

  ClassDefNV(NA6PVerTelReconstruction, 1);
};

#endif
