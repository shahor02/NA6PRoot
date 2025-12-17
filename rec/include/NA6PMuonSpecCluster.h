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

#ifndef NA6P_MUONSPEC_CLUSTER_H
#define NA6P_MUONSPEC_CLUSTER_H

#include "NA6PBaseCluster.h"

// Muon Spectrometer cluster class

class NA6PMuonSpecCluster : public NA6PBaseCluster
{
 public:

  NA6PMuonSpecCluster() = default;
  NA6PMuonSpecCluster(float x, float y, float z, int clusiz);
  NA6PMuonSpecCluster(const NA6PMuonSpecCluster&) = default;
  NA6PMuonSpecCluster& operator=(const NA6PMuonSpecCluster&) = default;
  virtual ~NA6PMuonSpecCluster() {}
  
  ClassDefNV(NA6PMuonSpecCluster, 1);
};

#endif
