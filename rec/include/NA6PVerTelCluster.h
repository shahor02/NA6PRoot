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

#ifndef NA6P_VERTEL_CLUSTER_H
#define NA6P_VERTEL_CLUSTER_H

#include "NA6PBaseCluster.h"

// Vertex Telescope cluster class

class NA6PVerTelCluster : public NA6PBaseCluster
{
 public:

  NA6PVerTelCluster() = default;
  NA6PVerTelCluster(float x, float y, float z, int clusiz, int nDet);
  NA6PVerTelCluster(const NA6PVerTelCluster&) = default;
  NA6PVerTelCluster& operator=(const NA6PVerTelCluster&) = default;

  ClassDefNV(NA6PVerTelCluster, 1);
};

#endif
