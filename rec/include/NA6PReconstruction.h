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

#ifndef NA6P_RECONSTRUCTION_H
#define NA6P_RECONSTRUCTION_H

#include <string>
#include <Rtypes.h>

// Class to steer the reconstruction

class NA6PReconstruction
{
 public:
  NA6PReconstruction(const std::string& name) : mName(name) {}
  virtual ~NA6PReconstruction() {}

  const std::string& getName() const { return mName; }

  virtual bool init(const char* filename, const char* geoname = "NA6P");

  // methods to steer cluster reconstruction
  virtual void createClustersOutput();
  virtual void clearClusters();
  virtual void writeClusters();
  virtual void closeClustersOutput();
  // methods to steer vertexing
  virtual void createVerticesOutput();
  virtual void clearVertices();
  virtual void writeVertices();
  virtual void closeVerticesOutput();
  // methods to steer tracking
  virtual void createTracksOutput();
  virtual void clearTracks();
  virtual void writeTracks();
  virtual void closeTracksOutput();

 protected:
  std::string mName{"MothClass"}; // detector name
  bool mIsInitialized = false;    // flag for initialization
  ClassDefNV(NA6PReconstruction, 1);
};

#endif
