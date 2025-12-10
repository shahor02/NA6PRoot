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
#include <TVector3.h>

// Class to steer the VT reconstruction

#include "NA6PBaseCluster.h"
#include "NA6PTrack.h"
#include "NA6PVertex.h"
#include "NA6PReconstruction.h"

class TFile;
class TTree;
class NA6PVerTelHit;
class NA6PVertexerTracklets;
class NA6PTrackerCA;

class NA6PVerTelReconstruction : public NA6PReconstruction
{
 public:
  NA6PVerTelReconstruction();
  ~NA6PVerTelReconstruction() override = default;

  bool init(const char* filename, const char* geoname = "NA6P") override;
  // methods to steer cluster reconstruction
  void createClustersOutput() override;
  void clearClusters() override { mClusters.clear(); }
  void writeClusters() override;
  void closeClustersOutput() override;
  // fast method to smear the hits bypassing digitization and cluster finder
  void setClusterSpaceResolution(double clures) { mCluRes = clures; }
  void hitsToRecPoints(const std::vector<NA6PVerTelHit>& vtHits);

  void setPrimaryVertexPosition(double x, double y, double z)
  {
    mPrimaryVertex.SetXYZ(x, y, z);
  }
  // methods to steer tracking
  void setClusters(const std::vector<NA6PBaseCluster>& clusters)
  {
    mClusters = clusters;
    hClusPtr = &mClusters;
  }
  void createVerticesOutput() override;
  void clearVertices() override { mVertices.clear(); }
  void writeVertices() override;
  void closeVerticesOutput() override;
  void runVertexerTracklets();
  const std::vector<NA6PVertex>& getVertices() const { return mVertices; }

  void createTracksOutput() override;
  void clearTracks() override { mTracks.clear(); }
  void writeTracks() override;
  void closeTracksOutput() override;
  void runTracking();

 private:
  std::vector<NA6PBaseCluster> mClusters, *hClusPtr = &mClusters; // vector of clusters
  TFile* mClusFile = nullptr;                                     // file with clusters
  TTree* mClusTree = nullptr;                                     // tree of clusters
  double mCluRes = 5.e-4;                                         // cluster resolution, cm (for fast simu)
  std::vector<NA6PVertex> mVertices, *hVerticesPtr = &mVertices;  // vector of vertices
  TFile* mVertexFile = nullptr;                                   // file with vertices
  TTree* mVertexTree = nullptr;                                   // tree of vertices
  NA6PVertexerTracklets* mVTTrackletVertexer = nullptr;           // vertexer
  TVector3 mPrimaryVertex{0.0, 0.0, 0.0};                         // primary vertex position
  std::vector<NA6PTrack> mTracks, *hTrackPtr = &mTracks;          // vector of tracks
  TFile* mTrackFile = nullptr;                                    // file with tracks
  TTree* mTrackTree = nullptr;                                    // tree of tracks
  NA6PTrackerCA* mVTTracker = nullptr;                            // tracker

  ClassDefNV(NA6PVerTelReconstruction, 1);
};

#endif
