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

#ifndef NA6P_TRACKER_CA_H
#define NA6P_TRACKER_CA_H

#include <string>
#include <Rtypes.h>
#include "NA6PTrack.h"
#include <TVector3.h>

// Cellular Automaton track finder

class NA6PFastTrackFitter;

// structures for temporary objects used in track finding
struct TrackletCandidate {
  int startingLayer;
  int firstClusterIndex;
  int secondClusterIndex;
  double tanL;
  double phi;
  double pxpz;
  double pypz;
};

struct CellCandidate {
  int startingLayer;
  int firstTrackletIndex;
  int secondTrackletIndex;
  std::array<int,3> cluIDs;
  NA6PTrack trackFitFast;
  //    genfit::Track trackFit;
};

struct TrackCandidate {
  int innerLayer;
  int innerCellIndex;
  int outerLayer;
  int outerCellIndex;
  std::vector<int> cluIDs;
};

struct TrackFitted{
  int innerLayer;
  int outerLayer;
  int nClus;
  std::vector<int> cluIDs;
  NA6PTrack trackFitFast;
  //    genfit::Track trackFit;
  double chi2ndf;
};

enum class ExtendDirection { kInward, kOutward };

class NA6PTrackerCA
{
public:
  static constexpr int maxIterationsCA = 10;

  NA6PTrackerCA();
  ~NA6PTrackerCA();

  // setters for configurable parameters
  void setNLayers(int n) {mNLayers = n;}
  void SetMaxNumberOfSharedClusters(int n) {mMaxSharedClusters = n;
  }
  bool loadGeometry(const char* filename, const char* geoname = "NA6P");
  void findTracks(std::vector<NA6PBaseCluster>& cluArr, TVector3 primVert);

protected:
  // methods used in tracking
  void sortClustersByLayerAndEta(std::vector<NA6PBaseCluster>& cluArr,
				 std::vector<int>& firstIndex,
				 std::vector<int>& lastIndex);
  void sortTrackletsByLayerAndIndex(std::vector<TrackletCandidate>& tracklets,
				    std::vector<int>& firstIndex,
				    std::vector<int>& lastIndex);
  void sortCellsByLayerAndIndex(std::vector<CellCandidate>& cells,
				std::vector<int>& firstIndex,
				std::vector<int>& lastIndex);
  void computeLayerTracklets(const std::vector<NA6PBaseCluster>& cluArr,
			     const std::vector<int>& firstIndex,
			     const std::vector<int>& lastIndex,
			     std::vector<TrackletCandidate>& tracklets,
			     double deltaThetaMax,
			     double deltaPhiMax);
  void computeLayerCells(const std::vector<TrackletCandidate>& tracklets,
			 const std::vector<int>& firstIndex,
			 const std::vector<int>& lastIndex,
			 const std::vector<NA6PBaseCluster>& cluArr,
			 std::vector<CellCandidate>& cells,
			 double deltaTanLMax,
			 double deltaPhiMax,
			 double deltaPxPzMax,
			 double deltaPyPzMax,
			 double maxChi2TrClu,
			 double maxChi2NDF);
  double computeTrackToClusterChi2(const NA6PTrack& track,
				   const NA6PBaseCluster& clu);
  bool fitTrackPointsFast(const std::vector<int>& cluIDs,
			  const std::vector<NA6PBaseCluster>& cluArr,
			  NA6PTrack& fitTrack,
			  double maxChi2TrClu,
			  double maxChi2NDF);
  void findCellsNeighbours(const std::vector<CellCandidate>& cells,
			   const std::vector<int>& firstIndex,
			   const std::vector<int>& lastIndex,
			   std::vector<std::pair<int,int>>& cneigh,
			   const std::vector<NA6PBaseCluster>& cluArr,
			   double maxChi2TrClu);
  std::vector<TrackCandidate> prolongSeed(const TrackCandidate& seed,
					  const std::vector<CellCandidate>& cells,
					  const std::vector<int>& firstIndex,
					  const std::vector<int>& lastIndex,
					  const std::vector<NA6PBaseCluster>& cluArr,
					  double maxChi2TrClu,
					  ExtendDirection dir);
  void findRoads(const std::vector<std::pair<int,int>>& cneigh,
		 const std::vector<CellCandidate>& cells,
		 const std::vector<int>& firstIndex,
		 const std::vector<int>& lastIndex,
		 const std::vector<TrackletCandidate>& tracklets,
		 const std::vector<NA6PBaseCluster>& cluArr,
		 std::vector<TrackCandidate>& trackCands,
		 double maxChi2TrClu);
  void fitAndSelectTracks(const std::vector<TrackCandidate>& trackCands,
			  const std::vector<NA6PBaseCluster>& cluArr,
			  std::vector<TrackFitted>& tracks,
			  double maxChi2TrClu,
			  int minNClu,
			  double maxChi2NDF);
  template<typename T>
  void printStats(const std::vector<T>& candidates,
		  const std::vector<NA6PBaseCluster>& cluArr,
		  const std::vector<CellCandidate>& cells,
		  const std::string& label,
		  int requiredClus = -1);




private:
  
  int    mNLayers = 5;
  double mPrimVertPos[3] = {};
  NA6PFastTrackFitter* mTrackFitter = nullptr;
  std::vector<bool> mIsClusterUsed = {};
  std::vector<TrackFitted> mFinalTracks = {};
  int    mMaxSharedClusters = 0;
  int    mNIterationsCA = 2;
  double maxDeltaThetaTrackletsCA[maxIterationsCA] = {0.04, 0.1, 0.15, 0.3, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  double maxDeltaPhiTrackletsCA[maxIterationsCA] = {0.1, 0.2, 0.25, 0.5, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  double maxDeltaTanLCellsCA[maxIterationsCA] = {4., 9., 18., 40., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  double maxDeltaPhiCellsCA[maxIterationsCA] = {0.4, 0.6, 1.2, 2.4, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  double maxDeltaPxPzCellsCA[maxIterationsCA] = {0.02, 0.05, 0.1, 0.2, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  double maxDeltaPyPzCellsCA[maxIterationsCA] = {2e-3, 7.5e-3, 1.5e-2, 3.e-2, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  double maxChi2TrClCellsCA[maxIterationsCA] = {5., 10., 10., 10.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  double maxChi2ndfCellsCA[maxIterationsCA] = {5., 10., 10., 10., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  double maxChi2ndfTracksCA[maxIterationsCA] = {5., 10., 10., 10., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  int    minNClusTracksCA[maxIterationsCA] = {5, 3, 0, 0, 0, 0, 0, 0, 0, 0};



 
  

  ClassDefNV(NA6PTrackerCA, 1);
};

#endif
