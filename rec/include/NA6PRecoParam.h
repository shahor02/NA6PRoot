// NA6PCCopyright

#ifndef NA6P_RECO_PARAM_H_
#define NA6P_RECO_PARAM_H_

#include "ConfigurableParam.h"
#include "ConfigurableParamHelper.h"

struct NA6PRecoParam : public na6p::conf::ConfigurableParamHelper<NA6PRecoParam> {

  static constexpr int MaxIterationsTrackerCA = 10;

  // VerTel reconstruction parameters
  int nLayers = 5;
  // tracklet vertexer
  int vertexerLayerToStart = 0;
  float vertexerMaxDeltaThetaTracklet = 0.6;
  float vertexerMaxDeltaPhiTracklet = 0.05;
  float vertexerMaxDeltaTanLamInOut = 1.;
  float vertexerMaxDeltaPhiInOut = 0.2;
  float vertexerMaxDeltaPxPzInOut = 99.;
  float vertexerMaxDeltaPyPzInOut = 0.005;
  std::string vertexerRecoType = "YZ";
  float vertexerMaxDCAxy = 0.25;
  std::string vertexerPeakMethod = "KDE";
  std::string vertexerWeightedMeanOption = "NoWeight";
  float vertexerZMin = -20.0;
  float vertexerZMax = 5.;
  float vertexerZWindowWidth = 1.25;
  int vertexerNBinsForPeakFind = 250;
  int vertexerPeakWidthBins = 3;
  int vertexerMinCountsInPeak = 3;
  std::string vertexerKDEOption = "Standard";
  int vertexerNGridKDE = 500;
  float vertexerKDEBandwidth = 0.5;
  float vertexerMaxPairDCA = 0.1;
  float vertexerMaxPairVertRadius = 0.5;
  float vertexerMinCandidateDistanceZ = 0.1;
  float vertexerMinCandidateDistance3D = 0.2;
  std::string vertexerMultiVertexMode = "Iterative";
  bool vertexerAllowSingleConstribClusters = false;
  // CA tracker
  int nIterationsTrackerCA = 2;
  float maxDeltaThetaTrackletsCA[MaxIterationsTrackerCA] = {0.04, 0.1, 0.15, 0.3, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float maxDeltaPhiTrackletsCA[MaxIterationsTrackerCA] = {0.1, 0.2, 0.25, 0.5, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float maxDeltaTanLCellsCA[MaxIterationsTrackerCA] = {4., 9., 18., 40., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float maxDeltaPhiCellsCA[MaxIterationsTrackerCA] = {0.4, 0.6, 1.2, 2.4, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float maxDeltaPxPzCellsCA[MaxIterationsTrackerCA] = {0.02, 0.05, 0.1, 0.2, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float maxDeltaPyPzCellsCA[MaxIterationsTrackerCA] = {2e-3, 7.5e-3, 1.5e-2, 3.e-2, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float maxChi2TrClCellsCA[MaxIterationsTrackerCA] = {100., 500., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float maxChi2ndfCellsCA[MaxIterationsTrackerCA] = {100., 500., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float maxChi2ndfTracksCA[MaxIterationsTrackerCA] = {100., 500., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  int minNClusTracksCA[MaxIterationsTrackerCA] = {5, 3, 0, 0, 0, 0, 0, 0, 0, 0};

  NA6PParamDef(NA6PRecoParam, "reco");
};

#endif
