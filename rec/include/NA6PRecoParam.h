// NA6PCCopyright

#ifndef NA6P_RECO_PARAM_H_
#define NA6P_RECO_PARAM_H_

#include "ConfigurableParam.h"
#include "ConfigurableParamHelper.h"

struct NA6PRecoParam : public na6p::conf::ConfigurableParamHelper<NA6PRecoParam> {

  static constexpr int MaxIterationsTrackerCA = 10;

  // VerTel reconstruction parameters
  int vtNLayers = 5;
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
  // VerTel CA tracker
  bool vtDoOutwardPropagation = false;
  float vtZOutProp = 40.;
  int vtNIterationsTrackerCA = 2;
  float vtMaxDeltaThetaTrackletsCA[MaxIterationsTrackerCA] = {0.04, 0.1, 0.15, 0.3, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float vtMaxDeltaPhiTrackletsCA[MaxIterationsTrackerCA] = {0.1, 0.2, 0.25, 0.5, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float vtMaxDeltaTanLCellsCA[MaxIterationsTrackerCA] = {4., 9., 18., 40., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float vtMaxDeltaPhiCellsCA[MaxIterationsTrackerCA] = {0.4, 0.6, 1.2, 2.4, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float vtMaxDeltaPxPzCellsCA[MaxIterationsTrackerCA] = {0.02, 0.05, 0.1, 0.2, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float vtMaxDeltaPyPzCellsCA[MaxIterationsTrackerCA] = {2e-3, 7.5e-3, 1.5e-2, 3.e-2, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float vtMaxChi2TrClCellsCA[MaxIterationsTrackerCA] = {100., 500., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float vtMaxChi2ndfCellsCA[MaxIterationsTrackerCA] = {100., 500., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float vtMaxChi2ndfTracksCA[MaxIterationsTrackerCA] = {100., 500., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  int vtMinNClusTracksCA[MaxIterationsTrackerCA] = {5, 3, 0, 0, 0, 0, 0, 0, 0, 0};
  // MuonSpec CA tracker
  int msNLayers = 6;
  int msNIterationsTrackerCA = 2;
  float msMaxDeltaThetaTrackletsCA[MaxIterationsTrackerCA] = {0.06, 0.1, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float msMaxDeltaPhiTrackletsCA[MaxIterationsTrackerCA] = {0.1, 0.6, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float msMaxDeltaTanLCellsCA[MaxIterationsTrackerCA] = {6., 9., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float msMaxDeltaPhiCellsCA[MaxIterationsTrackerCA] = {0.6, 0.8, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float msMaxDeltaPxPzCellsCA[MaxIterationsTrackerCA] = {0.05, 0.08, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float msMaxDeltaPyPzCellsCA[MaxIterationsTrackerCA] = {0.05, 0.08, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float msMaxChi2TrClCellsCA[MaxIterationsTrackerCA] = {5., 100., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float msMaxChi2ndfCellsCA[MaxIterationsTrackerCA] = {5., 100., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  float msMaxChi2ndfTracksCA[MaxIterationsTrackerCA] = {5., 100., 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0, 999.0};
  int msMinNClusTracksCA[MaxIterationsTrackerCA] = {6, 6, 0, 0, 0, 0, 0, 0, 0};
  // Matching parameters
  float zMatching = 40;
  bool isZMatchingSet = false;
  bool mcMatching = true;
  int minVTHits = 5;
  int minMSHits = 6;
  float minTrackP = 2.0; // GeV/c
  float maxChi2Match = std::numeric_limits<float>::max();
  float maxChi2Refit = std::numeric_limits<float>::max();
  float pMatchWindow = 3; // GeV/c

  NA6PParamDef(NA6PRecoParam, "reco");
};

#endif
