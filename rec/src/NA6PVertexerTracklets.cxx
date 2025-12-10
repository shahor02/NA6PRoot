// NA6PCCopyright

#include <fmt/format.h>
#include <fairlogger/Logger.h>
#include "NA6PBaseCluster.h"
#include "NA6PVertex.h"
#include "NA6PVertexerTracklets.h"

ClassImp(NA6PVertexerTracklets)

  NA6PVertexerTracklets::NA6PVertexerTracklets()
{
  configurePeakFinding(mZMin, mZMax, mNBinsForPeakFind);
}

//______________________________________________________________________

void NA6PVertexerTracklets::sortClustersByLayerAndEta(std::vector<NA6PBaseCluster>& cluArr,
                                                      std::vector<int>& firstIndex,
                                                      std::vector<int>& lastIndex)
{
  // count hits per layer
  std::vector<int> count(mNLayersVT, 0);
  for (const auto& clu : cluArr) {
    int jDet = clu.getDetectorID();
    int jLay = jDet / 4;
    if (jLay >= 0 && jLay < mNLayersVT) {
      count[jLay]++;
    }
  }
  // starting offset for each layer
  firstIndex.resize(mNLayersVT);
  lastIndex.resize(mNLayersVT);
  firstIndex[0] = 0;
  lastIndex[0] = count[0];
  for (int jLay = 1; jLay < mNLayersVT; jLay++) {
    firstIndex[jLay] = firstIndex[jLay - 1] + count[jLay - 1];
    lastIndex[jLay] = firstIndex[jLay] + count[jLay];
  }
  // Create a reordered vector based on layer grouping
  std::vector<int> countReord = firstIndex;
  std::vector<NA6PBaseCluster> reordered(cluArr.size());
  for (const auto& clu : cluArr) {
    int jLay = clu.getDetectorID() / 4;
    if (jLay >= 0 && jLay < mNLayersVT)
      reordered[countReord[jLay]++] = clu;
  }
  cluArr = std::move(reordered);
  // sort by theta within each layer (use z^2/r^2 as a proxy of theta to avoid sqrt and atan)
  for (int jLay = 0; jLay < mNLayersVT; jLay++) {
    auto first = cluArr.begin() + firstIndex[jLay];
    auto last = cluArr.begin() + lastIndex[jLay];
    std::sort(first, last, [](const NA6PBaseCluster& a, const NA6PBaseCluster& b) {
      double xa = a.getX();
      double ya = a.getY();
      double za = a.getZ();
      double xb = b.getX();
      double yb = b.getY();
      double zb = b.getZ();
      double r2a = xa * xa + ya * ya;
      double r2b = xb * xb + yb * yb;
      return za * za * r2b < zb * zb * r2a;
    });
  }
}

//______________________________________________________________________

void NA6PVertexerTracklets::sortTrackletsByLayerAndIndex(std::vector<TrackletForVertex>& tracklets,
                                                         std::vector<int>& firstIndex,
                                                         std::vector<int>& lastIndex)
{
  // count clusters per layer
  std::vector<int> count(mNLayersVT - 1, 0);
  for (const auto& trkl : tracklets) {
    int jLay = trkl.startingLayer;
    if (jLay >= 0 && jLay < mNLayersVT - 1) {
      count[jLay]++;
    }
  }
  // starting offset for each layer
  firstIndex.resize(mNLayersVT - 1);
  lastIndex.resize(mNLayersVT - 1);
  firstIndex[0] = 0;
  lastIndex[0] = count[0];
  for (int jLay = 1; jLay < mNLayersVT - 1; jLay++) {
    firstIndex[jLay] = firstIndex[jLay - 1] + count[jLay - 1];
    lastIndex[jLay] = firstIndex[jLay] + count[jLay];
  }
  // Create a reordered vector based on layer grouping
  std::vector<int> countReord = firstIndex;
  std::vector<TrackletForVertex> reordered(tracklets.size());
  for (const auto& trkl : tracklets) {
    int jLay = trkl.startingLayer;
    if (jLay >= 0 && jLay < mNLayersVT - 1)
      reordered[countReord[jLay]++] = trkl;
  }
  tracklets = std::move(reordered);
  // sort by cluster index within each layer
  for (int jLay = 0; jLay < mNLayersVT - 1; jLay++) {
    auto first = tracklets.begin() + firstIndex[jLay];
    auto last = tracklets.begin() + lastIndex[jLay];
    std::sort(first, last, [](const TrackletForVertex& a, const TrackletForVertex& b) {
      return a.firstClusterIndex < b.firstClusterIndex || (a.firstClusterIndex == b.firstClusterIndex && a.secondClusterIndex < b.secondClusterIndex);
    });
  }
}

//______________________________________________________________________

void NA6PVertexerTracklets::computeLayerTracklets(const std::vector<NA6PBaseCluster>& cluArr,
                                                  const std::vector<int>& firstIndex,
                                                  const std::vector<int>& lastIndex,
                                                  std::vector<TrackletForVertex>& tracklets)
{

  tracklets.clear();

  for (int iLayer = mLayerToStart; iLayer < mLayerToStart + 2; ++iLayer) {
    auto layerBegin = cluArr.begin() + firstIndex[iLayer + 1];
    auto layerEnd = cluArr.begin() + lastIndex[iLayer + 1];
    for (int jClu1 = firstIndex[iLayer]; jClu1 < lastIndex[iLayer]; ++jClu1) {
      const NA6PBaseCluster& clu1 = cluArr[jClu1];
      double x1 = clu1.getX();
      double y1 = clu1.getY();
      double z1 = clu1.getZ();
      double r1 = std::sqrt(x1 * x1 + y1 * y1);
      double theta1 = std::atan2(z1, r1);
      double phi1 = std::atan2(y1, x1);
      double tanth2Min = std::max(0., std::tan(theta1 - 1.2 * mMaxDeltaThetaTracklet)); // 1.2 is a safety margin
      double tanth2Max = std::tan(std::min(M_PI / 2.001, theta1 + 1.2 * mMaxDeltaThetaTracklet));
      auto lower = std::partition_point(layerBegin, layerEnd,
                                        [&](const NA6PBaseCluster& clu) {
                                          double x = clu.getX();
                                          double y = clu.getY();
                                          double z = clu.getZ();
                                          double r2 = x * x + y * y;
                                          double tan2 = z * z / r2;
                                          return tan2 < tanth2Min * tanth2Min;
                                        });

      auto upper = std::partition_point(layerBegin, layerEnd,
                                        [&](const NA6PBaseCluster& clu) {
                                          double x = clu.getX();
                                          double y = clu.getY();
                                          double z = clu.getZ();
                                          double r2 = x * x + y * y;
                                          double tan2 = z * z / r2;
                                          return tan2 <= tanth2Max * tanth2Max;
                                        });
      int lowerIdx = std::distance(cluArr.begin(), lower);
      int upperIdx = std::distance(cluArr.begin(), upper);
      for (int jClu2 = lowerIdx; jClu2 < upperIdx; ++jClu2) {
        const NA6PBaseCluster& clu2 = cluArr[jClu2];
        double x2 = clu2.getX();
        double y2 = clu2.getY();
        double z2 = clu2.getZ();
        double r2 = std::sqrt(x2 * x2 + y2 * y2);
        double theta2 = std::atan2(z2, r2);
        double phi2 = std::atan2(y2, x2);
        double dphi = phi2 - phi1;
        if (dphi > M_PI)
          dphi -= 2 * M_PI;
        else if (dphi < -M_PI)
          dphi += 2 * M_PI;
        if (std::abs(theta2 - theta1) < mMaxDeltaThetaTracklet && std::abs(dphi) < mMaxDeltaPhiTracklet) {
          double phi = std::atan2(y2 - y1, x2 - x1);
          double tanL = (z2 - z1) / (r2 - r1);
          double pxpz = (x2 - x1) / (z2 - z1);
          double pypz = (y2 - y1) / (z2 - z1);
          bool signal = false;
          if (clu1.getParticleID() == clu2.getParticleID())
            signal = true;
          tracklets.emplace_back(iLayer, jClu1, jClu2, tanL, phi, pxpz, pypz, signal);
        }
      }
    }
  }
}

//______________________________________________________________________

void NA6PVertexerTracklets::selectTracklets(const std::vector<TrackletForVertex>& tracklets,
                                            const std::vector<int>& firstIndex,
                                            const std::vector<int>& lastIndex,
                                            const std::vector<NA6PBaseCluster>& cluArr,
                                            std::vector<TrackletForVertex>& selTracklets)
{

  selTracklets.clear();
  int iLayer = mLayerToStart;
  auto layer1Begin = tracklets.begin() + firstIndex[iLayer + 1];
  auto layer1End = tracklets.begin() + lastIndex[iLayer + 1];
  std::vector<bool> isTrackletSelected;
  isTrackletSelected.reserve(lastIndex[iLayer] - firstIndex[iLayer]);
  for (int jTrkl0 = firstIndex[iLayer]; jTrkl0 < lastIndex[iLayer]; ++jTrkl0) {
    isTrackletSelected[jTrkl0] = false;
    const TrackletForVertex& trkl01 = tracklets[jTrkl0];
    const int nextLayerClusterIndex = trkl01.secondClusterIndex;
    auto lower = std::partition_point(layer1Begin, layer1End,
                                      [&](const TrackletForVertex& trkl) {
                                        return trkl.firstClusterIndex < nextLayerClusterIndex;
                                      });
    auto upper = std::partition_point(layer1Begin, layer1End,
                                      [&](const TrackletForVertex& trkl) {
                                        return trkl.firstClusterIndex <= nextLayerClusterIndex;
                                      });
    for (auto it = lower; it != upper; ++it) {
      const TrackletForVertex& trkl12 = *it;
      if (trkl12.firstClusterIndex != nextLayerClusterIndex)
        continue;
      const double deltaTanLambda = std::abs(trkl12.tanL - trkl01.tanL);
      double dphi = trkl12.phi - trkl01.phi;
      if (dphi > M_PI)
        dphi -= 2 * M_PI;
      else if (dphi < -M_PI)
        dphi += 2 * M_PI;
      double deltapxpz = std::abs(trkl12.pxpz - trkl01.pxpz);
      double deltapypz = std::abs(trkl12.pypz - trkl01.pypz);
      if (deltapypz < mMaxDeltaPyPzInOut && deltapxpz < mMaxDeltaPxPzInOut && deltaTanLambda < mMaxDeltaTanLamInOut && std::abs(dphi) < mMaxDeltaPhiInOut) {
        // if(trkl01.isSignal) isTrackletSelected[jTrkl0] = true;
        isTrackletSelected[jTrkl0] = true;
      }
    }
    if (isTrackletSelected[jTrkl0])
      selTracklets.push_back(std::move(trkl01));
  }
}

//______________________________________________________________________

void NA6PVertexerTracklets::filterOutUsedTracklets(const std::vector<TrackletForVertex>& tracklets,
                                                   std::vector<TrackletForVertex>& usableTracklets)
{
  usableTracklets.clear();
  for (auto trkl : tracklets) {
    if (mIsClusterUsed[trkl.firstClusterIndex] || mIsClusterUsed[trkl.secondClusterIndex])
      continue;
    usableTracklets.push_back(std::move(trkl));
  }
}

//______________________________________________________________________

void NA6PVertexerTracklets::computeIntersections(const std::vector<TrackletForVertex>& selTracklets,
                                                 const std::vector<NA6PBaseCluster>& cluArr,
                                                 std::vector<TracklIntersection>& zIntersec)
{
  zIntersec.clear();
  for (auto& trkl : selTracklets) {
    auto clu0 = cluArr[trkl.firstClusterIndex];
    auto clu1 = cluArr[trkl.secondClusterIndex];
    double x0 = clu0.getX();
    double y0 = clu0.getY();
    double z0 = clu0.getZ();
    double x1 = clu1.getX();
    double y1 = clu1.getY();
    double z1 = clu1.getZ();
    double sigmayClu0 = std::sqrt(clu0.getSigYY());
    double sigmayClu1 = std::sqrt(clu1.getSigYY());
    double dx = x1 - x0;
    double dy = y1 - y0;
    double dz = z1 - z0;
    double zi = -99999.;
    double sigmazi = 0.1;
    switch (mRecoType) {
      case kYZ: {
        // work in the yz plane (non bending) to use straight line approximation of traks
        // slo = dy/dz  y = slo*z + q --> q = y0-slo*z0 --> zv = -q/slo = z0 - y0/slo
        double slo = dy / dz;
        zi = z0 - y0 / slo;
        double sigmaSlo = std::sqrt(sigmayClu0 * sigmayClu0 + sigmayClu1 * sigmayClu1) / std::abs(dz);
        double invSlo = 1.0 / slo;
        double term1 = (sigmayClu0 * invSlo) * (sigmayClu0 * invSlo);
        double term2 = (y0 * sigmaSlo / (slo * slo)) * (y0 * sigmaSlo / (slo * slo));
        sigmazi = std::sqrt(term1 + term2);
        break;
      }
      case kXZ: {
        // work in the xz plane (bending) just for debug purposes
        double slo = dx / dz;
        zi = z0 - x0 / slo;
        break;
      }
      case kRZ: {
        // solution in the r/z plane (as in AliRoot)
        double r0 = std::sqrt(x0 * x0 + y0 * y0);
        double r1 = std::sqrt(x1 * x1 + y1 * y1);
        double slo = (r1 - r0) / (z1 - z0);
        // r = slo * z + q --> q = r0 - slo * z0 --> slo * zv + q = = --> zv = -q/slo = -(r0 - slo*z0)/slo = z0 - r0/slo
        zi = z0 - r0 / slo;
        break;
      }
      case k3D: {
        // 3D solution for DCA between lines
        double denomXY = dx * dx + dy * dy;
        if (denomXY < 1e-8)
          continue;
        double tBeam = -((x0 - mBeamX) * dx + (y0 - mBeamY) * dy) / denomXY;
        double xi = x0 + tBeam * dx;
        double yi = y0 + tBeam * dy;
        double ri2 = xi * xi + yi * yi;
        if (ri2 > mMaxDCAxy * mMaxDCAxy)
          continue;
        zi = z0 + tBeam * dz;
        double sigmaDCA = sigmayClu0 * std::sqrt(2.0) * (1.0 + std::abs(tBeam));
        sigmazi = sigmaDCA * std::abs(dz) / std::sqrt(denomXY);
        break;
      }
      default:
        LOGP(error, "wrong option for vertex reco");
        continue;
    }

    if (zi >= mZMin && zi < mZMax)
      zIntersec.emplace_back(zi, sigmazi, trkl.tanL, trkl.firstClusterIndex, trkl.secondClusterIndex);
  }
}

//______________________________________________________________________

bool NA6PVertexerTracklets::findVertexHistoPeak(std::vector<TracklIntersection>& zIntersec,
                                                std::vector<NA6PVertex>& vertices)
{
  std::fill(mHistIntersec.begin(), mHistIntersec.end(), 0.0);
  for (const auto& zInt : zIntersec) {
    double zi = zInt.zeta;
    if (zi >= mZMin && zi < mZMax) {
      int bin = int((zi - mZMin) / mBinWidth);
      mHistIntersec[bin]++;
    }
  }
  // find peak (case of multiple peaks with same height is treated)
  int maxHeight = 0;
  for (int jBin = 0; jBin < mNBinsForPeakFind; ++jBin) {
    if (mHistIntersec[jBin] > maxHeight)
      maxHeight = mHistIntersec[jBin];
  }
  if (maxHeight == 0 || maxHeight < mMinCountsInPeak)
    return false;

  std::vector<std::pair<int, int>> peaks;
  for (int jBin = 0; jBin < mNBinsForPeakFind; ++jBin) {
    if (mHistIntersec[jBin] == maxHeight) {
      int lowBin = std::max(jBin - mPeakWidthBins, 0);
      int highBin = std::min(jBin + mPeakWidthBins, mNBinsForPeakFind - 1);
      double integral = 0;
      for (int jPeak = lowBin; jPeak <= highBin; ++jPeak)
        integral += mHistIntersec[jPeak];
      peaks.push_back(std::make_pair(jBin, integral));
    }
  }
  int nPeaks = peaks.size();
  if (nPeaks == 0)
    return false;
  std::sort(peaks.begin(), peaks.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
  int nPeaksSameIntegral = 0;
  int jMainPeak = -1;
  if (nPeaks > 1) {
    // check integrals around peak to define the "main" peak
    int integralMainPeak = 0;
    for (int jp = 0; jp < nPeaks; ++jp) {
      if (peaks[jp].second > integralMainPeak) {
        integralMainPeak = peaks[jp].second;
        jMainPeak = jp;
      }
    }
    for (int jp = 0; jp < nPeaks; ++jp) {
      if (peaks[jp].second == integralMainPeak)
        ++nPeaksSameIntegral;
    }
  } else if (nPeaks == 1) {
    jMainPeak = 0;
    nPeaksSameIntegral = 1;
  }
  double zPeak, zWinMin, zWinMax;
  bool peakFound = false;
  if (jMainPeak >= 0 && (nPeaks == 1 || (nPeaks > 1 && nPeaksSameIntegral == 1))) {
    // compute the z position for the peak
    zPeak = mZMin + (peaks[jMainPeak].first + 0.5) * mBinWidth;
    zWinMin = zPeak - mZWindowWidth;
    zWinMax = zPeak + mZWindowWidth;
    peakFound = true;
  } else {
    //  merge adjacent peaks
    struct Cluster {
      int integral = 0;
      int firstBin = -1;
      int lastBin = -1;
    };
    std::vector<Cluster> clusters;
    Cluster cur;
    cur.integral = peaks[0].second;
    cur.firstBin = peaks[0].first;
    cur.lastBin = peaks[0].first;
    for (int jPk = 1; jPk < nPeaks; ++jPk) {
      int prev = peaks[jPk - 1].first;
      int curr = peaks[jPk].first;
      if (std::abs(curr - prev) < mPeakWidthBins) {
        cur.integral += peaks[jPk].second;
        cur.lastBin = curr;
      } else {
        clusters.push_back(cur);
        cur.integral = peaks[jPk].second;
        cur.firstBin = peaks[jPk].first;
        cur.lastBin = peaks[jPk].first;
      }
    }
    clusters.push_back(cur);
    int jLargestClu = 0;
    int maxIntegral = clusters[0].integral;
    bool ambiguous = false;
    for (int jClu = 1; jClu < (int)clusters.size(); ++jClu) {
      if (clusters[jClu].integral > maxIntegral) {
        maxIntegral = clusters[jClu].integral;
        jLargestClu = jClu;
        ambiguous = false;
      } else if (clusters[jClu].integral == maxIntegral) {
        ambiguous = true;
      }
    }
    if (!ambiguous) {
      double zFirst = mZMin + (clusters[jLargestClu].firstBin + 0.5) * mBinWidth;
      double zLast = mZMin + (clusters[jLargestClu].lastBin + 0.5) * mBinWidth;
      zPeak = 0.5 * (zFirst + zLast);
      zWinMin = zFirst - mZWindowWidth;
      zWinMax = zLast + mZWindowWidth;
      peakFound = true;
    }
  }
  if (!peakFound) {
    // resort to use the lowest z peak
    zPeak = mZMin + (peaks[0].first + 0.5) * mBinWidth;
    zWinMin = zPeak - mZWindowWidth;
    zWinMax = zPeak + mZWindowWidth;
  }
  std::sort(zIntersec.begin(), zIntersec.end(),
            [](const TracklIntersection& a, const TracklIntersection& b) {
              return a.zeta < b.zeta;
            });
  // std::sort(zIntersec.begin(),zIntersec.end());
  auto lower = std::lower_bound(zIntersec.begin(), zIntersec.end(), zWinMin,
                                [](const TracklIntersection& a, double value) {
                                  return a.zeta < value;
                                });

  auto upper = std::upper_bound(zIntersec.begin(), zIntersec.end(), zWinMax,
                                [](double value, const TracklIntersection& a) {
                                  return value < a.zeta;
                                });
  int nContrib = 0;
  double sumZ = 0.;
  double sumZW = 0.;
  double sumW = 0.;
  for (auto it = lower; it != upper; ++it) {
    double zp = it->zeta;
    double w = 1.0 / (1.0 + it->tanl * it->tanl);
    if (mWeightedMeanOption == kSigma)
      w = 1.0 / (it->sigmazeta * it->sigmazeta);
    sumZ += zp;
    sumZW += zp * w;
    sumW += w;
    ++nContrib;
  }
  double zMean;
  if (mWeightedMeanOption == kNoWeight) {
    zMean = (nContrib > 0) ? sumZ / (double)nContrib : 0.0;
  } else {
    zMean = (sumW > 0) ? sumZW / sumW : 0.0;
  }
  // second average with better centered window
  zWinMin = zMean - mZWindowWidth;
  zWinMax = zMean + mZWindowWidth;
  lower = std::lower_bound(zIntersec.begin(), zIntersec.end(), zWinMin,
                           [](const TracklIntersection& a, double value) {
                             return a.zeta < value;
                           });

  upper = std::upper_bound(zIntersec.begin(), zIntersec.end(), zWinMax,
                           [](double value, const TracklIntersection& a) {
                             return value < a.zeta;
                           });
  nContrib = 0;
  sumZ = 0.;
  sumZW = 0.;
  sumW = 0.;
  for (auto it = lower; it != upper; ++it) {
    double zp = it->zeta;
    double w = 1.0 / (1.0 + it->tanl * it->tanl);
    if (mWeightedMeanOption == kSigma)
      w = 1.0 / (it->sigmazeta * it->sigmazeta);
    sumZ += zp;
    sumZW += zp * w;
    sumW += w;
    // assign the clusters as used (needed for pileup searches)
    mIsClusterUsed[it->firstClusterIndex] = true;
    mIsClusterUsed[it->secondClusterIndex] = true;
    ++nContrib;
  }
  if (mWeightedMeanOption == kNoWeight) {
    zMean = (nContrib > 0) ? sumZ / (double)nContrib : 0.0;
  } else {
    zMean = (sumW > 0) ? sumZW / sumW : 0.0;
  }
  double pos[3] = {mBeamX, mBeamY, zMean};
  NA6PVertex vert(pos, nContrib);
  vert.setVertexType(NA6PVertex::kTrackletPrimaryVertex);
  vertices.push_back(vert);
  return true;
}

//______________________________________________________________________

bool NA6PVertexerTracklets::findVertexKDE(const std::vector<TracklIntersection>& zIntersec,
                                          std::vector<NA6PVertex>& vertices)
{
  // KDE-based vertex finder
  if (zIntersec.empty())
    return false;

  // Prepare grid
  std::vector<double> gridZ(mNGridKDE);
  double dz = (mZMax - mZMin) / (mNGridKDE - 1);
  for (int jg = 0; jg < mNGridKDE; ++jg)
    gridZ[jg] = mZMin + jg * dz;
  // Evaluate KDE at each grid point
  std::vector<double> fkde(mNGridKDE, 0.0);
  for (int jg = 0; jg < mNGridKDE; ++jg) {
    double z = gridZ[jg];
    double sum = 0.0;
    double sumW = 0.0;
    double sumGW = 0.0;
    for (const auto& zInt : zIntersec) {
      double zi = zInt.zeta;
      double sigma = mKDEWidth;
      if (mKDEOption == kAdaptiveKDE)
        sigma = std::sqrt(mKDEWidth * mKDEWidth + zInt.sigmazeta * zInt.sigmazeta);
      double u = (z - zi) / sigma;
      double g = gaussKernel(u);
      double w = 1.0 / (1.0 + zInt.tanl * zInt.tanl);
      if (mKDEOption == kStandardKDE && mWeightedMeanOption == kSigma)
        w = 1.0 / (zInt.sigmazeta * zInt.sigmazeta);
      sum += g;
      sumGW += w * g;
      sumW += w;
    }
    if (mKDEOption == kAdaptiveKDE || mWeightedMeanOption == kNoWeight) {
      fkde[jg] = sum;
    } else {
      fkde[jg] = (sumW > 0 ? sumGW / sumW : 0.0);
    }
  }
  // Find global maximum on the grid
  double fMax = -1.0;
  int jMax = -1;
  for (int jg = 0; jg < mNGridKDE; ++jg) {
    if (fkde[jg] > fMax) {
      fMax = fkde[jg];
      jMax = jg;
    }
  }
  if (jMax < 0)
    return false;
  // Quadratic interpolation for sub-grid precision
  double zPeak = gridZ[jMax];
  if (jMax > 0 && jMax < mNGridKDE - 1) {
    double f0 = fkde[jMax - 1];
    double f1 = fkde[jMax];
    double f2 = fkde[jMax + 1];
    double denom = 2. * (f0 - 2. * f1 + f2);
    if (std::abs(denom) > 1e-12) {
      double delta = (f0 - f2) / denom;
      zPeak += delta * dz;
    }
  }
  // check peak height in small bin
  int nInPeak = 0;
  for (const auto& zInt : zIntersec) {
    double zi = zInt.zeta;
    if (zi > zPeak - 0.5 * mBinWidth && zi < zPeak + 0.5 * mBinWidth)
      ++nInPeak;
  }
  if (nInPeak < mMinCountsInPeak)
    return false;

  // count contributors and mark used hits
  double zWinMin = zPeak - mZWindowWidth;
  double zWinMax = zPeak + mZWindowWidth;
  int nContrib = 0;
  for (const auto& zInt : zIntersec) {
    double zi = zInt.zeta;
    if (zi > zWinMin && zi < zWinMax) {
      ++nContrib;
      mIsClusterUsed[zInt.firstClusterIndex] = true;
      mIsClusterUsed[zInt.secondClusterIndex] = true;
    }
  }
  double pos[3] = {mBeamX, mBeamY, zPeak};
  NA6PVertex vert(pos, nContrib);
  vert.setVertexType(NA6PVertex::kTrackletPrimaryVertex);
  vertices.push_back(vert);
  return true;
}

//______________________________________________________________________

void NA6PVertexerTracklets::findVertices(std::vector<NA6PBaseCluster>& cluArr,
                                         std::vector<NA6PVertex>& vertices)
{

  uint nClus = cluArr.size();
  LOGP(info, "Process event with nClusters {}", nClus);
  resetClusters(nClus);
  std::vector<TracklIntersection> zIntersec;
  std::vector<int> firstCluPerLay;
  std::vector<int> lastCluPerLay;
  sortClustersByLayerAndEta(cluArr, firstCluPerLay, lastCluPerLay);
  std::vector<TrackletForVertex> allTracklets;
  std::vector<TrackletForVertex> selTracklets;
  computeLayerTracklets(cluArr, firstCluPerLay, lastCluPerLay, allTracklets);
  std::vector<int> firstTrklPerLay;
  std::vector<int> lastTrklPerLay;
  sortTrackletsByLayerAndIndex(allTracklets, firstTrklPerLay, lastTrklPerLay);
  if (mVerbose)
    printStats(allTracklets, cluArr, "tracklets");
  selectTracklets(allTracklets, firstTrklPerLay, lastTrklPerLay, cluArr, selTracklets);
  if (mVerbose)
    printStats(selTracklets, cluArr, "selected tracklets");
  computeIntersections(selTracklets, cluArr, zIntersec);
  bool retCode;
  if (mMethod == kKDE)
    retCode = findVertexKDE(zIntersec, vertices);
  else
    retCode = findVertexHistoPeak(zIntersec, vertices);

  // pileup detection
  std::vector<TrackletForVertex> remainingTracklets;
  for (int jPil = 0; jPil < kMaxPileupVertices; jPil++) {
    filterOutUsedTracklets(selTracklets, remainingTracklets);
    if (mVerbose)
      printStats(remainingTracklets, cluArr, "remaining tracklets");
    computeIntersections(remainingTracklets, cluArr, zIntersec);
    if (mMethod == kKDE)
      retCode = findVertexKDE(zIntersec, vertices);
    else
      retCode = findVertexHistoPeak(zIntersec, vertices);
    if (!retCode)
      break;
    selTracklets.swap(remainingTracklets);
  }
  if (mVerbose) {
    int nVertices = vertices.size();
    LOGP(info, "Number of reconstructed vertices = {}", nVertices);
    int jv = 0;
    for (auto vert : vertices)
      LOGP(info, "Vertex %d, z = %f nContrib = %d\n", jv++, vert.getZ(), vert.getNContributors());
  }
}

//______________________________________________________________________
void NA6PVertexerTracklets::printStats(const std::vector<TrackletForVertex>& candidates,
                                       const std::vector<NA6PBaseCluster>& cluArr,
                                       const std::string& label)
{

  int nFound = candidates.size();
  LOGP(info, "Number of {} = {}", label.c_str(), nFound);
  int nGood = 0;
  int nSelected = 0;
  for (int j = 0; j < nFound; ++j) {
    const auto& tr = candidates[j];
    int nClus = 2;
    int idClus[2] = {tr.firstClusterIndex, tr.secondClusterIndex};
    int startLay = tr.startingLayer;
    int idPartTrack = -9999999;
    int nTrueClus = 0;
    for (int jClu = 0; jClu < nClus; jClu++) {
      int cluID = idClus[jClu];
      if (cluID >= 0) {
        ++nTrueClus;
        const auto& clu = cluArr[cluID];
        int jLay = clu.getDetectorID() / 4;
        if (jClu == 0 && jLay != startLay)
          LOGP(error, "mismatch in {} layers: {} {}", label.c_str(), jLay, startLay);
        int idPartClu = clu.getParticleID();
        if (jClu == 0)
          idPartTrack = idPartClu;
        else if (idPartClu != idPartTrack && idPartTrack > 0)
          idPartTrack *= (-1);
      }
    }
    nSelected++;
    if (idPartTrack > 0)
      nGood++;
  }
  if (nFound > 0)
    LOGP(info, "Fraction of good {} = {} / {} = {}",
         label.c_str(), nGood, nFound, (double)nGood / (double)nFound);
}
