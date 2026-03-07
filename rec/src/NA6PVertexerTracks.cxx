// NA6PCCopyright
#include <NA6PVertexerTracks.h>

ClassImp(NA6PVertexerTracks)

NA6PVertexerTracks::NA6PVertexerTracks()
{
  configurePeakFinding(mZMin, mZMax, mNBinsForPeakFind);
}

void NA6PVertexerTracks::createTracksPool(const std::vector<NA6PTrack>& tracks)
{
  mTracksPool.clear();
  auto ntGlo = tracks.size();
  mTracksPool.reserve(ntGlo);
  for (uint32_t i = 0; i < ntGlo; i++) {
    NA6PTrack trc = tracks[i];
    if (!trc.propagateToDCABeamAxis(mBeamX, mBeamY, mMaxDCA))
      continue;
    auto& tvf = mTracksPool.emplace_back(trc);
    if (!tvf.isValid()) {
      mTracksPool.pop_back(); // discard bad track
      continue;
    }
  }
}

void NA6PVertexerTracks::buildAndFillHistoZ()
{
  std::fill(mHistZ.begin(), mHistZ.end(), 0.f);
  mFilledBinsZ.clear();
  for (const auto& tvf : mTracksPool) {
    float z = tvf.mLine.mOriginPoint[2];
    if (z >= mZMin && z < mZMax) {
      int bin = int((z - mZMin) / mZBinWidth);
      if (mHistZ[bin] == 0.f)
        mFilledBinsZ.push_back(bin);
      mHistZ[bin] += tvf.mWeightHisto;  // weighted fill      
    }
  }
}

int NA6PVertexerTracks::findPeakBin()
{
  if (mFilledBinsZ.empty())
    return -1;
  
  int maxBin = -1, ib = mFilledBinsZ.size(), last = ib;
  float maxv = 0.f;
  while (ib--) {
    auto bin = mFilledBinsZ[ib];
    auto v = mHistZ[bin];
    if (v > maxv) {
      maxv = v;
      maxBin = bin;
    } else if (v <= 0.f) {                      // bin was emptied
      mFilledBinsZ[ib] = mFilledBinsZ[--last]; // move last non-empty bin in place of emptied one
    }
  }
  mFilledBinsZ.resize(last);
  return maxBin;
}

int NA6PVertexerTracks::findVertices()
{
  int nfound = 0, ntr = 0;
  buildAndFillHistoZ();
  int nTrials = 0;
  while (nfound < mMaxVerticesPerCluster && nTrials < mMaxTrialsPerCluster) {
    int peakBin = findPeakBin();
    if (peakBin < 0) {
      break;
    }
    float zv = mZMin + (peakBin + 0.5) * mZBinWidth;
  }
  return nfound;
}
