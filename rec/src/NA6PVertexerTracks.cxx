// NA6PCCopyright
#include <fairlogger/Logger.h>
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include <NA6PVertexerTracks.h>

ClassImp(NA6PVertexerTracks)

NA6PVertexerTracks::NA6PVertexerTracks()
{
  configurePeakFinding(mZMin, mZMax, mNBinsForPeakFind);
}

//___________________________________________________________________
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
      mHistZ[bin] += tvf.mSig2ZI;  // weighted fill      
    }
  }
}

//___________________________________________________________________
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

//___________________________________________________________________
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
    float zSeed = mZMin + (peakBin + 0.5f) * mZBinWidth;
    NA6PVertex vtx;
    if (fitVertex(zSeed, vtx)) {
    }
  }
  return nfound;
}

//___________________________________________________________________
bool NA6PVertexerTracks::fitVertex(float zSeed, NA6PVertex& vtx)
{
  VertexSeed vtxSeed;
  vtxSeed.x = mBeamX;
  vtxSeed.y = mBeamY;
  vtxSeed.z = zSeed;
  vtxSeed.setScale(mInitScaleSigma2, mTukey2I);
  vtxSeed.scaleSigma2Prev = mInitScaleSigma2;
  vtx.setChi2(1.e30);
  
  FitStatus result = FitStatus::IterateFurther;
  while (result == FitStatus::IterateFurther) {
    vtxSeed.resetForNewIteration();
    vtxSeed.nIterations++;
    result = fitIteration(vtxSeed);

    if (result == FitStatus::OK) {
      result = evalIterations(vtxSeed);
    } else if (result == FitStatus::NotEnoughTracks) {
      if (vtxSeed.nIterations <= mMaxIterations && upscaleSigma(vtxSeed)) {
        result = FitStatus::IterateFurther;
        continue; // redo with stronger rescaling
      } else {
        break;
      }
    } else if (result == FitStatus::PoolEmpty || result == FitStatus::Failure) {
      break;
    }
  }

  if (result != FitStatus::OK) {
    vtx.setChi2(vtxSeed.maxScaleSigma2Tested);
    return false;
  }
  vtx.setXYZ(vtxSeed.x, vtxSeed.y, vtxSeed.z);
  vtx.setNContributors(vtxSeed.nContributors);
  vtx.setChi2(vtxSeed.getChi2());
  vtx.setCov(vtxSeed.cxx, vtxSeed.cxy, vtxSeed.cxz, vtxSeed.cyy, vtxSeed.cyz, vtxSeed.czz);

  return true;
}

//___________________________________________________________________
NA6PVertexerTracks::FitStatus NA6PVertexerTracks::fitIteration(VertexSeed& vtxSeed)
{
  int nTested = 0;
  for (auto& tvf : mTracksPool) {
    if (tvf.canUse()) {
      accountTrack(tvf, vtxSeed);
      nTested++;
    }
  }

  vtxSeed.maxScaleSigma2Tested = vtxSeed.scaleSigma2;
  if (vtxSeed.getNContributors() < mMinTracksPerVtx) {
    return nTested < mMinTracksPerVtx ? FitStatus::PoolEmpty : FitStatus::NotEnoughTracks;
  }
  if (!solveVertex(vtxSeed)) {
    return FitStatus::Failure;
  }
  return FitStatus::OK;
}

//___________________________________________________________________
void NA6PVertexerTracks::accountTrack(TrackVF& trc, VertexSeed& vtxSeed) const
{
  // deltas defined as track - vertex
  float vtxPos[3] = {vtxSeed.x, vtxSeed.y, vtxSeed.z};
  auto res = trc.getResiduals(vtxPos); 
  float dx = res[0], dy = res[1], dz = res[2];
  auto chi2T = trc.evalChi2ToVertex(dx, dy);
  float wghT = (1.f - chi2T * vtxSeed.scaleSig2ITuk2I); // weighted distance to vertex
  if (wghT < kAlmost0F) {
    trc.wgh = 0.f;
    return;
  }
  wghT *= wghT;
  float sxxI = trc.mSig2XI * wghT;
  float syyI = trc.mSig2YI * wghT;
  float sxyI = trc.mSigXYI * wghT;
  trc.wgh = wghT;
  vtxSeed.wghSum += wghT;
  vtxSeed.wghChi2 += wghT * chi2T;
  vtxSeed.cxx += sxxI;
  vtxSeed.cyy += syyI;
  vtxSeed.cxy += sxyI;
  vtxSeed.cx0 += sxxI * dx + sxyI * dy;
  vtxSeed.cy0 += sxyI * dx + syyI * dy;
  vtxSeed.czz += trc.mSig2ZI * wghT;
  vtxSeed.cxz += trc.mSigXZI * wghT;
  vtxSeed.cyz += trc.mSigYZI * wghT;
  vtxSeed.cz0 += wghT * (trc.mSig2ZI * dz + trc.mSigXZI * dx + trc.mSigYZI * dy);

  vtxSeed.nContributors++;
}

//___________________________________________________________________
bool NA6PVertexerTracks::solveVertex(VertexSeed& vtxSeed) const
{
  ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> mat;
  mat(0, 0) = vtxSeed.cxx;
  mat(0, 1) = vtxSeed.cxy;
  mat(0, 2) = vtxSeed.cxz;
  mat(1, 1) = vtxSeed.cyy;
  mat(1, 2) = vtxSeed.cyz;
  mat(2, 2) = vtxSeed.czz;
  if (!mat.InvertFast()) {
    LOGP(error, "Failed to invert matrix");
    return false;
  }
  ROOT::Math::SVector<double, 3> rhs(vtxSeed.cx0, vtxSeed.cy0, vtxSeed.cz0);
  auto sol = mat * rhs;
  vtxSeed.x = sol(0);
  vtxSeed.y = sol(1);
  vtxSeed.z = sol(2);
  vtxSeed.cxx = mat(0,0);  vtxSeed.cxy = mat(1,0);  vtxSeed.cyy = mat(1,1);
  vtxSeed.cxz = mat(2,0);  vtxSeed.cyz = mat(2,1);  vtxSeed.czz = mat(2,2);

  vtxSeed.setChi2((vtxSeed.getNContributors() - vtxSeed.wghSum) / vtxSeed.scaleSig2ITuk2I); // calculate chi^2
  auto newScale = vtxSeed.wghSum > 0. ? vtxSeed.wghChi2 / vtxSeed.wghSum : mInitScaleSigma2;
  vtxSeed.setScale(newScale < mMinScale2 ? mMinScale2 : newScale, mTukey2I);
  return true;
}

//___________________________________________________________________
NA6PVertexerTracks::FitStatus NA6PVertexerTracks::evalIterations(VertexSeed& vtxSeed) const
{
  // decide if new iteration should be done, prepare next one if needed
  // if scaleSigma2 reached its lower limit stop
  FitStatus result = FitStatus::IterateFurther;

  if (vtxSeed.nIterations > mMaxIterations) {
    result = FitStatus::Failure;
    return result;
  } else if (vtxSeed.scaleSigma2Prev <= mMinScale2 + kAlmost0F) {
    result = FitStatus::OK;
  }

  if (result == FitStatus::OK) {
    auto chi2Mean = vtxSeed.getChi2() / vtxSeed.getNContributors();
    if (chi2Mean > mMaxChi2Mean) {
      result = FitStatus::Failure;
    } else {
      return result;
    }
  }

  if (vtxSeed.scaleSigma2 > vtxSeed.scaleSigma2Prev) {
    if (++vtxSeed.nScaleIncrease > mMaxNScaleIncreased) {
      result = FitStatus::Failure;
    }
  } else if (vtxSeed.scaleSigma2 > mSlowConvergenceFactor * vtxSeed.scaleSigma2Prev) {
    if (++vtxSeed.nScaleSlowConvergence > mMaxNScaleSlowConvergence) {
      if (vtxSeed.scaleSigma2 < mAcceptableScale2) {
        vtxSeed.setScale(mMinScale2, mTukey2I);
        result = FitStatus::IterateFurther;
      } else {
        result = FitStatus::Failure;
      }
    }
  } else {
    vtxSeed.nScaleSlowConvergence = 0;
  }
  return result;
}

//___________________________________________________________________
bool NA6PVertexerTracks::upscaleSigma(VertexSeed& vtxSeed) const
{
  // scale upward the scaleSigma2 if needed
  if (vtxSeed.scaleSigma2 < mMaxScale2) {
    auto s = vtxSeed.scaleSigma2 * mUpscaleFactor;
    vtxSeed.setScale(s > mMaxScale2 ? mMaxScale2 : s, mTukey2I);
    return true;
  }
  return false;
}

