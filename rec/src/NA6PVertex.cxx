// NA6PCCopyright

#include <fmt/format.h>
#include <fairlogger/Logger.h>
#include "NA6PVertex.h"

ClassImp(NA6PVertex)

//_______________________________________________________________________
NA6PVertex::NA6PVertex(const ROOT::Math::XYZPointF& pos, const std::array<float, kNCov>& cov, int nCont, float chi2) :
  mPos{pos},
  mCov{cov},
  mChi2{chi2},
  mNContributors{nCont},
  mVertexType{kBaseVertex}
{
}
//_______________________________________________________________________
NA6PVertex::NA6PVertex(const double *xyz, int nCont) :
  mPos{},
  mCov{},
  mChi2{0.0},
  mNContributors{nCont},
  mVertexType{kBaseVertex}
{
  setXYZ(xyz[0],xyz[1],xyz[2]);  
}
//_______________________________________________________________________
void NA6PVertex::init(const double *xyz, const double *cov, int nCont, float chi2)
{
  setXYZ(xyz[0],xyz[1],xyz[2]);
  setCov(cov[kCovXX],cov[kCovXY],cov[kCovYY],cov[kCovXZ],cov[kCovYZ],cov[kCovZZ]);
  mChi2 = chi2;
  mNContributors = nCont;
}

//_______________________________________________________________________
std::string NA6PVertex::asString() const
{
  return fmt::format("Vtx {{{:+.4e},{:+.4e},{:+.4e}}} Cov.:{{{{{:.3e}..}},{{{:.3e},{:.3e}..}},{{{:.3e},{:.3e},{:.3e}}}}} Contributors {}  Chi2 {}",
                     mPos.X(), mPos.Y(), mPos.Z(), mCov[0], mCov[1], mCov[2], mCov[3], mCov[4], mCov[5],mNContributors,mChi2);
}

//_______________________________________________________________________
void NA6PVertex::print() const
{
  LOGP(info, "{}", asString());
}

