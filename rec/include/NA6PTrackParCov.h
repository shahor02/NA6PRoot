// NA6PCCopyright
#ifndef NA6P_TRACKPARCOV_H
#define NA6P_TRACKPARCOV_H

#include "NA6PTrackPar.h"
#include "NA6PBaseCluster.h"

/*
  Add covariance matrix and related methods to NA6PTrackPar
*/

class NA6PTrackParCov : public NA6PTrackPar
{
 public:
  using NA6PTrackPar::NA6PTrackPar;

  NA6PTrackParCov() = default;
  NA6PTrackParCov(const NA6PTrackParCov& src) = default;
  NA6PTrackParCov(float z, const std::array<float, 5>& par, std::array<float, 15> cov) : NA6PTrackPar(z, par), mC{cov} {}
  NA6PTrackParCov(const float* xyz, const float* pxyz, int sign, float errLoose = -2);
  ~NA6PTrackParCov() = default;
  NA6PTrackParCov& operator=(const NA6PTrackParCov&) = default;
  void init(const float* xyz, const float* pxyz, int sign, float errLoose = -2);

  void setCov(const std::array<float, 15>& c) { mC = c; }
  const std::array<float, 15>& getCov() const { return mC; }
  std::array<float, 15>& getCov() { return mC; }
  float getSigmaX2() const { return mC[kXX]; }
  float getSigmaYX() const { return mC[kYX]; }
  float getSigmaY2() const { return mC[kYY]; }
  float getSigmaTxX() const { return mC[kTxX]; }
  float getSigmaTxY() const { return mC[kTxY]; }
  float getSigmaTx2() const { return mC[kTxTx]; }
  float getSigmaTyX() const { return mC[kTyX]; }
  float getSigmaTyY() const { return mC[kTyY]; }
  float getSigmaTyTx() const { return mC[kTyTx]; }
  float getSigmaTy2() const { return mC[kTyTy]; }
  float getSigmaQ2PX() const { return mC[kQ2PX]; }
  float getSigmaQ2PY() const { return mC[kQ2PY]; }
  float getSigmaQ2PTx() const { return mC[kQ2PTx]; }
  float getSigmaQ2PTy() const { return mC[kQ2PTy]; }
  float getSigmaQ2P2() const { return mC[kQ2PQ2P]; }
  float getCovMatElem(int i, int j) const { return mC[CovarMap[i][j]]; }
  void setCovMatElem(int i, int j, float v) { mC[CovarMap[i][j]] = v; }

  // --- propagation in dipole field B = (0, By, 0) ---

  bool propagateToZ(float z, float by);
  bool propagateToZ(float z, const float* b) { return propagateToZ(z, b[1]); }                   // TODO
  bool propagateToZ(float z, const std::array<float, 3> b) { return propagateToZ(z, b.data()); } // TODO
  bool propagateToDCA(float xv, float yv, float zv, float by, float epsZ = 1e-4, float epsDCA = 1e-5, int maxIt = 30);

  float getPredictedChi2(float xm, float ym, float sx2, float syx, float sy2) const;
  float getPredictedChi2(const std::array<float, 2>& meas, const std::array<float, 3>& cov) const;
  float getPredictedChi2(const float meas[2], const float cov[3]) const;
  float getPredictedChi2(const NA6PBaseCluster& cl) const;

  bool update(float xm, float ym, float sx2, float sxy, float sy2);
  bool update(const std::array<float, 2>& meas, const std::array<float, 3>& cov);
  bool update(const float meas[2], const float cov[3]);
  bool update(const NA6PBaseCluster& cl);

  bool correctForMeanMaterial(float xOverX0, float xTimesRho);

  void resetCovariance(float s2 = -1.);

  std::string asString() const;

  // access to covariance matrix by row and column
  constexpr static int CovarMap[5][5] = {{kXX, kYX, kTxX, kTyX, kQ2PX},
                                         {kYX, kYY, kTxY, kTyY, kQ2PY},
                                         {kTxX, kTxY, kTxTx, kTyTx, kQ2PTx},
                                         {kTyX, kTyY, kTyTx, kTyTy, kQ2PTy},
                                         {kQ2PX, kQ2PY, kQ2PTx, kQ2PTy, kQ2PQ2P}};

 protected:
  void propagateCov(double f02, double f03, double f04, double f13, double f22, double f23, double f24);

  std::array<float, 15> mC{}; // covariance matrix in lower triangle representation

  ClassDefNV(NA6PTrackParCov, 1);
};

inline bool NA6PTrackParCov::update(const std::array<float, 2>& meas, const std::array<float, 3>& cov)
{
  return update(meas[kX], meas[kY], cov[kXX], cov[kYX], cov[kYY]);
}

inline bool NA6PTrackParCov::update(const float meas[2], const float cov[3])
{
  return update(meas[kX], meas[kY], cov[kXX], cov[kYX], cov[kYY]);
}

inline bool NA6PTrackParCov::update(const NA6PBaseCluster& cl)
{
  return update(cl.getX(), cl.getY(), cl.getSigXX(), cl.getSigYX(), cl.getSigYY());
}

inline float NA6PTrackParCov::getPredictedChi2(const NA6PBaseCluster& cl) const
{
  return getPredictedChi2(cl.getX(), cl.getY(), cl.getSigXX(), cl.getSigYX(), cl.getSigYY());
}

inline float NA6PTrackParCov::getPredictedChi2(const std::array<float, 2>& meas, const std::array<float, 3>& cov) const
{
  return getPredictedChi2(meas[0], meas[1], cov[0], cov[1], cov[2]);
}

inline float NA6PTrackParCov::getPredictedChi2(const float meas[2], const float cov[3]) const
{
  return getPredictedChi2(meas[0], meas[1], cov[0], cov[1], cov[2]);
}

inline bool NA6PTrackParCov::propagateToDCA(float xv, float yv, float zv, float by, float epsZ, float epsDCA, int maxIt)
{
  NA6PTrackPar tTmp(*this);
  return tTmp.propagateParamToDCA(xv, yv, zv, by, epsZ, epsDCA, maxIt) && propagateToZ(tTmp.getZ(), by);
}

#endif // NA6P_TRACKPARCOV_H
