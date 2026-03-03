#ifndef NA6P_TRACKPARCOV_H
#define NA6P_TRACKPARCOV_H

#include "NA6PTrackPar.h"
#include "NA6PBaseCluster.h"

// Keep the same correlation checks knobs as your original
#define _CHECK_BAD_CORRELATIONS_
//#define _FIX_BAD_CORRELATIONS_x
#define _PRINT_BAD_CORRELATIONS_

/*
  Add covariance matrix and related methods to NA6PTrackPar
  State: {x, y, tx=px/pxz, ty=py/pxz, q/pxz}
*/

class NA6PTrackParCov : public NA6PTrackPar
{
 public:
  using prec_t = double;
  enum CovLables : int {
    kXX,
    kYX, kYY,
    kTxX, kTxY, kTxTx,
    kTyX, kTyY, kTyTx, kTyTy,
    kQ2PxzX, kQ2PxzY, kQ2PxzTx, kQ2PxzTy, kQ2PxzQ2Pxz };
  
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
  float getSigmaQ2PxzX() const { return mC[kQ2PxzX]; }
  float getSigmaQ2PxzY() const { return mC[kQ2PxzY]; }
  float getSigmaQ2PxzTx() const { return mC[kQ2PxzTx]; }
  float getSigmaQ2PzTy() const { return mC[kQ2PxzTy]; }
  float getSigmaQ2Pxz2() const { return mC[kQ2PxzQ2Pxz]; }

  float getCovMatElem(int i, int j) const { return mC[CovarMap[i][j]]; }
  void setCovMatElem(int i, int j, float v) { mC[CovarMap[i][j]] = v; }

  // --- propagation: state + Jacobian + covariance transport ---
  bool propagateToZ(float z, float by);
  bool propagateToZ(float z, float by, NA6PTrackPar& linRef);
  bool propagateToZ(float z, const float* b);
  bool propagateToZ(float z, const float* b, NA6PTrackPar& linRef);
  bool propagateToZ(float z, const std::array<float, 3> b) { return propagateToZ(z, b.data()); }
  bool propagateToZ(float z, const std::array<float, 3> b, NA6PTrackPar& linRef) { return propagateToZ(z, b.data(), linRef); }

  float getPredictedChi2(float xm, float ym, float sx2, float syx, float sy2) const;
  float getPredictedChi2(const std::array<float, 2>& meas, const std::array<float, 3>& cov) const { return getPredictedChi2(meas[0], meas[1], cov[0], cov[1], cov[2]); }
  float getPredictedChi2(const float meas[2], const float cov[3]) const { return getPredictedChi2(meas[0], meas[1], cov[0], cov[1], cov[2]); }
  float getPredictedChi2(const NA6PBaseCluster& cl) const { return getPredictedChi2(cl.getX(), cl.getY(), cl.getSigXX(), cl.getSigYX(), cl.getSigYY()); }

  bool update(float xm, float ym, float sx2, float sxy, float sy2);
  bool update(const std::array<float, 2>& meas, const std::array<float, 3>& cov) { return update(meas[kX], meas[kY], cov[kXX], cov[kYX], cov[kYY]); }
  bool update(const float meas[2], const float cov[3]) { return update(meas[kX], meas[kY], cov[kXX], cov[kYX], cov[kYY]); }
  bool update(const NA6PBaseCluster& cl) { return update(cl.getX(), cl.getY(), cl.getSigXX(), cl.getSigYX(), cl.getSigYY()); }

  bool correctForMaterial(float x2x0, float xrho, bool anglecorr = false);
  bool correctForMaterial(float x2x0, float xrho, NA6PTrackPar& linRef, bool anglecorr = false);

  
  void resetCovariance(float s2 = -1.);
  void checkCorrelations();
  void fixCorrelations();
  std::string asString() const;
  void printCorr() const;

  // covariance access map (lower triangle)
  constexpr static int CovarMap[5][5] = {{kXX,   kYX,   kTxX,  kTyX,  kQ2PxzX},
                                         {kYX,   kYY,   kTxY,  kTyY,  kQ2PxzY},
                                         {kTxX,  kTxY,  kTxTx, kTyTx, kQ2PxzTx},
                                         {kTyX,  kTyY,  kTyTx, kTyTy, kQ2PxzTy},
                                         {kQ2PxzX, kQ2PxzY, kQ2PxzTx,kQ2PxzTy,kQ2PxzQ2Pxz}};

 protected:
  std::array<float, 15> mC{}; // lower triangle representation

  void transportCovariance(prec_t f02, prec_t f04, prec_t f12, prec_t f13, prec_t f14, prec_t f24);

  ClassDefNV(NA6PTrackParCov, 1);
};


#endif
