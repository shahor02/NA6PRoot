// NA6PCCopyright
#ifndef NA6P_TRACKPAR_H
#define NA6P_TRACKPAR_H

#include <Rtypes.h>
#include <array>
#include <cmath>
#include "PID.h"

// uncomment this to account for the Bethe-Bloch varaion over the correction range
// #define _BB_NONCONST_CORR_

/*
  Units: cm, kGauss, GeV.
  The dominant B-field in Y dimension

  Propagation of the track parameters a = { x, y, tx=px/pz, ty=py/pz, q/p }
  from plane z0 to plane z1 in a uniform dipole field B = (0, By, 0).

  Invariants in a pure magnetic field (no energy loss):
  |p|    = const
  py     = const
  |p_xz| = sqrt(px^2 + pz^2) = const
  q/p   = const

  The motion in the x–z plane is circular with curvature
  kappa = q * By / |p_xz| = kB2C * By * (q/p_xz), where
  p_xz / pz = sqrt(1 + tx^2) = D,
  p   / pz  = sqrt(1 + tx^2 + ty^2) = A,
  q/p_xz = (q/p) * (A/D).

  => kappa = kB2C*By * q/p * A/D = kB2C*By * q/p * sqrt(1 + ty^2/(1+tx^2))

  Let psi be the angle of the momentum projection in the x–z plane:
  tan(psi) = tx  = px/pz
  sin(psi) = tx / D
  cos(psi) = 1  / D.

  We make a step of dz in z direction (from z0 to z1)

  From the geometry: sinPsi1 - sinPsi0 = dz * kappa
  => sin(Psi0+dPsi) - sinPsi0 \approx dPsi * cosPsi0
  => dPsi \approx kappa * dz / cosPsi0

  The exact z–geometry of the circular motion is used to determine
  the new angle psi1 at z1 such that the z-step is satisfied.

  x1 = x0 + (cosPsi0 - cosPsi1) / kappa                   (1)

  From psi1 we get the new slope:
  tx1 = tan(psi1)                                         (2)

  The y motion is a drift along the direction component
  n_y = py / |p_xz| = ty0 / D  and the displacement in y is linear in the arc-length in the x–z plane:
  s_perp = (psi1 - psi0) / kappa
  dy     = n_y * s_perp = (ty0 / D) * (psi1 - psi0) / kappa.
  => y1 = y0 + ty0/D * (psi1 - psi0) / kappa              (3)

  py and q/p and pxz  do not change
  => ty1 = py/pz1 = ty0*pz0/pz1 = ty0*D1/D0 = ty0*cosPsi0/cosPsi1 (4)

  For the covariance matrix calculation we use simplified (up to dz^2) state propagation
  x1 = x0 + tx0 * dz + kappa/2 * dz^2                     (1a)
  tx1 = tx0 + kappa * dz                                  (2a)
  y1 = y0 + ty0 * dz                                      (3a)
  ty1 = ty0                                               (4a)

  //-------- cov. matrix transformation ------------//
  Since the curvature kappa depends on q/p, tx and ty, we need its derivatives:

  dkappa/d(q2pt) = kB2C*By*A/D  or kappa / (q/pt) for q!=0
  dkappa/dtx = kappa * tx * (1/A^2 - 1/D^2)  = -kappa tx * ty^2(A^2 *D^2)
  dkappa/dty = kappa * ty /A^2

  Hence, the Jacobian of the state transformation F = dPar(new) / dPar(old) has the following non-0 elements
  F00 = F11 = F33 = F44 = 1
  F02 = dz + 0.5*dz^2 kappa * tx * (1/A^2 - 1/D^2)
  F03 = 0.5*dz^2 * kappa * ty / A^2
  F04 = 0.5*dz^2 * kB2C*By* A/D  | or 0.5*dz^2 kappa/(q/p) for q/p!=0
  F13 = dz
  F22 = 1 + dz * kappa * tx * (1/A^2 - 1/D^2)
  F23 = dz * kappa * ty / A^2
  F24 = dz * kB2C*By* A/D    |   or dz * kapa / (q/p)  for q/p != 0
*/

class NA6PTrackPar
{
 public:
  enum ParLabels : int { kX,
                         kY,
                         kTx,
                         kTy,
                         kQ2P }; // x,y coordinates, px/pz, py/pz, q/P
  enum CovLables : int { kXX,
                         kYX,
                         kYY,
                         kTxX,
                         kTxY,
                         kTxTx,
                         kTyX,
                         kTyY,
                         kTyTx,
                         kTyTy,
                         kQ2PX,
                         kQ2PY,
                         kQ2PTx,
                         kQ2PTy,
                         kQ2PQ2P };
  static constexpr float kB2C = -0.299792458e-3f; // GeV/(kG*cm)
  static constexpr float kTiny = 1e-12f;
  static constexpr float kPI = 3.14159265358979323846f;
  static constexpr float ELoss2EKinThreshInv = 1. / 0.025; // do not allow E.Loss correction step with dE/Ekin above the inverse of this value
  static constexpr float kMSConst2 = 0.0136f * 0.0136f;
  static constexpr float kMinP = 0.01f;   // kill below this momentum
  static constexpr int MaxELossIter = 20; // max number of iteration for the ELoss to account for BB dependence on beta*gamma

  NA6PTrackPar() = default;
  NA6PTrackPar(const NA6PTrackPar& src) = default;
  NA6PTrackPar(float z, const std::array<float, 5>& par) : mZ{z}, mP{par} {}
  NA6PTrackPar(const std::array<float, 3>& xyz, const std::array<float, 3>& pxyz, int q);
  NA6PTrackPar(const float* xyz, const float* pxyz, int sign) { initParam(xyz, pxyz, sign); }
  ~NA6PTrackPar() = default;
  NA6PTrackPar& operator=(const NA6PTrackPar&) = default;
  void initParam(const float* xyz, const float* pxyz, int sign);

  float getZ() const { return mZ; }
  float getX() const { return mP[kX]; }
  float getY() const { return mP[kY]; }
  float getTx() const { return mP[kTx]; }
  float getTy() const { return mP[kTy]; }
  float getQ2P() const { return mP[kQ2P]; }
  float getParam(ParLabels l) const { return mP[l]; }
  float getP() const { return mP[kQ2P] != 0.f ? 1. / std::abs(mP[kQ2P]) : 0.; }
  PID getPID() const { return mPID; }
  int getSign() const { return mP[kQ2P] < 0.f ? -1 : 1; }
  void setPID(PID v) { mPID = v; }

  void setParam(ParLabels l, float v) { mP[l] = v; }
  void incParam(ParLabels l, float v) { mP[l] += v; }

  const std::array<float, 5>& getParam() const { return mP; }
  std::array<float, 5>& getParam() { return mP; }

  void setZ(float v) { mZ = v; }
  void setX(float v) { mP[kX] = v; }
  void setY(float v) { mP[kY] = v; }
  void setTx(float v) { mP[kTx] = v; }
  void setTy(float v) { mP[kTy] = v; }
  void setQ2P(float v) { mP[kQ2P] = v; }

  // --- position helpers ---
  template <typename T = float>
  void getXYZ(std::array<T, 3>& xyz) const;
  template <typename T = float>
  void getXYZ(T xyz[3]) const;
  std::array<float, 3> getXYZ() const { return {getX(), getY(), getZ()}; }

  // --- momentum helpers ---
  float getPx() const { return getTx() * getPz(); }
  float getPy() const { return getTy() * getPz(); }
  float getPz() const;

  std::array<float, 3> getPXYZ() const;
  template <typename T = float>
  void getPXYZ(std::array<T, 3>& pxyz) const;
  template <typename T = float>
  void getPXYZ(T pxyz[3]) const;

  template <typename T = float>
  T getCurvature(float by) const
  {
    return kB2C * by * getQ2P() * std::sqrt(T(1) + T(getTy()) * T(getTy()) / getPxz2Pz2<T>());
  }

  // --- propagation in dipole field B = (0, By, 0) ---
  bool propagateParamToZ(float z, float by);
  bool propagateParamToZ(float z, const float* bxyz) { return propagateParamToZ(z, bxyz[1]); }                   // TODO
  bool propagateParamToZ(float z, const std::array<float, 3> bxyz) { return propagateParamToZ(z, bxyz.data()); } // TODO
  bool propagateToZ(float z, float by) { return propagateParamToZ(z, by); }
  bool propagateToZ(float z, const float* bxyz) { return propagateParamToZ(z, bxyz); }
  bool propagateToZ(float z, const std::array<float, 3> bxyz) { return propagateParamToZ(z, bxyz); }

  bool propagateParamToDCA(float xv, float yv, float zv, float by, float epsZ = 1e-4, float epsDCA = 1e-5, int maxIt = 30);
  bool correctForMeanMaterial(float, float xTimesRho);

  template <typename T = float>
  std::array<T, 3> getCircleParams(float by) const;

  std::string asString() const;

 protected:
  float BetheBlochSolidOpt(float bg) const;
  float BetheBlochSolidDerivative(float dedx, float bg) const;

  // helpers
  template <typename T = float>
  T getPxz2Pz2() const
  {
    return T(1) + T(getTx()) * T(getTx());
  } // 1 + tx^2 = pxz^2 / pz^2
  template <typename T = float>
  T getP2Pz2() const
  {
    return getPxz2Pz2<T>() + T(getTy()) * T(getTy());
  } // 1 + tx^2 + + ty^2 = (p/pz)^2
  template <typename T = float>
  T getP2Pz() const
  {
    return std::sqrt(getP2Pz2<T>());
  }
  template <typename T = float>
  T getCosXZ() const
  {
    return T(1) / std::sqrt(getPxz2Pz2());
  }
  template <typename T>
  void getSinCosXZ(T& s, T& c) const;

  float mZ{};                // evaluation z
  std::array<float, 5> mP{}; // parameters: {x,y,px/pz,py/pz,q/p}
  PID mPID{PID::Pion};       // 8 bit PID

  ClassDefNV(NA6PTrackPar, 1);
};

template <typename T>
inline void NA6PTrackPar::getXYZ(std::array<T, 3>& xyz) const
{
  xyz[0] = getX();
  xyz[1] = getY();
  xyz[2] = getZ();
}

template <typename T>
inline void NA6PTrackPar::getXYZ(T xyz[3]) const
{
  xyz[0] = getX();
  xyz[1] = getY();
  xyz[2] = getZ();
}

inline float NA6PTrackPar::getPz() const
{
  auto pzinv = getP2Pz() * std::abs(getQ2P());
  return std::abs(pzinv) > kTiny ? 1.f / pzinv : 1e6f;
}

inline std::array<float, 3> NA6PTrackPar::getPXYZ() const
{
  auto pz = getPz();
  return {getTx() * pz, getTy() * pz, pz};
}

template <typename T>
inline void NA6PTrackPar::getPXYZ(std::array<T, 3>& pxyz) const
{
  pxyz[2] = getPz();
  pxyz[0] = getTx() * pxyz[2];
  pxyz[1] = getTy() * pxyz[2];
}

template <typename T>
inline void NA6PTrackPar::getPXYZ(T pxyz[3]) const
{
  pxyz[2] = getPz();
  pxyz[0] = getTx() * pxyz[2];
  pxyz[1] = getTy() * pxyz[2];
}

template <typename T>
inline void NA6PTrackPar::getSinCosXZ(T& s, T& c) const
{
  // get sin and cos of track in the bending plane
  c = getCosXZ<T>();
  s = getTx() * c;
}

template <typename T>
std::array<T, 3> NA6PTrackPar::getCircleParams(float by) const
{
  // get circle params as {R, Xcenter, Zcenter}, for straight line just set to current coordinates
  constexpr T ZeroCurv = 1e-8;
  auto crv = getCurvature<T>(by);
  if (std::abs(crv) > ZeroCurv) {
    auto crvi = T(1) / crv;
    T sn, cs;
    getSinCosXZ(sn, cs);
    return {std::abs(crvi), getX() + cs * crvi, getZ() - sn * crvi};
  } else {
    return {T(0), T(getX()), T(getZ())};
  }
}

#endif // NA6P_TRACKPAR_H
