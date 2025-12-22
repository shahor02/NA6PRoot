// NA6PCCopyright

#include "NA6PTrackParCov.h"
#include <stdexcept>
#include <algorithm>
#include <fairlogger/Logger.h>

NA6PTrackParCov::NA6PTrackParCov(const float* xyz, const float* pxyz, int sign, float errLoose) : NA6PTrackPar(xyz, pxyz, sign)
{
  setCov({1e-6f, 0.f, 1e-6f, 0.f, 0.f, 1e-6f, 0.f, 0.f, 0.f, 1e-6f, 0.f, 0.f, 0.f, 0.f, 1e-2f});
  if (errLoose >= -2)
    resetCovariance(errLoose);
}

void NA6PTrackParCov::init(const float* xyz, const float* pxyz, int sign, float errLoose)
{
  initParam(xyz, pxyz, sign);
  setCov({1e-6f, 0.f, 1e-6f, 0.f, 0.f, 1e-6f, 0.f, 0.f, 0.f, 1e-6f, 0.f, 0.f, 0.f, 0.f, 1e-2f});
  if (errLoose >= -2)
    resetCovariance(errLoose);
}

std::string NA6PTrackParCov::asString() const
{
  return fmt::format("{} Cov=[{:+.3e}, {:+.3e},{:+.3e}, {:+.3e},{:+.3e},{:+.3e}, {:+.3e},{:+.3e},{:+.3e},{:+.3e}, {:+.3e},{:+.3e},{:+.3e},{:+.3e},{:+.3e}]",
                     NA6PTrackPar::asString(), mC[0], mC[1], mC[2], mC[3], mC[4], mC[5], mC[6], mC[7], mC[8], mC[9], mC[10], mC[11], mC[12], mC[13], mC[14]);
}

bool NA6PTrackParCov::propagateToZ(float z, float by)
{
  using prec_t = double;
  const float dz = z - getZ();
  if (std::abs(dz) < 1e-6f) {
    return true;
  }
  if (std::abs(by) < kTiny || std::abs(getQ2P()) < kTiny) { // Straight line propagation in absence of bending
    incParam(kX, getTx() * dz);
    incParam(kY, getTy() * dz);
    double f23{0.}, f02{dz}, f03{0.}, f04{0.}, f13{dz}, f22{1.}, f24{0.};
    propagateCov(f02, f03, f04, f13, f22, f23, f24);
    mZ = z;
    return true;
  }
  const auto pxz2pz2 = getPxz2Pz2<prec_t>();                // D^2 of the into
  const auto pxz2pz = std::sqrt(pxz2pz2);                   // D of the into
  const auto p2pz = std::sqrt(pxz2pz2 + getTy() * getTy()); // A of the intro
  //   tan(psi0) = tx0 = p_x/p_z
  const auto cosPsi0 = prec_t(1) / pxz2pz;
  const auto sinPsi0 = getTx() * cosPsi0;
  const auto kappa = kB2C * by * getQ2P() * (p2pz * cosPsi0); // curvature
  // dz = (1/kappa) [ sin(psi1) - sin(psi0) ]
  // => sin(psi1) = sinPsi0 + kappa * dz
  const auto sinPsi1 = sinPsi0 + kappa * dz;
  auto cosPsi1 = (prec_t(1) - sinPsi1) * (prec_t(1) + sinPsi1);
  if (cosPsi1 < kTiny) {
    return false;
  }
  cosPsi1 = std::sqrt(cosPsi1);
  if (cosPsi0 < 0.f) {
    cosPsi1 = -cosPsi1;
  }
  auto dPsi = cosPsi0 * sinPsi1 - cosPsi1 * sinPsi0; // could be simpler for large curvature/small step : kappa*dz/cosPsi0
  if (std::abs(dPsi) > prec_t(1) - kTiny) {
    return false;
  }
  dPsi = std::asin(dPsi); // dPsi = (std::atan2(sinPsi1,cosPsi1) - std::atan2(sinPsi0,cosPsi0));
  if ((sinPsi0 * sinPsi1 < 0.f) && (sinPsi0 * sinPsi0 + sinPsi1 * sinPsi1 > prec_t(1))) {
    dPsi = sinPsi1 > 0.f ? (kPI - dPsi) : (-kPI - dPsi);
  }
  auto cosPsi1Inv = prec_t(1) / cosPsi1;
  float xUpd = (cosPsi0 - cosPsi1) / kappa;
  float yUpd = getTy() / (pxz2pz * kappa) * dPsi;
  float txNew = sinPsi1 * cosPsi1Inv;             // tx1 = tan(psi1) = sin(psi1)/cos(psi1)
  float tyNew = getTy() * (cosPsi0 * cosPsi1Inv); // ty = p_y/p_z => ty1 = ty0 * (cosPsi0 / cosPsi1)   simplify this
  if (std::abs(kappa) < kTiny) {
    xUpd = getTx() * dz;
    yUpd = getTy() * dz;
  }
  // ----------------------------- precalculate for the covariance transport (double) ------------------------------
  // non-0 elements of the Jacobian
  auto dzh = 0.5 * dz;
  double tmpA2 = 1. / (p2pz * p2pz), tmpD2 = 1. / (pxz2pz * pxz2pz), f24{dz * kB2C * by * p2pz / pxz2pz};
  double f23{dz * kappa * getTy() * tmpA2}, tmpXA2D2 = dz * kappa * getTx() * (tmpA2 - tmpD2);
  double f02{dz + dzh * tmpXA2D2}, f03{dzh * f23}, f04{dzh * f24}, f13{dz}, f22{1. + tmpXA2D2};
  // ---------------------------------------------------------------------------------------------------------------
  incParam(kX, xUpd);
  incParam(kY, yUpd);
  setTx(txNew);
  setTy(tyNew);
  mZ = z;
  propagateCov(f02, f03, f04, f13, f22, f23, f24);
  return true;
}

void NA6PTrackParCov::propagateCov(double f02, double f03, double f04, double f13, double f22, double f23, double f24)
{
  const double Cxx = mC[kXX], Cyx = mC[kYX], Cyy = mC[kYY], CTxX = mC[kTxX], CTxY = mC[kTxY], CTxTx = mC[kTxTx], CTyX = mC[kTyX], CTyY = mC[kTyY];
  const double CTyTx = mC[kTyTx], CTyTy = mC[kTyTy], CQpx = mC[kQ2PX], CQpy = mC[kQ2PY], CQpTx = mC[kQ2PTx], CQpTy = mC[kQ2PTy], CQpQp = mC[kQ2PQ2P];

  // --- frequently reused products with f22,f23,f24 ---
  const double f22TxTx = f22 * CTxTx;                       // f22 * C(2,2)
  const double f23TyTx = f23 * CTyTx;                       // f23 * C(3,2)
  const double f23TyTy = f23 * CTyTy;                       // f23 * C(3,3)
  const double f24Q2PTy = f24 * CQpTy;                      // f24 * C(4,3)
  const double f24Q2PQ2P = f24 * CQpQp;                     // f24 * C(4,4)
  const double f22Q2PTx = f22 * CQpTx;                      // f22 * C(4,2)
  const double f23Q2PTy = f23 * CQpTy;                      // f23 * C(4,3)
  const double f22TyTx = f22 * CTyTx;                       // f22 * C(3,2)
  const double sumTy = f22TyTx + f23TyTy + f24Q2PTy;        // common combination: f22*CTyTx + f23*CTyTy + f24*CQpTy
  const double sumTxQpTx = f22TxTx + f23TyTx + f24 * CQpTx; // common combination: f22*CTxTx + f23*CTyTx + f24*CQpTx
  const double sumQp = f22Q2PTx + f23Q2PTy + f24Q2PQ2P;     // common combination: f22*CQpTx + f23*CQpTy + f24*CQpQp

  // --- frequently reused products with f02,f03,f04 ---
  const double f02TyTx = f02 * CTyTx;   // used in c01, c03
  const double f03TyTy = f03 * CTyTy;   // used in c00, c03
  const double f04Q2PTy = f04 * CQpTy;  // used in c00, c03
  const double f02Q2PTx = f02 * CQpTx;  // used in c00, c04
  const double f03Q2PTy = f03 * CQpTy;  // used in c02, c04
  const double f04Q2PQ2P = f04 * CQpQp; // used in c00, c04

  // --- frequently reused products with f13 ---
  const double f13TyTy = f13 * CTyTy;  // used in c11, c01, c13
  const double f13TyTx = f13 * CTyTx;  // used in c01, c12
  const double f13Q2PTy = f13 * CQpTy; // used in c01, c12, c14

  double c22 = f22 * (f22TxTx + 2. * f23TyTx) + f23 * (f23TyTy + 2. * f24Q2PTy) + f24 * (f24Q2PQ2P + 2. * f22Q2PTx);
  double c23 = sumTy;
  double c24 = sumQp;
  double c00 = Cxx + f02 * (f02 * CTxTx + 2. * (CTxX + f03 * CTyTx)) + f03 * (f03 * CTyTy + 2. * (CTyX + f04Q2PTy)) + f04 * (f04Q2PQ2P + 2. * (CQpx + f02Q2PTx));
  double c11 = Cyy + f13 * (f13TyTy + 2. * CTyY);
  double c01 = Cyx + f02 * CTxY + f13 * (CTyX + f02TyTx) + f03 * (CTyY + f13TyTy) + f04 * (CQpy + f13Q2PTy);
  double c13 = CTyY + f13TyTy;
  double c14 = CQpy + f13Q2PTy;
  double c03 = CTyX + f02TyTx + f03TyTy + f04Q2PTy;
  double c04 = CQpx + f02Q2PTx + f03Q2PTy + f04Q2PQ2P;
  double c02 = f22 * CTxX + f23 * CTyX + f24 * CQpx + f02 * sumTxQpTx + f03 * sumTy + f04 * sumQp;
  double c12 = f22 * CTxY + f23 * CTyY + f24 * CQpy + f13 * sumTy;

  mC[kXX] = c00;
  mC[kYX] = c01;
  mC[kYY] = c11;
  mC[kTxX] = c02;
  mC[kTxY] = c12;
  mC[kTxTx] = c22;
  mC[kTyX] = c03;
  mC[kTyY] = c13;
  mC[kTyTx] = c23;
  // mC[kTyTy]  unchanged
  mC[kQ2PX] = c04;
  mC[kQ2PY] = c14;
  mC[kQ2PTx] = c24;
  // mC[kQ2PTy] unchanged
  // mC[kQ2PQ2P] unchanged
}

// calculate chi2 between the track and the cluster
float NA6PTrackParCov::getPredictedChi2(float xm, float ym, float sx2, float syx, float sy2) const
{
  auto exx = static_cast<double>(getSigmaX2()) + static_cast<double>(sx2);
  auto eyx = static_cast<double>(getSigmaYX()) + static_cast<double>(syx);
  auto eyy = static_cast<double>(getSigmaY2()) + static_cast<double>(sy2);
  auto det = exx * eyy - eyx * eyx;

  if (std::abs(det) < 1e-16) {
    return 1.e16;
  }

  float dx = getX() - xm;
  float dy = getY() - ym;
  float chi2 = (dx * (eyy * dx - eyx * dy) + dy * (exx * dy - dx * eyx)) / det;
  if (chi2 < 0.) {
    LOGP(warning, "Negative chi2={}, Cluster Sig2: {} {} {} Dy:{} Dz:{} | exx:{} eyx:{} eyy:{} det:{}", chi2, sx2, syx, sy2, dx, dy, exx, eyx, eyy, det);
    LOGP(warning, "Track: {}", asString());
  }
  return chi2;
}

// Kalman update with 2D measurement (xm, ym), and covariance (sx2, sxy, sy2)
bool NA6PTrackParCov::update(float xm, float ym, float sx2, float sxy, float sy2)
{
  double dx = xm - getX(), dy = ym - getY(); // innovation (residual)
  const double Cxx = mC[kXX], Cyx = mC[kYX], Cyy = mC[kYY], CTxX = mC[kTxX], CTxY = mC[kTxY], CTxTx = mC[kTxTx], CTyX = mC[kTyX], CTyY = mC[kTyY];
  const double CTyTx = mC[kTyTx], CTyTy = mC[kTyTy], CQpx = mC[kQ2PX], CQpy = mC[kQ2PY], CQpTx = mC[kQ2PTx], CQpTy = mC[kQ2PTy], CQpQp = mC[kQ2PQ2P];

  // innovation covariance S = H C H^T + R (2x2) ---
  const double S00 = Cxx + sx2, S01 = Cyx + sxy, S11 = Cyy + sy2, detS = S00 * S11 - S01 * S01;
  if (std::abs(detS) < kTiny) { // Singular or ill-conditioned S: skip update
    return false;
  }
  const double invDetS = 1. / detS;
  const double invS00 = S11 * invDetS;
  const double invS11 = S00 * invDetS;
  const double invS01 = -S01 * invDetS;

  // Kalman gain K (5x2)
  const double K0_0 = Cxx * invS00 + Cyx * invS01;
  const double K0_1 = Cxx * invS01 + Cyx * invS11;
  const double K1_0 = Cyx * invS00 + Cyy * invS01;
  const double K1_1 = Cyx * invS01 + Cyy * invS11;
  const double K2_0 = CTxX * invS00 + CTxY * invS01;
  const double K2_1 = CTxX * invS01 + CTxY * invS11;
  const double K3_0 = CTyX * invS00 + CTyY * invS01;
  const double K3_1 = CTyX * invS01 + CTyY * invS11;
  const double K4_0 = CQpx * invS00 + CQpy * invS01;
  const double K4_1 = CQpx * invS01 + CQpy * invS11;

  // State update: a += K r ---
  incParam(kX, K0_0 * dx + K0_1 * dy);
  incParam(kY, K1_0 * dx + K1_1 * dy);
  incParam(kTx, K2_0 * dx + K2_1 * dy);
  incParam(kTy, K3_0 * dx + K3_1 * dy);
  incParam(kQ2P, K4_0 * dx + K4_1 * dy);

  // Covariance update: C' = C - K S K^T ---
  // Precompute for each j:
  //   t0[j] = S00*K[j,0] + S01*K[j,1]
  //   t1[j] = S01*K[j,0] + S11*K[j,1]
  // then dP_ij = K[i,0]*t0[j] + K[i,1]*t1[j]

  const double t0_0 = S00 * K0_0 + S01 * K0_1;
  const double t0_1 = S00 * K1_0 + S01 * K1_1;
  const double t0_2 = S00 * K2_0 + S01 * K2_1;
  const double t0_3 = S00 * K3_0 + S01 * K3_1;
  const double t0_4 = S00 * K4_0 + S01 * K4_1;

  const double t1_0 = S01 * K0_0 + S11 * K0_1;
  const double t1_1 = S01 * K1_0 + S11 * K1_1;
  const double t1_2 = S01 * K2_0 + S11 * K2_1;
  const double t1_3 = S01 * K3_0 + S11 * K3_1;
  const double t1_4 = S01 * K4_0 + S11 * K4_1;

  mC[kXX] -= (K0_0 * t0_0 + K0_1 * t1_0);
  mC[kYX] -= (K1_0 * t0_0 + K1_1 * t1_0);
  mC[kYY] -= (K1_0 * t0_1 + K1_1 * t1_1);
  mC[kTxX] -= (K2_0 * t0_0 + K2_1 * t1_0);
  mC[kTxY] -= (K2_0 * t0_1 + K2_1 * t1_1);
  mC[kTxTx] -= (K2_0 * t0_2 + K2_1 * t1_2);
  mC[kTyX] -= (K3_0 * t0_0 + K3_1 * t1_0);
  mC[kTyY] -= (K3_0 * t0_1 + K3_1 * t1_1);
  mC[kTyTx] -= (K3_0 * t0_2 + K3_1 * t1_2);
  mC[kTyTy] -= (K3_0 * t0_3 + K3_1 * t1_3);
  mC[kQ2PX] -= (K4_0 * t0_0 + K4_1 * t1_0);
  mC[kQ2PY] -= (K4_0 * t0_1 + K4_1 * t1_1);
  mC[kQ2PTx] -= (K4_0 * t0_2 + K4_1 * t1_2);
  mC[kQ2PTy] -= (K4_0 * t0_3 + K4_1 * t1_3);
  mC[kQ2PQ2P] -= (K4_0 * t0_4 + K4_1 * t1_4);
  return true;
}

void NA6PTrackParCov::resetCovariance(float s2)
{
  constexpr float
    kCX2max = 100 * 100,    // SigmaX<=100cm
    kCY2max = 100 * 100,    // SigmaY<=100cm
    kCTx2max = 1 * 1,       // tgX<=1
    kCTy2max = 1 * 1,       // tgY<=1
    kC1P2max = 100 * 100,   // Sigma1/P<=100 1/GeV
    kMostProbablePt = 0.6f; // Most Probable Pt (GeV), for running with Bz=0

  // Reset the covarince matrix to "something big"
  float d0(kCX2max), d1(kCY2max), d2(kCTx2max), d3(kCTy2max), d4(kC1P2max);
  if (s2 > 1e-16f) {
    d0 = getSigmaX2() * s2;
    d1 = getSigmaY2() * s2;
    d2 = getSigmaTx2() * s2;
    d3 = getSigmaTy2() * s2;
    d4 = getSigmaQ2P2() * s2;
    if (d0 > kCX2max) {
      d0 = kCX2max;
    }
    if (d1 > kCY2max) {
      d1 = kCY2max;
    }
    if (d2 > kCTx2max) {
      d2 = kCTx2max;
    }
    if (d3 > kCTy2max) {
      d3 = kCTy2max;
    }
    if (d4 > kC1P2max) {
      d4 = kC1P2max;
    }
  }
  for (int i = 0; i < 15; i++) {
    mC[i] = 0;
  }
  mC[kXX] = d0;
  mC[kYY] = d1;
  mC[kTxTx] = d2;
  mC[kTyTy] = d3;
  mC[kQ2PQ2P] = d4;
}

bool NA6PTrackParCov::correctForMeanMaterial(float xOverX0, float xTimesRho)
{
  // multiple scattering in Rossi param, no log term

  auto p = getP(), p0 = p, p2 = p * p;
  auto en2 = p * p + mPID.getMass2();
  auto beta2 = p2 / en2;
  if (xOverX0 != 0.f) {
    auto theta2 = kMSConst2 * xOverX0 / (beta2 * p2);
    mC[kTxTx] += theta2;
    mC[kTyTy] += theta2;
  }
  if (xTimesRho != 0.f) {
    auto m = mPID.getMass();
    auto en = std::sqrt(en2), en0 = en, ekin = en - m, betagamma = p * mPID.getMassInv(), dedx = BetheBlochSolidOpt(betagamma);
#ifdef _BB_NONCONST_CORR_
    auto dedxDer = 0.f, dedx1 = dedx;
#endif
    // if (charge2 != 1) dedx *= charge2;
    auto dE = dedx * xTimesRho, dETot = 0.f;
    int na = 1 + int(std::abs(dE) / ekin * ELoss2EKinThreshInv);
    if (na > MaxELossIter) {
      na = MaxELossIter;
    }
    if (na > 1) {
      dE /= na;
      xTimesRho /= na;
#ifdef _BB_NONCONST_CORR_
      dedxDer = BetheBlochSolidDerivative(dedx1, betagamma); // require correction for non-constantness of dedx vs betagamma
                                                             // if (charge2 != 1) dedxDer *= charge2;
#endif
    }
    while (na--) {
#ifdef _BB_NONCONST_CORR_
      if (dedxDer != 0.f) { // correction for non-constantness of dedx vs beta*gamma (in linear approximation): for a single step dE -> dE * [(exp(dedxDer) - 1)/dedxDer]
        if (xTimesRho < 0.f) {
          dedxDer = -dedxDer; // E.loss ( -> positive derivative)
        }
        auto corrC = (std::exp(dedxDer) - 1.f) / dedxDer;
        dE *= corrC;
      }
#endif
      en += dE;
      if (en > m) { // stopped
        p = std::sqrt(en * en - mPID.getMass2());
      } else {
        return false;
      }
      if (na) {
        betagamma = p * mPID.getMassInv();
        dedx = BetheBlochSolidOpt(betagamma);
#ifdef _BB_NONCONST_CORR_
        dedxDer = BetheBlochSolidDerivative(dedx, betagamma);
#endif
        // if (charge2 != 1) {
        //   dedx *= charge2;
        // #ifdef _BB_NONCONST_CORR_
        //   dedxDer *= charge2;
        // #endif
        // }
        dE = dedx * xTimesRho;
      }
    }
    if (p < kMinP) {
      return false;
    }
    dETot = en - en0;
    // Linearised transform u' = f(u) for u = q/p:
    //   u0 = q/p0, p0 = q/u0
    //   u1 = q/(p0 - dE) = f(u0)
    //   f'(u0) = (p0 / p1)^2
    auto q2pscale = p0 / p, der = q2pscale * q2pscale;
    mC[kQ2PX] = der;
    mC[kQ2PY] = der;
    mC[kQ2PTx] = der;
    mC[kQ2PTy] = der;
    mC[kQ2PQ2P] = der * der;
    mP[kQ2P] = getSign() / p;
  }
  return true;
}
