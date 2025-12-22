// NA6PCCopyright

#include "NA6PTrackPar.h"
#include <stdexcept>
#include <algorithm>
#include <fairlogger/Logger.h>

std::string NA6PTrackPar::asString() const
{
  return fmt::format("Z={:+.3e} Par=[{:+.3e},{:+.3e},{:+.3e},{:+.3e},{:+.3e}]", getZ(), getX(), getY(), getTx(), getTy(), getQ2P());
}

NA6PTrackPar::NA6PTrackPar(const std::array<float, 3>& xyz, const std::array<float, 3>& pxyz, int q)
{
  if (pxyz[2] != 0) {
    mP[kQ2P] = q / std::sqrt(pxyz[0] * pxyz[0] + pxyz[1] * pxyz[1] + pxyz[2] * pxyz[2]);
    auto pzinv = 1. / pxyz[2];
    mP[kTx] = pxyz[0] * pzinv;
    mP[kTy] = pxyz[1] * pzinv;
  } else {
    throw std::invalid_argument("forbidden pz = 0 is provided");
  }
  mP[kX] = xyz[0];
  mP[kY] = xyz[1];
  mZ = xyz[2];
}

void NA6PTrackPar::initParam(const float* xyz, const float* pxyz, int sign)
{
  mZ = xyz[2];
  mP[kX] = xyz[0];
  mP[kY] = xyz[1];
  float p2 = pxyz[0] * pxyz[0] + pxyz[1] * pxyz[1] + pxyz[2] * pxyz[2];
  if (pxyz[2] < 1e-6 || p2 < 1e-6) {
    LOGP(error, "Invalid momentum Pxyz={{{},{},{}}}", pxyz[0], pxyz[1], pxyz[2]);
    throw std::invalid_argument(fmt::format("Invalid momentum Pxyz={{{},{},{}}}", pxyz[0], pxyz[1], pxyz[2]));
  }
  auto pzinv = 1.f / pxyz[2], pinv = 1.f / std::sqrt(p2);
  mP[kTx] = pxyz[0] * pzinv;
  mP[kTy] = pxyz[1] * pzinv;
  mP[kQ2P] = sign > 0 ? pinv : -pinv;
}

bool NA6PTrackPar::propagateParamToZ(float z, float by)
{
  using prec_t = double;
  const float dz = z - getZ();
  if (std::abs(dz) < 1e-6f) {
    return true;
  }
  if (std::abs(by) < kTiny || std::abs(getQ2P()) < kTiny) { // Straight line propagation in absence of bending
    incParam(kX, getTx() * dz);
    incParam(kY, getTy() * dz);
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
  if (cosPsi0 < 0) {
    cosPsi1 = -cosPsi1;
  }
  auto dPsi = cosPsi0 * sinPsi1 - cosPsi1 * sinPsi0; // could be simpler for large curvature/small step : kappa*dz/cosPsi0
  if (std::abs(dPsi) > prec_t(1) - kTiny) {
    return false;
  }
  dPsi = std::asin(dPsi); // dPsi = (std::atan2(sinPsi1,cosPsi1) - std::atan2(sinPsi0,cosPsi0));
  if ((sinPsi0 * sinPsi1 < prec_t(0)) && (sinPsi0 * sinPsi0 + sinPsi1 * sinPsi1 > prec_t(1))) {
    dPsi = sinPsi1 > prec_t(0) ? (kPI - dPsi) : (-kPI - dPsi);
  }
  float xUpd = (cosPsi0 - cosPsi1) / kappa;
  LOGP(debug, "dz:{} kappa:{} c0 {} -> c1 {} | s0 {} -> s1 {} | dphi: {} |  dx {} x {}", dz, kappa, cosPsi0, cosPsi1, sinPsi0, sinPsi1, dPsi, xUpd, getX() + xUpd);

  auto cosPsi1Inv = prec_t(1) / cosPsi1;
  float yUpd = getTy() / (pxz2pz * kappa) * dPsi;
  float txNew = sinPsi1 * cosPsi1Inv;             // tx1 = tan(psi1) = sin(psi1)/cos(psi1)
  float tyNew = getTy() * (cosPsi0 * cosPsi1Inv); // ty = p_y/p_z => ty1 = ty0 * (cosPsi0 / cosPsi1)   simplify this
  if (std::abs(kappa) < kTiny) {
    xUpd = getTx() * dz;
    yUpd = getTy() * dz;
  }
  incParam(kX, xUpd);
  incParam(kY, yUpd);
  setTx(txNew);
  setTy(tyNew);
  mZ = z;
  return true;
}

bool NA6PTrackPar::propagateParamToDCA(float xv, float yv, float zv, float by, float epsZ, float epsDCA, int maxIt)
{
  // propagate track to the point of closest aproach to the vertex
  using prec_t = double;
  // 1) PCA in the bending XZ plane: the point should be on the parametric eq. x=xC+t*dvx, z=cZ+t*dvz
  constexpr float ZeroCurv = 1e-8;
  prec_t crv = getCurvature<prec_t>(by);
  float zdcaXZ = zv;
  if (std::abs(crv) > ZeroCurv) { // use straight line approximation
    prec_t sn, cs;
    getSinCosXZ(sn, cs);
    prec_t crvi = 1. / crv, xC = getX() + cs * crvi, zC = getZ() - sn * crvi; // center of the track circle
    float dvx = xv - xC, dvz = zv - zC, dv2 = dvx * dvx + dvz * dvz, t = std::sqrt(crvi * crvi / dv2);
    zdcaXZ = zC + t * dvz;
  }
  auto tTmp(*this);
  if (!tTmp.propagateParamToZ(zdcaXZ, by)) {
    return false;
  }
  // 2) from this point, perform Newton-Raphson minimization in the parabolic approximation
  int nit = 0;
  double dz = 0., dzc = 0., dxpca = tTmp.getX() - double(xv), dypca = tTmp.getY() - double(yv), dca2 = 0.5 * (dxpca * dxpca + dypca * dypca);

  auto getD2 = [&](double ddz) {
    double dzc = ddz * crv;
    double dxpca = tTmp.getX() - xv + dz * tTmp.getTx() + 0.5 * dz * dzc;
    double dypca = tTmp.getY() - yv + dz * tTmp.getTy();
    double dca2New = 0.5 * (dxpca * dxpca + dypca * dypca);
    return dca2New;
  };
  auto getDD2dz = [&](double ddz) {
    auto d2P = getD2(ddz + 0.01), d2M = getD2(ddz - 0.01);
    auto der = (d2P - d2M) / 0.02;
    return der;
  };

  do {
    auto dca2dz = dxpca * (tTmp.getTx() + dzc) + dypca * tTmp.getTy();
    auto dzCorr = -dca2 / dca2dz;

    auto valD2H = getD2(dz);
    auto valD2HD = getDD2dz(dz);
    LOGP(debug, "Iter{} dz={}, corr={} | D2={} D2D={} | Num: D2={}, D2D={}", nit, tTmp.getZ() + dz, dzCorr, dca2, dca2dz, valD2H, valD2HD);

    dz += dzCorr;
    if (std::abs(dzCorr) < epsZ) {
      break;
    }
    dzc = dz * crv;
    dxpca = tTmp.getX() - xv + dz * tTmp.getTx() + 0.5 * dz * dzc;
    dypca = tTmp.getY() - yv + dz * tTmp.getTy();
    double dca2New = 0.5 * (dxpca * dxpca + dypca * dypca);
    if (dca2New > dca2) {
      dz -= dzCorr;
      break;
    }
    if (dca2New + dca2 - 2 * std::sqrt(dca2New * dca2) < epsDCA * epsDCA) {
      break;
    }
    dca2 = dca2New;
  } while (++nit < maxIt);
  if (!tTmp.propagateParamToZ(tTmp.getZ() + dz, by)) {
    return false;
  }
  (*this) = tTmp;
  return true;
}

float NA6PTrackPar::BetheBlochSolidOpt(float bg) const
{
  //
  // This is the parameterization of the Bethe-Bloch formula inspired by Geant with hardcoded constants and better optimization
  //
  // bg  - beta*gamma
  // rho - density [g/cm^3]
  // kp1 - density effect first junction point
  // kp2 - density effect second junction point
  // meanI - mean excitation energy [GeV]
  // meanZA - mean Z/A
  //
  // The default values for the kp* parameters are for silicon.
  // The returned value is in [GeV/(g/cm^2)].
  //
  //  constexpr value_T rho = 2.33;
  //  constexpr value_T meanI = 173e-9;
  //  constexpr value_T me = 0.511e-3;    // [GeV/c^2]

  constexpr float mK = 0.307075e-3f; // [GeV*cm^2/g]
  constexpr float kp1 = 0.20f * 2.303f;
  constexpr float kp2 = 3.00f * 2.303f;
  constexpr float meanZA = 0.49848f;
  constexpr float lhwI = -1.7175226f;         // Log(28.816 * 1e-9 * Sqrt(rho * meanZA) / meanI);
  constexpr float log2muTomeanI = 8.6839805f; // Log( 2. * me / meanI);

  float bg2 = bg * bg, beta2 = bg2 / (1. + bg2);

  //*** Density effect
  float d2 = 0.;
  const float x = std::log(bg);
  if (x > kp2) {
    d2 = lhwI - 0.5f + x;
  } else if (x > kp1) {
    float r = (kp2 - x) / (kp2 - kp1);
    d2 = lhwI - 0.5 + x + (0.5 - lhwI - kp1) * r * r * r;
  }
  auto dedx = mK * meanZA / beta2 * (log2muTomeanI + x + x - beta2 - d2);
  return dedx > 0. ? dedx : 0.;
}

bool NA6PTrackPar::correctForMeanMaterial(float, float xTimesRho)
{
  if (xTimesRho != 0.f) {
    auto p = getP(), p0 = p, p2 = p * p;
    auto en2 = p * p + mPID.getMass2();
    auto beta2 = p2 / en2;
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
    mP[kQ2P] = getSign() / p;
  }
  return true;
}

float NA6PTrackPar::BetheBlochSolidDerivative(float dedx, float bg) const
{
  //
  // This is approximate derivative of the BB over betagamm, NO check for the consistency of the provided dedx and bg is done
  // Charge 1 particle is assumed for the provied dedx. For charge > 1 particles dedx/q^2 should be provided and obtained value must be scaled by q^2
  // The call should be usually done as
  // auto dedx = BetheBlochSolidOpt(bg);
  // // if derivative needed
  // auto ddedx = BetheBlochSolidDerivative(dedx, bg, bg*bg)
  //
  // dedx - precalculate dedx for bg
  // bg  - beta*gamma
  //
  constexpr float mK = 0.307075e-3f; // [GeV*cm^2/g]
  constexpr float meanZA = 0.49848f;
  auto bg2 = bg * bg;
  auto t1 = 1.f + bg2;
  //  auto derH = (mK * meanZA * (t1+bg2) - dedx*bg2)/(bg*t1);
  auto derH = (mK * meanZA * (t1 + 1.f / bg2) - dedx) / (bg * t1);
  return derH + derH;
}
