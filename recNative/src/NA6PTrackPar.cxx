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

bool NA6PTrackPar::propagateParam(float z, float by)
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
  if (!tTmp.propagateParam(zdcaXZ, by)) {
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
  if (!tTmp.propagateParam(tTmp.getZ() + dz, by)) {
    return false;
  }
  (*this) = tTmp;
  return true;
}
