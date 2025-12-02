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
  const auto pxz2pz2 = getPxz2Pz2();                        // D^2 of the into
  const auto pxz2pz = std::sqrt(pxz2pz2);                   // D of the into
  const auto p2pz = std::sqrt(pxz2pz2 + getTy() * getTy()); // A of the intro
  //   tan(psi0) = tx0 = p_x/p_z
  const auto cosPsi0 = 1.f / pxz2pz;
  const auto sinPsi0 = getTx() * cosPsi0;
  const auto kappa = kB2C * by * getQ2P() * (p2pz * cosPsi0); // curvature
  // dz = (1/kappa) [ sin(psi1) - sin(psi0) ]
  // => sin(psi1) = sinPsi0 + kappa * dz
  const auto sinPsi1 = sinPsi0 + kappa * dz;
  auto cosPsi1 = (1.f - sinPsi1) * (1.f + sinPsi1);
  if (cosPsi1 < kTiny) {
    return false;
  }
  cosPsi1 = std::sqrt(cosPsi1);
  if (cosPsi0 < 0.f) {
    cosPsi1 = -cosPsi1;
  }
  auto dPsi = cosPsi0 * sinPsi1 - cosPsi1 * sinPsi0; // could be simpler for large curvature/small step : kappa*dz/cosPsi0
  if (std::abs(dPsi) > 1.f - kTiny) {
    return false;
  }
  dPsi = std::asin(dPsi); // dPsi = (std::atan2(sinPsi1,cosPsi1) - std::atan2(sinPsi0,cosPsi0));
  if ((sinPsi0 * sinPsi1 < 0.f) && (sinPsi0 * sinPsi0 + sinPsi1 * sinPsi1 > 1.f)) {
    dPsi = sinPsi1 > 0.f ? (kPI - dPsi) : (-kPI - dPsi);
  }
  auto xUpd = (cosPsi0 - cosPsi1) / kappa;
  auto yUpd = getTy() / (pxz2pz * kappa) * dPsi;
  auto cosPsi1Inv = 1.f / cosPsi1;
  auto txNew = sinPsi1 * cosPsi1Inv;             // tx1 = tan(psi1) = sin(psi1)/cos(psi1)
  auto tyNew = getTy() * (cosPsi0 * cosPsi1Inv); // ty = p_y/p_z => ty1 = ty0 * (cosPsi0 / cosPsi1)   simplify this
  if (std::abs(kappa) < kTiny) {
    xUpd = getTx() * dz;
    yUpd = getTy() * dz;
  }
  auto dzh = 0.5 * dz;
  incParam(kX, xUpd);
  incParam(kY, yUpd);
  setTx(txNew);
  setTy(tyNew);
  mZ = z;
  return true;
}
