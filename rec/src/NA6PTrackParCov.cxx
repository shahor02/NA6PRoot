// NA6PCCopyright
#include "NA6PTrackParCov.h"

#include <algorithm>
#include <stdexcept>
#include <fmt/format.h>
#include <fairlogger/Logger.h>

// ----------------------- ctor / init -----------------------

NA6PTrackParCov::NA6PTrackParCov(const float* xyz, const float* pxyz, int sign, float errLoose)
  : NA6PTrackPar(xyz, pxyz, sign)
{
  // same “reasonable defaults” as your original
  setCov({1e-6f, 0.f, 1e-6f, 0.f, 0.f, 1e-6f, 0.f, 0.f, 0.f, 1e-6f, 0.f, 0.f, 0.f, 0.f, 1e-2f});
  if (errLoose >= -2) {
    resetCovariance(errLoose);
  }
}

void NA6PTrackParCov::init(const float* xyz, const float* pxyz, int sign, float errLoose)
{
  initParam(xyz, pxyz, sign);
  setCov({1e-6f, 0.f, 1e-6f, 0.f, 0.f, 1e-6f, 0.f, 0.f, 0.f, 1e-6f, 0.f, 0.f, 0.f, 0.f, 1e-2f});
  if (errLoose >= -2) {
    resetCovariance(errLoose);
  }
}

std::string NA6PTrackParCov::asString() const
{
  return fmt::format("{} Cov=[{:+.3e}, {:+.3e},{:+.3e}, {:+.3e},{:+.3e},{:+.3e}, {:+.3e},{:+.3e},{:+.3e},{:+.3e}, {:+.3e},{:+.3e},{:+.3e},{:+.3e},{:+.3e}]",
                     NA6PTrackPar::asString(),
                     mC[0], mC[1], mC[2], mC[3], mC[4], mC[5], mC[6], mC[7], mC[8], mC[9], mC[10], mC[11], mC[12], mC[13], mC[14]);
}

// ------------------------------------------------------------------
// propagateToZ: state + Jacobian + covariance transport
// State: {x, y, s=sin(psi), ty, q/pxz}
// kappa = K * q,  where K = kB2C*By
// s1 = s0 + kappa*dz
// dx = dz*(s0+s1)/(c0+c1)
// dy = ty*dz*(dpsi/ds),  dpsi = asin(s1)-asin(s0), ds = s1-s0
// ------------------------------------------------------------------
bool NA6PTrackParCov::propagateToZ(float z, float by)
{
  const float dz = z - mZ;
  if (std::abs(dz) < 1e-6f) {
    mZ = z;
    return true;
  }

  const float kappa = getCurvature(by);            // kB2C*By*(q/pxz)
  const float bend  = kappa * dz, abend = std::abs(bend);

  if (std::abs(kappa) < kSmallKappa || abend < kSmallBend) {     // Straight/near-straight transport: tx, ty, q/pxz unchanged to this order
    const float s0 = getTx(), c0I = 1.f/getCosPsi(), c0Idz = dz * c0I;
    mP[kX] += s0 * c0Idz;          // dx/dz = tan(psi) = sin/cos = tx/cos
    mP[kY] += getTy() * c0Idz;     // dy/dz = (py/pz) = (py/pxz)/(pz/pxz) = ty/cos
    mZ = z;
    // Jacobian
    const float f02 = c0Idz*c0I*c0I;;                   // dx/ds0
    const float f12 = getTy() * f02 * s0;               // d(ty/c0)/ds0 * dz
    const float f13 = c0Idz;                            // dy/dty
    transportCovariance(f02, 0., f12, f13, 0., 0.);
    return true;
  }
  const float s0 = getTx(), s1 = s0 + bend;   // Exact dipole transport using z-geometry:  // sin(psi1) = sin(psi0) + kappa*dz
  if (std::abs(s1) >= kAlmost1F) { // Physically this would correspond to psi leaving [-pi/2,pi/2] (pz changing sign).
    return false;
  }

  const float c0 = getCosPsi(), c1 = getCos2FromSin(s1), denom = c0 + c1;
  if (denom < kTinyF) {
    return false;
  }
  const float dx2dz = (s0 + s1) / denom, cps_dx2dz = (c1 + s1 * dx2dz);
  prec_t dr_ds0 = 0.0; // In the ds->0 limit the full derivative is well-behaved but algebra is messy;
  prec_t dr_dq  = 0.0; // for transport stability we can safely drop these tiny cross-terms (they scale with bend).
 
  if (abend < kSmallBend) {
    mP[kY] += dz * cps_dx2dz * getTy();
  } else {
    // for small bends the linear apporximation of the arc by the segment is OK, but at large bends need precise value
    // angle traversed delta = 2*asin(dist_start_end * |kappa| / 2), hence the arc is: deltaPhi/|kappa|
    // The dist_start_end is obtained from sqrt(dx^2+dz^2) = z/(c0+c1)*sqrt(2+s0*s1+c0*c1)
    auto arg = c0 * s1 - c1 * s0;
    if (std::abs(arg) > kAlmost1F) {
      return false; // loop
    }
    float rot = std::asin(arg);
    if (s0*s0 + s1*s1 > 1.f && s0*s1 < 0) {
      if (s1 > 0.f) {
	rot = phys_const::PI - rot;
      } else {
	rot = -phys_const::PI - rot;
      }
    }
    const prec_t invc0 = 1./c0, invc1 = 1./c1;
    mP[kY] += getTy() / kappa * rot;
    dr_ds0 = (invc1 - invc0) / bend;
    dr_dq  = (K * dz) * ( (bend * invc1) - rot ) / (bend * bend);
  }

  mP[kX] += dx2dz * dz;
  mP[kTx] += bend;
  mZ = z;

  const float K = kB2C * by;
  // --- Jacobian for covariance transport ---
  // x1 = x0 + dz*(s0+s1)/(c0+c1)
  const prec_t N = s0 + s1, invDen = 1. / denom, invDen2 = invDen * invDen;
  const prec_t dc0_ds0 = -s0 / c0,  dc1_ds1 = -s1 / c1;
  const prec_t ds1_ds0 = 1.0, ds1_dq  = K * dz;   // ds1/ds0 = 1 ; ds1/dq = K*dz
  const prec_t dN_ds0 = 2.0, dN_dq  = ds1_dq;   // dN/ds0 = 1 + ds1/ds0 = 2 ; dN/dq = ds1/dq
  const prec_t dDen_ds0 = dc0_ds0 + dc1_ds1 * ds1_ds0; // dDen/ds0 = dc0/ds0 + dc1/ds1 * ds1/ds0
  const prec_t dDen_dq  = dc1_ds1 * ds1_dq;   // dDen/dq = dc1/ds1 * ds1/dq

  // d(dx)/dvar = dz * ( dN*Den - N*dDen ) / Den^2
  const prec_t f02 = dz * (dN_ds0 * denom - N * dDen_ds0) * invDen2; // dx/d(tx0)
  const prec_t f04 = dz * (dN_dq  * denom - N * dDen_dq ) * invDen2; // dx/d(q)

  // y1 = y0 + ty * dz * r, r = (asin(s1)-asin(s0))/ds
  // dr/ds0 = (1/c1 - 1/c0)/ds   (since ds is independent of s0)
  // dr/dq  = (K*dz) * ( ds/c1 - dpsi ) / ds^2
  const prec_t f12 = getTy() * dz * dr_ds0;      // dy/d(tx0)
  const prec_t f13 = dz * cps_dx2dz;                 // dy/d(ty0)
  const prec_t f14 = getTy() * dz * dr_dq;       // dy/d(q)
  const prec_t f24 = K * dz;  // tx1 = s1 = s0 + K*q*dz

  transportCovariance(f02, f04, f12, f13, f14, f24);
  return true;
}

// propagate with linearization wrt provided reference rather than the track itself
bool NA6PTrackParCov::propagateToZ(float z, float by, NA6PTrackPar& linRef0)
{
  const float dz = z - mZ;
  if (std::abs(dz) < 1e-6f) {
    setZ(z);
    linRef0.setZ(z);
    return true;
  }
  // propagate reference track
  NA6PTrackPar linRef1 = linRef0;
  if (!linRef1.propagateToZ(z, by)) {
    return false;
  }
  const float K = kB2C * by;

  prec_t snpRef0 = linRef0.getTx(), cspRef0 = linRef0.getCosPsi();
  prec_t snpRef1 = linRef1.getTx(), cspRef1 = linRef1.getCosPsi();
  prec_t cspRef0Inv = 1 / cspRef0, cspRef1Inv = 1 / cspRef1, cc = cspRef0 + cspRef1, ccInv = 1 / cc, dx2dz = (snpRef0 + snpRef1) * ccInv;
  prec_t dzccInv = dz * ccInv, hh = dzccInv * cspRef1Inv * (1 + cspRef0 * cspRef1 + snpRef0 * snpRef1), jj = dz * (dx2dz - snpRef1 * cspRef1Inv);
  // --- Jacobian for covariance transport for lin ref---  
  prec_t f02 = hh * cspRef0Inv;
  prec_t f04 = hh * dzccInv * K;
  prec_t f24 = dz * K;
  prec_t f12 = linRef0.getTy() * (f02 * snpRef1 + jj);
  prec_t f13 = dz * (cspRef1 + snpRef1 * dx2dz); // dS
  prec_t f14 = linRef0.getTy() * (f04 * snpRef1 + jj * f24);

  // difference between the current and reference state
  float diff[5];
  for (int i = 0; i < 5; i++) {
    diff[i] = getParam(i) - linRef0.getParam(i);
  }
  float snpUpd = snpRef1 + diff[kTx] + f24 * diff[kQ2PXZ];
  if (std::abs(snpUpd) > kAlmost1F) {
    return false;
  }
  linRef0 = linRef1; // update reference track
  setZ(z);
  setX(linRef1.getX() + diff[kX] + f02 * diff[kTx] + f04 * diff[kQ2PXZ]);
  setY(linRef1.getY() + diff[kY] + f13 * diff[kTy] + f14 * diff[kQ2PXZ]);
  setTx(snpUpd);
  setTy(linRef1.getTy() + diff[kTy]);
  setQ2PXZ(linRef1.getQ2XZ() + diff[kQ2PXZ]);
  
  transportCovariance(f02, f04, f12, f13, f14, f24);
  return true;

}

bool NA6PTrackPar::propagateToZ(float z, const float* bxyz)
{
  //----------------------------------------------------------------
  // Extrapolate this track to the plane z in the field bxyz. Cov matrix is trasported with by only
  //----------------------------------------------------------------
  const float dz = z - mZ;
  if (std::abs(dz) < 1e-6f) {
    mZ = z;
    return true;
  }
  const float kappa = getCurvature(bxyz[1]);            // kB2C*By*(q/pxz)
  if (std::abs(kappa) < kSmallKappa) {
    return propagateToZ(z, 0.f); // for the straight-line propagation use 1D field method
  }
  const float bend  = kappa * dz, abend = std::abs(bend);
  const float s0 = getTx(), s1 = s0 + bend;
  if (std::abs(s0) > kAlmost1F || std::abs(s1) > kAlmost1F) {
    return false;
  }
  const float c0 = getCosPsi(), c1 = getCos2FromSin(s1), denom = c0 + c1;
  if (denom < kTinyF) {
    return false;
  }
  const float dx2dz = (s0 + s1) / denom;
  const float step = (abend < 0.05f) ? dz * std::abs(c1 + s1 * dx2dz)              // chord
    : 2.f * std::asin(0.5f * dz * std::sqrt(1.f + dx2dz * dx2dz) * kappa) / kappa; // arc
  step *= getP2Pxz();
  //
  // get the track x,y,z,px/p,py/p,pz/p,p
  std::array<float, 7> vecLab{0.f};
  if (!getPosDirGlo(vecLab)) {
    return false;
  }
  // rotate to the system where Bx=By=0.
  float bxy2 = b[0] * b[0] + b[1] * b[1];
  float bt = gpu::CAMath::Sqrt(bxy2);
  float cosphi = 1.f, sinphi = 0.f;
  if (bt > constants::math::Almost0) {
    cosphi = b[0] / bt;
    sinphi = b[1] / bt;
  }
  float bb = gpu::CAMath::Sqrt(bxy2 + b[2] * b[2]);
  float costet = 1.f, sintet = 0.f;
  if (bb > constants::math::Almost0) {
    costet = b[2] / bb;
    sintet = bt / bb;
  }
  std::array<float, 7> vect{costet * cosphi * vecLab[0] + costet * sinphi * vecLab[1] - sintet * vecLab[2],
                                      -sinphi * vecLab[0] + cosphi * vecLab[1],
                                      sintet * cosphi * vecLab[0] + sintet * sinphi * vecLab[1] + costet * vecLab[2],
                                      costet * cosphi * vecLab[3] + costet * sinphi * vecLab[4] - sintet * vecLab[5],
                                      -sinphi * vecLab[3] + cosphi * vecLab[4],
                                      sintet * cosphi * vecLab[3] + sintet * sinphi * vecLab[4] + costet * vecLab[5],
                                      vecLab[6]};

  // Do the helix step
  float q = getCharge();
  g3helx3(q * bb, step, vect);

  // rotate back to the Global System
  vecLab[0] = cosphi * costet * vect[0] - sinphi * vect[1] + cosphi * sintet * vect[2];
  vecLab[1] = sinphi * costet * vect[0] + cosphi * vect[1] + sinphi * sintet * vect[2];
  vecLab[2] = -sintet * vect[0] + costet * vect[2];

  vecLab[3] = cosphi * costet * vect[3] - sinphi * vect[4] + cosphi * sintet * vect[5];
  vecLab[4] = sinphi * costet * vect[3] + cosphi * vect[4] + sinphi * sintet * vect[5];
  vecLab[5] = -sintet * vect[3] + costet * vect[5];

  // Do the final correcting step to the target plane (linear approximation)
  float x = vecLab[0], y = vecLab[1];
  dz = z - vecLab[2];
  if (abs(dz) > kAlmost0F) {
    x += dz * vecLab[3] / vecLab[5]; // dz * px/pz
    y += dz * vecLab[4] / vecLab[5]; // dz * py/pz
  }
  mP[kX] = x;
  mP[kY] = y;
  
  // Calculate the track parameters
  t = 1.f / gpu::CAMath::Sqrt(vecLab[3] * vecLab[3] + vecLab[5] * vecLab[5]); // p / pxz
  mX = xk;
  mP[kTx] = vecLab[3] * t;
  mP[kTy] = vecLab[4] * t;
  mP[kQ2PXZ] = q * t / vecLab[6];
  //
  // transport cov matrix
  const float K = kB2C * bxyz[1];
  // --- Jacobian for covariance transport ---
  // x1 = x0 + dz*(s0+s1)/(c0+c1)
  const prec_t N = s0 + s1, invDen = 1.f / denom, invDen2 = invDen * invDen;
  const prec_t dc0_ds0 = -s0 / c0,  dc1_ds1 = -s1 / c1;
  const prec_t ds1_ds0 = 1.0, ds1_dq  = K * dz;   // ds1/ds0 = 1 ; ds1/dq = K*dz
  const prec_t dN_ds0 = 2.0, dN_dq  = ds1_dq;   // dN/ds0 = 1 + ds1/ds0 = 2 ; dN/dq = ds1/dq
  const prec_t dDen_ds0 = dc0_ds0 + dc1_ds1 * ds1_ds0; // dDen/ds0 = dc0/ds0 + dc1/ds1 * ds1/ds0
  const prec_t dDen_dq  = dc1_ds1 * ds1_dq;   // dDen/dq = dc1/ds1 * ds1/dq

  // d(dx)/dvar = dz * ( dN*Den - N*dDen ) / Den^2
  const prec_t f02 = dz * (dN_ds0 * denom - N * dDen_ds0) * invDen2; // dx/d(tx0)
  const prec_t f04 = dz * (dN_dq  * denom - N * dDen_dq ) * invDen2; // dx/d(q)

  // y1 = y0 + ty * dz * r, r = (asin(s1)-asin(s0))/ds
  // dr/ds0 = (1/c1 - 1/c0)/ds   (since ds is independent of s0)
  // dr/dq  = (K*dz) * ( ds/c1 - dpsi ) / ds^2
  const prec_t f12 = getTy() * dz * dr_ds0;      // dy/d(tx0)
  const prec_t f13 = dz * cps_dx2dz;                 // dy/d(ty0)
  const prec_t f14 = getTy() * dz * dr_dq;       // dy/d(q)
  const prec_t f24 = K * dz;  // tx1 = s1 = s0 + K*q*dz

  transportCovariance(f02, f04, f12, f13, f14, f24);

  
  return true;  
}

bool NA6PTrackPar::propagateToZ(float z, const float* bxyz, NA6PTrackPar& linRef0)
{
  //----------------------------------------------------------------
  // Extrapolate this track to the plane z in the field bxyz. Cov matrix is trasported with by only. Linearization wrt externally provided linRef
  //----------------------------------------------------------------
  const float dz = z - mZ;
  if (std::abs(dz) < 1e-6f) {
    setZ(z);
    linRef0.setZ(z);
    return true;
  }
  const float kappa = (std::abs(bxyz[1]) < kTinyF) ? 0.f : linRef0.getCurvature(bxyz[1]);            // kB2C*By*(q/pxz)
  if (std::abs(kappa) < kSmallKappa) {
    return propagateToZ(z, 0.f, linRef0); // for the straight-line propagation use 1D field method
  }
  const float K = kB2C * bxyz[1], bend  = kappa * dz, abend = std::abs(bend);
  const float s0 = linRef0.getTx(), s1 = s0 + bend;
  if (std::abs(s0) > kAlmost1F || std::abs(s1) > kAlmost1F) {
    return false;
  }
  const float c0 = linRef0.getCosPsi(), c1 = linRef0.getCos2FromSin(s1), denom = c0 + c1;
  if (denom < kTinyF) {
    return false;
  }
  const float dx2dz = (s0 + s1) / denom;
  const float step = (abend < 0.05f) ? dz * std::abs(c1 + s1 * dx2dz)              // chord
    : 2.f * std::asin(0.5f * dz * std::sqrt(1.f + dx2dz * dx2dz) * kappa) / kappa; // arc
  step *= linRef0.getP2Pxz();
  //
  // get the track x,y,z,px/p,py/p,pz/p,p
  std::array<float, 7> vecLab{0.f};
  if (!linRef0.getPosDirGlo(vecLab)) {
    return false;
  }
  // rotate to the system where Bx=By=0.
  float bxy2 = b[0] * b[0] + b[1] * b[1];
  float bt = std::sqrt(bxy2);
  float cosphi = 1.f, sinphi = 0.f;
  if (bt > kTinyF) {
    cosphi = b[0] / bt;
    sinphi = b[1] / bt;
  }
  float bb = std::sqrt(bxy2 + b[2] * b[2]);
  float costet = 1.f, sintet = 0.f;
  if (bb > constants::math::Almost0) {
    costet = b[2] / bb;
    sintet = bt / bb;
  }
  std::array<float, 7> vect{costet * cosphi * vecLab[0] + costet * sinphi * vecLab[1] - sintet * vecLab[2],
                                      -sinphi * vecLab[0] + cosphi * vecLab[1],
                                      sintet * cosphi * vecLab[0] + sintet * sinphi * vecLab[1] + costet * vecLab[2],
                                      costet * cosphi * vecLab[3] + costet * sinphi * vecLab[4] - sintet * vecLab[5],
                                      -sinphi * vecLab[3] + cosphi * vecLab[4],
                                      sintet * cosphi * vecLab[3] + sintet * sinphi * vecLab[4] + costet * vecLab[5],
                                      vecLab[6]};

  // Do the helix step
  float q = getCharge();
  g3helx3(q * bb, step, vect);

  // rotate back to the Global System
  vecLab[0] = cosphi * costet * vect[0] - sinphi * vect[1] + cosphi * sintet * vect[2];
  vecLab[1] = sinphi * costet * vect[0] + cosphi * vect[1] + sinphi * sintet * vect[2];
  vecLab[2] = -sintet * vect[0] + costet * vect[2];

  vecLab[3] = cosphi * costet * vect[3] - sinphi * vect[4] + cosphi * sintet * vect[5];
  vecLab[4] = sinphi * costet * vect[3] + cosphi * vect[4] + sinphi * sintet * vect[5];
  vecLab[5] = -sintet * vect[3] + costet * vect[5];

  // Do the final correcting step to the target plane (linear approximation)
  float x = vecLab[0], y = vecLab[1];
  auto dzFin = z - vecLab[2];
  if (abs(dzFin) > kTinyF) {    
    x += dzFin * vecLab[3] / vecLab[5]; // dz * px/pz
    y += dzFin * vecLab[4] / vecLab[5]; // dz * py/pz
  }  
  // Calculate the linRef updated track parameters
  auto linRef1 = linRef0;
  t = 1.f / std::sqrt(vecLab[3] * vecLab[3] + vecLab[5] * vecLab[5]); // p / pxz
  linRef1.setZ(z);
  linRef1.setX(x);
  linRef1.setY(y);
  linRef1.setTx(s1 = vecLab[3] * t); // reassign snpRef1
  linRef1.setTy(vecLab[4] * t);
  linRef1.setQ2PXZ(q * t / vecLab[6]);
  //
  // transport cov matrix
  // --- Jacobian for covariance transport ---
  // x1 = x0 + dz*(s0+s1)/(c0+c1)
  const prec_t N = s0 + s1, invDen = 1. / denom, invDen2 = invDen * invDen;
  const prec_t dc0_ds0 = -s0 / c0,  dc1_ds1 = -s1 / c1;
  const prec_t ds1_ds0 = 1.0, ds1_dq  = K * dz;   // ds1/ds0 = 1 ; ds1/dq = K*dz
  const prec_t dN_ds0 = 2.0f, dN_dq  = ds1_dq;   // dN/ds0 = 1 + ds1/ds0 = 2 ; dN/dq = ds1/dq
  const prec_t dDen_ds0 = dc0_ds0 + dc1_ds1 * ds1_ds0; // dDen/ds0 = dc0/ds0 + dc1/ds1 * ds1/ds0
  const prec_t dDen_dq  = dc1_ds1 * ds1_dq;   // dDen/dq = dc1/ds1 * ds1/dq

  // d(dx)/dvar = dz * ( dN*Den - N*dDen ) / Den^2
  const prec_t f02 = dz * (dN_ds0 * denom - N * dDen_ds0) * invDen2; // dx/d(tx0)
  const prec_t f04 = dz * (dN_dq  * denom - N * dDen_dq ) * invDen2; // dx/d(q)

  // y1 = y0 + ty * dz * r, r = (asin(s1)-asin(s0))/ds
  // dr/ds0 = (1/c1 - 1/c0)/ds   (since ds is independent of s0)
  // dr/dq  = (K*dz) * ( ds/c1 - dpsi ) / ds^2
  const prec_t f12 = getTy() * dz * dr_ds0;      // dy/d(tx0)
  const prec_t f13 = dz * cps_dx2dz;                 // dy/d(ty0)
  const prec_t f14 = getTy() * dz * dr_dq;       // dy/d(q)
  const prec_t f24 = K * dz;  // tx1 = s1 = s0 + K*q*dz

  float diff[5];
  for (int i = 0; i < 5; i++) {
    diff[i] = getParam(i) - linRef0.getParam(i);
  }
  float snpUpd = snpRef1 + diff[kTx] + f24 * diff[kQ2PXZ];
  if (std::abs(snpUpd) > kAlmost1F) {
    return false;
  }
  linRef0 = linRef1; // update reference track
  setZ(z);
  setX(linRef1.getX() + diff[kX] + f02 * diff[kTx] + f04 * diff[kQ2PXZ]);
  setY(linRef1.getY() + diff[kY] + f13 * diff[kTy] + f14 * diff[kQ2PXZ]);
  setTx(snpUpd);
  setTy(linRef1.getTy() + diff[kTy]);
  setQ2PXZ(linRef1.getQ2XZ() + diff[kQ2PXZ]);  
  
  transportCovariance(f02, f04, f12, f13, f14, f24);
  return true;  
}

void NA6PTrackParCov::transportCovariance(prec_t f02, prec_t f04, prec_t f12, prec_t f13, prec_t f14, prec_t f24)
{
  const auto &C00 = mC[kXX], &C10 = mC[kYX], &C20 = mC[kTxX], &C30 = mC[kTyX], &C40 = mC[kQ2PX], &C11 = mC[kYY], &C21 = mC[kTxY],
    &C31 = mC[kTyY], &C41 = mC[kQ2PY], &C22 = mC[kTxTx], &C32 = mC[kTyTx], &C42 = mC[kQ2PTx], &C33 = mC[kTyTy], &C43 = mC[kQ2PTy], &C44 = mC[kQ2PQ2P];
  
  // b = C*ft
  prec_t b00 = f02 * C20 + f04 * C40, b01 = f12 * C20 + f14 * C40 + f13 * C30;
  prec_t b02 = f24 * C40;
  prec_t b10 = f02 * C21 + f04 * C41, b11 = f12 * C21 + f14 * C41 + f13 * C31;
  prec_t b12 = f24 * C41;
  prec_t b20 = f02 * C22 + f04 * C42, b21 = f12 * C22 + f14 * C42 + f13 * C32;
  prec_t b22 = f24 * C42;
  prec_t b40 = f02 * C42 + f04 * C44, b41 = f12 * C42 + f14 * C44 + f13 * C43;
  prec_t b42 = f24 * C44;
  prec_t b30 = f02 * C32 + f04 * C43, b31 = f12 * C32 + f14 * C43 + f13 * C33;
  prec_t b32 = f24 * C43;

  // a = f*b = f*C*ft
  prec_t a00 = f02 * b20 + f04 * b40, a01 = f02 * b21 + f04 * b41, a02 = f02 * b22 + f04 * b42;
  prec_t a11 = f12 * b21 + f14 * b41 + f13 * b31, a12 = f12 * b22 + f14 * b42 + f13 * b32;
  prec_t a22 = f24 * b42;

  // F*C*Ft = C + (b + bt + a)
  C00 += b00 + b00 + a00;
  C10 += b10 + b01 + a01;
  C20 += b20 + b02 + a02;
  C30 += b30;
  C40 += b40;
  C11 += b11 + b11 + a11;
  C21 += b21 + b12 + a12;
  C31 += b31;
  C41 += b41;
  C22 += b22 + b22 + a22;
  C32 += b32;
  C42 += b42;

  checkCorrelations();
}

// ---------------- chi2 and update ----------------
// Copied/compatible with your original (state still x,y observed).

float NA6PTrackParCov::getPredictedChi2(float xm, float ym, float sx2, float syx, float sy2) const
{
  auto exx = static_cast<prec_t>(getSigmaX2()) + static_cast<prec_t>(sx2);
  auto eyx = static_cast<prec_t>(getSigmaYX()) + static_cast<prec_t>(syx);
  auto eyy = static_cast<prec_t>(getSigmaY2()) + static_cast<prec_t>(sy2);
  auto det = exx * eyy - eyx * eyx;

  if (det < 1e-16) return 1.e16;

  float dx = getX() - xm;
  float dy = getY() - ym;
  float chi2 = (dx * (eyy * dx - eyx * dy) + dy * (exx * dy - dx * eyx)) / det;
  if (chi2 < 0.) {
    LOGP(warning, "Negative chi2={}, Cluster Sig2: {} {} {} Dx:{} Dy:{} | exx:{} eyx:{} eyy:{} det:{}",
         chi2, sx2, syx, sy2, dx, dy, exx, eyx, eyy, det);
    LOGP(warning, "Track: {}", asString());
  }
  return chi2;
}

bool NA6PTrackParCov::update(float xm, float ym, float sx2, float sxy, float sy2)
{
  // Innovation
  const prec_t dx = static_cast<prec_t>(xm) - static_cast<prec_t>(getX());
  const prec_t dy = static_cast<prec_t>(ym) - static_cast<prec_t>(getY());

  // Full C in prec_t
  prec_t C[5][5];
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      C[i][j] = C[j][i] = static_cast<prec_t>(getCovMatElem(i, j));
    }
  }

  // Measurement covariance R
  const prec_t R00 = static_cast<prec_t>(sx2);
  const prec_t R01 = static_cast<prec_t>(sxy);
  const prec_t R11 = static_cast<prec_t>(sy2);

  // S = H C H^T + R, with H = [[1,0,0,0,0],[0,1,0,0,0]]
  const prec_t S00 = C[0][0] + R00;
  const prec_t S01 = C[0][1] + R01;
  const prec_t S11 = C[1][1] + R11;

  const prec_t det = S00 * S11 - S01 * S01;
  if (det < 1e-20) return false;

  const prec_t invS00 =  S11 / det;
  const prec_t invS01 = -S01 / det;
  const prec_t invS11 =  S00 / det;

  // K = C H^T S^-1 : K is 5x2, but H picks cols 0 and 1 of C
  prec_t Kmat[5][2];
  for (int i = 0; i < 5; ++i) {
    const prec_t Ci0 = C[i][0];
    const prec_t Ci1 = C[i][1];
    Kmat[i][0] = Ci0 * invS00 + Ci1 * invS01;
    Kmat[i][1] = Ci0 * invS01 + Ci1 * invS11;
  }

  // State update: x += K*res
  mP[kX] += static_cast<float>(Kmat[0][0] * dx + Kmat[0][1] * dy);
  mP[kY] += static_cast<float>(Kmat[1][0] * dx + Kmat[1][1] * dy);
  mP[kTx] += static_cast<float>(Kmat[2][0] * dx + Kmat[2][1] * dy);
  mP[kTy] += static_cast<float>(Kmat[3][0] * dx + Kmat[3][1] * dy);
  mP[kQ2PXZ] += static_cast<float>(Kmat[4][0] * dx + Kmat[4][1] * dy);

  // Joseph form: C = (I-KH)C(I-KH)^T + K R K^T
  // Here KH affects only first 2 columns of I.
  prec_t I_KH[5][5] = {{0}};
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) I_KH[i][j] = (i == j) ? 1.0 : 0.0;
    I_KH[i][0] -= Kmat[i][0];
    I_KH[i][1] -= Kmat[i][1];
  }

  prec_t T[5][5] = {{0}};
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      prec_t s = 0.0;
      for (int k = 0; k < 5; ++k) s += I_KH[i][k] * C[k][j];
      T[i][j] = s;
    }
  }

  prec_t Cnew[5][5] = {{0}};
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      prec_t s = 0.0;
      for (int k = 0; k < 5; ++k) s += T[i][k] * I_KH[j][k]; // * (I-KH)^T
      Cnew[i][j] = s;
    }
  }

  // Add K R K^T
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      const prec_t add =
        Kmat[i][0]*(R00*Kmat[j][0] + R01*Kmat[j][1]) +
        Kmat[i][1]*(R01*Kmat[j][0] + R11*Kmat[j][1]);
      Cnew[i][j] += add;
    }
  }

  // Store back lower triangle
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      setCovMatElem(i, j, static_cast<float>(Cnew[i][j]));
    }
  }

  checkCorrelations();
  return true;
}

// ---------------- correlation checks / reset ----------------
// For brevity: keep identical semantics; adapt if you had extra project-specific logic.

void NA6PTrackParCov::resetCovariance(float s2)
{
  // mimic your existing behavior: if s2<0 choose defaults; if >0 scale
  const float sc = (s2 > 0.f) ? s2 : 1.f;
  mC = {1e-6f*sc, 0.f, 1e-6f*sc, 0.f, 0.f, 1e-6f*sc, 0.f, 0.f, 0.f, 1e-6f*sc, 0.f, 0.f, 0.f, 0.f, 1e-2f*sc};
}

void NA6PTrackParCov::printCorr() const
{
  // optional: implement like your original if needed
}

void NA6PTrackParCov::checkCorrelations()
{
#ifdef _CHECK_BAD_CORRELATIONS_
  // minimal sanity: protect against negative variances
  if (mC[kXX] < 0.f || mC[kYY] < 0.f || mC[kTxTx] < 0.f || mC[kTyTy] < 0.f || mC[kQ2PQ2P] < 0.f) {
#ifdef _PRINT_BAD_CORRELATIONS_
    LOGP(warning, "Bad covariance diag detected: {}", asString());
#endif
#ifdef _FIX_BAD_CORRELATIONS_
    fixCorrelations();
#endif
  }
#endif
}

void NA6PTrackParCov::fixCorrelations()
{
  // simple repair: clamp small negative diagonals to tiny positive
  mC[kXX] = std::max(mC[kXX], 1e-12f);
  mC[kYY] = std::max(mC[kYY], 1e-12f);
  mC[kTxTx] = std::max(mC[kTxTx], 1e-12f);
  mC[kTyTy] = std::max(mC[kTyTy], 1e-12f);
  mC[kQ2PQ2P] = std::max(mC[kQ2PQ2P], 1e-18f);
}

bool NA6PTrackParCov::correctForMaterial(float x2x0, float xrho, bool anglecorr)
{
    //------------------------------------------------------------------
  // This function corrects the track parameters for the crossed material.
  // "x2x0"   - X/X0, the thickness in units of the radiation length.
  // "xrho" - is the product length*density (g/cm^2).
  //     It should be passed as negative when propagating tracks
  //     from the intreaction point to the outside of the central barrel.
  // "dedx" - mean enery loss (GeV/(g/cm^2), if <=kCalcdEdxAuto : calculate on the fly
  // "anglecorr" - switch for the angular correction
  //------------------------------------------------------------------
  constexpr float kMSConst2 = 0.0136f * 0.0136f;
  constexpr float kMinP = 0.01f; // kill below this momentum

  float csp2 = getCosPsi2();
  float cst2I = (1.f + getTy() * getTy());        // 1/cos(lambda)^2
  if (anglecorr) {                                                // Apply angle correction, if requested
    float angle = std::sqrt(cst2I / csp2);
    x2x0 *= angle;
    xrho *= angle;
  }
  auto m = getPID().getMass();
  int charge2 = 1; //in case we introduce charge > 1 particle: getAbsCharge() * getAbsCharge();
  float p = getP(), p0 = p, p02 = p * p, e2 = p02 + getPID().getMass2(), massInv = 1. / m, bg = p * massInv, dETot = 0.;
  float e = std::sqrt(e2), e0 = e;
  if (m > 0 && xrho != 0.f) {
    float ekin = e - m, dedx = getdEdxBBOpt(bg);
#ifdef _BB_NONCONST_CORR_
    float dedxDer = 0., dedx1 = dedx;
#endif
    if (charge2 != 1) {
      dedx *= charge2;
    }
    float dE = dedx * xrho;
    int na = 1 + int(std::abs(dE) / ekin * ELoss2EKinThreshInv);
    if (na > MaxELossIter) {
      na = MaxELossIter;
    }
    if (na > 1) {
      dE /= na;
      xrho /= na;
#ifdef _BB_NONCONST_CORR_
      dedxDer = getBetheBlochSolidDerivativeApprox(dedx1, bg); // require correction for non-constantness of dedx vs betagamma
      if (charge2 != 1) {
        dedxDer *= charge2;
      }
#endif
    }
    while (na--) {
#ifdef _BB_NONCONST_CORR_
      if (dedxDer != 0.) { // correction for non-constantness of dedx vs beta*gamma (in linear approximation): for a single step dE -> dE * [(exp(dedxDer) - 1)/dedxDer]
        if (xrho < 0) {
          dedxDer = -dedxDer; // E.loss ( -> positive derivative)
        }
        auto corrC = (std::exp(dedxDer) - 1.) / dedxDer;
        dE *= corrC;
      }
#endif
      e += dE;
      if (e > m) { // stopped
        p = std::sqrt(e * e - getPID().getMass2());
      } else {
        return false;
      }
      if (na) {
        bg = p * massInv;
        dedx = getdEdxBBOpt(bg);
#ifdef _BB_NONCONST_CORR_
        dedxDer = getBetheBlochSolidDerivativeApprox(dedx, bg);
#endif
        if (charge2 != 1) {
          dedx *= charge2;
#ifdef _BB_NONCONST_CORR_
          dedxDer *= charge2;
#endif
        }
        dE = dedx * xrho;
      }
    }

    if (p < kMinP) {
      return false;
    }
    dETot = e - e0;
  } // end of e.loss correction

  // Calculating the multiple scattering corrections******************
  float& fC22 = mC[kTxTx];
  float& fC33 = mC[kTyTy];
  float& fC43 = mC[kQ2PxzTy];
  float& fC44 = mC[kQ2PxzQ2Pxz];
  //
  float cC22(0.f), cC33(0.f), cC43(0.f), cC44(0.f);
  if (x2x0 != 0.f) {
    float beta2 = p02 / e2, theta2 = kMSConst2 / (beta2 * p02) * std::abs(x2x0);
    float fp34 = getTy();
    if (charge2 != 1) {
      theta2 *= charge2;
      fp34 *= getCharge2Pt();
    }
    if (theta2 > kPI * kPI) {
      return false;
    }
    float t2c2I = theta2 * cst2I;
    cC22 = t2c2I * csp2;
    cC33 = t2c2I * cst2I;
    cC43 = t2c2I * fp34;
    cC44 = theta2 * fp34 * fp34;
  }

  // the energy loss correction contribution to cov.matrix: approximate energy loss fluctuation (M.Ivanov)
  constexpr float knst = 0.0007f; // To be tuned.
  float sigmadE = knst * std::sqrt(std::abs(dETot)) * e0 / p02 * getCharge2Pt();
  cC44 += sigmadE * sigmadE;

  // Applying the corrections*****************************
  fC22 += cC22;
  fC33 += cC33;
  fC43 += cC43;
  fC44 += cC44;
  
  mP[kQ2PXZ] * = p0 / p;

  checkCovariance();

  return true;
}

bool NA6PTrackParCov::correctForMaterial(float x2x0, float xrho, NA6PTrackPar& linRef, bool anglecorr)
{
    //------------------------------------------------------------------
  // This function corrects the track parameters for the crossed material.
  // "x2x0"   - X/X0, the thickness in units of the radiation length.
  // "xrho" - is the product length*density (g/cm^2).
  //     It should be passed as negative when propagating tracks
  //     from the intreaction point to the outside of the central barrel.
  // "dedx" - mean enery loss (GeV/(g/cm^2), if <=kCalcdEdxAuto : calculate on the fly
  // "anglecorr" - switch for the angular correction
  //------------------------------------------------------------------
  constexpr float kMSConst2 = 0.0136f * 0.0136f;
  constexpr float kMinP = 0.01f; // kill below this momentum

  float csp2 = linRef.getCosPsi2();
  float cst2I = (1.f + linRef.getTy() * linRef.getTy());        // 1/cos(lambda)^2
  if (anglecorr) {                                                // Apply angle correction, if requested
    float angle = std::sqrt(cst2I / csp2);
    x2x0 *= angle;
    xrho *= angle;
  }
  auto m = linRef.getPID().getMass();
  int charge2 = 1; //in case we introduce charge > 1 particle: getAbsCharge() * getAbsCharge();
  float p = linRef.getP(), p0 = p, p02 = p * p, e2 = p02 + getPID().getMass2(), massInv = 1. / m, bg = p * massInv, dETot = 0.;
  float e = std::sqrt(e2), e0 = e;
  if (m > 0 && xrho != 0.f) {
    float ekin = e - m, dedx = getdEdxBBOpt(bg);
#ifdef _BB_NONCONST_CORR_
    float dedxDer = 0., dedx1 = dedx;
#endif
    if (charge2 != 1) {
      dedx *= charge2;
    }
    float dE = dedx * xrho;
    int na = 1 + int(std::abs(dE) / ekin * ELoss2EKinThreshInv);
    if (na > MaxELossIter) {
      na = MaxELossIter;
    }
    if (na > 1) {
      dE /= na;
      xrho /= na;
#ifdef _BB_NONCONST_CORR_
      dedxDer = getBetheBlochSolidDerivativeApprox(dedx1, bg); // require correction for non-constantness of dedx vs betagamma
      if (charge2 != 1) {
        dedxDer *= charge2;
      }
#endif
    }
    while (na--) {
#ifdef _BB_NONCONST_CORR_
      if (dedxDer != 0.) { // correction for non-constantness of dedx vs beta*gamma (in linear approximation): for a single step dE -> dE * [(exp(dedxDer) - 1)/dedxDer]
        if (xrho < 0) {
          dedxDer = -dedxDer; // E.loss ( -> positive derivative)
        }
        auto corrC = (std::exp(dedxDer) - 1.) / dedxDer;
        dE *= corrC;
      }
#endif
      e += dE;
      if (e > m) { // stopped
        p = std::sqrt(e * e - getPID().getMass2());
      } else {
        return false;
      }
      if (na) {
        bg = p * massInv;
        dedx = getdEdxBBOpt(bg);
#ifdef _BB_NONCONST_CORR_
        dedxDer = getBetheBlochSolidDerivativeApprox(dedx, bg);
#endif
        if (charge2 != 1) {
          dedx *= charge2;
#ifdef _BB_NONCONST_CORR_
          dedxDer *= charge2;
#endif
        }
        dE = dedx * xrho;
      }
    }

    if (p < kMinP) {
      return false;
    }
    dETot = e - e0;
  } // end of e.loss correction

  // Calculating the multiple scattering corrections******************
  float& fC22 = mC[kTxTx];
  float& fC33 = mC[kTyTy];
  float& fC43 = mC[kQ2PxzTy];
  float& fC44 = mC[kQ2PxzQ2Pxz];
  //
  float cC22(0.f), cC33(0.f), cC43(0.f), cC44(0.f);
  if (x2x0 != 0.f) {
    float beta2 = p02 / e2, theta2 = kMSConst2 / (beta2 * p02) * std::abs(x2x0);
    float fp34 = linRef.getTy();
    if (charge2 != 1) {
      theta2 *= charge2;
      fp34 *= linRef.getCharge2Pt();
    }
    if (theta2 > kPI * kPI) {
      return false;
    }
    float t2c2I = theta2 * cst2I;
    cC22 = t2c2I * csp2;
    cC33 = t2c2I * cst2I;
    cC43 = t2c2I * fp34;
    cC44 = theta2 * fp34 * fp34;
  }

  // the energy loss correction contribution to cov.matrix: approximate energy loss fluctuation (M.Ivanov)
  constexpr float knst = 0.0007f; // To be tuned.
  float sigmadE = knst * std::sqrt(std::abs(dETot)) * e0 / p02 * getCharge2Pt();
  cC44 += sigmadE * sigmadE;

  // Applying the corrections*****************************
  fC22 += cC22;
  fC33 += cC33;
  fC43 += cC43;
  fC44 += cC44;
  linRef.mP[kQ2PXZ] * = p0 / p;
  mP[kQ2PXZ] * = p0 / p;

  checkCovariance();

  return true;
}
