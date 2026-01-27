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
  auto sav = *this; // RS TOREM

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
  const auto pxz2pz2 = getPxz2Pz2<prec_t>();                // D^2 of the intro
  const auto pxz2pz = std::sqrt(pxz2pz2);                   // D of the intro
  const auto p2pz = std::sqrt(pxz2pz2 + getTy() * getTy()); // A of the intro
  //   tan(psi0) = tx0 = p_x/p_z
  const auto cosPsi0 = prec_t(1) / pxz2pz;
  const auto sinPsi0 = getTx() * cosPsi0;
  const auto K = kB2C * by, K2QP = K * getQ2P();
  const auto kappa = K2QP * (p2pz * cosPsi0); // curvature
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
  /* // The case of small kappa is already treated above.
  if (std::abs(kappa) < kTiny) {
    xUpd = getTx() * dz;
    yUpd = getTy() * dz;
  }
  */
  // ----------------------------- precalculate for the covariance transport (double) ------------------------------
  // non-0 elements of the Jacobian
  /*
  auto dzh = 0.5 * dz;
  double tmpA2 = 1. / (p2pz * p2pz), tmpD2 = 1. / (pxz2pz * pxz2pz), f24{dz * kB2C * by * p2pz / pxz2pz};
  double f23{dz * kappa * getTy() * tmpA2}, tmpXA2D2 = dz * kappa * getTx() * (tmpA2 - tmpD2);
  double f02{dz + dzh * tmpXA2D2}, f03{dzh * f23}, f04{dzh * f24}, f13{dz}, f22{1. + tmpXA2D2};
  */
  // ALT
  auto invA = 1. / p2pz;
  auto invC13 = cosPsi1Inv * cosPsi1Inv * cosPsi1Inv;
  auto invD = cosPsi0, invD3 = 1. / (pxz2pz2 * pxz2pz);

  // ------------------------------------------------------------------
  // Derivatives of kappa wrt (tx, ty, qop)
  // kappa = K*qop*(A/D)
  // ------------------------------------------------------------------
  // dA/dtx = tx/A, dA/dty = ty/A, dD/dtx = tx/D
  // dk/dtx = K*qop * tx*( 1/(A*D) - A/(D^3) )
  // dk/dty = K*qop * ty/(A*D)
  // dk/dq  = K * (A/D)
  auto dk_dtx = K2QP * getTx() * (invA * invD - p2pz * invD3);
  auto dk_dty = K2QP * getTy() * (invA * invD);
  auto dk_dq = K * (p2pz * invD);
  // ------------------------------------------------------------------
  // Derivatives of sinPsi1 wrt (tx, ty, qop)
  // s1 = s0 + kappa*dz, with ds0/dtx = 1/D^3
  // ------------------------------------------------------------------
  auto ds0_dtx = invD3;
  auto ds1_dtx = ds0_dtx + dz * dk_dtx;
  auto ds1_dty = dz * dk_dty;
  auto ds1_dq = dz * dk_dq;
  // ------------------------------------------------------------------
  // f22,f23,f24 from tx1 = tan(psi1) = s1/c1
  // d(tx1)/du = (1/c1^3) * ds1/du
  // ------------------------------------------------------------------
  double f22 = invC13 * ds1_dtx;
  double f23 = invC13 * ds1_dty;
  double f24 = invC13 * ds1_dq;
  // ------------------------------------------------------------------
  // f02,f03,f04 from x1 = x0 + (c0 - c1)/kappa
  // dx/du = ((dc0/du - dc1/du)/kappa) - (c0 - c1)/kappa^2 * dk/du
  // dc0/dtx = -tx/D^3 ; dc0/dty=0 ; dc0/dq=0
  // dc1/du  = -(s1/c1) * ds1/du
  // ------------------------------------------------------------------
  auto invK = 1. / kappa;
  auto invK2 = invK * invK;

  auto dc0_dtx = -getTx() * invD3;

  auto s1_over_c1 = sinPsi1 * cosPsi1Inv;
  auto dc1_dtx = -s1_over_c1 * ds1_dtx;
  auto dc1_dty = -s1_over_c1 * ds1_dty;
  auto dc1_dq = -s1_over_c1 * ds1_dq;

  auto c0_minus_c1 = cosPsi0 - cosPsi1; // simplify

  double f02 = ((dc0_dtx - dc1_dtx) * invK) - (c0_minus_c1 * invK2) * dk_dtx;
  double f03 = (-dc1_dty * invK) - (c0_minus_c1 * invK2) * dk_dty;
  double f04 = (-dc1_dq * invK) - (c0_minus_c1 * invK2) * dk_dq;

  // ------------------------------------------------------------------
  // f13 from y1 = y0 + (ty/D) * (dpsi/kappa)
  // here we keep the sparse structure: only dy/dty is retained.
  //
  // Let ny = ty/D (depends on ty only via ty), D depends on tx only.
  // s_perp = dpsi/kappa
  //
  // f13 = (1/D)*s_perp + (ty/D) * d(s_perp)/dty
  // d(dpsi)/dty = d(psi1)/dty = (1/c1) * ds1_dty
  // d(s_perp)/dty = ( (d(dpsi)/dty)*kappa - dpsi*dk_dty ) / kappa^2
  // ------------------------------------------------------------------
  auto s_perp = dPsi * invK;             // dpsi/kappa
  auto ddPsi_dty = cosPsi1Inv * ds1_dty; // (1/c1)*ds1_dty
  auto dsperp_dty = (ddPsi_dty * kappa - dPsi * dk_dty) * invK2;
  double f13 = invD * s_perp + (getTy() * invD) * dsperp_dty;

  // --- additional Jacobian terms for Option B (ty evolves) ---

  // d(psi1)/du = ds1_du / cosPsi1 ; d(psi0)/dtx = ds0_dtx / cosPsi0
  double ddPsi_dtx = cosPsi1Inv * ds1_dtx - (ds0_dtx / cosPsi0);
  ////  double ddPsi_dty = cosPsi1Inv * ds1_dty;     // you already used this
  double ddPsi_dq  = cosPsi1Inv * ds1_dq;

  // ds_perp/du where s_perp = dPsi/kappa
  double dsperp_dtx = (ddPsi_dtx * kappa - dPsi * dk_dtx) * invK2;
  ////  double dsperp_dty = (ddPsi_dty * kappa - dPsi * dk_dty) * invK2; // same as your current
  double dsperp_dq  = (ddPsi_dq  * kappa - dPsi * dk_dq ) * invK2;

  // ny = ty/D = ty*cosPsi0  (since cosPsi0 = 1/D)
  double ny0 = getTy() * cosPsi0;

  // y1 = y0 + ny0 * s_perp
  // f13 is existing dy/dty; now add dy/dtx and dy/d(q/p)
  double dny_dtx = getTy() * dc0_dtx;     // dc0_dtx = d(1/D)/dtx
  double f12 = dny_dtx * s_perp + ny0 * dsperp_dtx;
  double f14 = ny0 * dsperp_dq;

  // ty1 = ty0 * (cosPsi0/cosPsi1)
  double R = cosPsi0 * cosPsi1Inv;
  double invC12 = cosPsi1Inv * cosPsi1Inv;

  // f33 = dty1/dty0
  // d(1/cosPsi1)/dty = (sinPsi1 * ds1_dty)/cosPsi1^3 = sinPsi1*ds1_dty*invC13
  double f33 = R + getTy() * cosPsi0 * sinPsi1 * ds1_dty * invC13;
  
  // f34 = dty1/d(q/p)
  double f34 = getTy() * cosPsi0 * sinPsi1 * ds1_dq * invC13;

  // f32 = dty1/dtx0
  // dR/dtx = (dc0_dtx*cosPsi1 - cosPsi0*dc1_dtx)/cosPsi1^2
  double f32 = getTy() * (dc0_dtx * cosPsi1 - cosPsi0 * dc1_dtx) * invC12;
  //
  //------------------------
  auto dzh = 0.5 * dz;
  double tmpA2 = 1. / (p2pz * p2pz), tmpD2 = 1. / (pxz2pz * pxz2pz), f24XX{dz * kB2C * by * p2pz / pxz2pz};
  double f23XX{dz * kappa * getTy() * tmpA2}, tmpXA2D2 = dz * kappa * getTx() * (tmpA2 - tmpD2);
  double f02XX{dz + dzh * tmpXA2D2}, f03XX{dzh * f23}, f04XX{dzh * f24}, f13XX{dz}, f22XX{1. + tmpXA2D2};
  //

  // ---------------------------------------------------------------------------------------------------------------
  incParam(kX, xUpd);
  incParam(kY, yUpd);
  setTx(txNew);
  setTy(tyNew);
  mZ = z;
  //  propagateCov(f02, f03, f04, f13, f22, f23, f24);
  propagateCovB(f02, f03, f04, f12, f13, f14, f22, f23, f24, f32, f33, f34);
  return true;
}

void NA6PTrackParCov::propagateCovB(double f02, double f03, double f04,
                                    double f12, double f13, double f14,
                                    double f22, double f23, double f24,
                                    double f32, double f33, double f34)
{
  // Load old covariance (lower triangle)
  const double Cxx   = mC[kXX];
  const double Cyx   = mC[kYX];
  const double Cyy   = mC[kYY];

  const double CTxX  = mC[kTxX];
  const double CTxY  = mC[kTxY];
  const double CTxTx = mC[kTxTx];

  const double CTyX  = mC[kTyX];
  const double CTyY  = mC[kTyY];
  const double CTyTx = mC[kTyTx];
  const double CTyTy = mC[kTyTy];

  const double CQpx  = mC[kQ2PX];
  const double CQpy  = mC[kQ2PY];
  const double CQpTx = mC[kQ2PTx];
  const double CQpTy = mC[kQ2PTy];
  const double CQpQp = mC[kQ2PQ2P];

  // --- useful linear combos for u = (tx,ty,q) block ---
  // Cuu entries: [ CTxTx  CTyTx  CQpTx
  //               CTyTx  CTyTy  CQpTy
  //               CQpTx  CQpTy  CQpQp ]

  // For tx' row coefficients (f22,f23,f24): compute cov(tx', tx), cov(tx', ty), cov(tx', q)
  const double txp_tx = f22*CTxTx + f23*CTyTx + f24*CQpTx;
  const double txp_ty = f22*CTyTx + f23*CTyTy + f24*CQpTy;
  const double txp_q  = f22*CQpTx + f23*CQpTy + f24*CQpQp;

  // For ty' row coefficients (f32,f33,f34)
  const double typ_tx = f32*CTxTx + f33*CTyTx + f34*CQpTx;
  const double typ_ty = f32*CTyTx + f33*CTyTy + f34*CQpTy;
  const double typ_q  = f32*CQpTx + f33*CQpTy + f34*CQpQp;

  // For x' coeffs (f02,f03,f04)
  const double xp_tx  = f02*CTxTx + f03*CTyTx + f04*CQpTx;
  const double xp_ty  = f02*CTyTx + f03*CTyTy + f04*CQpTy;
  const double xp_q   = f02*CQpTx + f03*CQpTy + f04*CQpQp;

  // For y' coeffs (f12,f13,f14)
  const double yp_tx  = f12*CTxTx + f13*CTyTx + f14*CQpTx;
  const double yp_ty  = f12*CTyTx + f13*CTyTy + f14*CQpTy;
  const double yp_q   = f12*CQpTx + f13*CQpTy + f14*CQpQp;

  // --- Updated variances for tx', ty' (no x,y contributions) ---
  const double c22 =
      f22*txp_tx + f23*txp_ty + f24*txp_q;

  const double c33 =
      f32*typ_tx + f33*typ_ty + f34*typ_q;

  const double c23 =
      f32*txp_tx + f33*txp_ty + f34*txp_q;  // cov(tx',ty')

  // --- Cov with q' (q stays unchanged) ---
  const double c24 = txp_q;
  const double c34 = typ_q;

  // --- x' variance and covariances ---
  // x' = x + f02*tx + f03*ty + f04*q
  const double c00 =
      Cxx
    + 2.0*(f02*CTxX + f03*CTyX + f04*CQpx)
    + (f02*xp_tx + f03*xp_ty + f04*xp_q);

  // cov(x', tx') = cov(x,tx') + [f02 f03 f04] * cov(u,tx')
  // cov(x,tx') = f22*CTxX + f23*CTyX + f24*CQpx
  const double c02 =
      (f22*CTxX + f23*CTyX + f24*CQpx)
    + (f02*txp_tx + f03*txp_ty + f04*txp_q);

  // cov(x', ty')
  // cov(x,ty') = f32*CTxX + f33*CTyX + f34*CQpx
  const double c03 =
      (f32*CTxX + f33*CTyX + f34*CQpx)
    + (f02*typ_tx + f03*typ_ty + f04*typ_q);

  // cov(x', q')
  const double c04 =
      CQpx + f02*CQpTx + f03*CQpTy + f04*CQpQp;

  // --- y' variance and covariances ---
  // y' = y + f12*tx + f13*ty + f14*q
  const double c11 =
      Cyy
    + 2.0*(f12*CTxY + f13*CTyY + f14*CQpy)
    + (f12*yp_tx + f13*yp_ty + f14*yp_q);

  // cov(y', tx')
  // cov(y,tx') = f22*CTxY + f23*CTyY + f24*CQpy
  const double c12 =
      (f22*CTxY + f23*CTyY + f24*CQpy)
    + (f12*txp_tx + f13*txp_ty + f14*txp_q);

  // cov(y', ty')
  // cov(y,ty') = f32*CTxY + f33*CTyY + f34*CQpy
  const double c13 =
      (f32*CTxY + f33*CTyY + f34*CQpy)
    + (f12*typ_tx + f13*typ_ty + f14*typ_q);

  // cov(y', q')
  const double c14 =
      CQpy + f12*CQpTx + f13*CQpTy + f14*CQpQp;

  // --- cov(x', y') ---
  // cov(x',y') = cov(x,y) + ax*cov(u,y) + ay*cov(u,x) + ax*Cuu*ay^T
  const double c01 =
      Cyx
    + (f02*CTxY + f03*CTyY + f04*CQpy)
    + (f12*CTxX + f13*CTyX + f14*CQpx)
    + (f02*yp_tx + f03*yp_ty + f04*yp_q);

  // --- Write back packed lower triangle ---
  mC[kXX]     = c00;

  mC[kYX]     = c01;
  mC[kYY]     = c11;

  mC[kTxX]    = c02;
  mC[kTxY]    = c12;
  mC[kTxTx]   = c22;

  mC[kTyX]    = c03;
  mC[kTyY]    = c13;
  mC[kTyTx]   = c23;
  mC[kTyTy]   = c33;

  mC[kQ2PX]   = c04;
  mC[kQ2PY]   = c14;
  mC[kQ2PTx]  = c24;
  mC[kQ2PTy]  = c34;
  // mC[kQ2PQ2P] unchanged

  checkCorrelations();
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
  checkCorrelations();
}

// calculate chi2 between the track and the cluster
float NA6PTrackParCov::getPredictedChi2(float xm, float ym, float sx2, float syx, float sy2) const
{
  auto exx = static_cast<double>(getSigmaX2()) + static_cast<double>(sx2);
  auto eyx = static_cast<double>(getSigmaYX()) + static_cast<double>(syx);
  auto eyy = static_cast<double>(getSigmaY2()) + static_cast<double>(sy2);
  auto det = exx * eyy - eyx * eyx;

  if (det < 1e-16) {
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
  auto sav = *this; // TOREM

  double dx = xm - getX(), dy = ym - getY(); // innovation (residual)
  const double Cxx = mC[kXX], Cyx = mC[kYX], Cyy = mC[kYY], CTxX = mC[kTxX], CTxY = mC[kTxY], CTxTx = mC[kTxTx], CTyX = mC[kTyX], CTyY = mC[kTyY];
  const double CTyTx = mC[kTyTx], CTyTy = mC[kTyTy], CQpx = mC[kQ2PX], CQpy = mC[kQ2PY], CQpTx = mC[kQ2PTx], CQpTy = mC[kQ2PTy], CQpQp = mC[kQ2PQ2P];

  // innovation covariance S = H C H^T + R (2x2) ---
  /*
  if (sx2 / Cxx < 1e-4) {
    sx2 = Cxx * 1e-4;
  }
  if (sy2 / Cyy < 1e-4) {
    sy2 = Cyy * 1e-4;
  }
  */
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
  checkCorrelations();
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
  auto sav = *this; // TOREM

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
    mC[kQ2PX] *= der;
    mC[kQ2PY] *= der;
    mC[kQ2PTx] *= der;
    mC[kQ2PTy] *= der;
    mC[kQ2PQ2P] *= der * der;
    mP[kQ2P] = getSign() / p;
  }
  checkCorrelations();
  return true;
}

void NA6PTrackParCov::checkCorrelations()
{
#ifdef _CHECK_BAD_CORRELATIONS_
  // This function forces the abs of correlation coefficients to be <1.
  constexpr float MaxCorr = 0.999;

  bool corrBad = false;
  for (int i = 5; i--;) {
    for (int j = i; j--;) {
      auto sig2 = getCovMatElem(i, i) * getCovMatElem(j, j);
      auto cov = getCovMatElem(i, j);
      if (cov * cov >= MaxCorr * sig2) { // constrain correlation
                                         // cov = gpu::CAMath::Sqrt(sig2) * (cov > 0. ? MaxCorr : -MaxCorr);
        LOGP(warn, "Correlation {} {} is too high: {} vs {} {}", i, j, cov, getCovMatElem(i, i), getCovMatElem(j, j));
        corrBad = true;
        break;
      }
      if (corrBad) {
        break;
      }
    }
  }
  if (corrBad) {
#ifdef _PRINT_BAD_CORRELATIONS_
    printCorr();
#endif

#ifdef _FIX_BAD_CORRELATIONS_
    fixCorrelations();
#endif
  }
#endif
}

void NA6PTrackParCov::fixCorrelations()
{
  // This function forces the abs of correlation coefficients to be <1.
  constexpr float MaxCorr = 0.99;
  for (int i = 1; i < 5; i++) {
    for (int j = 0; j < i; j++) {
      auto sig2 = getCovMatElem(i, i) * getCovMatElem(j, j);
      auto cov = getCovMatElem(i, j);
      if (cov * cov >= MaxCorr * sig2) { // constrain correlation
        setCovMatElem(i, j, std::sqrt(sig2) * (cov > 0 ? MaxCorr : -MaxCorr));
      }
    }
  }
}

void NA6PTrackParCov::printCorr() const
{
  LOGP(info, "... {}", ((NA6PTrackPar*)this)->asString().c_str());
  for (int i = 0; i < 5; i++) {
    std::string ss = "... ";
    for (int j = 0; j <= i; j++) {
      auto sig2 = getCovMatElem(i, i) * getCovMatElem(j, j);
      auto cov = getCovMatElem(i, j);
      cov = i == j ? getCovMatElem(i, i) : (sig2 > 0 ? cov / std::sqrt(sig2) : -999);
      ss += fmt::format("{:+.3e} ", cov);
    }
    LOGP(info, "{}", ss);
  }
};
