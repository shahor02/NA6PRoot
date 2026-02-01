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



// ------------------------------------------------------------------
// propagateToZ: state + full Jacobian + covariance transport
// Includes:
//  - (2) straight-line branch for negligible curvature (0 field or q/p≈0)
//  - (3) double precision internal math
// ------------------------------------------------------------------
bool NA6PTrackParCov::propagateToZ(float znew, float b_y)
{
  constexpr double kB2C = -0.299792458e-3; // GeV/(kG*cm)
  constexpr double kTiny  = 1e-15;
  constexpr double kSmallBend = 1e-6;      // threshold on |kappa*dz|
  constexpr double kSmallDs   = 1e-10;     // threshold on |s1-s0|
  constexpr double kSmallKappa = 1e-14;    // threshold on |kappa| for "0 field"

  const double dz = static_cast<double>(znew) - static_cast<double>(mZ);
  if (dz == 0.0) return true;

  const double K = kB2C * static_cast<double>(b_y);

  // State in double
  double x  = static_cast<double>(mP[0]);
  double y  = static_cast<double>(mP[1]);
  const double tx0 = static_cast<double>(mP[2]);
  const double ty0 = static_cast<double>(mP[3]);
  const double q0  = static_cast<double>(mP[4]);

  // Straight-line / negligible curvature branch
  // Triggered by very small |kappa| (i.e. small K or small q/p).
  {
    const double D  = std::sqrt(1.0 + tx0*tx0);
    const double A  = std::sqrt(1.0 + tx0*tx0 + ty0*ty0);
    const double kappa = K * q0 * (A / D);

    if (std::fabs(kappa) < kSmallKappa) {
      // State update: x += tx*dz ; y += ty*dz ; tx,ty,q unchanged
      x += tx0 * dz;
      y += ty0 * dz;

      mP[0] = static_cast<float>(x);
      mP[1] = static_cast<float>(y);
      // tx, ty, q/p unchanged
      mZ = znew;

      // Jacobian for straight line:
      // x1 = x0 + dz*tx0
      // y1 = y0 + dz*ty0
      // tx1=tx0, ty1=ty0, q1=q0
      const double f02 = dz, f03 = 0.0, f04 = 0.0;
      const double f12 = 0.0, f13 = dz,  f14 = 0.0;
      const double f22 = 1.0, f23 = 0.0, f24 = 0.0;
      const double f32 = 0.0, f33 = 1.0, f34 = 0.0;

      transportCovarianceFullJac(f02,f03,f04, f12,f13,f14, f22,f23,f24, f32,f33,f34);
      return true;
    }
  }

  // ---- Dipole propagation (Option B) ----
  const double tx2 = tx0*tx0;
  const double ty2 = ty0*ty0;
  const double D   = std::sqrt(1.0 + tx2);         // sqrt(1+tx^2)
  const double A   = std::sqrt(1.0 + tx2 + ty2);   // sqrt(1+tx^2+ty^2)
  const double invD  = 1.0 / D;
  const double invD2 = invD*invD;
  const double invD3 = invD2*invD;
  const double invA  = 1.0 / A;

  const double s0 = tx0 * invD; // sin(psi0)
  const double c0 = invD;       // cos(psi0) = 1/D

  const double kappa = K * q0 * (A * invD);

  double s1 = s0 + kappa * dz;
  s1 = clamp_pm1(s1);

  double c1sq = 1.0 - s1*s1;
  if (c1sq < 0.0) c1sq = 0.0;
  double c1 = std::sqrt(c1sq);
  if (c1 < kTiny) c1 = kTiny;

  const double invC1  = 1.0 / c1;
  const double invC12 = invC1 * invC1;
  const double invC13 = invC12 * invC1;

  // robust dpsi using atan2(s,c)
  const double psi0 = std::atan2(s0, c0);
  const double psi1 = std::atan2(s1, c1);
  double dpsi = psi1 - psi0;

  const double ds = s1 - s0;
  double ratio_dpsi_ds;
  if (std::fabs(ds) > kSmallDs) {
    ratio_dpsi_ds = dpsi / ds;
  } else {
    // dpsi/ds ≈ 1/c0 + (s0*ds)/(2*c0^3)
    const double invc0  = 1.0 / c0;
    const double invc03 = invc0*invc0*invc0;
    ratio_dpsi_ds = invc0 + 0.5*s0*ds*invc03;
    dpsi = ds * ratio_dpsi_ds;
  }

  // s_perp = dpsi/kappa = dz*(dpsi/ds) (stable)
  const double s_perp = dz * ratio_dpsi_ds;

  // State update (stable)
  double denom_x = (c0 + c1);
  if (std::fabs(denom_x) < kTiny) denom_x = (denom_x >= 0.0 ? kTiny : -kTiny);

  x += dz * (s0 + s1) / denom_x;

  const double ny0 = ty0 * c0;
  y += ny0 * s_perp;

  const double tx1 = s1 * invC1;
  const double R   = c0 * invC1;
  const double ty1 = ty0 * R;

  // Store state
  mP[0] = static_cast<float>(x);
  mP[1] = static_cast<float>(y);
  mP[2] = static_cast<float>(tx1);
  mP[3] = static_cast<float>(ty1);
  // q/p unchanged
  mZ = znew;

  // ---- Jacobian ----
  // dkappa/du
  const double dk_dq  = K * (A * invD);
  const double dk_dty = K * q0 * (ty0 * (invA * invD));
  const double dk_dtx = K * q0 * (tx0 * (invA*invD - A*invD3));

  const double ds0_dtx = invD3;

  const double ds1_dtx = ds0_dtx + dz * dk_dtx;
  const double ds1_dty =          dz * dk_dty;
  const double ds1_dq  =          dz * dk_dq;

  const double dc0_dtx = -tx0 * invD3;

  const double s1_over_c1 = s1 * invC1;
  const double dc1_dtx = -s1_over_c1 * ds1_dtx;
  const double dc1_dty = -s1_over_c1 * ds1_dty;
  const double dc1_dq  = -s1_over_c1 * ds1_dq;

  // tx1 derivatives
  const double f22 = invC13 * ds1_dtx;
  const double f23 = invC13 * ds1_dty;
  const double f24 = invC13 * ds1_dq;

  // x1 derivatives (use series if small bend)
  double f02, f03, f04;
  const double kd = std::fabs(kappa * dz);
  if (kd < kSmallBend) {
    const double half_dz2 = 0.5 * dz * dz;
    f02 = dz + half_dz2 * dk_dtx;
    f03 =      half_dz2 * dk_dty;
    f04 =      half_dz2 * dk_dq;
  } else {
    const double invK  = 1.0 / kappa;
    const double invK2 = invK * invK;
    const double c0mc1 = (c0 - c1);
    f02 = ((dc0_dtx - dc1_dtx) * invK) - (c0mc1 * invK2) * dk_dtx;
    f03 = ((0.0     - dc1_dty) * invK) - (c0mc1 * invK2) * dk_dty;
    f04 = ((0.0     - dc1_dq ) * invK) - (c0mc1 * invK2) * dk_dq;
  }

  // dpsi derivatives
  const double dpsi_dtx = invC1 * ds1_dtx - (1.0 / c0) * ds0_dtx;
  const double dpsi_dty = invC1 * ds1_dty;
  const double dpsi_dq  = invC1 * ds1_dq;

  // y derivatives (small bend: straight-line leading order)
  double f12, f13, f14;
  if (kd < kSmallBend) {
    f12 = 0.0;
    f13 = dz;
    f14 = 0.0;
  } else {
    const double invK  = 1.0 / kappa;
    const double invK2 = invK * invK;

    const double dsperp_dtx = (dpsi_dtx * kappa - dpsi * dk_dtx) * invK2;
    const double dsperp_dty = (dpsi_dty * kappa - dpsi * dk_dty) * invK2;
    const double dsperp_dq  = (dpsi_dq  * kappa - dpsi * dk_dq ) * invK2;

    f13 = c0 * s_perp + ny0 * dsperp_dty;
    f12 = (ty0 * dc0_dtx) * s_perp + ny0 * dsperp_dtx;
    f14 = ny0 * dsperp_dq;
  }

  // Option B: ty1 = ty0*(c0/c1)
  const double f33 = R + ty0 * c0 * s1 * ds1_dty * invC13;
  const double f34 = ty0 * c0 * s1 * ds1_dq  * invC13;
  const double f32 = ty0 * (dc0_dtx * c1 - c0 * dc1_dtx) * invC12;

  // ---- Covariance transport ----
  transportCovarianceFullJac(f02,f03,f04, f12,f13,f14, f22,f23,f24, f32,f33,f34);
  return true;
}



// ------------------------------------------------------------------
// 1) Separate covariance transport: C <- F C F^T
// Jacobian sparse structure (Option B):
// row0: [1,0,f02,f03,f04]
// row1: [0,1,f12,f13,f14]
// row2: [0,0,f22,f23,f24]
// row3: [0,0,f32,f33,f34]
// row4: [0,0, 0,  0, 1 ]
// ------------------------------------------------------------------
void NA6PTrackParCov::transportCovarianceFullJac(double f02, double f03, double f04,
                                                 double f12, double f13, double f14,
                                                 double f22, double f23, double f24,
                                                 double f32, double f33, double f34)
{
  // Load old covariance to a full symmetric double matrix (5x5).
  // This avoids overwrite hazards and improves numerical stability.
  double C[5][5];
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      const double v = static_cast<double>( getCovMatElem(i,j) );
      C[i][j] = v;
      C[j][i] = v;
    }
  }

  // Extract u-subblock (tx,ty,q) where indices are 2,3,4
  const double C22 = C[2][2], C23 = C[2][3], C24 = C[2][4];
  const double C33 = C[3][3], C34 = C[3][4];
  const double C44 = C[4][4];

  auto dot_u = [](double a2, double a3, double a4, double b2, double b3, double b4) -> double {
    return a2*b2 + a3*b3 + a4*b4;
  };

  auto aT_Cuu = [&](double a2, double a3, double a4, double& v2, double& v3, double& v4) {
    // v = a^T * Cuu  (components correspond to columns tx,ty,q)
    v2 = a2*C22 + a3*C23 + a4*C24;
    v3 = a2*C23 + a3*C33 + a4*C34;
    v4 = a2*C24 + a3*C34 + a4*C44;
  };

  double Cp[5][5] = {}; // new covariance

  // ---- block for rows/cols (2,3,4): depends only on Cuu ----
  {
    double v2,v3,v4;

    // C'22
    aT_Cuu(f22,f23,f24, v2,v3,v4);
    Cp[2][2] = dot_u(f22,f23,f24, v2,v3,v4);

    // C'33
    aT_Cuu(f32,f33,f34, v2,v3,v4);
    Cp[3][3] = dot_u(f32,f33,f34, v2,v3,v4);

    // C'23
    double w2,w3,w4;
    aT_Cuu(f22,f23,f24, w2,w3,w4);
    Cp[2][3] = dot_u(f32,f33,f34, w2,w3,w4);

    // C'24, C'34, C'44
    Cp[2][4] = f22*C24 + f23*C34 + f24*C44;
    Cp[3][4] = f32*C24 + f33*C34 + f34*C44;
    Cp[4][4] = C44;
  }

  // ---- x' row/col ----
  {
    const double C00 = C[0][0];
    const double C02 = C[0][2], C03 = C[0][3], C04 = C[0][4];

    double v2,v3,v4;
    aT_Cuu(f02,f03,f04, v2,v3,v4);

    // C'00
    Cp[0][0] =
        C00
      + 2.0*(f02*C02 + f03*C03 + f04*C04)
      + dot_u(f02,f03,f04, v2,v3,v4);

    // C'02
    const double cov_x_txp = f22*C02 + f23*C03 + f24*C04;
    Cp[0][2] =
        cov_x_txp
      + (f02*(f22*C22 + f23*C23 + f24*C24)
       + f03*(f22*C23 + f23*C33 + f24*C34)
       + f04*(f22*C24 + f23*C34 + f24*C44));

    // C'03
    const double cov_x_typ = f32*C02 + f33*C03 + f34*C04;
    Cp[0][3] =
        cov_x_typ
      + (f02*(f32*C22 + f33*C23 + f34*C24)
       + f03*(f32*C23 + f33*C33 + f34*C34)
       + f04*(f32*C24 + f33*C34 + f34*C44));

    // C'04
    Cp[0][4] = C04 + (f02*C24 + f03*C34 + f04*C44);
  }

  // ---- y' row/col ----
  {
    const double C11 = C[1][1];
    const double C12 = C[1][2], C13 = C[1][3], C14 = C[1][4];

    double v2,v3,v4;
    aT_Cuu(f12,f13,f14, v2,v3,v4);

    // C'11
    Cp[1][1] =
        C11
      + 2.0*(f12*C12 + f13*C13 + f14*C14)
      + dot_u(f12,f13,f14, v2,v3,v4);

    // C'12
    const double cov_y_txp = f22*C12 + f23*C13 + f24*C14;
    Cp[1][2] =
        cov_y_txp
      + (f12*(f22*C22 + f23*C23 + f24*C24)
       + f13*(f22*C23 + f23*C33 + f24*C34)
       + f14*(f22*C24 + f23*C34 + f24*C44));

    // C'13
    const double cov_y_typ = f32*C12 + f33*C13 + f34*C14;
    Cp[1][3] =
        cov_y_typ
      + (f12*(f32*C22 + f33*C23 + f34*C24)
       + f13*(f32*C23 + f33*C33 + f34*C34)
       + f14*(f32*C24 + f33*C34 + f34*C44));

    // C'14
    Cp[1][4] = C14 + (f12*C24 + f13*C34 + f14*C44);
  }

  // ---- C'01 ----
  {
    const double C01 = C[0][1];
    const double C02 = C[0][2], C03 = C[0][3], C04 = C[0][4];
    const double C12 = C[1][2], C13 = C[1][3], C14 = C[1][4];

    double v2,v3,v4;
    aT_Cuu(f02,f03,f04, v2,v3,v4);

    Cp[0][1] =
        C01
      + (f02*C12 + f03*C13 + f04*C14)
      + (f12*C02 + f13*C03 + f14*C04)
      + dot_u(f12,f13,f14, v2,v3,v4);
  }
  
  // Symmetrize for storage convenience
  Cp[1][0]=Cp[0][1];
  Cp[2][0]=Cp[0][2]; Cp[2][1]=Cp[1][2];
  Cp[3][0]=Cp[0][3]; Cp[3][1]=Cp[1][3]; Cp[3][2]=Cp[2][3];
  Cp[4][0]=Cp[0][4]; Cp[4][1]=Cp[1][4]; Cp[4][2]=Cp[2][4]; Cp[4][3]=Cp[3][4];

  // Store back packed lower triangle
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      setCovMatElem(i,j, static_cast<float>(Cp[i][j]));
    }
  }
  checkCorrelations();

}

/*
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
*/

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

/* // old update
// Kalman update with 2D measurement (xm, ym), and covariance (sx2, sxy, sy2)
bool NA6PTrackParCov::update(float xm, float ym, float sx2, float sxy, float sy2)
{
  auto sav = *this; // TOREM

  double dx = xm - getX(), dy = ym - getY(); // innovation (residual)
  const double Cxx = mC[kXX], Cyx = mC[kYX], Cyy = mC[kYY], CTxX = mC[kTxX], CTxY = mC[kTxY], CTxTx = mC[kTxTx], CTyX = mC[kTyX], CTyY = mC[kTyY];
  const double CTyTx = mC[kTyTx], CTyTy = mC[kTyTy], CQpx = mC[kQ2PX], CQpy = mC[kQ2PY], CQpTx = mC[kQ2PTx], CQpTy = mC[kQ2PTy], CQpQp = mC[kQ2PQ2P];

  // innovation covariance S = H C H^T + R (2x2) ---
  //  if (sx2 / Cxx < 1e-4) {
  //    sx2 = Cxx * 1e-4;
  //  }
  //  if (sy2 / Cyy < 1e-4) {
  //    sy2 = Cyy * 1e-4;
  //  }
  //  const double S00 = Cxx + sx2, S01 = Cyx + sxy, S11 = Cyy + sy2, detS = S00 * S11 - S01 * S01;
  //  if (std::abs(detS) < kTiny) { // Singular or ill-conditioned S: skip update
  //    return false;
  //  }
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
   
 */


// Kalman update with 2D measurement (xm, ym), and covariance (sx2, sxy, sy2)
bool NA6PTrackParCov::update(float xm, float ym, float sx2, float sxy, float sy2)
{
  auto sav = *this; // TOREM
  // Innovation (residual)
  const double dx = static_cast<double>(xm) - static_cast<double>(getX());
  const double dy = static_cast<double>(ym) - static_cast<double>(getY());

  // ------------------------------------------------------------------
  // Build full covariance matrix in double and (if needed) repair it to
  // be numerically SPD by adding a tiny diagonal "jitter".
  //
  // Motivation: even very small negative eigenvalues (from transport
  // roundoff) can yield negative diagonal elements after an extremely
  // strong measurement update (tiny R).
  // ------------------------------------------------------------------
  double C[5][5];
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      C[i][j] = C[j][i] = static_cast<double>(getCovMatElem(i, j));
    }
  }

  auto make_spd_by_jitter = [&](double M[5][5]) {
    // Try a simple Cholesky test. If it fails, add diagonal jitter and retry.
    // This is cheap for 5x5 and avoids costly eigen-decomposition.
    const double maxDiag = std::max({std::fabs(M[0][0]), std::fabs(M[1][1]), std::fabs(M[2][2]), std::fabs(M[3][3]), std::fabs(M[4][4]), 1.0});
    const double eps = 1e-12 * maxDiag; // pivot threshold
    double jitter = 0.0;

    for (int it = 0; it < 6; ++it) {
      // Cholesky attempt: M + jitter*I
      double L[5][5] = {{0.0}};
      bool ok = true;
      for (int i = 0; i < 5 && ok; ++i) {
        for (int j = 0; j <= i; ++j) {
          double s = M[i][j];
          if (i == j) {
            s += jitter;
          }
          for (int k = 0; k < j; ++k) {
            s -= L[i][k] * L[j][k];
          }
          if (i == j) {
            if (s <= eps) {
              ok = false;
              break;
            }
            L[i][j] = std::sqrt(s);
          } else {
            L[i][j] = s / L[j][j];
          }
        }
      }
      if (ok) {
        // Apply jitter to the matrix in-place if it was needed.
        if (jitter != 0.0) {
          for (int i = 0; i < 5; ++i) {
            M[i][i] += jitter;
          }
        }
        return;
      }
      // Increase jitter (geometric growth). Start very small relative to scale.
      jitter = (jitter == 0.0) ? (1e-10 * maxDiag) : (jitter * 10.0);
    }

    // Last resort: ensure strictly positive diagonals.
    for (int i = 0; i < 5; ++i) {
      if (M[i][i] < eps) {
        M[i][i] = eps;
      }
    }
  };

  //  make_spd_by_jitter(C);

  // Pull the needed covariance elements from the repaired matrix.
  const double Cxx   = C[0][0];
  const double Cyx   = C[1][0];
  const double Cyy   = C[1][1];
  const double CTxX  = C[2][0];
  const double CTxY  = C[2][1];
  const double CTyX  = C[3][0];
  const double CTyY  = C[3][1];
  const double CQpx  = C[4][0];
  const double CQpy  = C[4][1];

  // Innovation covariance S = H C H^T + R (2x2)
  const double R00 = static_cast<double>(sx2);
  const double R01 = static_cast<double>(sxy);
  const double R11 = static_cast<double>(sy2);

  const double S00 = Cxx + R00;
  const double S01 = Cyx + R01;
  const double S11 = Cyy + R11;
  const double detS = S00 * S11 - S01 * S01;
  if (std::abs(detS) < kTiny) { // Singular or ill-conditioned S: skip update
    return false;
  }
  const double invDetS = 1. / detS;
  const double invS00 = S11 * invDetS;
  const double invS11 = S00 * invDetS;
  const double invS01 = -S01 * invDetS;

  // Kalman gain K (5x2). H selects (x,y) so K = C H^T S^{-1}.
  const double K0_0 = Cxx  * invS00 + Cyx  * invS01;
  const double K0_1 = Cxx  * invS01 + Cyx  * invS11;
  const double K1_0 = Cyx  * invS00 + Cyy  * invS01;
  const double K1_1 = Cyx  * invS01 + Cyy  * invS11;
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

  // ------------------------------------------------------------------
  // Joseph-stabilized covariance update (PSD-safe in finite precision):
  //   C+ = (I - K H) C (I - K H)^T + K R K^T
  // Here H selects (x,y), so KH only modifies columns 0 and 1.
  // We compute in double and store back to float.
  // ------------------------------------------------------------------

  // Note: we reuse the already repaired C[5][5] built above.

  // A = I - K H, only columns 0 and 1 differ from identity
  double A[5][5] = {0};
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      A[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }
  A[0][0] -= K0_0;  A[0][1] -= K0_1;
  A[1][0] -= K1_0;  A[1][1] -= K1_1;
  A[2][0] -= K2_0;  A[2][1] -= K2_1;
  A[3][0] -= K3_0;  A[3][1] -= K3_1;
  A[4][0] -= K4_0;  A[4][1] -= K4_1;

  // Temp = A*C
  double T[5][5];
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      T[i][j] = 0.0;
      for (int k = 0; k < 5; ++k) {
        T[i][j] += A[i][k] * C[k][j];
      }
    }
  }

  // Cp = Temp * A^T
  double Cp[5][5];
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      Cp[i][j] = 0.0;
      for (int k = 0; k < 5; ++k) {
        Cp[i][j] += T[i][k] * A[j][k]; // A^T[k][j] = A[j][k]
      }
    }
  }

  // Add K R K^T
  // KR rows: (K_i0*K_i1) times R
  auto add_KRKt = [&](int i, int j, double Ki0, double Ki1, double Kj0, double Kj1) {
    return Ki0 * (R00 * Kj0 + R01 * Kj1) + Ki1 * (R01 * Kj0 + R11 * Kj1);
  };

  const double K0[2] = {K0_0, K0_1};
  const double K1[2] = {K1_0, K1_1};
  const double K2[2] = {K2_0, K2_1};
  const double K3[2] = {K3_0, K3_1};
  const double K4[2] = {K4_0, K4_1};
  const double* Krows[5] = {K0, K1, K2, K3, K4};

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      Cp[i][j] += add_KRKt(i, j, Krows[i][0], Krows[i][1], Krows[j][0], Krows[j][1]);
    }
  }

  // Enforce symmetry and write back packed lower triangle
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      const double sym = 0.5 * (Cp[i][j] + Cp[j][i]);
      setCovMatElem(i,j, static_cast<float>(sym));
    }
  }
  checkCorrelations();
  return true;
}

void NA6PTrackParCov::resetCovariance(float s2)
{
  constexpr float
    kCX2max = 5 * 5,        // SigmaX<=5cm
    kCY2max = 5 * 5,        // SigmaY<=5cm
    kCTx2max = 0.3 * 0.3,   // tgX<=1
    kCTy2max = 0.3 * 0.3,   // tgY<=1
    kC1P2max = 30 * 30,   // Sigma1/P<=30 1/GeV
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
