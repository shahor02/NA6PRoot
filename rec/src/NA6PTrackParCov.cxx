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
  const auto dz = znew - mZ;
  if (std::abs(dz) == 1e-6f) {
    return true;
  }
  const float K = kB2C * b_y;

  // State in prec_t
  const prec_t tx0 = mP[2];
  const prec_t ty0 = mP[3];
  
  const auto D2 = getPxz2Pz2<prec_t>();                // D^2 of the into == pxz^2/pz^2 == 1.0 + tx0*tx0
  const auto D = std::sqrt(D2);                        // D of the into 
  const auto A = std::sqrt(D2 + getTy()*getTy());      // A of the intro = p/pz = sqrt(1 + tx0*tx0 + ty0*ty0)
  const prec_t KQ0 = K * getQ2P();
  const prec_t kappa = KQ0 * (A / D);

  // Straight-line / negligible curvature branch
  // Triggered by very small |kappa| (i.e. small K or small q/p).
  if (std::abs(kappa) < kSmallKappa) {
    // State update: x += tx*dz ; y += ty*dz ; tx,ty,q unchanged
    mP[0] += mP[2] * dz;
    mP[1] += mP[3] * dz;
    mZ = znew;    
    // Jacobian for straight line:
    const prec_t f02 = dz, f03 = 0.0, f04 = 0.0;
    const prec_t f12 = 0.0, f13 = dz,  f14 = 0.0;
    const prec_t f22 = 1.0, f23 = 0.0, f24 = 0.0;
    const prec_t f32 = 0.0, f33 = 1.0, f34 = 0.0;    
    transportCovarianceFullJac(f02,f03,f04, f12,f13,f14, f22,f23,f24, f32,f33,f34);
    return true;
  }

  // ---- Dipole propagation ----
  const prec_t bend = kappa * dz;
  const prec_t invD  = 1.0 / D, invD2 = invD*invD, invD3 = invD2*invD, invA  = 1.0 / A;
  const prec_t s0 = tx0 * invD, c0 = invD; // sin(psi0) and cos(psi0) = 1/D
  prec_t s1 = s0 + bend; // sin(psi1)
  if (std::abs(s1) >= 1.) {
    return false;
  }
  const prec_t c1 = std::sqrt((1.0 - s1)*(1.0 + s1)); // cos(psi1)
  const prec_t invC1  = 1.0 / c1, invC12 = invC1 * invC1, invC13 = invC12 * invC1;

  // robust dpsi using atan2(s,c)
  const prec_t ds = s1 - s0;
  prec_t dpsi, dpsi2ds;
  if (std::abs(ds) > kSmallDs) {
    const prec_t psi0 = std::atan2(s0, c0);
    const prec_t psi1 = std::atan2(s1, c1);
    dpsi = psi1 - psi0;
    dpsi2ds = dpsi / ds;
  } else {  // dpsi/ds ≈ 1/c0 + (s0*ds)/(2*c0^3)
    const prec_t invc0  = 1.0 / c0;
    const prec_t invc03 = invc0*invc0*invc0;
    dpsi2ds = invc0 + 0.5*s0*ds*invc03;
    dpsi = ds * dpsi2ds;
  }
  const prec_t s_perp = dz * dpsi2ds;  // s_perp = dpsi/kappa = dz*(dpsi/ds) (stable)

  // State update (stable)
  prec_t denom_x = (c0 + c1);
  if (std::abs(denom_x) < kTiny) {
    denom_x = (denom_x >= 0.0 ? kTiny : -kTiny);
  }
  const prec_t ny0 = ty0 * c0, R = c0 * invC1;
  mP[0] += dz * (s0 + s1) / denom_x;
  mP[1] += ny0 * s_perp;
  mP[2] = s1 * invC1;
  mP[3] = ty0 * R;
  // q/p unchanged
  mZ = znew;

  // ---- Jacobian ----
  // dkappa/du
  const prec_t dk_dq  = K * (A * invD);
  const prec_t dk_dty = KQ0 * (ty0 * (invA * invD));
  const prec_t dk_dtx = KQ0 * (tx0 * (invA*invD - A*invD3));

  const prec_t ds0_dtx = invD3;

  const prec_t ds1_dtx = ds0_dtx + dz * dk_dtx;
  const prec_t ds1_dty =          dz * dk_dty;
  const prec_t ds1_dq  =          dz * dk_dq;

  const prec_t dc0_dtx = -tx0 * invD3;

  const prec_t s1_over_c1 = s1 * invC1;
  const prec_t dc1_dtx = -s1_over_c1 * ds1_dtx;
  const prec_t dc1_dty = -s1_over_c1 * ds1_dty;
  const prec_t dc1_dq  = -s1_over_c1 * ds1_dq;

  // tx1 derivatives
  const prec_t f22 = invC13 * ds1_dtx;
  const prec_t f23 = invC13 * ds1_dty;
  const prec_t f24 = invC13 * ds1_dq;

  // dpsi derivatives
  const prec_t dpsi_dtx = invC1 * ds1_dtx - (1.0 / c0) * ds0_dtx;
  const prec_t dpsi_dty = invC1 * ds1_dty;
  const prec_t dpsi_dq  = invC1 * ds1_dq;
  
  prec_t f02, f03, f04;  // x1 derivatives (use series if small bend)
  prec_t f12, f13, f14;  // y derivatives (small bend: straight-line leading order)
  bool smallBend = std::abs(bend) < static_cast<prec_t>(kSmallBend);
  if (smallBend) {
    const prec_t half_dz2 = 0.5 * dz * dz;
    f02 = dz + half_dz2 * dk_dtx;
    f03 =      half_dz2 * dk_dty;
    f04 =      half_dz2 * dk_dq;

    f12 = 0.0;
    f13 = dz;
    f14 = 0.0;    
  } else {
    const prec_t invK  = 1.0 / kappa;
    const prec_t invK2 = invK * invK;
    const prec_t c0mc1toK2 = (c0 - c1) * invK2;
    f02 = ((dc0_dtx - dc1_dtx) * invK) - c0mc1toK2 * dk_dtx;
    f03 = -( dc1_dty * invK + c0mc1toK2 * dk_dty);
    f04 = -( dc1_dq  * invK + c0mc1toK2 * dk_dq);

    const prec_t dsperp_dtx = (dpsi_dtx * kappa - dpsi * dk_dtx) * invK2;
    const prec_t dsperp_dty = (dpsi_dty * kappa - dpsi * dk_dty) * invK2;
    const prec_t dsperp_dq  = (dpsi_dq  * kappa - dpsi * dk_dq ) * invK2;
    f13 = c0 * s_perp + ny0 * dsperp_dty;
    f12 = (ty0 * dc0_dtx) * s_perp + ny0 * dsperp_dtx;
    f14 = ny0 * dsperp_dq;  
  }

  // ty1 = ty0*(c0/c1)
  const prec_t f33 = R + ty0 * c0 * s1 * ds1_dty * invC13;
  const prec_t f34 = ty0 * c0 * s1 * ds1_dq  * invC13;
  const prec_t f32 = ty0 * (dc0_dtx * c1 - c0 * dc1_dtx) * invC12;

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
void NA6PTrackParCov::transportCovarianceFullJac(prec_t f02, prec_t f03, prec_t f04,
                                                 prec_t f12, prec_t f13, prec_t f14,
                                                 prec_t f22, prec_t f23, prec_t f24,
                                                 prec_t f32, prec_t f33, prec_t f34)
{
  // Load old covariance to a full symmetric prec_t matrix (5x5).
  // This avoids overwrite hazards and improves numerical stability.

  const prec_t C00 = mC[kXX], C01 = mC[kYX], C02 = mC[kTxX], C03 = mC[kTyX], C04 = mC[kQ2PX];
  const prec_t C11 = mC[kYY], C12 = mC[kTxY], C13 = mC[kTyY], C14 = mC[kQ2PY];
  const prec_t C22 = mC[kTxTx], C23 = mC[kTyTx], C24 = mC[kQ2PTx];
  const prec_t C33 = mC[kTyTy], C34 = mC[kQ2PTy], C44 = mC[kQ2PQ2P];
  
  auto dot_u = [](prec_t a2, prec_t a3, prec_t a4, prec_t b2, prec_t b3, prec_t b4) -> prec_t {
    return a2*b2 + a3*b3 + a4*b4;
  };

  auto aT_Cuu = [&](prec_t a2, prec_t a3, prec_t a4, prec_t& v2, prec_t& v3, prec_t& v4) {    // v = a^T * Cuu  (components correspond to columns tx,ty,q)
    v2 = a2*C22 + a3*C23 + a4*C24;
    v3 = a2*C23 + a3*C33 + a4*C34;
    v4 = a2*C24 + a3*C34 + a4*C44;
  };

  // ---- block for rows/cols (2,3,4): depends only on Cuu ----
  prec_t v2,v3,v4;  
  aT_Cuu(f22,f23,f24, v2,v3,v4);
  mC[kTxTx] = dot_u(f22,f23,f24, v2,v3,v4);

  aT_Cuu(f32,f33,f34, v2,v3,v4);
  mC[kTyTy] = dot_u(f32,f33,f34, v2,v3,v4);

  aT_Cuu(f22,f23,f24, v2,v3,v4);
  mC[kTyTx] = dot_u(f32,f33,f34, v2,v3,v4);
  mC[kQ2PTx] = f22*C24 + f23*C34 + f24*C44;
  mC[kQ2PTy] = f32*C24 + f33*C34 + f34*C44;
  
  // ---- x' row/col ----
  aT_Cuu(f02,f03,f04, v2,v3,v4);
  mC[kXX] += 2.0*(f02*C02 + f03*C03 + f04*C04) + dot_u(f02,f03,f04, v2,v3,v4);

  const prec_t cov_x_txp = f22*C02 + f23*C03 + f24*C04;
  mC[kTxX] = cov_x_txp + (f02*(f22*C22 + f23*C23 + f24*C24) + f03*(f22*C23 + f23*C33 + f24*C34) + f04*(f22*C24 + f23*C34 + f24*C44));

  const prec_t cov_x_typ = f32*C02 + f33*C03 + f34*C04;
  mC[kTyX] = cov_x_typ + (f02*(f32*C22 + f33*C23 + f34*C24) + f03*(f32*C23 + f33*C33 + f34*C34) + f04*(f32*C24 + f33*C34 + f34*C44));

  mC[kQ2PX] +=  (f02*C24 + f03*C34 + f04*C44);

  // ---- y' row/col ----
  aT_Cuu(f12,f13,f14, v2,v3,v4);
  mC[kYY] += 2.0*(f12*C12 + f13*C13 + f14*C14)  + dot_u(f12,f13,f14, v2,v3,v4);

  const prec_t cov_y_txp = f22*C12 + f23*C13 + f24*C14;
  mC[kTxY] = cov_y_txp + (f12*(f22*C22 + f23*C23 + f24*C24) + f13*(f22*C23 + f23*C33 + f24*C34) + f14*(f22*C24 + f23*C34 + f24*C44));

  const prec_t cov_y_typ = f32*C12 + f33*C13 + f34*C14;
  mC[kTyY] = cov_y_typ + (f12*(f32*C22 + f33*C23 + f34*C24) + f13*(f32*C23 + f33*C33 + f34*C34) + f14*(f32*C24 + f33*C34 + f34*C44));
  
  mC[kQ2PY] += (f12*C24 + f13*C34 + f14*C44);

  aT_Cuu(f02,f03,f04, v2,v3,v4);
  mC[kYX] += (f02*C12 + f03*C13 + f04*C14) + (f12*C02 + f13*C03 + f14*C04) + dot_u(f12,f13,f14, v2,v3,v4);
  checkCorrelations();

}

// calculate chi2 between the track and the cluster
float NA6PTrackParCov::getPredictedChi2(float xm, float ym, float sx2, float syx, float sy2) const
{
  auto exx = static_cast<prec_t>(getSigmaX2()) + static_cast<prec_t>(sx2);
  auto eyx = static_cast<prec_t>(getSigmaYX()) + static_cast<prec_t>(syx);
  auto eyy = static_cast<prec_t>(getSigmaY2()) + static_cast<prec_t>(sy2);
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
  // Innovation (residual)
  const prec_t dx = static_cast<prec_t>(xm) - static_cast<prec_t>(getX());
  const prec_t dy = static_cast<prec_t>(ym) - static_cast<prec_t>(getY());

  // ------------------------------------------------------------------
  // Build full covariance matrix in prec_t and (if needed) repair it to
  // be numerically SPD by adding a tiny diagonal "jitter".
  //
  // Motivation: even very small negative eigenvalues (from transport
  // roundoff) can yield negative diagonal elements after an extremely
  // strong measurement update (tiny R).
  // ------------------------------------------------------------------
  prec_t C[5][5];
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      C[i][j] = C[j][i] = static_cast<prec_t>(getCovMatElem(i, j));
    }
  }

  auto make_spd_by_jitter = [&](prec_t M[5][5]) {
    // Try a simple Cholesky test. If it fails, add diagonal jitter and retry.
    // This is cheap for 5x5 and avoids costly eigen-decomposition.
    const prec_t maxDiag = std::max({std::abs(M[0][0]), std::abs(M[1][1]), std::abs(M[2][2]), std::abs(M[3][3]), std::abs(M[4][4]), 1.0});
    const prec_t eps = 1e-12 * maxDiag; // pivot threshold
    prec_t jitter = 0.0;

    for (int it = 0; it < 6; ++it) {
      // Cholesky attempt: M + jitter*I
      prec_t L[5][5] = {{0.0}};
      bool ok = true;
      for (int i = 0; i < 5 && ok; ++i) {
        for (int j = 0; j <= i; ++j) {
          prec_t s = M[i][j];
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
  const prec_t Cxx   = C[0][0];
  const prec_t Cyx   = C[1][0];
  const prec_t Cyy   = C[1][1];
  const prec_t CTxX  = C[2][0];
  const prec_t CTxY  = C[2][1];
  const prec_t CTyX  = C[3][0];
  const prec_t CTyY  = C[3][1];
  const prec_t CQpx  = C[4][0];
  const prec_t CQpy  = C[4][1];

  // Innovation covariance S = H C H^T + R (2x2)
  const prec_t R00 = static_cast<prec_t>(sx2);
  const prec_t R01 = static_cast<prec_t>(sxy);
  const prec_t R11 = static_cast<prec_t>(sy2);

  const prec_t S00 = Cxx + R00;
  const prec_t S01 = Cyx + R01;
  const prec_t S11 = Cyy + R11;
  const prec_t detS = S00 * S11 - S01 * S01;
  if (std::abs(detS) < static_cast<prec_t>(kTiny)) { // Singular or ill-conditioned S: skip update
    return false;
  }
  const prec_t invDetS = 1. / detS;
  const prec_t invS00 = S11 * invDetS;
  const prec_t invS11 = S00 * invDetS;
  const prec_t invS01 = -S01 * invDetS;

  // Kalman gain K (5x2). H selects (x,y) so K = C H^T S^{-1}.
  const prec_t K0_0 = Cxx  * invS00 + Cyx  * invS01;
  const prec_t K0_1 = Cxx  * invS01 + Cyx  * invS11;
  const prec_t K1_0 = Cyx  * invS00 + Cyy  * invS01;
  const prec_t K1_1 = Cyx  * invS01 + Cyy  * invS11;
  const prec_t K2_0 = CTxX * invS00 + CTxY * invS01;
  const prec_t K2_1 = CTxX * invS01 + CTxY * invS11;
  const prec_t K3_0 = CTyX * invS00 + CTyY * invS01;
  const prec_t K3_1 = CTyX * invS01 + CTyY * invS11;
  const prec_t K4_0 = CQpx * invS00 + CQpy * invS01;
  const prec_t K4_1 = CQpx * invS01 + CQpy * invS11;

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
  // We compute in prec_t and store back to float.
  // ------------------------------------------------------------------

  // Note: we reuse the already repaired C[5][5] built above.

  // A = I - K H, only columns 0 and 1 differ from identity
  prec_t A[5][5] = {0};
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
  prec_t T[5][5];
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      T[i][j] = 0.0;
      for (int k = 0; k < 5; ++k) {
        T[i][j] += A[i][k] * C[k][j];
      }
    }
  }

  // Cp = Temp * A^T
  prec_t Cp[5][5];
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
  auto add_KRKt = [&](int i, int j, prec_t Ki0, prec_t Ki1, prec_t Kj0, prec_t Kj1) {
    return Ki0 * (R00 * Kj0 + R01 * Kj1) + Ki1 * (R01 * Kj0 + R11 * Kj1);
  };

  const prec_t K0[2] = {K0_0, K0_1};
  const prec_t K1[2] = {K1_0, K1_1};
  const prec_t K2[2] = {K2_0, K2_1};
  const prec_t K3[2] = {K3_0, K3_1};
  const prec_t K4[2] = {K4_0, K4_1};
  const prec_t* Krows[5] = {K0, K1, K2, K3, K4};

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      Cp[i][j] += add_KRKt(i, j, Krows[i][0], Krows[i][1], Krows[j][0], Krows[j][1]);
    }
  }

  // Enforce symmetry and write back packed lower triangle
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      const prec_t sym = 0.5 * (Cp[i][j] + Cp[j][i]);
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
