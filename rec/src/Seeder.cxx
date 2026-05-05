// NA6PCCopyright

#include "Seeder.h"
#include "NA6PLine.h"

Seeder::CircularSeed Seeder::buildCircularSeed(bool reverse, const dimf3& p0, const dimf3& df01, const dimf3& df02)
{
  const double rhs1 = df01[0] * double(p0[0] + p0[0] + df01[0]) + df01[2] * double(p0[2] + p0[2] + df01[2]);
  const double rhs2 = df02[0] * double(p0[0] + p0[0] + df02[0]) + df02[2] * double(p0[2] + p0[2] + df02[2]);
  const double det = double(df01[0]) * df02[2] - double(df01[2]) * df02[0];
  CircularSeed c;
  if (std::abs(det) > 1.e-14f) {
    const double detI = 0.5 / det;
    c.xc = (df02[2] * rhs1 - df01[2] * rhs2) * detI;
    c.zc = (-df02[0] * rhs1 + df01[0] * rhs2) * detI;
    const auto dx0 = p0[0] - c.xc, dz0 = p0[2] - c.zc;
    c.r2 = dx0 * dx0 + dz0 * dz0;
    c.xSign = (dx0 >= 0.f) ? -1 : 1;
  }
  float sgn = 1; // reverse ? 1.f : -1.f;
  c.tgl = (reverse ? 0.5 : -0.5) * (df01[1] / std::hypot(df01[0], df01[2]) + (df02[1] - df01[1]) / std::hypot(df02[0] - df01[0], df02[2] - df01[2]));
  c.dydz = df02[1] / df02[2];
  c.x0 = p0[0];
  c.y0 = p0[1];
  c.z0 = p0[2];
  return c;
}

void Seeder::addResidualJacobian(const Seeder::QuadNodes& nodes, float u, float beta, float& R, float& Ju, float& Jb)
{
  for (const auto& n : nodes.node) {
    const float wv = u + beta * n.A;
    const float oneMinus = std::max(1.e-12f, (1.f - wv) * (1.f + wv));
    const float invSqrt = 1.f / std::sqrt(oneMinus);

    const float f = wv * invSqrt;       // dx/dz
    const float g = invSqrt / oneMinus; // d(dx/dz)/dw

    R += n.w * f;
    Ju += n.w * g;
    Jb += n.w * n.A * g;
  }
}

void Seeder::addResidualJacobianR(const Seeder::QuadNodes& nodes, float u, float beta, float& R)
{
  for (const auto& n : nodes.node) {
    const float wv = u + beta * n.A;
    const float oneMinus = std::max(1.e-12f, (1.f - wv) * (1.f + wv));
    const float f = wv / std::sqrt(oneMinus); // dx/dz
    R += n.w * f;
  }
}

bool Seeder::create(NA6PTrackPar& seed, bool reverse, const dimf3& p0, const dimf3& p1, const dimf3& p2, float improvePrecThreshod)
{
  // try to improve the precision if circular seed estimate is worse than improvePrecThreshod
  constexpr int maxNewtonIter = 2;
  const auto df01 = NA6PLine::getDiff(p0, p1), df02 = NA6PLine::getDiff(p0, p2);
  CircularSeed cSeed = buildCircularSeed(reverse, p0, df01, df02);
  seed.setXYZ(p0);
  seed.setTy(cSeed.tgl);
  if (cSeed.r2 == 0 && std::abs(Propagator::Instance()->getBy(p1)) < 0.001f) { // straight line in XZ
    const float tgp = df02[0] / df02[2];
    seed.setTx(tgp / std::sqrt(1.f + tgp * tgp));
    seed.setQ2Pxz(1.f); // 1./kMostProbablePt RSTODO
    return true;
  }
  cSeed.byAv = cSeed.integrateBRange(p2[2]) / df02[2]; // Average By on the same circular-X / polyline-Y reference path.
  float rInv = std::sqrt(1. / cSeed.r2);
  seed.setTx(cSeed.xSign * rInv * (p0[2] - cSeed.zc));
  seed.setQ2Pxz(rInv * cSeed.xSign / (cSeed.byAv * NA6PTrackPar::kB2C));
  std::cout << "built seed " << seed.asString() << "\n";
  if (improvePrecThreshod < 0.f) { // check and improvements are not requested
    return true;
  }
  float beta = seed.getQ2Pxz() * NA6PTrackPar::kB2C, u = seed.getTx();
  auto nodes01 = makeNodes(cSeed, p0[2], p1[2]);
  auto nodes12 = makeNodes(cSeed, p1[2], p2[2]);
  {
    //    std::cout << "Before:\n";
    //    printResid(seed, p1);
    //    printResid(seed, p2);
  }

  // Exact-model residual at ITS seed:  Rj = x_model(zj; uITS, betaITS) - dxj
  float R1 = -df01[0], R2 = -df02[0];
  addResidualJacobianR(nodes01, u, beta, R1);
  addResidualJacobianR(nodes01, u, beta, R2);
  addResidualJacobianR(nodes12, u, beta, R2);
  // Bending lever arms: F1 = int_{z0}^{z1} A(z) dz and  F2 = int_{z0}^{z2} A(z) dz
  float F1 = nodes01.integral(), F2 = F1 + nodes12.integral();
  // S has the same dimension as beta
  const float denF = std::sqrt(F1 * F1 + F2 * F2);
  const float numR = std::sqrt(R1 * R1 + R2 * R2);
  const float S = numR / std::max(denF, 1.e-12f);
  const float epsBeta = S / std::max(std::abs(beta), 1.e-12f);
  /*
  printf("ITS residual consistency:\n");
  printf("  R1      = %+e\n", R1);
  printf("  R2      = %+e\n", R2);
  printf("  F1      = %+e\n", F1);
  printf("  F2      = %+e\n", F2);
  printf("  S       = %+e\n", S);
  printf("  epsBeta = %+e\n", epsBeta);
  */
  if (epsBeta < improvePrecThreshod) {
    // printf("ITS curvature likely good to ~1%%: skip Newton\n");
    return true;
  } else {
    printf("ITS curvature inconsistency %f >1%%: Newton may help | %s\n", epsBeta, seed.asString().c_str());
  }

  for (int it = 0; it < maxNewtonIter; ++it) {
    float R1 = -df01[0], Ju1 = 0.f, Jb1 = 0.f;
    addResidualJacobian(nodes01, u, beta, R1, Ju1, Jb1);
    float R2 = -df02[0], Ju2 = 0.f, Jb2 = 0.f;
    addResidualJacobian(nodes01, u, beta, R2, Ju2, Jb2);
    addResidualJacobian(nodes12, u, beta, R2, Ju2, Jb2);

    const float norm = R1 * R1 + R2 * R2;
    if (norm < 1.e-10f) {
      break;
    }
    const float Jdet = Ju1 * Jb2 - Ju2 * Jb1;
    const float Jscale = std::abs(Ju1 * Jb2) + std::abs(Ju2 * Jb1) + 1.e-30f;
    if (std::abs(Jdet) < 1.e-7f * Jscale) {
      break;
    }
    const float du = (-R1 * Jb2 + R2 * Jb1) / Jdet;
    const float db = (-Ju1 * R2 + Ju2 * R1) / Jdet;

    bool accepted = false;
    for (float lambda = 1.f; lambda >= 1.f / 64.f; lambda *= 0.5f) {
      const float uTry = u + lambda * du;
      const float bTry = beta + lambda * db;
      // if (!validState(uTry, bTry, nodes01, nodes12)) continue;
      float tR1 = -df01[0], tJu1 = 0.f, tJb1 = 0.f;
      addResidualJacobian(nodes01, uTry, bTry, tR1, tJu1, tJb1);
      float tR2 = -df02[0], tJu2 = 0.f, tJb2 = 0.f;
      addResidualJacobian(nodes01, uTry, bTry, tR2, tJu2, tJb2);
      addResidualJacobian(nodes12, uTry, bTry, tR2, tJu2, tJb2);
      const float tryNorm = tR1 * tR1 + tR2 * tR2;
      if (tryNorm <= norm) {
        u = uTry;
        beta = bTry;
        accepted = true;
        break;
      }
    }
    if (!accepted) {
      break;
    }
  }
  if (std::abs(u) < 1.f) {
    seed.setTx(u);
    seed.setQ2Pxz(beta / NA6PTrackPar::kB2C);
    //    std::cout << "After:\n";
    //    printResid(seed, p1);
    //    printResid(seed, p2);
  } else {
    seed.invalidate();
    return false;
  }
  return true;
}
