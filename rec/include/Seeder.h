// NA6PCCopyright

#ifndef NA6P_SEEDER_H
#define NA6P_SEEDER_H

#include "NA6PTrackPar.h"
#include "Propagator.h"
#include <array>

class Seeder
{
 public:
  using dimf3 = std::array<float, 3>;

  static bool create(NA6PTrackPar& seed, bool reverse, const dimf3& p0, const dimf3& p1, const dimf3& p2, float improvePrecThreshod = 0.01f);
  static NA6PTrackPar create(bool reverse, const dimf3& p0, const dimf3& p1, const dimf3& p2, float improvePrecThreshod = 0.01f)
  {
    NA6PTrackPar seed;
    create(seed, reverse, p0, p1, p2, improvePrecThreshod);
    return seed;
  }

  static void printResid(const NA6PTrackPar& seed, const dimf3& pos);

 private:
  struct CircularSeed {
    double xc = 0.;
    double zc = 0.;
    double r2 = 0.;
    float tgl = 0.f;
    float dydz = 0.f;
    float y0 = 0.f;
    float x0 = 0.f;
    float z0 = 0.f;
    float byAv = 0.f;
    int xSign = 1;

    // point along the seed at Z (almost)
    dimf3 evalFast(float z) const { return {float(xc - xSign * std::sqrt(r2 - (z - zc) * (z - zc))), y0 + dydz * (z - z0), z}; }

    inline float integrateBRange(float z) const
    {
      const float mid = 0.5f * (z0 + z), half = 0.5f * (z - z0);
      float s = 0.f;
      for (int i = 0; i < 3; ++i) {
        s += WG3[i] * Propagator::Instance()->getBy(evalFast(mid + half * XG3[i]));
      }
      return half * s;
    }
  };

  struct QuadNode {
    float w = 0.f;
    float A = 0.f; // A(z) = int_z0^z By(ref(z')) dz'
  };

  struct QuadNodes {
    static constexpr int N = 3;
    std::array<QuadNode, N> node{};
    auto integral() const { return node[0].w * node[0].A + node[1].w * node[1].A + node[2].w * node[2].A; }
  };

  static CircularSeed buildCircularSeed(bool reverse, const dimf3& p0, const dimf3& df01, const dimf3& df02); // p0, p1-p0, p2-p0

  static QuadNodes makeNodes(CircularSeed& cSeed, float zA, float zB)
  {
    const float mid = 0.5f * (zA + zB), half = 0.5f * (zB - zA);
    auto qnode = [mid, half, &cSeed](int i) { return QuadNode{half * WG3[i], cSeed.integrateBRange(mid + half * XG3[i])}; };
    return {qnode(0), qnode(1), qnode(2)};
  }

  static void addResidualJacobian(const QuadNodes& nodes, float u, float beta, float& R, float& Ju, float& Jb);
  static void addResidualJacobianR(const QuadNodes& nodes, float u, float beta, float& R);

  static constexpr float XG3[3] = {-0.7745966692414834f, 0.f, 0.7745966692414834f};
  static constexpr float WG3[3] = {0.5555555555555556f, 0.8888888888888888f, 0.5555555555555556f};

  ClassDefNV(Seeder, 1);
};

#endif
