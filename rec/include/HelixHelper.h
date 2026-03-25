// NA6PCCopyright

/// \file HelixHelper.h
/// \brief Helper classes for helical tracks (ZX-plane bending) manipulations (based on the O2 HelixHelper)
/// \author ruben.shahoyan@cern.ch

#ifndef _NA6P_HELIX_HELPER_
#define _NA6P_HELIX_HELPER_

#include "CommonConstants/MathConstants.h"
#include "MathUtils/Utils.h"
#include "MathUtils/Primitive2D.h"

//__________________________________________________________
//< crossing coordinates of 2 circles
template <typename T = float>
struct CrossInfo {
  using CInfo = std::array<T, 3>; // track XZ circle params: R, Xc, Zc
  enum { kR,
         kXC,
         kZC };
  static constexpr T MaxDistZXDef = T(10);
  T xDCA[2] = {};
  T zDCA[2] = {};
  int nDCA = 0;

  int circlesCrossInfo(const CInfo& trax0, const CInfo& trax1, T maxDistZX = MaxDistZXDef, bool isCollinear = false)
  {
    const auto& trcA = trax0[kR] > trax1[kR] ? trax0 : trax1; // designate the largest circle as A
    const auto& trcB = trax0[kR] > trax1[kR] ? trax1 : trax0;
    nDCA = 0;
    auto zDist = trcB[kZC] - trcA[kZC], xDist = trcB[kXC] - trcA[kXC];
    auto dist2 = zDist * zDist + xDist * xDist, dist = std::sqrt(dist2), rsum = trcA[kR] + trcB[kR];
    if (dist < 1e-12) {
      return nDCA; // circles are concentric?
    }
    if (dist > rsum) { // circles don't touch, chose a point in between
      // the parametric equation of lines connecting the centers is
      // x = x0 + t/dist * (x1-x0), z = z0 + t/dist * (z1-z0)
      if (dist - rsum > maxDistZX) { // too large distance
        return nDCA;
      }
      notTouchingZX(dist, zDist, xDist, trcA, trcB[kR], isCollinear);
    } else if (auto dfr = dist + trcB[kR] - trcA[kR]; dfr < 0.) { // the small circle is nestled into large one w/o touching
      if (dfr > -maxDistZX) {
        // select the point of closest approach of 2 circles
        notTouchingZX(dist, zDist, xDist, trcA, -trcB[kR], isCollinear);
      } else {
        return nDCA;
      }
    } else { // 2 intersection points
      if (isCollinear) {
        /// collinear tracks, e.g. electrons from photon conversion
        /// if there are 2 crossings of the circle it is better to take
        /// a weighted average of the crossing points as a radius
        T r2r = trcA[kR] + trcB[kR], r1_r = trcA[kR] / r2r, r2_r = trcB[kR] / r2r;
        zDCA[0] = r2_r * trcA[kZC] + r1_r * trcB[kZC];
        xDCA[0] = r2_r * trcA[kXC] + r1_r * trcB[kXC];
        nDCA = 1;
      } else if (o2::gpu::GPUCommonMath::Abs(zDist) < o2::gpu::GPUCommonMath::Abs(xDist)) {
        // to simplify calculations, we move to new frame x->x+Xc0, y->y+Yc0, so that
        // the 1st one is centered in origin
        T a = (trcA[kR] * trcA[kR] - trcB[kR] * trcB[kR] + dist2) / (T(2) * xDist), b = -zDist / xDist, ab = a * b, bb = b * b;
        T det = ab * ab - (T(1) + bb) * (a * a - trcA[kR] * trcA[kR]);
        if (det > 0.) {
          det = std::sqrt(det);
          zDCA[0] = (-ab + det) / (T(1) + b * b);
          xDCA[0] = a + b * zDCA[0] + trcA[kXC];
          zDCA[0] += trcA[kZC];
          zDCA[1] = (-ab - det) / (T(1) + b * b);
          xDCA[1] = a + b * zDCA[1] + trcA[kXC];
          zDCA[1] += trcA[kZC];
          nDCA = 2;
        } else { // due to the finite precision the det<=0, i.e. the circles are barely touching, fall back to this special case
          notTouchingZX(dist, zDist, xDist, trcA, trcB[kR]);
        }
      } else {
        T a = (trcA[kR] * trcA[kR] - trcB[kR] * trcB[kR] + dist2) / (T(2) * zDist), b = -xDist / zDist, ab = a * b, bb = b * b;
        T det = ab * ab - (T(1) + bb) * (a * a - trcA[kR] * trcA[kR]);
        if (det > 0.) {
          det = std::sqrt(det);
          xDCA[0] = (-ab + det) / (T(1) + bb);
          zDCA[0] = a + b * xDCA[0] + trcA[kZC];
          xDCA[0] += trcA[kXC];
          xDCA[1] = (-ab - det) / (T(1) + bb);
          zDCA[1] = a + b * xDCA[1] + trcA[kZC];
          xDCA[1] += trcA[kXC];
          nDCA = 2;
        } else { // due to the finite precision the det<=0, i.e. the circles are barely touching, fall back to this special case
          notTouchingZX(dist, zDist, xDist, trcA, trcB[kR]);
        }
      }
    }
    return nDCA;
  }

  void notTouchingZX(T dist, T zDist, T xDist, const CInfo& trcA, T rBSign, bool isCollinear = false)
  {
    if (isCollinear) {
      /// for collinear tracks it is better to take
      /// a weighted average of the crossing points as a radius
      T r2r = trcA[kR] + rBSign, r1_r = trcA[kR] / r2r, r2_r = rBSign / r2r;
      zDCA[0] = r2_r * trcA[kZC] + r1_r * (zDist + trcA[kZC]);
      xDCA[0] = r2_r * trcA[kXC] + r1_r * (xDist + trcA[kXC]);
    } else {
      // fast method to calculate DCA between 2 circles, assuming that they don't touch each outer:
      // the parametric equation of lines connecting the centers is x = xA + t/dist * zDist, z = zA + t/dist * zDist
      // with zA,xA being the center of the circle A ( = trcA[kZC], trcA[kXC] ), zDist = trcB[kZC] - trcA[kZC] ...
      // There are 2 special cases:
      // (a) small circle is inside the large one: provide rBSign as -trcB[kR]
      // (b) circle are side by side: provide rBSign as trcB[kR]
      auto t2d = (dist + trcA[kR] - rBSign) / dist;
      zDCA[0] = trcA[kZC] + T(0.5) * (zDist * t2d);
      xDCA[0] = trcA[kXC] + T(0.5) * (xDist * t2d);
    }
    nDCA = 1;
  }

  template <typename T>
  GPUd() int linesCrossInfo(const CInfo& trax0, const T& tr0,
                            const CInfo& trax1, const T& tr1, T maxDistZX = MaxDistZXDef)
  {
    /// closest approach of 2 straight lines
    ///  TrackParam propagation can be parameterized in lab in a form
    ///  z(t) = z + t
    ///  x(t) = y + t*snp/csp
    ///  y(t) = y + t*tgl/csp
    ///  where t is the z-step, x,y,z are reference track coordinates filled by CInfo for straight line tracks.
    ///  snp = pX/pZX, csp = cos(asin(snp)), tgl = pY/pZX
    ///
    ///  Therefore, for the parametric track equation in lab 3D we have (wrt tracking-X increment t)
    ///  x(t) = x + t Kx;  Kx = 1
    ///  y(t) = y + t Ky;  Ky = snp/csp
    ///  z(t) = z + t Kz;  Kz = tgl / csp
    ///  Note that Kx^2 + Ky^2 + Kz^2 = (1+tgl^2) / csp^2
    nDCA = 0;
    T dz = trax1[kZC] - trax0[kZC]; // for straight line CInfo stores lab coordinates at referene point!!!
    T dx = trax1[kXC] - trax0[kXC]; //
    T dy = tr1.getY() - tr0.getY();
    const auto csp0i2 = T(1) / tr0.getCosPsi2(), csp0i = std::sqrt(csp0i2), tgp0 = tr0.getSinPsi() * csp0i;
    const auto csp1i2 = T(1) / tr1.getCosPsi2(), csp1i = std::sqrt(csp1i2), tgp1 = tr1.getSinPsi() * csp1i;
    const T kZ0 = T(1), kX0 = tgp0, kY0 = tr0.getTy() * csp0i;
    const T kZ1 = T(1), kX1 = tgp1, kY1 = tr1.getTy() * csp1i;
    /// Minimize |vecL1 - vecL0|^2 wrt t0 and t1: point of closest approach
    /// Leads to system
    /// A Dx = B with Dx = {dx0, dx1}
    /// with A =
    ///  |      kZ0^2+kX0^2+kY0^2     -(kZ0*kx1+kX0*ky1+kY0*kz1) | =  (1+tgl0^2) / csp0^2           ....
    ///  | -(kZ0*kx1+kX0*ky1+kY0*kz1)     kZ0^2+kX0^2+kY0^2      |     .....                   (1+tgl1^2) / csp1^2
    /// and B = {(dx KZ0 + dy KX0 + dz KY0), -(dx Kx1 + dy Ky1 + dz Kz1) }
    ///
    T a00 = (T(1) + tr0.getTy() * tr0.getTy()) * csp0i2, a11 = (T(1) + tr1.getTy() * tr1.getTy()) * csp1i2, a01 = -(kZ0 * kx1 + kX0 * ky1 + kY0 * kz1);
    T b0 = dx * kZ0 + dy * kX0 + dz * kY0, b1 = -(dx * kx1 + dy * ky1 + dz * kz1);
    T det = a00 * a11 - a01 * a01, det0 = b0 * a11 - b1 * a01, det1 = a00 * b1 - a01 * b0;
    if (std::sqrt(det) > o2::constants::math::Almost0) {
      auto detI = T(1) / det;
      auto t0 = det0 * detI;
      auto t1 = det1 * detI;
      T addx0 = kZ0 * t0, addy0 = kX0 * t0, addx1 = kx1 * t1, addy1 = ky1 * t1;
      dx += addx1 - addx0; // recalculate ZX distance at DCA
      dy += addy1 - addy0;
      if (dx * dx + dy * dy > maxDistZX * maxDistZX) {
        return nDCA;
      }
      zDCA[0] = (trax0[kZC] + addx0 + trax1[kZC] + addx1) * T(0.5);
      xDCA[0] = (trax0[kXC] + addy0 + trax1[kXC] + addy1) * T(0.5);
      nDCA = 1;
    }
    return nDCA;
  }

  template <typename T>
  GPUd() int circleLineCrossInfo(const CInfo& trax0, const T& tr0,
                                 const CInfo& trax1, const T& tr1, T maxDistZX = MaxDistZXDef)
  {
    /// closest approach of line and circle
    ///  TrackParam propagation can be parameterized in lab in a form
    ///  xLab(t) = (x*cosAlp - y*sinAlp) + t*(cosAlp - sinAlp* snp/csp) = xLab0 + t*(cosAlp - sinAlp* snp/csp)
    ///  yLab(t) = (x*sinAlp + y*cosAlp) + t*(sinAlp + cosAlp* snp/csp) = yLab0 + t*(sinAlp + cosAlp* snp/csp)
    ///  zLab(t) = z + t * tgl / csp = zLab0 + t * tgl / csp
    ///  where t is the x-step in the track alpha-frame, xLab,yLab,zLab are reference track coordinates in lab
    ///  frame (filled by CInfo for straight line tracks).
    ///
    ///  Therefore, for the parametric track equation in lab 3D we have (wrt tracking-X increment t)
    ///  xL(t) = xL + t Kx;  Kx = (cosAlp - sinAlp* snp/csp)
    ///  yL(t) = yL + t Ky;  Ky = (sinAlp + cosAlp* snp/csp)
    ///  zL(t) = zL + t Kz;  Kz = tgl / csp
    ///  Note that Kx^2 + Ky^2  = 1 / csp^2

    const auto& traxH = trax0[kR] > trax1[kR] ? trax0 : trax1; // circle (for the line rC is set to 0)
    const auto& traxL = trax0[kR] > trax1[kR] ? trax1 : trax0; // line
    const auto& trcL = trax0[kR] > trax1[kR] ? tr1 : tr0;      // track of the line

    // solve quadratic equation of line crossing the circle
    T dx = traxL[kZC] - traxH[kZC]; // X distance between the line lab reference and circle center
    T dy = traxL[kXC] - traxH[kXC]; // Y...
    // t^2(kx^2+ky^2) + 2t(dx*kx+dy*ky) + dx^2 + dy^2 - r^2 = 0
    auto cspi2 = T(1) / trcL.getCosPsi2(); // 1 / csp^2 == kx^2 +  ky^2
    auto cspi = std::sqrt(cspi2);
    auto tgp = trcL.getSinPsi() * cspi;
    T kx = traxL.c - traxL.s * tgp;
    T ky = traxL.s + traxL.c * tgp;
    double dk = dx * kx + dy * ky;
    double det = dk * dk - cspi2 * (dx * dx + dy * dy - traxH[kR] * traxH[kR]);
    if (det > 0) { // 2 crossings
      det = std::sqrt(det);
      T t0 = (-dk + det) * cspi2;
      T t1 = (-dk - det) * cspi2;
      zDCA[0] = traxL[kZC] + kx * t0;
      xDCA[0] = traxL[kXC] + ky * t0;
      zDCA[1] = traxL[kZC] + kx * t1;
      xDCA[1] = traxL[kXC] + ky * t1;
      nDCA = 2;
    } else {
      // there is no crossing, find the point of the closest approach on the line which is closest to the circle center
      T t = -dk * cspi2;
      T xL = traxL[kZC] + kx * t, yL = traxL[kXC] + ky * t; // point on the line, need to average with point on the circle
      T dxc = xL - traxH[kZC], dyc = yL - traxH[kXC], dist = std::sqrt(dxc * dxc + dyc * dyc);
      if (dist - traxH[kR] > maxDistZX) { // too large distance
        return nDCA;
      }
      T drcf = traxH[kR] / dist; // radius / distance to circle center
      T xH = traxH[kZC] + dxc * drcf, yH = traxH[kXC] + dyc * drcf;
      zDCA[0] = (xL + xH) * T(0.5);
      xDCA[0] = (yL + yH) * T(0.5);
      nDCA = 1;
    }
    return nDCA;
  }

  template <typename T>
  GPUd() int set(const CInfo& trax0, const T& tr0, const CInfo& trax1, const T& tr1, T maxDistZX = MaxDistZXDef, bool isCollinear = false)
  {
    // calculate up to 2 crossings between 2 circles
    nDCA = 0;
    if (trax0[kR] > o2::constants::math::Almost0 && trax1[kR] > o2::constants::math::Almost0) { // both are not straight lines
      nDCA = circlesCrossInfo(trax0, trax1, maxDistZX, isCollinear);
    } else if (trax0[kR] < o2::constants::math::Almost0 && trax1[kR] < o2::constants::math::Almost0) { // both are straigt lines
      nDCA = linesCrossInfo(trax0, tr0, trax1, tr1, maxDistZX);
    } else {
      nDCA = circleLineCrossInfo(trax0, tr0, trax1, tr1, maxDistZX);
    }
    //
    return nDCA;
  }

  GPUdDefault() CrossInfo() = default;

  template <typename T>
  GPUd() CrossInfo(const CInfo& trax0, const T& tr0, const CInfo& trax1, const T& tr1, T maxDistZX = MaxDistZXDef, bool isCollinear = false)
  {
    set(trax0, tr0, trax1, tr1, maxDistZX, isCollinear);
  }
  ClassDefNV(CrossInfo, 1);
};

} // namespace track
} // namespace o2

#endif
