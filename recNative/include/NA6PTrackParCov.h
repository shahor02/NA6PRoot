// NA6PCCopyright
#ifndef NA6P_TRACKPARCOV_H
#define NA6P_TRACKPARCOV_H

#include "NA6PTrackPar.h"

/*
  Add covariance matrix and related methods to NA6PTrackPar
*/

class NA6PTrackParCov : public NA6PTrackPar
{
 public:
  using NA6PTrackPar::NA6PTrackPar;

  NA6PTrackParCov() = default;
  NA6PTrackParCov(const NA6PTrackParCov& src) = default;
  NA6PTrackParCov(float z, const std::array<float, 5>& par, std::array<float, 15> cov) : NA6PTrackPar(z, par), mC{cov} {}
  ~NA6PTrackParCov() = default;
  NA6PTrackParCov& operator=(const NA6PTrackParCov&) = default;

  void setCov(const std::array<float, 15>& c) { mC = c; }
  const std::array<float, 15>& getCov() const { return mC; }
  std::array<float, 15>& getCov() { return mC; }

  // --- propagation in dipole field B = (0, By, 0) ---
  bool propagate(float z, float by);

  bool update(const float xm, const float ym, const float sx2, const float sxy, const float sy2);
  bool update(const std::array<float, 2>& meas, const std::array<float, 3>& cov);

  std::string asString() const;

 protected:
  void propagateCov(double f02, double f03, double f04, double f13, double f22, double f23, double f24);

  std::array<float, 15> mC{}; // covariance matrix in lower triangle representation

  ClassDefNV(NA6PTrackParCov, 1);
};

inline bool NA6PTrackParCov::update(const std::array<float, 2>& meas, const std::array<float, 3>& cov)
{
  return update(meas[kX], meas[kY], cov[kXX], cov[kYX], cov[kYY]);
}

#endif // NA6P_TRACKPARCOV_H
