// NA6PCCopyright

/// \file Propagator
/// \brief Singleton class for track propagation routines
/// \author ruben.shahoyan@cern.ch

#ifndef NA6P_PROPAGATOR_
#define NA6P_PROPAGATOR_

#include "NA6PTrackParCov.h"
#include "MagneticField.h"
#include <array>
#include <string>

struct MatBudget {
  float meanXRho = 0.f; ///< mean density, g/cm^3 time length
  float meanX2X0 = 0.f; ///< fraction of radiaton lenght
  ClassDefNV(MatBudget, 1);
};

class Propagator
{
 public:
  static Propagator* Instance(bool uninitialized = false)
  {
    static Propagator instance(uninitialized);
    return &instance;
  }

  static constexpr float MAX_STEP = 2.0f;
  enum class MatCorrType : int {
    USEMatCorrNONE, // flag to not use material corrections
    USEMatCorrTGeo, // flag to use TGeo for material queries
    //    USEMatCorrLUT // TODO
  }; // flag to use LUT for material queries (user must provide a pointer

  struct PropOpt {
    MatCorrType matCorr{MatCorrType::USEMatCorrTGeo};
    float maxStep{MAX_STEP};
    bool fixCorrelations{false};
  };

  bool propagateToZ(NA6PTrackParCov& track, float z, const PropOpt& opt = {MatCorrType::USEMatCorrTGeo, MAX_STEP, false}) const;
  bool propagateToZ(NA6PTrackPar& track, float z, const PropOpt& opt = {MatCorrType::USEMatCorrTGeo, MAX_STEP, false}) const;

  Propagator(Propagator const&) = delete;
  Propagator(Propagator&&) = delete;
  Propagator& operator=(Propagator const&) = delete;
  Propagator& operator=(Propagator&&) = delete;

  MatBudget getMeanMaterial(const float* start, const float* end) const;
  MatBudget getMeanMaterial(const std::array<float, 3>& start, const std::array<float, 3>& end) const { return getMeanMaterial(start.data(), end.data()); }

  void getMeanMaterialBudgetFromGeom(const float* start, const float* end, float* mparam) const;
  void getMeanMaterialBudgetFromGeom(const std::array<float, 3>& start, const std::array<float, 3>& end, std::array<float, 7>& mparam) const { getMeanMaterialBudgetFromGeom(start.data(), end.data(), mparam.data()); }

  void updateField();
  bool hasMagFieldSet() const { return mField != nullptr; }
  bool hasGeometryLoaded() const;

  template <typename T = float>
  void getFieldXYZ(const T* xyz, T* bxyz) const
  {
    mField->getField(xyz, bxyz);
  }
  template <typename T = float>
  void getFieldXYZ(const std::array<T, 3> xyz, std::array<T, 3>& bxyz) const
  {
    getFieldXYZ(xyz.data(), bxyz.data());
  }
  template <typename T = float>
  float getBy(const T* xyz) const
  {
    T b[3];
    mField->getField(xyz, b);
    return b[1];
  }
  template <typename T = float>
  float getBy(const std::array<T, 3> xyz) const
  {
    return getBy(xyz.data());
  }

 private:
  static constexpr float Epsilon = 0.00001; // precision of propagation to Z

  Propagator(bool uninitialized = false);
  ~Propagator() = default;

  MagneticField* mField = nullptr; ///< External nominal field map

  ClassDefNV(Propagator, 0);
};

inline MatBudget Propagator::getMeanMaterial(const float* start, const float* end) const
{
  float mat[7];
  getMeanMaterialBudgetFromGeom(start, end, mat);
  return {mat[0] * mat[4], mat[1]};
}

#endif
