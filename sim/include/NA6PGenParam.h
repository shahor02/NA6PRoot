// NA6PCCopyright

#ifndef NA6P_GENPARAM_H
#define NA6P_GENPARAM_H

#include <TF1.h>
#include "NA6PGenerator.h"

class NA6PGenParam : public NA6PGenerator
{
 public:
  NA6PGenParam() = default;
  NA6PGenParam(const std::string& name, int pdg, float mult, const std::string& parTrans, const std::string& parLong,
               float ptmin, float ptmax, int npartr, float lmin, float lmax, int nparl, bool istransPt = true, bool islongY = true, bool isPoisson = true);
  ~NA6PGenParam() override = default;
  void generate() override;

  void init() override;

  auto getPDGCode() const { return mPDGCode; }
  auto getFunTrans() const { return mFunTrans.get(); }
  auto getFunLong() const { return mFunLong.get(); }
  float getMultiplicity() const { return mMult; }
  bool isLongY() const { return mLongIsY; }
  bool isTransPt() const { return mTransIsPt; }
  bool isPoisson() const { return mPoisson; }

  void setParametersTrans(const std::vector<float>& v);
  void setParametersLong(const std::vector<float>& v);
  void setPDGCode(int c) { mPDGCode = c; }
  void setFunTrans(const std::string& fname, float vmin, float vmax, int np);
  void setFunLong(const std::string& fname, float vmin, float vmax, int np);
  void setFunTrans(const std::string& name, const std::string& formula, float vmin, float vmax);
  void setFunLong(const std::string& name, const std::string& formula, float vmin, float vmax);
  void setMultiplicity(float m) { mMult = m; }
  void setPoisson(bool v) { mPoisson = v; }
  void setLongIsY(bool v) { mLongIsY = v; }
  void setTransIsPt(bool v) { mTransIsPt = v; }

 protected:
  static constexpr int ParamsSetBit = 0x1 << 14;
  std::unique_ptr<TF1> mFunTrans; // holder of transverse distribution (pT or mT)
  std::unique_ptr<TF1> mFunLong;  // holder or longitudinal distribuition (Y or eta)
  int mPDGCode = 0;               // particle type
  float mMult = 0;                // mean mult. of particles to generate
  bool mPoisson = true;
  bool mLongIsY = true;
  bool mTransIsPt = true;

  ClassDefOverride(NA6PGenParam, 1); // Square box (in variables of NA6PGenCutParam) random generator
};

#endif