// NA6PCCopyright

#ifndef NA6P_GENHISTO_H
#define NA6P_GENHISTO_H

#include <TH1.h>
#include <TH2.h>
#include <memory>
#include <string>

#include "NA6PGenerator.h"

class NA6PGenHisto : public NA6PGenerator
{
 public:
  using NA6PGenerator::NA6PGenerator;
  NA6PGenHisto() = default;
  NA6PGenHisto(const std::string& name, int pdg, float mult, bool isPoisson = true, TH2* ptYHisto = nullptr, float ycm = 0.f);
  NA6PGenHisto(const std::string& name, int pdg, float mult, bool isPoisson = true, TH1* ptHisto = nullptr, TH1* yHisto = nullptr, float ycm = 0.f);
  ~NA6PGenHisto() override = default;

  void init() override;
  void generate() override;

  void setPDGCode(int c) { mPDGCode = c; }
  int getPDGCode() const { return mPDGCode; }

  void setMultiplicity(float m) { mMult = m; }
  float getMultiplicity() const { return mMult; }

  void setPoisson(bool v) { mPoisson = v; }
  bool isPoisson() const { return mPoisson; }

  void setPtYHistogram(TH2* histo);
  void setPtHistogram(TH1* histo);
  void setYHistogram(TH1* histo);
  const TH2* getPtYHistogram() const { return mPtYHisto; }
  const TH1* getPtHistogram() const { return mPtHisto; }
  const TH1* getYHistogram() const { return mYHisto; }
 protected:
  int mPDGCode = 0;
  float mMult = -1.f;
  bool mPoisson = true;
  float mYCM = 0.f;

  std::string mSourceFileName;
  std::string mSourceHistoName;
  TH2* mPtYHisto = nullptr;
  TH1* mPtHisto = nullptr;
  TH1* mYHisto = nullptr;

  ClassDefOverride(NA6PGenHisto, 1);
};

#endif