// NA6PCCopyright

#include "NA6PGenHisto.h"

#include "NA6PGenCutParam.h"
#include "NA6PMCStack.h"

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom.h>
#include <Math/Vector4D.h>
#include <Math/Boost.h>

#include <cmath>
#include <algorithm>

NA6PGenHisto::NA6PGenHisto(const std::string& name, int pdg, float mult, bool isPoisson, TH3* mPtYHisto, float ycm)
  : NA6PGenerator(name), mPDGCode(pdg), mMult(mult), mPoisson(isPoisson), mMPtYHisto(mPtYHisto), mYCM(ycm)
{
}

NA6PGenHisto::NA6PGenHisto(const std::string& name, int pdg, float mult, bool isPoisson, TH2* ptYHisto, float ycm)
  : NA6PGenerator(name), mPDGCode(pdg), mMult(mult), mPoisson(isPoisson), mPtYHisto(ptYHisto), mYCM(ycm)
{
}

NA6PGenHisto::NA6PGenHisto(const std::string& name, int pdg, float mult, bool isPoisson, TH1* ptHisto, TH1* yHisto, float ycm)
  : NA6PGenerator(name), mPDGCode(pdg), mMult(mult), mPoisson(isPoisson), mPtHisto(ptHisto), mYHisto(yHisto), mYCM(ycm)
{
}

void NA6PGenHisto::setMPtYHistogram(TH3* histo)
{
  mMPtYHisto = histo;
}

void NA6PGenHisto::setPtYHistogram(TH2* histo)
{
  mPtYHisto = histo;
}

void NA6PGenHisto::setPtHistogram(TH1* histo)
{
  mPtHisto = histo;
}

void NA6PGenHisto::setYHistogram(TH1* histo)
{
  mYHisto = histo;
}

void NA6PGenHisto::init()
{
  if (isInitDone()) {
    return;
  }
  if (!TDatabasePDG::Instance()->GetParticle(mPDGCode)) {
    LOGP(fatal, "Particle PDG={} is not defined", mPDGCode);
  }
  if (mMult <= 0.f) {
    LOGP(fatal, "Multiplicity for particle PDG={} is {}", mPDGCode, mMult);
  }

  if (mPrimaryDimuonDecay && (mPDGCode != 23 || !mMPtYHisto)) {
    LOGP(fatal, "Primary dimuon decay requires PDG 23 and a mass-pT-y histogram");
  }

  NA6PGenerator::init();
}

void NA6PGenHisto::generate()
{
  generatePrimaryVertex();

  if (!mMPtYHisto && !mPtYHisto && !(mPtHisto && mYHisto)) {
    LOGP(fatal, "Generator {} is not initialized", getName());
  }

  const auto& cutParam = NA6PGenCutParam::Instance();

  const int nPart = mPoisson ? gRandom->Poisson(mMult) : TMath::Nint(mMult);

  for (int i = 0; i < nPart; ++i) {
    double mass = TDatabasePDG::Instance()->GetParticle(mPDGCode)->Mass();
    double pt = 0;
    double y = 0;
    if (mMPtYHisto) {
      mMPtYHisto->GetRandom3(mass, pt, y);
    } else if (mPtYHisto) {
      mPtYHisto->GetRandom2(pt, y);
    } else {
      pt = mPtHisto->GetRandom();
      y = mYHisto->GetRandom();
    }
    y += mYCM;
    const double phi = gRandom->Uniform(0., TMath::TwoPi());
    const double px = pt * std::cos(phi);
    const double py = pt * std::sin(phi);
    const double mt = std::sqrt(pt * pt + mass * mass);
    const double pz = mt * std::sinh(y);
    const double en = std::sqrt(px * px + py * py + pz * pz + mass * mass);

    int dummy = 0;
    if (mPrimaryDimuonDecay) {
      const double muMass = TDatabasePDG::Instance()->GetParticle(13)->Mass();
      if (!std::isfinite(mass) || mass <= 2.0 * muMass || !std::isfinite(en)) {
        LOGP(fatal, "Invalid virtual photon mass/energy for dimuon decay: M={}, E={}", mass, en);
      }
      // Parent is generator truth only: never transport/decay it in VMC.
      int parentID = -1;
      getStack()->PushTrack(0, -1, 23, px, py, pz, en,
                           getOriginX(), getOriginY(), getOriginZ(),
                           0., 0., 0., 0., TMCProcess::kPPrimary, parentID, 1., 2);

      const double eStar = 0.5 * mass;
      const double pStar = std::sqrt((eStar - muMass) * (eStar + muMass));
      const double cosTheta = gRandom->Uniform(-1., 1.);
      const double sinTheta = std::sqrt(std::max(0., 1. - cosTheta * cosTheta));
      const double phiStar = gRandom->Uniform(0., TMath::TwoPi());
      ROOT::Math::PxPyPzEVector muMinus;
      muMinus.SetPxPyPzE(pStar * sinTheta * std::cos(phiStar), pStar * sinTheta * std::sin(phiStar), pStar * cosTheta, eStar);
      ROOT::Math::PxPyPzEVector muPlus;
      muPlus.SetPxPyPzE(-muMinus.Px(), -muMinus.Py(), -muMinus.Pz(), eStar);
      const ROOT::Math::Boost boostToLab(px / en, py / en, pz / en);
      muMinus = boostToLab(muMinus);
      muPlus = boostToLab(muPlus);
      for (int j = 0; j < 2; ++j) {
        const auto& mu = j == 0 ? muMinus : muPlus;
        getStack()->PushTrack(1, parentID, j == 0 ? 13 : -13,
                             mu.Px(), mu.Py(), mu.Pz(), mu.E(),
                             getOriginX(), getOriginY(), getOriginZ(),
                             0., 0., 0., 0., TMCProcess::kPPrimary, dummy, 1., 1);
      }
    } else {
      getStack()->PushTrack(true, -1, mPDGCode,
                            px, py, pz, en,
                            getOriginX(), getOriginY(), getOriginZ(),
                            0., 0., 0., 0.,
                            TMCProcess::kPPrimary, dummy, 1., 0);
    }
  }

  auto mcHead = getStack()->getEventHeader();
  std::string info = getName() + "_pdg" + std::to_string(mPDGCode) + "x" + std::to_string(getMultiplicity()) + "_TH";
  const int nPrimaries = mPrimaryDimuonDecay ? 3 * nPart : nPart;
  if (mPrimaryDimuonDecay) {
    info += "_primaryDimuon_isotropic";
  }
  mcHead->getGenHeaders().emplace_back(nPrimaries, 0, mcHead->getNPrimaries(), 0, info);
  mcHead->incNPrimaries(nPrimaries);
}
