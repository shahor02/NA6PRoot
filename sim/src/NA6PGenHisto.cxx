// NA6PCCopyright

#include "NA6PGenHisto.h"

#include "NA6PGenCutParam.h"
#include "NA6PMCStack.h"

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom.h>

#include <cmath>

NA6PGenHisto::NA6PGenHisto(const std::string& name, int pdg, float mult, bool isPoisson, TH2* ptYHisto, float ycm)
  : NA6PGenerator(name), mPDGCode(pdg), mMult(mult), mPoisson(isPoisson), mPtYHisto(ptYHisto), mYCM(ycm)
{
}

NA6PGenHisto::NA6PGenHisto(const std::string& name, int pdg, float mult, bool isPoisson, TH1* ptHisto, TH1* yHisto, float ycm)
  : NA6PGenerator(name), mPDGCode(pdg), mMult(mult), mPoisson(isPoisson), mPtHisto(ptHisto), mYHisto(yHisto), mYCM(ycm)
{
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

  NA6PGenerator::init();
}

void NA6PGenHisto::generate()
{
  generatePrimaryVertex();

  if (!mPtYHisto && !(mPtHisto && mYHisto)) {
    LOGP(fatal, "Generator {} is not initialized", getName());
  }

  const auto& cutParam = NA6PGenCutParam::Instance();
  const float mass = TDatabasePDG::Instance()->GetParticle(mPDGCode)->Mass();

  const int nPart = mPoisson ? gRandom->Poisson(mMult) : TMath::Nint(mMult);

  for (int i = 0; i < nPart; ++i) {
    double pt = 0;
    double y = 0;
    if (mPtYHisto)
    {
      mPtYHisto->GetRandom2(pt, y);
    }
    else {
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
    getStack()->PushTrack(true, -1, mPDGCode,
                          px, py, pz, en,
                          getOriginX(), getOriginY(), getOriginZ(),
                          0., 0., 0., 0.,
                          TMCProcess::kPPrimary, dummy, 1., 0);
  }

  auto mcHead = getStack()->getEventHeader();
  std::string info = getName() + "_pdg" + std::to_string(mPDGCode) + "x" + std::to_string(getMultiplicity()) + "_TH";
  mcHead->getGenHeaders().emplace_back(nPart, 0, mcHead->getNPrimaries(), 0, info);
  mcHead->incNPrimaries(nPart);
}