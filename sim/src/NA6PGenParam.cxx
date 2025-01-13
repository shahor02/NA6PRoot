// NA6PCCopyright

#include "NA6PGenParam.h"
#include "NA6PGenCutParam.h"
#include "NA6PMCStack.h"
#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TMath.h>
#include <cmath>

NA6PGenParam::NA6PGenParam(const std::string& name, int pdg, float mult, const std::string& parTrans, const std::string& parLong,
                           float ptmin, float ptmax, int npartr, float lmin, float lmax, int nparl, bool istransPt, bool islongY, bool isPoisson) : NA6PGenerator(name), mPDGCode(pdg), mMult(mult), mPoisson(isPoisson), mTransIsPt(istransPt), mLongIsY(islongY)
{
  setFunTrans(parTrans, ptmin, ptmax, npartr);
  setFunLong(parLong, lmin, lmax, nparl);
}

void NA6PGenParam::setFunTrans(const std::string& fnameTrans, float ptmin, float ptmax, int npartr)
{
  // requires name of interpreted function
  if (!mTransIsPt) { // convert to Mt
    if (!TDatabasePDG::Instance()->GetParticle(mPDGCode)) {
      LOGP(fatal, "Particle PDG={} is not defined", mPDGCode);
    }
    auto m = TDatabasePDG::Instance()->GetParticle(mPDGCode)->Mass();
    ptmin = std::sqrt(m * m + ptmin * ptmin);
    ptmax = std::sqrt(m * m + ptmax * ptmax);
  }
  mFunTrans = std::make_unique<TF1>(fnameTrans.c_str(), ptmin, ptmax, npartr);
}

void NA6PGenParam::setFunLong(const std::string& fnameLong, float lmin, float lmax, int nparl)
{
  // requires name of interpreted function
  mFunLong = std::make_unique<TF1>(fnameLong.c_str(), lmin, lmax, nparl);
}

void NA6PGenParam::setFunTrans(const std::string& fname, const std::string& formulaTrans, float ptmin, float ptmax)
{
  // requires formula!
  if (!mTransIsPt) { // convert to Mt
    if (!TDatabasePDG::Instance()->GetParticle(mPDGCode)) {
      LOGP(fatal, "Particle PDG={} is not defined", mPDGCode);
    }
    auto m = TDatabasePDG::Instance()->GetParticle(mPDGCode)->Mass();
    ptmin = std::sqrt(m * m + ptmin * ptmin);
    ptmax = std::sqrt(m * m + ptmax * ptmax);
  }
  mFunTrans = std::make_unique<TF1>(fname.c_str(), formulaTrans.c_str(), ptmin, ptmax);
}

void NA6PGenParam::setFunLong(const std::string& fname, const std::string& formulaLong, float lmin, float lmax)
{
  // requires formula!
  mFunLong = std::make_unique<TF1>(fname.c_str(), formulaLong.c_str(), lmin, lmax);
}

void NA6PGenParam::init()
{
  if (isInitDone()) {
    return;
  }
  if (!TDatabasePDG::Instance()->GetParticle(mPDGCode)) {
    LOGP(fatal, "Particle PDG={} is not defined", mPDGCode);
  }
  if (mMult <= 0) {
    LOGP(fatal, "Multiplicity for particle PDG={} is {}", mPDGCode, mMult);
  }
  auto getLbl = [this](int i) { return i == 0 ? (mTransIsPt ? "Pt" : "Mt") : (mLongIsY ? "Y" : "Eta"); };
  TF1* ff[2] = {mFunTrans.get(), mFunLong.get()};
  for (int i = 0; i < 2; i++) {
    auto f = ff[i];
    if (!f) {
      LOGP(fatal, "{}-parametrization for PDG={} is not set", getLbl(i), mPDGCode);
    }
    if (f->GetNpar() > 0 && !f->TestBit(ParamsSetBit)) {
      LOGP(fatal, "{}-parametrization {} for PDG={} needs {} parameters which are not set", getLbl(i), f->GetName(), mPDGCode, f->GetNpar());
    }
  }

  NA6PGenerator::init();
}

void NA6PGenParam::setParametersTrans(const std::vector<float>& v)
{
  int npf = mFunTrans->GetNpar();
  if (npf != int(v.size())) {
    LOGP(fatal, "{}-parametrization {} needs {} parameters, {} provided", mTransIsPt ? "Pt" : "Mt", mFunTrans->GetName(), npf, v.size());
  }
  for (int i = 0; i < npf; i++) {
    mFunTrans->SetParameter(i, v[i]);
  }
  mFunTrans->SetBit(ParamsSetBit);
}

void NA6PGenParam::setParametersLong(const std::vector<float>& v)
{
  int npf = mFunLong->GetNpar();
  if (npf != int(v.size())) {
    LOGP(fatal, "{}-parametrization {} needs {} parameters, {} provided", mLongIsY ? "Y" : "Eta", mFunTrans->GetName(), npf, v.size());
  }
  for (int i = 0; i < npf; i++) {
    mFunLong->SetParameter(i, v[i]);
  }
  mFunLong->SetBit(ParamsSetBit);
}

void NA6PGenParam::generate()
{
  generatePrimaryVertex(); // will generate if not generated yet, otherwise, will set the origin from the MCHeader of the stack.
  const auto& param = NA6PGenCutParam::Instance();
  auto m = TDatabasePDG::Instance()->GetParticle(mPDGCode)->Mass();
  int nPart = mPoisson ? gRandom->Poisson(mMult) : TMath::Nint(mMult);
  for (int i = 0; i < nPart; i++) {
    int ntrial = 0;
    do {
      double vTrans, vLong, phi, pt2, pt, mt2, mt, pZ, en;
      try {
        vTrans = mFunTrans->GetRandom();
        vLong = mFunLong->GetRandom();
        phi = gRandom->Uniform(0, TMath::TwoPi());
        pt2 = vTrans * vTrans;
        mt2 = pt2;
        pt = vTrans;
        mt = pt;
        if (!mTransIsPt) { // mT
          pt2 -= m * m;
          pt = std::sqrt(pt2);
        } else {
          mt2 += m * m;
          mt = std::sqrt(mt2);
        }
        double pZ = 0.;
        if (mLongIsY) {
          pZ = mt * std::sinh(vLong);
        } else { // eta
          auto tht = 2. * std::atan(std::exp(-vLong));
          pZ = pt / std::tan(tht);
        }
        en = std::sqrt(mt2 + pZ * pZ);
      } catch (std::exception const& e) {
        LOGP(info, "Exception at trial {} of generating particle PDG={}, VTrans={}, vLong={}, error: {}", ntrial, mPDGCode, vTrans, vLong, e.what());
        if (ntrial++ > param.maxTrailsPerParticle) {
          LOGP(fatal, "Failed to generate track PDG:{} after max allowed {} trials, check generator", mPDGCode, ntrial);
        }
        continue;
      }
      ntrial = 0;
      int dummy = 0;
      float sn = std::sin(phi), cs = std::cos(phi);
      getStack()->PushTrack(true, -1, mPDGCode, pt * cs, pt * sn, pZ, en, getOriginX(), getOriginY(), getOriginZ(),
                            0., 0., 0., 0., TMCProcess::kPPrimary, dummy, 1., 0);
    } while (ntrial != 0);
  }
  // register generatot header in the MCHeader
  auto mcHead = getStack()->getEventHeader();
  static std::string info = fmt::format("{}_pdg{}x{}", getName(), mPDGCode, nPart);
  mcHead->getGenHeaders().emplace_back(nPart, 0, mcHead->getNPrimaries(), 0, info);
}
