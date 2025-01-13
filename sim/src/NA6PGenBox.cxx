// NA6PCCopyright
#include "NA6PGenBox.h"
#include "NA6PGenCutParam.h"
#include "NA6PMCStack.h"
#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TMath.h>
#include <cmath>

void NA6PGenBox::init()
{
  if (!TDatabasePDG::Instance()->GetParticle(mPDGCode)) {
    LOGP(fatal, "Particle PDG={} is not defined", mPDGCode);
  }
  NA6PGenerator::init();
}

void NA6PGenBox::generate()
{
  generatePrimaryVertex(); // will generate if not generated yet, otherwise, will set the origin from the MCHeader of the stack.
  const auto& param = NA6PGenCutParam::Instance();
  float rnd[3 * mNTracks];
  gRandom->RndmArray(3 * mNTracks, rnd);

  float m = TDatabasePDG::Instance()->GetParticle(mPDGCode)->Mass();
  int nTrials = 0;
  for (int i = 0; i < mNTracks; i++) {
    if (nTrials > param.maxTrailsPerParticle) {
      LOGP(fatal, "Failed to generate track {} with PDG {} after max allowed {} trials, check generator", i, mPDGCode, nTrials);
    }
    int offs = i * 3;
    float eta = param.etaMin + rnd[offs] * (param.etaMax - param.etaMin);
    if (eta <= 0.) {
      i--;
      nTrials++;
      continue;
    }
    float pt = param.ptMin + rnd[offs + 1] * (param.ptMax - param.ptMin);
    float phi = param.phiMin + rnd[offs + 2] * (param.phiMax - param.phiMin);
    float sn = std::sin(phi), cs = std::cos(phi), tht = 2. * std::atan(std::exp(-eta));
    if (tht <= 1e-6 || tht > TMath::Pi() - 1e-6) {
      i--;
      nTrials++;
      continue;
    }
    int dummy = 0;
    float p = pt / std::sin(tht), pz = pt / std::tan(tht), e = std::sqrt(p * p + m * m);
    getStack()->PushTrack(true, -1, mPDGCode, pt * cs, pt * sn, pz, e, getOriginX(), getOriginY(), getOriginZ(),
                          0., 0., 0., 0., TMCProcess::kPPrimary, dummy, 1., 0);
    nTrials = 0;
  }
  // register generatot header in the MCHeader
  auto mcHead = getStack()->getEventHeader();
  static std::string info = fmt::format("{}_pdg{}x{}", getName(), mPDGCode, mNTracks);
  mcHead->getGenHeaders().emplace_back(mNTracks, 0, mcHead->getNPrimaries(), 0, info);
}
