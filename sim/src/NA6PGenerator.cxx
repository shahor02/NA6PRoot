// NA6PCCopyright

#include "NA6PGenerator.h"
#include "NA6PMCStack.h"
#include "NA6PDetector.h"
#include "NA6PTarget.h"
#include <TMethodCall.h>

#include <fairlogger/Logger.h>

void NA6PGenerator::init()
{
  if (isInitDone()) {
    LOGP(warn, "Generator {} was already initialized", getName());
  }
  if (!mStack) {
    LOGP(fatal, "MC Stack is not initialized");
  }
}

void NA6PGenerator::generatePrimaryVertex(int maxTrials)
{
  // will generate if not generated yet, otherwise, will set the origin from the MCHeader of the stack
  auto mcH = mStack->getEventHeader();
  if (mOriginSet) {
    mcH->setVX(mOrigin[0]);
    mcH->setVY(mOrigin[1]);
    mcH->setVZ(mOrigin[2]);
    mStack->setPVGenerated(true);
    if (mVerbosity) {
      LOGP(info, "Use vertex imposed via setOrigin {:.4f},{:4f},{:4f}", mcH->getVX(), mcH->getVY(), mcH->getVZ());
    }
    return;
  }
  if (getStack()->isPVGenerated()) {
    setOrigin(mcH->getVX(), mcH->getVY(), mcH->getVZ(), false);
    if (mVerbosity) {
      LOGP(warn, "{}: the primary vertex was already generated, refusing!", getName());
      LOGP(info, "Use existing vertex at {:.4f},{:4f},{:4f}", mcH->getVX(), mcH->getVY(), mcH->getVZ());
    }
    return;
  }
  float x, y, z;
  if (getUserVertexMethod()) {
    const void* args[3] = {&x, &y, &z};
    int ret = 0;
    getUserVertexMethod()->Execute(nullptr, args, 3, &ret);
  } else {
    static NA6PTarget* tgt = (NA6PTarget*)NA6PDetector::instance().getModule("Target");
    if (!tgt->generateVertex(x, y, z, maxTrials)) {
      LOGP(fatal, "Failed to generate primary vertex");
    }
  }
  setOrigin(x, y, z, false);
  mcH->setVX(x);
  mcH->setVY(y);
  mcH->setVZ(z);
  mStack->setPVGenerated(true);
  if (mVerbosity) {
    LOGP(info, "Generated primary vertex at {:.4f},{:4f},{:4f}", x, y, z);
  }
}

void NA6PGenerator::clear()
{
  mOriginSet = false;
}
