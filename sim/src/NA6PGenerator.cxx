// NA6PCCopyright

#include "NA6PGenerator.h"
#include "NA6PMCStack.h"
#include "NA6PDetector.h"
#include "NA6PTarget.h"

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
  if (getStack()->isPVGenerated()) {
    setOrigin(mcH->getVX(), mcH->getVY(), mcH->getVZ());
    if (mVerbosity) {
      LOGP(warn, "{}: the primary vertex was already generated, refusing!", getName());
      return;
    }
  }
  static NA6PTarget* tgt = (NA6PTarget*)NA6PDetector::instance().getModule("Target");
  float x, y, z;
  if (!tgt->generateVertex(x, y, z, maxTrials)) {
    LOGP(fatal, "Failed to generate primary vertex");
  }
  setOrigin(x, y, z);
  mcH->setVX(x);
  mcH->setVY(y);
  mcH->setVZ(z);
  if (mVerbosity) {
    LOGP(info, "Generated primary vertex at {:.4f},{:4f},{:4f}", x, y, z);
  }
}
