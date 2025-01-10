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
  static NA6PTarget* tgt = (NA6PTarget*)NA6PDetector::instance().getModule("Target");
  if (getStack()->isPVGenerated()) {
    LOGP(warn, "The primary vertex was already generated, refusing!");
    return;
  }
  float x,y,z;
  if (!tgt->generateVertex(x,y,z, maxTrials)) {
    LOGP(fatal, "Failed to generate primary vertex");
  }
  setOrigin(x,y,z);
  auto mcH = mStack->getEventHeader();
  mcH->setVX(x);
  mcH->setVY(y);
  mcH->setVZ(z);
  if (mVerbosity) {
    LOGP(info, "Generated primary vertex at {:.4f},{:4f},{:4f}", x, y, z);
  }  
}
