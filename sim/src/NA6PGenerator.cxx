// NA6PCCopyright

#include "NA6PGenerator.h"
#include "NA6PMCStack.h"
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
