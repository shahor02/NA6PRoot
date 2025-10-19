// NA6PCCopyright

#include "NA6PGenCocktail.h"
#include "NA6PGenCutParam.h"
#include "NA6PMCStack.h"
#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TMath.h>
#include <cmath>

void NA6PGenCocktail::init()
{
  if (isInitDone()) {
    LOGP(warn, "Generator {} was already initialized", getName());
  }
  NA6PGenerator::init();
  for (auto g : mGenerators) {
    ((NA6PGenerator*)g)->init();
  }
}

void NA6PGenCocktail::setStack(NA6PMCStack* stack)
{
  NA6PGenerator::setStack(stack);
  for (auto g : mGenerators) {
    ((NA6PGenerator*)g)->setStack(stack);
  }
}

void NA6PGenCocktail::generate()
{
  generatePrimaryVertex(); // will generate if not generated yet, otherwise, will set the origin from the MCHeader of the stack.
  // register generatot header in the MCHeader: since this is just a wrapper for other generators, just register the name,
  // the subGenerators will register their own headers with NPrimaries and respective offsets
  auto mcHead = getStack()->getEventHeader();
  mcHead->getGenHeaders().emplace_back(0, 0, 0, 0, getName());
  for (auto g : mGenerators) {
    ((NA6PGenerator*)g)->generate();
  }
}

long NA6PGenCocktail::canGenerateMaxEvents() const
{
  long n = -1;
  for (auto g : mGenerators) {
    auto ng = ((NA6PGenerator*)g)->canGenerateMaxEvents();
    if (ng >= 0 && (n < 0 || n > ng)) {
      n = ng;
    }
  }
  return n;
}
