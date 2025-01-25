// NA6PCCopyright

#ifndef NA6P_GENCOCKTAIL_H
#define NA6P_GENCOCKTAIL_H

#include <TObjArray.h>
#include "NA6PGenerator.h"

class NA6PGenCocktail : public NA6PGenerator
{
 public:
  using NA6PGenerator::NA6PGenerator;
  NA6PGenCocktail() = default;
  ~NA6PGenCocktail() override = default;
  void generate() override;
  void init() override;
  void setStack(NA6PMCStack* stack) override;
  void setOwner(bool v) { mGenerators.SetOwner(v); }
  void addGenerator(NA6PGenerator* g) { mGenerators.AddLast(g); }

  int getNGenerators() const { return mGenerators.GetEntriesFast(); }
  auto& getGenerators() const { return mGenerators; }
  auto getGenerator(int i) const { return mGenerators[i]; }

 protected:
  TObjArray mGenerators;

  ClassDefOverride(NA6PGenCocktail, 1); // Set of used-added generators
};

#endif
