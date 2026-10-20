#if !defined(__CINT__) || defined(__MAKECINT__)
#include "NA6PGenBeam.h"
#endif

NA6PGenerator* genBeam()
{
  auto gen = new NA6PGenBeam();
  return gen;
}
