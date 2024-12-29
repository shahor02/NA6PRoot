#if !defined(__CINT__) || defined(__MAKECINT__)
#include "NA6PGenBox.h"
#endif

NA6PGenerator* genbox()
{
  auto gen = new NA6PGenBox();
  gen->setNTracks(5);
  return gen;
}
