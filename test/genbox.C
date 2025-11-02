#if !defined(__CINT__) || defined(__MAKECINT__)
#include "NA6PGenBox.h"
#endif

NA6PGenerator* genbox(int npart=5)
{
  auto gen = new NA6PGenBox();
  gen->setNTracks(npart);
  return gen;
}
