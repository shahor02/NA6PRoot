#include "NA6PLayoutParam.h"
#include <TRandom.h>
#include <TMath.h>

void vtxInDump(float& x, float& y, float& z)
{
  static bool first = true;
  static float lAbs = 0, zAbsStart = 0;
  if (first) {
    const auto& layoutPar = NA6PLayoutParam::Instance();
    for (int i = 0; i < layoutPar.nAbsorberSlices - layoutPar.nAbsorberSlicesWall; i++) {
      lAbs += layoutPar.thicknessAbsorber[i];
    }
    zAbsStart = layoutPar.posZStartAbsorber;
    first = false;
  }
  float lambda = 5; // int.length in cm
  do {
    z = -TMath::Log(1. - gRandom->Rndm()) * lambda;
  } while (z > lAbs);
  z += zAbsStart;
  x = y = 0.f;
}
