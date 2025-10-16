#if !defined(__CINT__) || defined(__MAKECINT__)
#include "NA6PGenHepMC.h"
#include <TMath.h>
#endif

NA6PGenerator* genHEPMC(const std::string& inpFileName, const std::string& genName = "genHEPMC", int opt = 0)
{
  if (inpFileName.empty()) {
    LOGP(fatal, "HEPMC input file name is not provided");
  }
  auto gen = new NA6PGenHepMC(genName, inpFileName);

  return gen;
}
