#if !defined(__CINT__) || defined(__MAKECINT__)
#include "NA6PGenHepMC.h"
#include <TMath.h>
#endif

NA6PGenerator* genHEPMC(const std::string& inpFileName, bool storeDecayed = true, const std::string& genName = "genHEPMC")
{
  if (inpFileName.empty()) {
    LOGP(fatal, "HEPMC input file name is not provided");
  }
  auto gen = new NA6PGenHepMC(genName, inpFileName, storeDecayed);

  return gen;
}
